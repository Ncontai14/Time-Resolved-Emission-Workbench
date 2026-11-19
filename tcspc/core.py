"""Exponential fitting for uniformly sampled emission decays.

All times are in ns. Use raw photon counts for count-based objectives.
Use an appropriate least-squares objective for other nonnegative signals.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
import csv
import json
from typing import Sequence

import numpy as np
from scipy.optimize import least_squares, nnls
from scipy.signal import fftconvolve
from scipy.special import xlogy
from scipy.sparse import csr_matrix, issparse
from threadpoolctl import threadpool_limits
from ._release import WORKBENCH_NAME, __version__
from .time_axis import channel_grid


# Dense trust-region factorization is fast for individual decays. Global TRES
# Jacobians have only shared columns plus a few local entries in each row.
_SPARSE_TRACE_THRESHOLD = 8
_DENSE_JACOBIAN_ENTRY_LIMIT = 2_000_000


@dataclass
class FitConfig:
    n_components: int = 1
    method: str = "reconvolution"
    statistic: str = "poisson"
    fit_start: float | None = None
    fit_end: float | None = None
    initial_lifetimes: Sequence[float] | None = None
    fixed_lifetimes: Sequence[float | None] | None = None
    fit_background: bool = True
    background: float = 0.0
    irf_background: float = 0.0
    irf_start: float | None = None
    irf_end: float | None = None
    fit_shift: bool = False
    shift: float = 0.0
    shift_bounds: tuple[float, float] = (-1.0, 1.0)
    repetition_period: float | None = None
    n_starts: int = 3
    max_nfev: int = 1000
    seed: int = 0
    blas_threads: int | None = 1


def _clean_json(value):
    if isinstance(value, dict):
        return {str(k): _clean_json(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_clean_json(v) for v in value]
    if isinstance(value, np.ndarray):
        return _clean_json(value.tolist())
    if isinstance(value, (np.floating, float)):
        return float(value) if np.isfinite(value) else None
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.bool_):
        return bool(value)
    return value


@dataclass
class FitResult:
    lifetimes: np.ndarray
    lifetime_stderr: np.ndarray
    amplitudes: np.ndarray
    backgrounds: np.ndarray
    shift: float
    time: np.ndarray
    observed: np.ndarray
    fitted: np.ndarray
    residuals: np.ndarray
    fit_mask: np.ndarray
    diagnostics: dict
    warnings: list[str]
    config: FitConfig

    def to_dict(self):
        a, tau = self.amplitudes, self.lifetimes
        fractions = a / np.maximum(a.sum(axis=1, keepdims=True), 1e-300)
        areas = a * tau
        area_fractions = areas / np.maximum(areas.sum(axis=1, keepdims=True), 1e-300)
        return _clean_json({
            "software": f"{WORKBENCH_NAME} {__version__}", "time_unit": "ns",
            "lifetimes_ns": tau, "lifetime_stderr_ns": self.lifetime_stderr,
            "amplitudes": a, "backgrounds_counts_per_bin": self.backgrounds,
            "irf_shift_ns": self.shift,
            "amplitude_fractions": fractions, "integrated_intensity_fractions": area_fractions,
            "amplitude_weighted_lifetime_ns": (fractions * tau).sum(axis=1),
            "intensity_weighted_lifetime_ns": (area_fractions * tau).sum(axis=1),
            "diagnostics": self.diagnostics, "warnings": self.warnings,
            "config": asdict(self.config),
            "amplitude_convention": (
                "Signal in input units at first included fit channel (tail)."
                if self.config.method == "tail" else
                "Underlying exponential amplitude at lag zero, before bin averaging; IRF has unit sum."
            ),
            "uncertainty_note": "Local Gaussian standard errors, conditional on the measured IRF; not confidence intervals.",
        })

    def export(self, directory):
        path = Path(directory)
        path.mkdir(parents=True, exist_ok=True)
        (path / "result.json").write_text(json.dumps(self.to_dict(), indent=2, allow_nan=False) + "\n")
        with (path / "curves.csv").open("w", newline="") as stream:
            writer = csv.writer(stream)
            header = ["time_ns", "included"]
            for j in range(self.observed.shape[1]):
                header.extend([f"counts_{j+1}", f"fit_{j+1}", f"residual_{j+1}"])
            writer.writerow(header)
            for i, t in enumerate(self.time):
                row = [t, int(self.fit_mask[i])]
                for j in range(self.observed.shape[1]):
                    row.extend([self.observed[i, j],
                                self.fitted[i, j] if np.isfinite(self.fitted[i, j]) else "",
                                self.residuals[i, j] if self.fit_mask[i] else ""])
                writer.writerow(row)
        with (path / "amplitudes.csv").open("w", newline="") as stream:
            writer = csv.writer(stream)
            writer.writerow(["trace", "component", "lifetime_ns", "amplitude", "A_times_tau"])
            for j, aa in enumerate(self.amplitudes):
                for k, a in enumerate(aa):
                    writer.writerow([j + 1, k + 1, self.lifetimes[k], a, a * self.lifetimes[k]])
        return path


def _time_grid(time):
    t, audit = channel_grid(time, min_bins=4)
    return t, audit["bin_width_ns"]


def _histogram(counts, n, name):
    y = np.asarray(counts, dtype=float)
    if y.ndim == 1:
        y = y[:, None]
    if y.ndim != 2 or y.shape[0] != n or y.shape[1] < 1:
        raise ValueError(f"{name} must have shape [time] or [time, trace].")
    if not np.all(np.isfinite(y)) or np.any(y < 0):
        raise ValueError(f"{name} must contain finite nonnegative signal values.")
    if np.any(y.sum(axis=0) <= 0):
        raise ValueError(f"Each {name} trace must contain a positive signal value.")
    return y


def _prepare_irf(t, dt, irf, irf_time, irf_background, m, period, irf_start=None, irf_end=None):
    if irf is None:
        raise ValueError("Reconvolution requires a measured IRF; choose tail mode to fit without it.")
    raw = np.asarray(irf, dtype=float)
    if raw.ndim == 1:
        raw = raw[:, None]
    if raw.ndim != 2 or raw.shape[1] not in (1, m) or not np.all(np.isfinite(raw)) or np.any(raw < 0):
        raise ValueError("IRF must contain nonnegative signal values, shared or one column per trace.")
    if not np.isfinite(irf_background) or irf_background < 0:
        raise ValueError("IRF background must be finite and nonnegative.")
    ti = t if irf_time is None else np.asarray(irf_time, dtype=float)
    ti, dti = _time_grid(ti)
    if len(ti) != len(raw):
        raise ValueError("IRF and IRF time have different lengths.")
    if not np.isclose(dti, dt, rtol=1e-5, atol=dt * 1e-7):
        raise ValueError("IRF and decay bin widths must match. Rebin to a common calibrated grid first.")
    if period is not None and (len(ti) != len(t) or not np.allclose(ti, t, atol=dt*1e-5, rtol=0)):
        raise ValueError("Periodic reconvolution needs the IRF and decay on the same full-period grid.")
    corrected = np.maximum(raw - irf_background, 0.0)
    ir = np.column_stack([np.interp(t, ti, corrected[:, j], left=0, right=0)
                          for j in range(raw.shape[1])])
    if irf_start is not None:
        if not np.isfinite(irf_start):
            raise ValueError("IRF start must be finite.")
        ir[t < irf_start] = 0
    if irf_end is not None:
        if not np.isfinite(irf_end):
            raise ValueError("IRF end must be finite.")
        ir[t > irf_end] = 0
    if np.any(ir.sum(axis=0) <= 0):
        raise ValueError("IRF has no positive area after background correction and grid alignment.")
    return ir / ir.sum(axis=0, keepdims=True)


def _shift_irf(t, ir, shift, period):
    if period is None:
        shifted = np.column_stack([np.interp(t - shift, t, ir[:, j], left=0, right=0)
                                   for j in range(ir.shape[1])])
    else:
        query = (t - t[0] - shift) % period + t[0]
        tp = np.r_[t, t[0] + period]
        shifted = np.column_stack([np.interp(query, tp, np.r_[ir[:, j], ir[0, j]])
                                   for j in range(ir.shape[1])])
    totals = shifted.sum(axis=0)
    if np.any(totals <= 1e-15):
        raise ValueError("The time shift moves the entire IRF out of the acquisition window.")
    return shifted / totals


def _basis(t, dt, tau, ir, shift, method, reference, period):
    if method == "tail":
        # The tail model has no defined rise before its reference channel.
        # Do not exponentiate backward extrapolations: they can overflow despite
        # being outside the fit window. Public outputs mark this domain NaN.
        lag = np.maximum(t - reference, 0)
        return np.exp(-lag[:, None]/tau)[:, None, :]
    h = _shift_irf(t, ir, shift, period)
    lag = np.arange(len(t), dtype=float) * dt
    # Average continuous exponential over each histogram channel.
    correction = -np.expm1(-dt / tau) * tau / dt
    kernel = np.exp(-lag[:, None] / tau) * correction
    if period is None:
        return np.stack([fftconvolve(h, kernel[:, k, None], mode="full", axes=0)[:len(t)]
                         for k in range(len(tau))], axis=2)
    kernel /= -np.expm1(-period / tau)
    transformed_ir = np.fft.rfft(h, axis=0)
    return np.stack([np.fft.irfft(transformed_ir * np.fft.rfft(kernel[:, k])[:, None],
                                 n=len(t), axis=0) for k in range(len(tau))], axis=2)


def predict(time, lifetimes, amplitudes, backgrounds=0, irf=None, irf_time=None,
            shift=0, irf_background=0, method="reconvolution", fit_start=None,
            repetition_period=None, irf_start=None, irf_end=None):
    """Evaluate the forward model. Output shape is [time, trace].

    Tail predictions are undefined (NaN) before the first included channel.
    """
    t, dt = _time_grid(time)
    tau = np.asarray(lifetimes, dtype=float)
    a = np.atleast_2d(np.asarray(amplitudes, dtype=float))
    if tau.ndim != 1 or len(tau) == 0 or not np.all(np.isfinite(tau)) or np.any(tau <= 0):
        raise ValueError("Lifetimes must be finite and positive.")
    if a.ndim != 2 or a.shape[1] != len(tau) or not np.all(np.isfinite(a)) or np.any(a < 0):
        raise ValueError("Amplitudes must be finite nonnegative values with shape [trace, component].")
    if method not in ("reconvolution", "tail"):
        raise ValueError("Method must be 'reconvolution' or 'tail'.")
    if not np.isfinite(shift):
        raise ValueError("IRF shift must be finite.")
    bg = np.broadcast_to(np.asarray(backgrounds, dtype=float), (a.shape[0],))
    if not np.all(np.isfinite(bg)) or np.any(bg < 0):
        raise ValueError("Decay backgrounds must be finite and nonnegative.")
    _check_period(t, dt, method, repetition_period)
    if method == "tail" and shift != 0:
        raise ValueError("A tail fit has no IRF time-shift parameter.")
    if fit_start is not None and (not np.isfinite(fit_start) or fit_start > t[-1]):
        raise ValueError("Fit start must be finite and intersect the time window.")
    reference = t[0] if fit_start is None else t[np.searchsorted(t, fit_start)]
    ir = None if method == "tail" else _prepare_irf(t, dt, irf, irf_time, irf_background, a.shape[0], repetition_period, irf_start, irf_end)
    basis = _basis(t, dt, tau, ir, shift, method, reference, repetition_period)
    prediction = np.maximum(np.einsum("njk,jk->nj", np.broadcast_to(basis, (len(t), a.shape[0], len(tau))), a) + bg, 1e-12)
    if method == "tail":
        prediction[t < reference] = np.nan
    return prediction


def _check_period(t, dt, method, period):
    if period is not None:
        if method != "reconvolution":
            raise ValueError("Repetition-period correction is supported only in reconvolution mode.")
        if not np.isfinite(period) or period <= 0 or not np.isclose(len(t) * dt, period, rtol=1e-5, atol=dt * 1e-5):
            raise ValueError("Periodic reconvolution requires a complete period: number of bins × bin width = repetition period (no duplicate endpoint).")


def _residual(mu, y, statistic):
    mu = np.maximum(mu, 1e-12)
    if statistic == "uniform":
        return y - mu, -np.ones_like(mu)
    if statistic == "neyman":
        weight = 1 / np.sqrt(np.maximum(y, 1))
        return (y - mu) * weight, -weight
    d2 = 2 * (mu - y + xlogy(y, y / mu))
    positive = y > 0
    q = np.zeros_like(y)
    q[positive] = (mu[positive] - y[positive]) / y[positive]
    close = positive & (np.abs(q) < 1e-4)
    z = q[close]
    d2[close] = 2 * y[close] * z*z * (0.5 - z/3 + z*z/4 - z*z*z/5)
    r = np.sign(y - mu) * np.sqrt(np.maximum(d2, 0))
    derivative = np.empty_like(mu)
    near = np.abs(y - mu) < 1e-6 * np.maximum(mu, 1e-12)
    derivative[near] = -1 / np.sqrt(mu[near])
    away = ~near
    derivative[away] = (1 - y[away] / mu[away]) / r[away]
    return r, derivative


def fit_decay(time_ns, counts, irf=None, irf_time=None, config=None):
    """Fit one histogram. For multiple histograms use fit_global."""
    y = np.asarray(counts)
    if y.ndim == 2 and y.shape[1] != 1:
        raise ValueError("fit_decay expects one trace. Use fit_global for a matrix.")
    return fit_global(time_ns, counts, irf=irf, irf_time=irf_time, config=config)


def fit_global(time_ns, counts, irf=None, irf_time=None, config=None):
    """Fit shared exponential lifetimes with trace-specific amplitudes/backgrounds."""
    cfg = FitConfig() if config is None else config
    if cfg.blas_threads is not None:
        if isinstance(cfg.blas_threads, bool) or not isinstance(cfg.blas_threads, (int, np.integer)) or cfg.blas_threads < 1:
            raise ValueError("blas_threads must be a positive integer or None.")
        # Sparse fits make many small BLAS calls; spawning dozens of threads per
        # call can be much slower than fitting. Scope the limit to this operation.
        with threadpool_limits(limits=int(cfg.blas_threads), user_api="blas"):
            return _fit_global_impl(time_ns, counts, irf, irf_time, cfg)
    return _fit_global_impl(time_ns, counts, irf, irf_time, cfg)


def _fit_global_impl(time_ns, counts, irf, irf_time, cfg):
    t, time_audit = channel_grid(time_ns, min_bins=4)
    dt = time_audit["bin_width_ns"]
    y = _histogram(counts, len(t), "Decay")
    m = y.shape[1]
    k = cfg.n_components
    if not isinstance(k, (int, np.integer)) or not 1 <= k <= 6:
        raise ValueError("Choose between one and six exponential components.")
    if cfg.method not in ("tail", "reconvolution") or cfg.statistic not in ("poisson", "neyman", "uniform"):
        raise ValueError("Unknown method or statistic.")
    if cfg.n_starts < 1 or cfg.max_nfev < 1:
        raise ValueError("n_starts and max_nfev must be positive.")
    if not np.isfinite(cfg.background) or cfg.background < 0:
        raise ValueError("Decay background must be finite and nonnegative.")
    if not np.isfinite(cfg.shift):
        raise ValueError("IRF shift must be finite.")
    if cfg.method == "tail" and (cfg.fit_shift or cfg.shift != 0):
        raise ValueError("Disable IRF shift when using a tail fit.")
    if cfg.fit_start is not None and not np.isfinite(cfg.fit_start):
        raise ValueError("Fit start must be finite.")
    if cfg.fit_end is not None and not np.isfinite(cfg.fit_end):
        raise ValueError("Fit end must be finite.")
    mask = np.ones(len(t), dtype=bool)
    if cfg.fit_start is not None:
        mask &= t >= cfg.fit_start
    if cfg.fit_end is not None:
        mask &= t <= cfg.fit_end
    if mask.sum() < 4:
        raise ValueError("The fitting window must include at least four time bins.")
    if np.any(y[mask].sum(axis=0) <= 0):
        raise ValueError("Each trace must contain counts inside the fit window.")
    _check_period(t, dt, cfg.method, cfg.repetition_period)
    ir = None if cfg.method == "tail" else _prepare_irf(t, dt, irf, irf_time, cfg.irf_background, m, cfg.repetition_period, cfg.irf_start, cfg.irf_end)
    reference = t[mask][0]
    span = len(t) * dt
    tau_bounds = (dt * 0.05, span * 100)
    fixed = np.full(k, np.nan)
    if cfg.fixed_lifetimes is not None:
        if len(cfg.fixed_lifetimes) != k:
            raise ValueError("fixed_lifetimes must have one value (or None) per component.")
        for i, val in enumerate(cfg.fixed_lifetimes):
            if val is not None:
                if not np.isfinite(val) or val <= 0:
                    raise ValueError("Fixed lifetimes must be finite and positive.")
                fixed[i] = val
    free = np.flatnonzero(np.isnan(fixed))
    nf = len(free)
    if cfg.initial_lifetimes is not None:
        tau0 = np.asarray(cfg.initial_lifetimes, dtype=float).copy()
        if tau0.shape != (k,) or np.any(tau0 <= 0) or not np.all(np.isfinite(tau0)):
            raise ValueError("initial_lifetimes must contain one finite positive value per component.")
        if np.any((tau0[free] <= tau_bounds[0]) | (tau0[free] >= tau_bounds[1])):
            raise ValueError(f"Initial free lifetimes must be inside {tau_bounds} ns.")
    else:
        summed = y.sum(axis=1)
        base = np.quantile(summed, 0.1)
        signal = np.maximum(summed - base, 0)
        if cfg.method == "tail":
            signal = signal * mask
            guess = np.sum((t - reference) * signal) / max(signal.sum(), 1)
        else:
            guess = np.sum(t * signal) / max(signal.sum(), 1) - np.mean(np.sum(t[:, None] * ir, axis=0)) - cfg.shift
        guess = float(np.clip(guess, dt * 2, span / 2))
        tau0 = np.array([guess]) if k == 1 else guess * np.geomspace(0.25, 2.5, k)
    tau0[np.isfinite(fixed)] = fixed[np.isfinite(fixed)]
    ns = int(cfg.fit_shift)
    amp_start = nf + ns
    bg_start = amp_start + m*k
    npar = bg_start + (m if cfg.fit_background else 0)
    dof = int(mask.sum() * m - npar)
    if dof <= 0:
        raise ValueError("The fit window has too few observations for the number of free parameters.")
    if ns and (len(cfg.shift_bounds) != 2 or not np.all(np.isfinite(cfg.shift_bounds)) or cfg.shift_bounds[0] >= cfg.shift_bounds[1] or not cfg.shift_bounds[0] <= cfg.shift <= cfg.shift_bounds[1]):
        raise ValueError("Set finite increasing shift bounds that contain the initial shift.")
    lower = np.r_[np.full(nf, np.log(tau_bounds[0])),
                  [cfg.shift_bounds[0]] if ns else [], np.zeros(m*k + (m if cfg.fit_background else 0))]
    upper = np.r_[np.full(nf, np.log(tau_bounds[1])),
                  [cfg.shift_bounds[1]] if ns else [], np.full(m*k + (m if cfg.fit_background else 0), np.inf)]
    yy = y[mask]
    nt = len(yy)
    sparse_jacobian = (m >= _SPARSE_TRACE_THRESHOLD
                       or nt*m*npar >= _DENSE_JACOBIAN_ENTRY_LIMIT)
    # CSR row structure is fixed throughout a fit. Each observation depends on
    # the global lifetimes/shift and only its own trace's amplitudes/background.
    if sparse_jacobian:
        n_shared = nf + ns
        n_entries_per_row = n_shared + k + int(cfg.fit_background)
        trace_columns = np.column_stack([
            np.broadcast_to(np.arange(n_shared), (m, n_shared)),
            amp_start + np.arange(m)[:, None]*k + np.arange(k),
            (bg_start + np.arange(m))[:, None]
            if cfg.fit_background else np.empty((m, 0), dtype=int),
        ])
        csr_indices = np.tile(trace_columns, (nt, 1)).ravel()
        csr_indptr = np.arange(nt*m + 1, dtype=np.int64)*n_entries_per_row

    def unpack(p):
        tau = fixed.copy()
        tau[free] = np.exp(p[:nf])
        shift = p[nf] if ns else cfg.shift
        a = p[amp_start:bg_start].reshape(m, k)
        bg = p[bg_start:] if cfg.fit_background else np.full(m, cfg.background)
        return tau, shift, a, bg

    def basis(tau, shift):
        bb = _basis(t, dt, tau, ir, shift, cfg.method, reference, cfg.repetition_period)
        return np.broadcast_to(bb, (len(t), m, k))

    cached_p = None
    cached = None
    def evaluate(p):
        nonlocal cached_p, cached
        if cached_p is not None and np.array_equal(p, cached_p):
            return cached
        tau, shift, a, bg = unpack(p)
        bb = basis(tau, shift)[mask]
        mu = np.maximum(np.einsum("nmk,mk->nm", bb, a) + bg, 1e-12)
        r, deriv = _residual(mu, yy, cfg.statistic)
        cached_p, cached = p.copy(), (tau, shift, a, bb, mu, r, deriv)
        return cached

    def fun(p):
        return evaluate(p)[5].ravel()

    def jac(p):
        tau, shift, a, bb, mu, r, deriv = evaluate(p)
        if sparse_jacobian:
            # Build only actual entries; never allocate [observations, npar].
            entries = np.empty((nt*m, n_entries_per_row))
        else:
            jj = np.zeros((nt*m, npar))
        eps = 1e-5
        for col, idx in enumerate(free):
            tp, tm = tau.copy(), tau.copy()
            tp[idx] *= np.exp(eps)
            tm[idx] *= np.exp(-eps)
            db = (basis(tp, shift)[mask, :, idx] - basis(tm, shift)[mask, :, idx]) / (2*eps)
            values = (deriv * db * a[:, idx]).ravel()
            if sparse_jacobian:
                entries[:, col] = values
            else:
                jj[:, col] = values
        if ns:
            h = max(dt * 1e-4, 1e-7)
            db = (basis(tau, shift+h)[mask] - basis(tau, shift-h)[mask]) / (2*h)
            values = (deriv * np.einsum("nmk,mk->nm", db, a)).ravel()
            if sparse_jacobian:
                entries[:, nf] = values
            else:
                jj[:, nf] = values
        if sparse_jacobian:
            entries[:, n_shared:n_shared+k] = (deriv[:, :, None]*bb).reshape(nt*m, k)
            if cfg.fit_background:
                entries[:, -1] = deriv.ravel()
            return csr_matrix((entries.ravel(), csr_indices, csr_indptr),
                              shape=(nt*m, npar))
        for j in range(m):
            rows = np.arange(nt)*m + j
            jj[rows, amp_start+j*k:amp_start+(j+1)*k] = deriv[:, j, None] * bb[:, j, :]
            if cfg.fit_background:
                jj[rows, bg_start+j] = deriv[:, j]
        return jj

    rng = np.random.default_rng(cfg.seed)
    best = None
    runs = []
    for trial in range(cfg.n_starts):
        tau = tau0.copy()
        if trial:
            tau[free] *= np.exp(rng.normal(0, 0.65, nf))
            tau[free] = np.clip(tau[free], tau_bounds[0]*1.01, tau_bounds[1]*0.99)
        bb = basis(tau, cfg.shift)[mask]
        aa, backgrounds = [], []
        for j in range(m):
            design = np.column_stack([bb[:, j], np.ones(nt)]) if cfg.fit_background else bb[:, j]
            target = yy[:, j] if cfg.fit_background else yy[:, j] - cfg.background
            weights = np.ones(nt) if cfg.statistic == "uniform" else 1 / np.sqrt(np.maximum(yy[:, j], 1))
            coefficients, _ = nnls(design * weights[:, None], target * weights, maxiter=max(100, 10*(k+1)))
            aa.extend(np.maximum(coefficients[:k], 1e-6))
            if cfg.fit_background:
                backgrounds.append(max(coefficients[-1], 1e-6))
        p0 = np.r_[np.log(tau[free]), [cfg.shift] if ns else [], aa, backgrounds]
        solver_options = ({"tr_solver": "lsmr",
                           "tr_options": {"atol": 1e-10, "btol": 1e-10,
                                          "maxiter": max(20, 2*npar)}}
                          if sparse_jacobian else {})
        run = least_squares(fun, p0, jac=jac, bounds=(lower, upper), x_scale="jac",
                            max_nfev=cfg.max_nfev, ftol=1e-9, xtol=1e-9, gtol=1e-8,
                            **solver_options)
        objective = float(np.dot(run.fun, run.fun))
        runs.append({"objective": objective, "success": bool(run.success), "nfev": run.nfev})
        if best is None or objective < 2*best.cost:
            best = run
    tau, shift, a, bg = unpack(best.x)
    bb = basis(tau, shift)
    fitted = np.maximum(np.einsum("nmk,mk->nm", bb, a) + bg, 1e-12)
    if cfg.method == "tail":
        fitted[t < reference] = np.nan
    residuals = np.full_like(y, np.nan)
    residuals[mask] = _residual(fitted[mask], yy, cfg.statistic)[0]
    notes = []
    if not best.success:
        notes.append("Optimizer did not converge: " + best.message)
    if cfg.statistic == "poisson" and not np.allclose(y, np.round(y), atol=1e-8, rtol=0):
        notes.append("Input has fractional counts. Poisson uncertainties assume raw, unscaled photon counts.")
    if cfg.statistic == "neyman":
        notes.append("Observed-count weighting can bias fits at low counts; compare with Poisson likelihood.")
    if cfg.method == "tail":
        peak_time = t[np.argmax(y.sum(axis=1))]
        if reference <= peak_time:
            notes.append("Tail window begins at or before the observed peak; move its start beyond the IRF-influenced region.")
    if cfg.method == "reconvolution":
        if cfg.repetition_period is None and float(ir[:2].sum() + ir[-2:].sum()) / ir.shape[1] > 0.01:
            notes.append("IRF has appreciable counts at acquisition boundaries; check background, clipping, or periodic excitation.")
        if cfg.irf_background == 0:
            notes.append("IRF background is set to zero. A measured IRF baseline should be supplied independently.")
        if cfg.repetition_period is None:
            shifted_mass = np.array([np.interp(t-shift, t, ir[:, j], left=0, right=0).sum() for j in range(ir.shape[1])])
            if np.any(shifted_mass < 0.99):
                notes.append("Time shifting truncates more than 1% of the IRF; acquire a wider time window.")
    if ns and (shift-cfg.shift_bounds[0] < dt*0.1 or cfg.shift_bounds[1]-shift < dt*0.1):
        notes.append("IRF shift is near a bound; review alignment and shift bounds.")
    sorted_indices = np.argsort(tau)
    if k > 1 and np.any(np.diff(tau[sorted_indices]) / tau[sorted_indices][:-1] < 0.15):
        notes.append("Lifetimes are closely spaced; individual components may not be identifiable.")
    # A spectral component can legitimately vanish at selected wavelengths.
    # Assess support over the complete global fit, rather than each local trace.
    total_amplitude = a.sum(axis=0)
    if np.any(total_amplitude < max(float(total_amplitude.max()), 1)*1e-5):
        notes.append("A component has negligible amplitude across the fitted traces; it may be unsupported by these data.")
    if np.any((tau[free] < tau_bounds[0]*1.05) | (tau[free] > tau_bounds[1]*0.95)):
        notes.append("A lifetime approaches the numerical bounds; check model identifiability and time calibration.")
    stderr = np.full(k, np.nan)
    stderr[np.isfinite(fixed)] = 0.0
    if issparse(best.jac):
        # The dense object here is only [npar, npar], never the full TRES
        # Jacobian. Gram matrices square the condition number, so reject
        # uncertain eigenvalues more conservatively than the dense SVD route.
        information = (best.jac.T @ best.jac).toarray()
        norms = np.sqrt(np.maximum(np.diag(information), 0))
        safe_norms = np.maximum(norms, 1e-300)
        information /= safe_norms[:, None]
        information /= safe_norms[None, :]
        eigenvalues, vectors = np.linalg.eigh((information + information.T)*0.5)
        maximum = max(float(eigenvalues[-1]), 1e-300)
        rank = int(np.sum(eigenvalues > maximum*1e-12))
        condition = float(np.sqrt(maximum/max(float(eigenvalues[0]), 1e-300)))
        reliable = rank == npar and condition < 1e6
        if reliable:
            variance = np.sum(vectors*vectors/eigenvalues, axis=1)/safe_norms**2
        covariance_method = "scaled sparse Gram matrix; condition cutoff 1e6"
    else:
        norms = np.linalg.norm(best.jac, axis=0)
        safe_norms = np.maximum(norms, 1e-300)
        normalized = best.jac / safe_norms
        _, singular, vt = np.linalg.svd(normalized, full_matrices=False)
        rank = int(np.sum(singular > singular[0]*1e-10))
        condition = float(singular[0] / max(singular[-1], 1e-300))
        reliable = rank == npar and condition < 1e8
        if reliable:
            variance = np.sum(vt.T*vt.T/singular**2, axis=1)/safe_norms**2
        covariance_method = "scaled dense Jacobian SVD; condition cutoff 1e8"
    if reliable:
        if cfg.statistic == "uniform":
            variance *= 2*best.cost/dof
        stderr[free] = tau[free] * np.sqrt(np.maximum(variance[:nf], 0))
    else:
        notes.append("Fit Jacobian is ill-conditioned; free-lifetime standard errors are unavailable.")
    for error, val in zip(stderr, tau):
        if np.isfinite(error) and error > 0.5 * val:
            notes.append("A lifetime has a large local uncertainty; compare simpler models and different fit windows.")
            break
    pearson = float(np.sum((yy-fitted[mask])**2 / fitted[mask]))
    dev = float(np.sum(_residual(fitted[mask], yy, "poisson")[0]**2))
    objective = float(2*best.cost)
    score = dev if cfg.statistic == "poisson" else objective
    acf = []
    for j in range(m):
        rr = residuals[mask, j] - residuals[mask, j].mean()
        denom = np.dot(rr, rr)
        acf.append([1.0] + [float(np.dot(rr[:-lag], rr[lag:])/denom) if denom > 0 else 0.0
                            for lag in range(1, min(50, len(rr)-1)+1)])
    diagnostics = {
        "time_axis": time_audit,
        "converged": bool(best.success), "message": best.message,
        "statistic": cfg.statistic, "objective": objective,
        "degrees_of_freedom": dof, "n_observations": int(yy.size), "n_free_parameters": npar,
        "pearson_chi_square": pearson, "reduced_pearson_chi_square": pearson/dof,
        "poisson_deviance": dev, "reduced_poisson_deviance": dev/dof,
        "aic_relative": (score + 2*npar) if cfg.statistic == "poisson" else None,
        "bic_relative": (score + npar*np.log(yy.size)) if cfg.statistic == "poisson" else None,
        "information_criterion_note": "Poisson likelihood relative constants omitted; compare only identical observations and windows.",
        "jacobian_rank": rank, "scaled_jacobian_condition": condition,
        "jacobian_storage": "sparse CSR" if sparse_jacobian else "dense",
        "covariance_method": covariance_method,
        "fit_start_ns": float(reference), "fit_end_ns": float(t[mask][-1]),
        "bin_width_ns": dt, "starts": runs,
        "residual_autocorrelation": acf,
        "uncertainty": "Local Gaussian approximation; fixed lifetimes have zero error. IRF noise and boundary effects excluded.",
    }
    return FitResult(tau[sorted_indices], stderr[sorted_indices], a[:, sorted_indices], bg,
                     float(shift), t, y, fitted, residuals, mask, diagnostics, notes, cfg)
