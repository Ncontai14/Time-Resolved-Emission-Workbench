"""Command-line interface for independent TCSPC lifetime fitting."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import re
import sys

import numpy as np

from . import FitConfig, fit_decay, fit_global
from .io import load_decay, load_tres
from .plotbox import load_plotbox


def _positive_float(value: str) -> float:
    number = float(value)
    if not np.isfinite(number) or not number > 0:
        raise argparse.ArgumentTypeError("must be finite and greater than zero")
    return number


def _json_default(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    return str(value)


def _is_plotbox(path):
    return Path(path).read_bytes().lstrip(b"\xef\xbb\xbf\r\n \t").startswith(b"PlotBox Export")


def parse_curve_indices(value, *, limit=None):
    """Parse explicit zero-based indices and inclusive ranges, e.g. 1,3,5-8."""
    text = value if isinstance(value, str) else " ".join(str(item) for item in value)
    tokens = [token for token in re.split(r"[,\s]+", text.strip()) if token]
    if not tokens:
        raise ValueError("Select at least one decay curve index.")
    selected = []
    for token in tokens:
        match = re.fullmatch(r"(\d+)(?:-(\d+))?", token)
        if match is None:
            raise ValueError("Use zero-based nonnegative indices or inclusive ranges, such as 1,3,5-8.")
        first, last = int(match.group(1)), int(match.group(2) or match.group(1))
        if last < first:
            raise ValueError("Curve-index ranges must run from smaller to larger indices.")
        if limit is not None and last >= limit:
            raise ValueError(f"Curve index {last} is outside this file's 0–{limit-1} range.")
        if last - first > 100000:
            raise ValueError("Curve-index range is too large.")
        selected.extend(range(first, last+1))
    if len(set(selected)) != len(selected):
        raise ValueError("Select each decay curve only once; ranges must not overlap.")
    return selected


def _add_options(parser: argparse.ArgumentParser, *, plotbox=False) -> None:
    parser.add_argument("--irf", type=Path, help="Separate IRF ASCII file: time, raw counts")
    parser.add_argument("--irf-file-curve", type=int, help="Zero-based IRF curve when --irf is a multi-curve PlotBox export")
    parser.add_argument("--irf-time-unit", choices=["ps", "ns", "us", "ms", "s"], help="Unit of a separate plain-text IRF; PlotBox units come from its header")
    parser.add_argument("--out", type=Path, required=True, help="Directory for JSON and CSV fit results")
    unit_help = "Optional confirmation of the PlotBox header unit; explicit mismatch is rejected" if plotbox else "Input unit for plain decay/TRES files (default ns); output settings are always in ns"
    parser.add_argument("--time-unit", choices=["ps", "ns", "us", "ms", "s"], default=None if plotbox else "ns", help=unit_help)
    parser.add_argument("--method", choices=["reconvolution", "tail"], default="reconvolution")
    parser.add_argument("--components", type=int, default=1, help="Number of exponential components (default: 1)")
    parser.add_argument("--statistic", choices=["poisson", "neyman", "uniform"], default="poisson", help="Poisson deviance, observed-count weighted least squares, or unweighted least squares")
    parser.add_argument("--fit-start", type=float, help="First fitted time, in ns after input unit conversion")
    parser.add_argument("--fit-end", type=float, help="Last fitted time, in ns after input unit conversion")
    lifetime_options = parser.add_mutually_exclusive_group()
    lifetime_options.add_argument("--tau", nargs="+", type=_positive_float, help="Initial lifetimes, in ns; supply one per component")
    lifetime_options.add_argument("--fixed-tau", nargs="+", type=_positive_float, help="Fix all lifetimes to these ns values; supply one per component")
    parser.add_argument("--background", type=float, default=0.0, help="Initial or fixed decay background, counts per bin")
    parser.add_argument("--fixed-background", action="store_true", help="Hold decay background at --background")
    parser.add_argument("--irf-background", type=float, default=0.0, help="Constant to subtract from IRF counts before normalization")
    parser.add_argument("--irf-start", type=float, help="First included IRF pulse time, in ns; excludes noisy baseline outside the pulse range")
    parser.add_argument("--irf-end", type=float, help="Last included IRF pulse time, in ns")
    parser.add_argument("--shift", type=float, default=0.0, help="Initial or fixed IRF delay, in ns (positive delays IRF)")
    parser.add_argument("--fit-shift", action="store_true", help="Optimize one shared IRF shift")
    parser.add_argument("--shift-bounds", type=float, nargs=2, default=(-1.0, 1.0), metavar=("LOW", "HIGH"), help="IRF shift bounds in ns")
    parser.add_argument("--period", type=_positive_float, help="Laser repetition period in ns, when periodic excitation is needed")
    parser.add_argument("--starts", type=int, default=3, help="Number of deterministic multistart attempts (default: 3)")
    parser.add_argument("--max-nfev", type=int, default=1000, help="Maximum optimizer evaluations per attempt")
    parser.add_argument("--seed", type=int, default=0, help="Random seed for repeatable starting points")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="python -m tcspc",
        description="Fit raw-count TCSPC decays and global series with shared exponential lifetimes. Internal time units are ns.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)
    decay_parser = subparsers.add_parser("fit", help="Fit one plain-text decay")
    decay_parser.add_argument("--decay", type=Path, required=True, help="ASCII file with time and raw-count columns")
    _add_options(decay_parser)
    global_parser = subparsers.add_parser("global", help="Fit a plain-text time-by-wavelength TRES matrix jointly")
    global_parser.add_argument("--tres", type=Path, required=True, help="ASCII matrix: first column time, following columns raw-count decays")
    global_parser.add_argument("--wavelengths", type=_positive_float, nargs="+", help="Wavelengths in nm, one per count column; alternatively supply a numeric wavelength header")
    _add_options(global_parser)
    inspect_parser = subparsers.add_parser("inspect", help="Inspect PlotBox curves, units and import audit without fitting")
    inspect_parser.add_argument("--input", type=Path, required=True)
    inspect_parser.add_argument("--out", type=Path, help="Optional JSON inventory output file")
    for name, plural in [("curve-fit", False), ("curve-global", True)]:
        bundle_parser = subparsers.add_parser(name, help="Fit explicitly selected curves from a PlotBox ASCII export")
        bundle_parser.add_argument("--input", type=Path, required=True)
        if plural:
            bundle_parser.add_argument("--decay-curves", nargs="+", required=True, help="Zero-based indices or inclusive ranges, e.g. 1-51 or 1,3,5-8; exclude the IRF")
            bundle_parser.add_argument("--wavelengths", type=_positive_float, nargs="+", help="Optional confirmed wavelengths in nm, in selected-curve order")
        else:
            bundle_parser.add_argument("--decay-curve", type=int, required=True, help="Zero-based decay curve index")
        bundle_parser.add_argument("--irf-curve", type=int, help="Zero-based shared IRF curve in the same file")
        _add_options(bundle_parser, plotbox=True)
    return parser


def _inspect(args):
    bundle = load_plotbox(args.input)
    curves = []
    for index, label in enumerate(bundle.labels):
        trace = bundle.counts[:, index]
        valid = np.isfinite(trace)
        peak = int(np.nanargmax(trace)) if np.any(valid) else None
        curves.append({"curve_index": index, "label": label,
                       "valid_bins": int(valid.sum()), "missing_bins": int((~valid).sum()),
                       "total_valid_counts": float(np.nansum(trace)),
                       "peak_counts": float(trace[peak]) if peak is not None else None,
                       "peak_time_ns": float(bundle.time_ns[peak]) if peak is not None else None})
    inventory = {"source_file": str(args.input), "source_sha256": hashlib.sha256(args.input.read_bytes()).hexdigest(), "source_time_unit": bundle.source_time_unit,
                 "time_unit": "ns", "time_bins": len(bundle.time_ns), "curves": curves,
                 "warnings": bundle.warnings, "import_metadata": bundle.metadata,
                 "interpretation": "Curve roles and wavelengths are not inferred from crv indices."}
    content = json.dumps(inventory, indent=2, default=_json_default, allow_nan=False)
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(content + "\n", encoding="utf-8")
    print(content)


def _validate_fit(parser, args):
    is_bundle = args.command.startswith("curve-")
    if is_bundle and args.irf is not None and args.irf_curve is not None:
        parser.error("choose either --irf-curve in the same file or a separate --irf")
    if args.irf is None and (args.irf_file_curve is not None or args.irf_time_unit is not None):
        parser.error("--irf-file-curve and --irf-time-unit require a separate --irf file")
    if args.method == "reconvolution" and args.irf is None and (not is_bundle or args.irf_curve is None):
        parser.error("reconvolution needs --irf or a same-file --irf-curve; use --method tail for a tail fit")
    if args.method == "tail" and args.fit_start is None:
        parser.error("tail fitting requires an explicit --fit-start after the IRF-affected rise")
    if args.components < 1:
        parser.error("--components must be positive")
    for name in ("tau", "fixed_tau"):
        values = getattr(args, name)
        if values is not None and len(values) != args.components:
            parser.error(f"--{name.replace('_', '-')} requires exactly {args.components} value(s)")
    if args.starts < 1 or args.max_nfev < 1:
        parser.error("--starts and --max-nfev must be positive")
    if args.shift_bounds[0] >= args.shift_bounds[1]:
        parser.error("--shift-bounds must satisfy LOW < HIGH")


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.command != "inspect":
        _validate_fit(parser, args)
    try:
        if args.command == "inspect":
            _inspect(args)
            return 0
        config = FitConfig(
            n_components=args.components, method=args.method, statistic=args.statistic,
            fit_start=args.fit_start, fit_end=args.fit_end,
            initial_lifetimes=args.tau, fixed_lifetimes=args.fixed_tau,
            fit_background=not args.fixed_background, background=args.background,
            irf_background=args.irf_background, irf_start=args.irf_start, irf_end=args.irf_end,
            fit_shift=args.fit_shift, shift=args.shift, shift_bounds=tuple(args.shift_bounds),
            repetition_period=args.period, n_starts=args.starts,
            max_nfev=args.max_nfev, seed=args.seed,
        )
        kwargs = {"config": config}
        separate_irf_metadata = None
        if args.irf:
            content = args.irf.read_bytes()
            if content.lstrip(b"\xef\xbb\xbf\r\n \t").startswith(b"PlotBox Export"):
                irf_bundle = load_plotbox(args.irf, time_unit=args.irf_time_unit)
                irf_index = args.irf_file_curve
                if irf_index is None:
                    if len(irf_bundle.labels) != 1:
                        raise ValueError("A multi-curve --irf PlotBox export needs --irf-file-curve INDEX.")
                    irf_index = 0
                it, ic = irf_bundle.select_curves([irf_index])
                kwargs.update(irf=ic[:, 0], irf_time=it)
                separate_irf_metadata = {"source_file": str(args.irf), "source_sha256": hashlib.sha256(content).hexdigest(), "bundle": irf_bundle.metadata,
                                         "selection": irf_bundle.selection_summary([irf_index])}
            else:
                if args.irf_file_curve is not None:
                    raise ValueError("--irf-file-curve applies to a PlotBox export; use a two-column file for a plain-text IRF.")
                irf = load_decay(args.irf, time_unit=args.irf_time_unit or args.time_unit or "ns")
                kwargs.update(irf=irf.counts, irf_time=irf.time_ns)
        wavelengths, import_metadata, curve_indices = None, None, None
        if args.command.startswith("curve-"):
            bundle = load_plotbox(args.input, time_unit=args.time_unit)
            curve_indices = [args.decay_curve] if args.command == "curve-fit" else parse_curve_indices(args.decay_curves, limit=len(bundle.labels))
            if len(set(curve_indices)) != len(curve_indices):
                raise ValueError("Select each decay curve only once.")
            if args.irf_curve in curve_indices:
                raise ValueError("The IRF curve must be excluded from the decay selection.")
            time, counts = bundle.select_curves(curve_indices)
            import_metadata = {"source_file": str(args.input), "source_sha256": hashlib.sha256(args.input.read_bytes()).hexdigest(), "bundle": bundle.metadata,
                               "decay_selection": bundle.selection_summary(curve_indices)}
            if len(time) < len(bundle.time_ns):
                print(f"Import warning: selected curves omit {len(bundle.time_ns)-len(time)} missing endpoint rows; see import_metadata.json.", file=sys.stderr)
            if args.irf_curve is not None:
                irf_time, irf_counts = bundle.select_curves([args.irf_curve])
                kwargs.update(irf=irf_counts[:, 0], irf_time=irf_time)
                import_metadata["irf_selection"] = bundle.selection_summary([args.irf_curve])
            if args.command == "curve-global" and args.wavelengths is not None:
                wavelengths = np.asarray(args.wavelengths)
                if len(wavelengths) != len(curve_indices) or len(np.unique(wavelengths)) != len(wavelengths):
                    raise ValueError("Provide one unique positive wavelength per selected decay curve, in selection order.")
            fitter = fit_decay if args.command == "curve-fit" else fit_global
            result = fitter(time, counts[:, 0] if args.command == "curve-fit" else counts, **kwargs)
            for warning in bundle.warnings:
                print(f"Import warning: {warning}", file=sys.stderr)
        elif args.command == "fit":
            if _is_plotbox(args.decay):
                raise ValueError("PlotBox export detected. Use curve-fit --input FILE --decay-curve INDEX and assign an IRF explicitly; its header supplies the time unit.")
            data = load_decay(args.decay, time_unit=args.time_unit)
            import_metadata = {"decay_time_axis": data.metadata}
            result = fit_decay(data.time_ns, data.counts, **kwargs)
        else:
            if _is_plotbox(args.tres):
                raise ValueError("PlotBox export detected. Use curve-global --input FILE --decay-curves INDICES --irf-curve INDEX; its header supplies the time unit.")
            data = load_tres(args.tres, time_unit=args.time_unit, wavelengths=args.wavelengths)
            import_metadata = {"decay_time_axis": data.metadata}
            wavelengths = data.wavelengths
            result = fit_global(data.time_ns, data.counts, **kwargs)
        if separate_irf_metadata is not None:
            import_metadata = {} if import_metadata is None else import_metadata
            import_metadata["separate_irf"] = separate_irf_metadata
        result.export(args.out)
        if wavelengths is not None:
            with (args.out / "wavelengths.csv").open("w", newline="", encoding="utf-8") as handle:
                writer = csv.writer(handle)
                writer.writerow(["trace", "wavelength_nm"])
                writer.writerows(enumerate(wavelengths, start=1))
        if import_metadata is not None:
            (args.out / "import_metadata.json").write_text(json.dumps(import_metadata, indent=2, default=_json_default, allow_nan=False) + "\n", encoding="utf-8")
        if curve_indices is not None:
            with (args.out / "curve_mapping.csv").open("w", newline="", encoding="utf-8") as handle:
                writer = csv.writer(handle)
                writer.writerow(["trace", "source_curve_index"])
                writer.writerows(enumerate(curve_indices, start=1))
        print("Lifetimes (ns): " + ", ".join(f"{tau:.8g}" for tau in result.lifetimes))
        diagnostics = result.to_dict().get("diagnostics", {})
        summary_keys = ["converged", "statistic", "degrees_of_freedom", "reduced_pearson_chi_square", "reduced_poisson_deviance"]
        print("Diagnostics: " + json.dumps({key: diagnostics[key] for key in summary_keys if key in diagnostics}, allow_nan=False))
        for message in result.warnings:
            print(f"Warning: {message}", file=sys.stderr)
        print(f"Results: {args.out.resolve()}")
        return 0
    except (ValueError, OSError, RuntimeError) as error:
        print(f"Error: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
