# Scientific methods and interpretation

Time-Resolved Emission Workbench implements standard exponential emission-decay models. It estimates the parameters of a stated forward model from photon histograms or compatible nonnegative intensity observations. Assess model suitability and acquisition calibration before interpreting a converged fit.

## Inputs and channel coordinates

Photon histograms are imported as raw nonnegative counts. Background is part of the fitted expectation rather than subtracted from the observed decay before a Poisson fit. Units must be explicit; internal times and reported lifetimes are in ns.

The imported channel coordinate can contain small numerical rounding even when the acquisition channels are equally spaced. The importer accepts agreement with one fixed-width sequence within 0.0001 of a bin over the full axis. It records original coordinates and the maximum adjustment. Photon counts stay in their original channels. This is not count interpolation, rebinning, smoothing, or normalization.

The importer and fitter require a uniform time grid. They reject missing internal channels and genuinely nonuniform acquisition coordinates, including for tail fitting. For supported multi-curve bundles, missing endpoints can be trimmed to a common complete interval and are audited.

Fluorescence upconversion and other pulsed photoluminescence measurements can use the exponential fitter when their sampled observations and kinetics meet these conditions. For upconversion reconvolution, use the measured time response of the gating experiment with matching time spacing and an appropriate time-zero calibration. A single shared IRF and shift may be inadequate when the time response varies across wavelengths. The Python API accepts per-trace response curves with a shared fitted shift.

## Tail model

For a fit beginning at the first included channel time \(t_a\), the expected counts are

\[
\mu_i=B+\sum_{k=1}^{K}A_k\exp[-(t_i-t_a)/\tau_k].
\]

Lifetimes are positive; amplitudes and decay background are nonnegative. Tail amplitudes refer to the first fitted channel. Changing the start changes the amplitude convention even if the lifetime is stable. Select a start beyond the IRF-influenced rise. Predictions before the tail origin are undefined and are exported as blank fitted cells.

## Measured-IRF reconvolution

The intrinsic isolated-pulse response is a positive exponential mixture. A measured IRF is corrected for a specified constant IRF background, restricted to a chosen pulse support if requested, clipped at zero, and normalized to unit sum. The model convolves that response with causal exponentials and adds a decay background.

The discrete exponential kernel uses a bin-average factor

\[
\frac{\tau_k}{\Delta t}\left[1-\exp(-\Delta t/\tau_k)\right].
\]

A reconvolution amplitude refers to the underlying exponential at zero lag, before this factor. Its scaling need not match other software's amplitude conventions. A positive IRF shift delays the response. Fractional shifts use interpolation of the response function.

Reconvolution is evaluated over the available full imported grid before selecting the decay channels included in the objective. Cropping the rise or relevant IRF tail before convolution can discard excitation history. A pulse-support limit is a distinct operation: it selects which measured IRF channels represent the response.

Background/support choices deserve inspection. Clipping a long noisy IRF baseline leaves a positive pedestal; an overly short support removes a real response tail. The measured IRF is treated as known during ordinary fitting, so its finite counting noise and calibration uncertainty are not included in the reported local lifetime errors.

## Objectives and residuals

For independent raw photon counts \(y_i\) and expectations \(\mu_i>0\), the Poisson deviance is

\[
D=2\sum_i[\mu_i-y_i+y_i\log(y_i/\mu_i)],
\]

where the logarithmic contribution is zero for an empty observed channel. The optimizer minimizes signed square-root deviance residuals. This has the same optimum as the Poisson negative log likelihood, up to constants independent of the model parameters. Evaluation near equality uses a numerically stable expression.

Alternative objectives are observed-count weighted least squares, with variance \(\max(y_i,1)\), and uniform least squares. Observed-count weighting can bias sparse fits. For non-counting intensity observations with approximately constant noise variance, uniform least squares is a suitable choice. Photon-count diagnostics such as reduced Poisson deviance require a count model and do not have that statistical interpretation on an arbitrary intensity scale.

Pearson chi-square uses expected-count variance, whereas Poisson deviance and observed-count weighted least squares are different quantities. Their reduced statistics must retain distinct labels. A reduced deviance is not generally expected to equal exactly one for very sparse data.

Degrees of freedom count included observations minus independently free parameters. Fixed lifetimes are not counted as free. A parameter shared globally is counted once. Inspect residual structure and its autocorrelation as well as objective values, particularly near the rise and in the late background.

## Global analysis and amplitude spectra

A global parallel-mixture model shares lifetimes across traces while allowing each trace its own amplitudes and background. The GUI and CLI use a shared measured IRF; the Python API can also accept one IRF column per trace. Timing shift is fixed or shared according to the selected configuration.

For \(M\) traces and \(K\) free shared lifetimes, with \(MK\) amplitudes, \(M\) backgrounds, and one fitted shared shift, there are \(K+MK+M+1\) free parameters before other fixing choices. There is no unique per-trace reduced chi-square allocation for shared parameters; per-trace residual summaries remain useful without asserting independent local degrees of freedom.

The curves \(A_k(\lambda)\) are positive decay-associated amplitude spectra for this model. Exponential areas are proportional to \(A_k(\lambda)\tau_k\), a different quantity. Amplitude fractions are \(A_k/\sum_jA_j\), while integrated-intensity fractions are \(A_k\tau_k/\sum_jA_j\tau_j\). The amplitude-weighted lifetime is \(\sum A_k\tau_k/\sum A_k\); the intrinsic mean photon delay is \(\sum A_k\tau_k^2/\sum A_k\tau_k\).

A component may vanish at some wavelengths while remaining supported by other traces. Aggregate amplitude support and identifiability matter more than one local zero. Raw amplitude spectra are not automatically corrected for detector response, bandwidth, unequal acquisition time, laser drift, or sample changes.

A positive parallel model describes independently decaying contributions. Rise components, energy transfer, solvation, and sequential populations may require signed amplitudes or an explicit kinetic target model; the current fitter does not implement those models.

## Steady-state-scaled TRES and TRANES reconstruction

A pre-exponential spectrum `A_k(lambda)` describes how one fitted component contributes across wavelengths. It is not, by itself, a time-resolved emission spectrum. The measured count surface also includes the IRF, background, count accumulation, and spectral acquisition factors. Spectral reconstruction uses the intrinsic exponential model after separating those effects within the selected fit model.

Let

\[
g(\lambda,t)=\sum_k A_k(\lambda)\exp(-t/\tau_k),\qquad
Q(\lambda)=\sum_k A_k(\lambda)\tau_k.
\]

For nonnegative amplitudes, positive lifetimes and an excitation-origin time coordinate, the isolated-pulse time integral of `g` is `Q`. If `S(lambda)` is the measured steady-state fluorescence intensity on the desired spectral scale, the reconstruction is

\[
I(\lambda,t)=S(\lambda)\frac{g(\lambda,t)}{Q(\lambda)},
\qquad \int_0^\infty I(\lambda,t)\,dt=S(\lambda).
\]

This is the operative formula in the supplied reference notebook. Multiplying raw fitted amplitudes by `S(lambda)` without the `Q(lambda)` denominator would retain arbitrary local count scaling. Multiplying all amplitudes at one wavelength by the same positive constant cancels between numerator and denominator; the normalized kinetic shape and reconstructed spectrum are unchanged. Normalizing the amplitude fractions first is therefore unnecessary, provided the same coefficients appear in both sums.

Steady-state spectra must belong to the same sample and relevant excitation/measurement conditions as the fitted kinetics. Wavelength-dependent response, bandwidth, acquisition time, drift, and blank/background corrections require experimental assessment; this reconstruction does not infer those calibrations. An arbitrary normalized fluorescence spectrum provides a relative spectral scale, not absolute photon intensity or proof of a mechanism.

### Sampling and scaling the steady-state spectrum

The linked workflow takes wavelengths, amplitudes, and paired lifetimes directly from the saved global result. The additional experimental input is its matching steady-state fluorescence intensity. An optional external fitted-table route follows the same calculation but requires explicit lifetime-column pairing.

With already-normalized values supplied at every fitted wavelength, the workbench matches whole wavelength/intensity rows and uses the supplied intensities unchanged. It requires one unique match per fitted wavelength and no extra rows. Row order may differ. Normalized values lie between 0 and 1; a subset of a larger normalized spectrum need not include its peak. This path performs no interpolation or further scaling.

With a full uploaded spectrum, normalization precedes linear interpolation at the global-fit wavelengths. The source grid may be denser or unevenly spaced but must cover every fitted wavelength. Extrapolation is rejected. The source spectrum, normalized values, normalization offset and divisor, matching method, and sampled values are recorded.

The two browser normalization choices are

\[
S_{\mathrm{peak}}(\lambda)=\frac{S(\lambda)}{\max S},\qquad
S_{\mathrm{minmax}}(\lambda)=\frac{S(\lambda)-\min S}{\max S-\min S}.
\]

Peak normalization requires nonnegative intensities and a positive maximum. It multiplies the entire reconstruction by one constant and leaves TRANES unchanged. Min–max normalization allows a finite negative raw offset but rejects a constant spectrum. Its subtraction changes the spectral shape unless the minimum was already zero, so it can change TRANES. This is an explicit selected scaling convention, not a physical blank correction inferred by the software.

The Python API also retains preserve-scale and wavelength-area normalization for compatible existing scripts. Area normalization uses the actual source wavelength spacing. These global spectrum normalizations differ from the per-wavelength kinetic-area denominator `Q(lambda)`.

### Area normalization over wavelength

At each selected delay, TRANES is

\[
T(\lambda,t)=\frac{I(\lambda,t)}{
\int_{\lambda_{\min}}^{\lambda_{\max}} I(\lambda',t)\,d\lambda'},
\qquad \int T(\lambda,t)\,d\lambda=1.
\]

The implementation uses trapezoidal integration on the actual sampled wavelength grid, including nonuniform wavelength spacing. The unit area applies only to the supplied wavelength interval. When wavelengths are in nm, TRANES is a spectral density in inverse nm. A peak-normalized steady-state input gives reconstructed TRES in relative intensity per ns; its absolute scale remains unspecified.

Area normalization over wavenumber is a different density convention and requires the appropriate transformation Jacobian. The current implementation operates over wavelength and does not perform that conversion. Interpreting a plotted curve after changing only its x-axis labels would not preserve the correct density or area.

### Connections between spectral samples

Plot connections are selectable: straight lines, shape-preserving piecewise cubic Hermite interpolation (PCHIP), or an interpolating cubic B-spline. PCHIP is the initial smooth connection. Interpolants pass through the calculated wavelength samples and remain inside the sampled wavelength interval. With fewer than four wavelengths, the B-spline uses the highest supported degree up to three.

The B-spline can overshoot or become negative between nonnegative samples. Such values are displayed without clipping so the selected interpolation remains visible. They are not negative measured fluorescence. Interpolated maxima and crossings depend on curve style and do not establish finer experimental spectral resolution. PCHIP avoids negative overshoot for nonnegative input samples but is still a display interpolation.

Changing the plot style does not alter the calculated sample values, reconstruction equations, or exported numerical matrices. TRANES always uses trapezoidal integration on the actual fitted wavelength grid. The continuous area under a smooth plotted interpolant may differ from that discrete area; smooth curves are not separately renormalized. The selected display convention is recorded in the browser reconstruction export.

### Amplitude reference time and extrapolation

Reconstruction time zero means excitation-origin lag zero of the intrinsic model. Reconvolution amplitudes in this fitter use that convention. A tail-fit amplitude is instead defined at the first included tail channel. If its known delay after excitation is `d`, the conversion is

\[
A_{k,0}(\lambda)=A_{k,d}(\lambda)\exp(d/\tau_k).
\]

A tail result needs that delay explicitly; the software does not assume that an absolute acquisition timestamp is a physical excitation delay. Back-extrapolating short lifetimes from a late tail can magnify uncertainty severely. External amplitude tables also require an explicit, correct origin and lifetime-column mapping. The contributed three-component table uses A1, A2, A3 with 114.15, 17.415 and 4.613 ns respectively; sorting only the lifetimes changes the reconstruction incorrectly.

Reconstruction at times beyond the measured post-excitation interval is a model prediction. Five times the longest lifetime leaves about 0.67 percent of that exponential's initial amplitude, but does not establish a detector noise floor or guarantee negligible total signal. The relative contributions, spectral scale and measurement sensitivity still matter. Uncertainty in lifetimes, amplitudes and steady-state measurements is not automatically propagated to a reconstructed-spectrum confidence interval.

A wavelength with positive steady-state intensity and zero fitted decay area is inconsistent and rejected. Zero spectral area cannot be normalized to TRANES. The implementation also rejects invalid lifetimes, negative delays, non-finite values, negative fitted amplitudes, duplicate wavelengths, or requests outside steady-state spectral coverage. Negative raw steady-state intensities are accepted only when explicit min–max normalization maps them into nonnegative values. The supported lifetime fitter uses nonnegative parallel components; signed rise components and kinetic target models need a different treatment.

## Repeated excitation

For a known period \(T\), the steady periodic contribution of one exponential is proportional to

\[
\frac{\exp(-u/\tau)}{1-\exp(-T/\tau)},\qquad 0\leq u<T.
\]

Periodic mode uses circular IRF convolution only when the imported uniform grid contains one full nonduplicated period, with \(N\Delta t=T\). Ordinary isolated-pulse mode uses padded linear convolution and does not introduce wrap-around. A cropped window, duplicate endpoint, or last timestamp alone cannot establish the required period.

## Optimization and uncertainty

Fitting uses bounded nonlinear least squares with repeatable multiple starting guesses. Larger global problems use sparse Jacobians; this changes resource use rather than the objective. Fitting scopes BLAS threading by default to avoid oversubscription. Advanced callers can set `FitConfig(blas_threads=None)` to retain their own thread policy.

Reported lifetime standard errors are local Gaussian approximations conditional on the model, settings, measured IRF, and calibrated axis. They are not profile-likelihood confidence intervals. Fixed lifetimes have zero reported error because they were held fixed. Ill-conditioned fits can make free-lifetime errors unavailable. Sparse information-matrix estimates use a conservative condition cutoff.

Similar lifetimes, weak components, broad response functions, short windows, and background correlations can give unstable decomposition. Multiple starts help explore alternatives but do not prove a global optimum. Compare model sizes only on the same observations, fitting range, and likelihood when using relative AIC/BIC. Bounds and zero amplitudes complicate formal statistical interpretation.

## Supported scope

Supported features are exponential tail fitting, measured-IRF reconvolution, positive-amplitude shared-lifetime global fits, fixed/initial lifetimes, fit ranges, specified IRF background/support, fitted/fixed decay background, shared timing shifts, full-period correction where applicable, inspectable residual/result exports, visual single-decay and global fitting-range controls, and steady-state-scaled TRES/TRANES reconstruction from compatible fitted amplitudes and a fluorescence spectrum.

Rotational-anisotropy fitting, signed-amplitude target analysis, stretched exponentials, lifetime distributions, prompt-scatter fitting, profile likelihood, IRF-aware uncertainty, acquisition control, and native binary workspace reading are not implemented. Polarized intensity ratios should not be treated as ordinary Poisson decay histograms. A calibrated G-factor and coupled polarized forward model are needed for quantitative rotational analysis.

The measured example fits and the contributed fitted-amplitude/steady-state table are demonstrations under explicit settings, not ground-truth validation. The three-component spectral example is a different dataset from the 51-wavelength measured global-fit example and supplies no missing steady-state spectrum for it. Independent synthetic tests with known parameters test numerical behavior. Neither replaces acquisition-specific assessment of pile-up, afterpulsing, timing response, drift, and detector calibration.

## Manual normalized intensities

Manual table entries follow the exact-wavelength, preserve-scale convention. The fit supplies the wavelength/amplitude/lifetime pairing; entered normalized intensities are used without another scaling or interpolation. A sampled subset can have a maximum below one. Missing entries are not replaced by zero. The exported steady-state values and source metadata record the manual input.
