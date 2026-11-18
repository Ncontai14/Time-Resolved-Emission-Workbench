# Independent review of the supplied spectral reconstruction notebook

The supplied notebook was inspected as reference code. Its widget callbacks and
file-writing commands were not executed. The notebook JSON parses as format 4.5,
with two cells (one populated code cell and one empty code cell). The saved output
is a widget display, not saved numerical TRES results. A full `nbformat` schema
validator is not installed in this review environment. The calculation below was
instead evaluated independently from the supplied numeric table.

## Meaning of the calculation

For wavelength λ, let the fitted intrinsic decay model be

    f(λ,t) = Σk A(λ,k) exp(−t/τk),  t ≥ 0.

Its continuous time integral is Σk A(λ,k)τk. If S(λ) is the steady-state
fluorescence spectrum measured under comparable conditions, the notebook forms

    TRES(λ,t) = S(λ) f(λ,t) / Σk A(λ,k)τk.

Consequently, the integral of reconstructed TRES over time is S(λ). This is the
critical normalization: multiplying the raw pre-exponentials by S(λ) alone would
leave a spurious, wavelength-dependent factor Σ Aτ in that integral. The ratio
also cancels any positive common amplitude scale applied to a particular
wavelength, such as a different count target or acquisition duration, provided
that scale affects all fitted components at that wavelength equally.

TRANES is the reconstructed TRES divided by its spectral area separately at each
delay:

    TRANES(λ,t) = TRES(λ,t) / ∫ TRES(λ,t) dλ.

The provided wavelengths are not equally spaced throughout. The notebook
correctly supplies the actual wavelength array to trapezoidal integration.
Normalization is over the selected wavelength interval, not the unmeasured full
emission band. It is not peak normalization. With wavelength in nm, TRANES has
units 1/nm. Peak-normalized steady-state input gives relative TRES intensity per
ns; this does not establish absolute radiometric calibration.

## Supplied three-component example

The attached data file is comma-separated and includes a header. It contains 13
wavelength rows, 498–653 nm, and the columns wavelength, A1, A2, A3, and
steady-state intensity. The matching lifetime order from the notebook is:

- A1: 114.15 ns.
- A2: 17.415 ns.
- A3: 4.613 ns.

Component columns and lifetime values must stay paired; sorting lifetimes alone
would produce a different and incorrect result. The steady-state values are
already peak-normalized to one in this particular supplied table. They are not a
steady-state spectrum for the separate 51-wavelength demonstration measurement.

Independent scalar evaluation at delays 0, 0.025, 1, 10, 50 and 100 ns generated:

- `notebook_reference_tres.csv`.
- `notebook_reference_tranes.csv`.
- `notebook_reference_summary.json`.

Every TRANES trapezoidal spectral area is one to about 2.3e-16 numerical error.
Independent adaptive quadrature of the reconstructed decay over time recovers
the supplied steady-state intensity at each wavelength, with maximum absolute
error 2.23e-16.

Comparison with the new engine gives maximum absolute differences 5.56e-17 for
TRES and 6.08e-18 for TRANES. Reversing wavelength order together with the aligned
amplitudes and steady-state values reproduces the same sorted result. Independent
positive amplitude factors spanning 1e-200 to 1e200 at different wavelengths also
leave the reconstructed spectra unchanged, as required by the formula. A known
7 ns tail-reference transformation recovers the excitation-origin reconstruction.
These checks are recorded in `independent_engine_comparison.json`.

## Notebook issues to correct in the application

1. **Input parsing:** the notebook hardcodes a tab separator and `header=None`;
   that does not match the actual comma-separated, header-containing attachment.
   The application must recognize the real format and require the right number
   of lifetime values.
2. **Zero areas:** the notebook silently assigns zero when a fitted time area
   is nonpositive, then divides by an unchecked spectral area for TRANES. A
   positive steady-state intensity with no fitted decay support should instead
   be an explicit error, as should a zero/invalid TRANES denominator.
3. **Negative delays and invalid lifetimes:** the notebook has no adequate bounds
   or component-count checks. Reconstruction requires finite positive lifetimes
   and finite nonnegative delays for this causal model.
4. **Cubic interpolation:** the notebook's cubic B-spline display can overshoot
   even when every input intensity is nonnegative. For the supplied data it
   creates a minimum near −1.276e-5 at 50 ns and −2.244e-6 at 100 ns. This is an
   interpolation artifact, not negative emission. Numerical exports should keep
   the sampled reconstructed values; straight lines or shape-preserving display
   interpolation are preferable.
5. **Peak precision:** at 50 ns the sampled maximum is 603 nm, whereas the
   notebook's cubic interpolant reports about 600.12 nm. At 100 ns it reports
   about 600.64 nm. Such sub-grid peak positions are interpolation-dependent and
   should not be presented as measured spectral resolution or a precise physical
   shift without a justified spectral model and uncertainty estimate.
6. **Suggested time horizon:** five times the longest lifetime is 570.75 ns in
   this example. A pure exponential still has exp(−5) ≈ 0.006738 of its initial
   amplitude there. Whether that is negligible depends on component weight,
   noise and dynamic range. It is a preview heuristic, not a guaranteed cutoff.
7. **Widget portability:** the upload code assumes the newer tuple-style widget
   value; earlier widget versions use a different structure. A native application
   input surface avoids copying this version-dependent callback logic.
8. **Export completeness:** the notebook exports the TRES matrix only. A clear
   application export should include TRES, TRANES, sampled steady-state values,
   lifetime/column pairing, normalization settings and reference-time assumptions.

## Conditions needed for physical interpretation

- A dense uploaded steady-state spectrum may be normalized, then sampled at the
  fitted wavelengths. Interpolation should remain inside its measured coverage;
  extrapolation must not silently supply missing spectral information.
- Steady-state and lifetime measurements should describe comparable sample and
  excitation conditions. Wavelength response, acquisition duration and spectral
  bandwidth corrections must be understood. A relative spectrum is not an
  absolute emission calibration.
- Fitted nonnegative parallel-exponential models are supported. A reconstructed
  spectrum does not establish solvation, energy-transfer or sequential kinetics;
  those may require signed amplitudes or target models.
- Reconvolution model amplitudes in this workbench refer to underlying
  excitation-origin exponentials before the discrete bin-average factor. Their
  intrinsic continuous area is therefore Σ Aτ. That bin-average factor must not
  be applied again during this reconstruction. External amplitudes with different
  definitions (e.g. integrated-intensity fractions or bin-integrated coefficients)
  require the corresponding conversion and cannot be assumed interchangeable.
- Tail-fit amplitudes refer to the first tail fitting channel. If its known
  elapsed delay after excitation is d, convert Aref to A0 = Aref exp(d/τ). The
  absolute displayed detector-time coordinate alone is not enough to determine
  d. Early-time reconstruction from a tail fit is model extrapolation.
- Lifetime and steady-state interpolation uncertainties are not automatically
  propagated into reconstructed spectra. Delays beyond the measured
  post-excitation window depend on extrapolating the fitted model.
