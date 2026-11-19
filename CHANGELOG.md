# Changelog

## 0.1.1 — 2026-11-19

- Adopted **Time-Resolved Emission Workbench** as the project name in the app, guides, screenshots and citation metadata.
- Documented when uniformly sampled fluorescence upconversion and other emission decay measurements can use the existing exponential models.
- Updated the software name and version in fit and spectral reconstruction exports.
- Included the existing `tcspc-workbench` distribution, `tcspc` Python module and command, and GitHub Actions configuration in the repository package.

## 0.1.0 — 2026-11-18

Initial public release as **TCSPC Workbench**.

- Single-decay IRF reconvolution and exponential tail fitting.
- Shared-lifetime global TCSPC analysis with separate trace amplitudes and backgrounds.
- Visual fit-window controls and kinetics previews for single and global fitting.
- Steady-state-scaled TRES and area-normalized TRANES reconstruction from global-fit results.
- Steady-state upload/normalization and manual normalized-intensity entry.
- Selectable spectral display connections and numerical exports with provenance.
- Experimental-user guides, example data, 69 numerical/import/CLI tests, and GitHub Actions configuration.

Known model and acquisition limitations are documented in SCIENTIFIC_NOTES.md and README.md.
