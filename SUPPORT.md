# Getting help

Begin with USER_GUIDE.md for operation, INSTALLATION.md for setup, DATA_FORMATS.md for inputs, and SCIENTIFIC_NOTES.md for model assumptions. The software requires a uniform time axis, matching IRF and decay spacing for reconvolution, and a known excitation delay for reconstruction from a tail fit.

For a bug, open a [GitHub issue](https://github.com/Ncontai14/Time-Resolved-Emission-Workbench/issues) using the template and include:

- Workbench version, operating system and Python/dependency versions.
- Which workflow and controls you used, with units and the selected time interval.
- Expected behavior, actual behavior and the complete relevant error text.
- A small synthetic or permission-cleared input that reproduces the issue, when possible.

Do not post private experimental data, passwords, tokens or full local configuration files. Remove unrelated personal paths from logs. For sensitive security issues, use the owner's private reporting/contact route rather than a public issue; the owner must configure that route after creating the repository.

This repository is research software. Interpret fits in the context of acquisition calibration, model adequacy and independent validation. A report of numerical convergence is not a substitute for those checks.
