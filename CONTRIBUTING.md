# Contributing

Bug reports, documentation improvements, reproducible examples and carefully validated numerical changes are welcome. Use the repository's Issues and Pull requests pages once it is hosted on GitHub.

## Development setup

1. Clone or extract the full repository and enter its root folder.
2. Activate a working environment or create one using [INSTALLATION.md](INSTALLATION.md).
3. Install dependencies with `python -m pip install -r requirements.txt` when using the pip route.
4. Run `python -m unittest discover -s tests -v`.
5. Launch `python -m streamlit run app.py` for interface changes.

Keep a branch for each focused change. Explain the problem, resulting behavior and relevant validation in the pull request. Do not commit local environments, credentials, private measurements or generated personal exports. The tracked example results are reference outputs; regenerate them only when the change warrants it.

## Numerical changes

State units, time-origin assumptions, model equations and statistical conventions. Add an independent reference or a numerical invariant that would catch an incorrect implementation. Explain effects on reconvolution, bin widths, background, amplitude interpretation and uncertainty where relevant. An improvement in visual smoothness alone is not evidence of a better physical model.

For interface changes, check the affected workflow through a real fit or reconstruction and verify exported numerical values. Document new controls in USER_GUIDE.md and keep HTML editions consistent.

## Data contributions

Provide only data you are entitled to redistribute, with an explicit data licence or permission statement and acquisition/provenance notes. Prefer small synthetic reproductions for bugs. Do not assume the project's software licence covers third-party measurements.

## Documentation

The repository README is a concise overview. INSTALLATION.md contains the detailed setup reference, and USER_GUIDE.md is the experimental-user operating manual. Maintain their offline HTML editions. Use relative links so GitHub and extracted packages remain navigable.
