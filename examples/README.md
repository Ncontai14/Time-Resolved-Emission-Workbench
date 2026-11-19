# Example data and redistribution

These examples were used to develop and validate the workbench. Their inclusion in a working package is not a blanket grant of rights to publish third-party measurements.

| Files | Provenance and intended use |
|---|---|
| `synthetic_TRES.csv`, `synthetic_IRF.csv`, `synthetic_truth.json` | Project-generated synthetic validation data, with known lifetime inputs recorded in the truth file. These synthetic examples are provided under the project's MIT licence. |
| `Decay.dat`, `IRF.dat` | Contributed demonstration decay/response supplied during development. Original rights remain with their owners; this repository does not establish redistribution permission. |
| `measured_TRES.csv`, `measured_IRF.csv` | Contributed measured demonstration matrix/IRF. No corresponding steady-state fluorescence spectrum was supplied. Original rights remain with their owners. |
| `curve_exports/*.dat` | Contributed multi-curve ASCII demonstration exports. Source-curve mappings and provenance hashes are recorded in EXAMPLE_DATA.json. Original rights remain with their owners. |
| `fitted_three_components.csv` | Contributor-supplied fitted amplitudes and steady-state values from a separate dataset. Its lifetime pairing and source hash are recorded in EXAMPLE_DATA.json. It is not the steady-state spectrum for measured_TRES.csv. Original rights remain with their owners. |

`EXAMPLE_DATA.json` identifies available provenance, input hashes, mappings and lifetime conventions. The measured examples' public redistribution permissions have not been established by this packaging task. The repository owner should confirm them before a public release of the full collection. No third-party licence is inferred from the project's MIT licence.

`example_results/` contains derived reference outputs and therefore also needs consideration alongside its corresponding source measurements when redistributing. Screenshots and notebook-review numerical references may display or derive from contributed examples.

The app and tests depend on several contributed examples. Any edition that substitutes synthetic datasets should also update the dependent demonstrations, tests and guides. Version 0.1.1 retains the existing example collection.
