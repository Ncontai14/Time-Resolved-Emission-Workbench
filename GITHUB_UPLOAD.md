# Put this package on GitHub

This is the repository-ready first public release of TCSPC Workbench, version 0.1.0. No repository has been created or published for you.

## 1. Extract and identify the repository root

Extract the ZIP completely. Open the inner `tcspc-workbench` folder until you see `README.md`, `app.py`, `tcspc`, `tests`, and `.github` together. That is the repository root. Do not upload the ZIP itself as the repository's only file, and do not accidentally place all these files inside another extra folder on GitHub.

Keep dotfiles and the `.github` directory: they provide ignore rules, text-file handling and test automation. File managers on macOS/Linux may hide them; show hidden files or use GitHub Desktop to include them.

## 2. Personalize the public credit and review the examples

- `CITATION.cff` currently credits **TCSPC Workbench contributors**. Replace that collective entry with your intended author information, optionally ORCID. Do not invent a DOI or paper citation. Add the real repository URL after creating it.
- The existing MIT licence covers original project code/documentation. Keep applicable notices. Contributed measurements have separate rights; read `examples/README.md` and confirm redistribution permission before making the complete example collection public.
- Keep your own experimental data outside the repository or under the ignored `data/` or `private_data/` folders. Do not add tokens, `.env` files or `.streamlit/secrets.toml`. The supplied `.streamlit/config.toml` is a shared interface configuration and should be included.

The supplied examples are preserved because the app and tests use them. The package does not establish permission to redistribute third-party measurements. Keeping a new repository private while checking those rights preserves the current working package. If permission cannot be established, make a separately validated public edition with substitute synthetic examples and updated tests/documentation before publication.

## 3. Upload using GitHub Desktop

This route handles hidden configuration files and folder structure.

1. Install GitHub Desktop from its official website and sign into your GitHub account.
2. Choose **File → Add local repository**, then select the extracted repository root.
3. If Desktop says it is not a Git repository, use its **create a repository here** link. Keep the location pointing to the existing folder. Do not replace the supplied README, licence or ignore file with generated defaults.
4. Review the Changes list. You should see the code, tests, examples, guides, `.gitignore`, `.gitattributes`, `CITATION.cff`, and `.github/workflows/tests.yml`. Private measurements and environments should not be included.
5. Commit with a summary such as `Prepare TCSPC Workbench v0.1.0`.
6. Select **Publish repository**. A suitable repository name is `tcspc-workbench`. Choose visibility deliberately; finish the author/data review before making the full collection public.
7. Open the repository on GitHub and check the README, screenshots, guides and **Actions** tab.

If Desktop creates a different folder from the one containing `app.py`, cancel and correct the location before committing.

## 4. Alternative: upload in the GitHub website

1. Create an empty repository named `tcspc-workbench` (or your preferred name).
2. Leave GitHub's automatic README, `.gitignore` and licence creation unchecked; these are already supplied.
3. Use **Add file → Upload files** and upload the contents of the local repository root, including subfolders and hidden configuration files. Upload in batches if the website requests it.
4. Commit the files. Confirm that `app.py` is at the top level, not inside an extra `tcspc-workbench` subfolder.
5. Check that `.github/workflows/tests.yml` and `.gitignore` really arrived. If your browser/file manager omits them, use GitHub Desktop.

Uploading a ZIP as a normal file does not unpack it into source files. The ZIP is useful as a release attachment instead.

## 5. Check automated tests and create a release

The supplied workflow runs numerical/import/CLI tests on Windows and Linux with Python 3.10 and 3.12. It does not publish a website, upload data to another service or require repository secrets. Hosted test results appear only after GitHub executes the workflow; the local checks do not imply a successful hosted run.

After checking those results, create a GitHub Release from **Releases → Draft a new release**:

- Tag: `v0.1.0`.
- Title: `TCSPC Workbench v0.1.0`.
- Description: summarize the available workflows and link the user guide, scientific notes, validation and limitations.
- Optionally attach the complete end-user ZIP as a release asset. GitHub also creates source archives from the tag.

The tag identifies this initial public release. Start the public changelog here and record subsequent public changes under their own version entries.

Suggested short repository description:

> Local TCSPC lifetime fitting, global analysis and steady-state-scaled TRES/TRANES reconstruction in Python.

Suggested topics: `tcspc`, `fluorescence`, `fluorescence-lifetime`, `global-analysis`, `spectroscopy`, `python`, `streamlit`.

## 6. Make it discoverable and citable

Link the stable repository or release in your profile and demonstration video. Ask experimental users to report reproducible issues. A tagged release can later be archived with Zenodo, at which point you can add the real DOI and repository URL to CITATION.cff. The citation currently contains neither because those identifiers do not yet exist.

## What is included in this repository

- A repository overview, detailed INSTALLATION.md and an experimental-user operating guide.
- Contribution/support guides, citation metadata and upload instructions.
- GitHub Actions tests plus issue and pull-request templates.
- Ignore rules for environments, credentials and designated personal data/output folders.
- Text attributes preserving histogram bytes/checksums across platforms.
- Data-provenance notes and updated offline HTML navigation.

This distribution has no Git history or account credentials. GitHub Desktop or your Git client initializes the repository locally. The public changelog starts with v0.1.0.
