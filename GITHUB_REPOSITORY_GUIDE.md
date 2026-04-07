# GitHub Repository Guide

This note is the practical upload guide for turning the thesis workspace into a clean GitHub repository.

## Recommended Principle

Keep the repository focused on:

- analysis code
- environment definitions
- concise workflow documentation
- small curated reference tables that are needed to understand the analysis

Keep local-only anything that contains participant-level data, raw sequencing files, generated outputs, or machine-specific runtimes.

## Recommended Public Repository Layout

For this thesis, the safest approach is to **keep the numbered scripts at the repository root**. The filenames already match the stage numbering used in the manuscript and in the launcher, so moving them into a new folder would make the thesis harder to cross-reference.

Recommended top-level contents:

- numbered analysis scripts: `01_...` through `26_...`
- shared helpers: `pipeline_utils.py`, `pipeline_validation.py`, `association_runtime.py`, `figure_style.py`, `objective2_*.py`
- launchers: `run_full_pipeline.sh`, `.run_full_pipeline_clean.sh`
- environment files: `environment.yml`, `requirements.txt`
- repository metadata: `README.md`, `CONTRIBUTING.md`, `.gitignore`
- sharable config template: `pipeline_config.example.toml`
- methods and workflow notes under `docs/notes/`
- thesis-facing narrative drafts under `docs/thesis_drafts/`
- small curated evidence tables that do not expose sensitive sample-level data

## Upload These

- source code files (`.py`, `.sh`, `.R`)
- `README.md`
- `.gitignore`
- `CONTRIBUTING.md`
- `environment.yml`
- `requirements.txt`
- `script_classification_table.tsv`
- `pipeline_config.example.toml`
- `pipeline_config.toml` only if you first remove private absolute paths and replace them with safe placeholders or relative paths
- `docs/notes/`
- `docs/thesis_drafts/`
- `docs/presentations/*.py` if you want to preserve slide-building code
- `docs/curated_snp_functional_evidence.tsv`

## Do Not Upload These

- `breast/` and `endometrium/` raw cohort directories
- `analysis_results/` generated outputs
- `archive/`
- `manifests/` if they contain project-specific sample identifiers you do not want public
- local environments: `.venv/`, `tfm_env/`, `miniconda3/`, `.vep/`
- bundled third-party software or installers: `ensembl-vep/`, `beagle.jar`, `setup/installers/`
- large or licensed reference resources: `ref/`, `ref_alt/`, `contigs/`, `setup/reference_helpers/`
- private clinical workbooks in `docs/source_workbooks/`
- derived workbooks that still contain participant-level information
- PowerPoint exports unless you specifically want them public

## Configuration Recommendation

The clean GitHub pattern for this project is:

1. Upload `pipeline_config.example.toml`.
2. Keep your personal machine-specific config as `pipeline_config.toml` locally.
3. If a collaborator wants to run the code, they copy the template and edit only the path values.
4. Optionally point to a custom config with `PIPELINE_CONFIG=/path/to/config.toml`.

This avoids exposing OneDrive, Downloads, or other private directories in the public repo.

## Suggested Folder Story For Examiners

If someone opens the repository for the first time, they should understand it in this order:

1. `README.md` for the overall workflow split and how to run it.
2. `run_full_pipeline.sh` for the actual execution order.
3. `script_classification_table.tsv` for the role of each script.
4. `docs/notes/PIPELINE.md` and `docs/notes/METHODS_ALIGNMENT.md` for thesis-facing interpretation.

## Final Pre-Push Checklist

- confirm `.gitignore` excludes raw data, outputs, environments, and private workbooks
- make sure no absolute personal paths remain in files you plan to upload
- remove temporary lock files such as `~$*.pptx`
- keep only the scripts and docs that are needed to explain or rerun the analysis
- add a repository description on GitHub that states this is the analysis code for a Master's thesis on GSDMB-focused targeted sequencing and RNA integration
