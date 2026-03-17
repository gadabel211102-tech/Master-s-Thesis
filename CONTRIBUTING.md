# Contributing Guide

This repository was developed as an academic analysis pipeline. If you or a collaborator extend it, the following practices will help keep it consistent and reviewable.

## General Principles

- keep outputs reproducible
- avoid hard-coding new paths in individual scripts
- prefer the shared helpers in `pipeline_utils.py` and `pipeline_validation.py`
- preserve the current tumour/control definitions unless the study design is intentionally being changed
- document any analytical assumption that affects interpretation

## Code Style

For Python scripts:

- add a short module header explaining purpose, inputs, outputs, and method
- use clear section comments for major analytical stages
- keep comments explanatory rather than decorative
- prefer shared config and utility functions over duplicated logic

For shell and R scripts:

- keep headers concise but academically clear
- explain workflow steps that would not be obvious to an external reader

## Pull Request Checklist

Before merging changes, try to confirm that:

- the script still runs syntactically
- any new required columns or files are validated
- new comparison plots use proportional metrics when sample-size bias would matter
- output labels remain consistent with the current `Control` terminology
- documentation is updated if the workflow changes

## Large Changes

If a change affects study design, statistical thresholds, or grouping strategy, document it explicitly in:

- the relevant script header
- [README.md](README.md)
- [REPRODUCIBILITY.md](REPRODUCIBILITY.md)

That makes it much easier for supervisors and collaborators to understand what changed and why.
