# Snakemake workflow: `<scPASU>`

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.3.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/<owner>/<repo>/workflows/Tests/badge.svg?branch=main)](https://github.com/<owner>/<repo>/actions?query=branch%3Amain+workflow%3ATests)


Installation requirements:
Most of the needed packages and tools are included in the Snakemake environment.yml files. The following tools need to be installed separately:
* [cellranger](https://www.10xgenomics.com/support/software/cell-ranger/latest/tutorials/cr-tutorial-in)
* [polyAfilter](https://github.com/MarekSvob/polyAfilter)
* [goldmine](https://github.com/jeffbhasin/goldmine)
Download the sample dataset FASTQs here: [LINK]
Run cellranger alignment on the eight samples in the dataset using the GRCh38_v2024 reference available from 10XGenomics. Make sure the --id argument matches the sample name as listed in the dataset.
Edit the config.yaml file to match your data. To run the sample dataset, the only variables that should be changed are ncores, work_dir, cellrangeroutputpath, cellrangerrefpath, polyAfilterpath, cellrangerpath, and goldminepath.  
It is recommended to run the pipeline on an HPC cluster: Snakemake Cluster Execution
Finally, since the pipeline is built using conda environments, make sure to specify the –use-conda argument when running. Example command:
snakemake run_module_1 --jobs 20 --cluster 'bsub -q e80medium -W 6:00 -u AEKrylova@mdanderson.org -n 32 -M 200 -R "rusage[mem=200]"' --use-conda
General Snakemake Command Line Reference



## Usage

The usage of this workflow is described in the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog/?usage=<owner>%2F<repo>).

If you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this (original) <repo>sitory and its DOI (see above).

# TODO

* Replace `<owner>` and `<repo>` everywhere in the template (also under .github/workflows) with the correct `<repo>` name and owning user or organization.
* Replace `<name>` with the workflow name (can be the same as `<repo>`).
* Replace `<description>` with a description of what the workflow does.
* The workflow will occur in the snakemake-workflow-catalog once it has been made public. Then the link under "Usage" will point to the usage instructions if `<owner>` and `<repo>` were correctly set.
