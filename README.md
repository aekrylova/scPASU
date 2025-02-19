# Snakemake workflow: scPASU

## Installation requirements:
Most of the needed packages and tools are included in the Snakemake environment.yml files. The following tools need to be installed separately:
* [cellranger](https://www.10xgenomics.com/support/software/cell-ranger/latest/tutorials/cr-tutorial-in)
* [polyAfilter](https://github.com/MarekSvob/polyAfilter)
* [goldmine](https://github.com/jeffbhasin/goldmine)
  * [Further goldmine documentation](https://www.jeffbio.com/goldmine/docs/)

Use the following command to clone this repository:
```
git clone https://github.com/aekrylova/scPASU.git
```
To run, navigate to the workflow folder within the directory structure and run the snakemake command from there. 

## Sample data setup instructions:
Download the sample dataset FASTQs and Seurat object here: [Preview dataset](https://data.mendeley.com/preview/6bf3w4wrj8?a=10479bdc-60fa-4632-84f8-89545ed45cf3) 

Run cellranger alignment on the eight samples in the dataset using the GRCh38_v2024 reference available from 10XGenomics. Make sure the --id argument matches the sample name as listed in the dataset.

Edit the config/config.yaml file to match your data. To run the sample dataset, the only variables that should be changed are ncores, work_dir, cellrangeroutputpath, cellrangerrefpath, polyAfilterpath, cellrangerpath, seurat_obj_file, polyAfilter_workdir, and goldminepath to match the locations and specifications for your system.  

It is recommended to run the pipeline on an HPC cluster: [Snakemake Cluster Execution](https://snakemake.readthedocs.io/en/v7.19.1/executing/cluster.html)

Finally, since the pipeline is built using conda environments, make sure to specify the –use-conda argument when running. 

If you are using snakemake version < 8.0, you can use the --cluster command line argument to submit to an HPC. Example command:
```bash
snakemake run_module_1 --jobs 20 --cluster 'bsub -q e80medium -W 6:00 -u AEKrylova@mdanderson.org -n 32 -M 200 -R "rusage[mem=200]"' --use-conda
```
If you are using snakemake version 8.01 or greater, the --cluster command is deprecated. You will need to install the [cluster-generic plugin](https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/cluster-generic.html) and use that in the command line. Alternatively, Snakemake also has lsf- or slurm-specific plugins you can install instead. Example command:
```bash
snakemake run_module_1 --jobs 20 --executor cluster-generic --cluster-generic-submit-cmd 'bsub -q e80medium -W 6:00 -u AEKrylova@mdanderson.org -n 32 -M 200 -R "rusage[mem=200]"' --use-conda
```

Additional guidance for using the snakemake tool: [General Snakemake Command Line Reference](https://snakemake.readthedocs.io/en/v7.19.1/executing/cli.html)

## DAG of pipeline workflow:

![dag](https://github.com/user-attachments/assets/4062fef2-eb3a-4a7e-9fb2-3fc82fa4e15f)

