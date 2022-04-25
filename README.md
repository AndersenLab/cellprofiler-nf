# cellprofiler-nf
A nextflow pipeline to run various CellProfiler pipelines on raw images and process output data.

## QUEST usage
```bash
module load python/anaconda3.6
module load singularity
source activate /projects/b1059/software/conda_envs/nf20_env/

# run with your project directory
nextflow run main.nf --project /projects/b1059/projects/Tim/cellprofiler-nf/CP_test
```