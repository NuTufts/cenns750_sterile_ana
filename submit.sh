#!/bin/bash

#SBATCH --job-name=sterilenu
#SBATCH --output=grid_cenns750_sterilenu.log
#SBATCH --mem-per-cpu=4000
#SBATCH --time=10:00
#SBATCH --array=10-121

CONTAINER=/cluster/tufts/wongjiradlab/twongj01/coherent/coherent_snowglobes_20200304.simg
WORKDIR=/cluster/tufts/wongjiradlab/twongj01/coherent/run_cenns750_sterile_jobs
OUTDIR=${WORKDIR}/grid_output
PARAM_LIST=${WORKDIR}/param_file_dm2_Ue4sq.dat

mkdir -p $OUTDIR

module load singularity
singularity exec $CONTAINER bash -c "cd $WORKDIR && source run_job.sh $PARAM_LIST $OUTDIR $WORKDIR"
# source run_job.sh $PWD/param_file_dm2_Ue4sq.dat $PWD/grid_output/ $PWD/
