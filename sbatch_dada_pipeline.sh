#!/bin/bash

#SBATCH --job-name="dada_array" #name of the job submitted

#SBATCH -p short #name of the queue you are submitting to

#SBATCH -N 1 #number of nodes in this job

#SBATCH -n 40 #number of cores/tasks in job

#SBATCH -t 48:00:00 #time allocated for job in hrs:min:sec

#SBATCH --mail-user=rmkepler@gmail.com #email to get messages

#SBATCH -o "stdout_dada_%j.%N" #standard out %j adds job number, %N adds node name

#SBATCH -e "stderr_dada_%j.%N" #optional file for standard error

#SBATCH -a 1-12

#SBATCH --workdir=/project/php-fungi #set the working directory

date # optional prints timestamp in stdout when job starts

module load perl/gcc/64/5.24.1
module load cutadapt/gcc/64/1.8.3
module load r/gcc/64/3.4.1

mkdir -p sequence_variants

perl dada_pipeline.pl --plate ${SLURM_ARRAY_TASK_ID} \
	--outfolder sequence_variants \
	--list fungi_dada_list.txt \
	--procs 40


date #optional prints timestamp when the job ends to standard out

#End of file
