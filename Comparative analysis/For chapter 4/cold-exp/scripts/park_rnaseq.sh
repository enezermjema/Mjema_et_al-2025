#!/bin/bash

#SBATCH -J park_rnaseq
#SBATCH --time=2-00:00
#SBATCH --mem-per-cpu=16G
#SBATCH --cpus-per-task=40
#SBATCH --chdir=/work/hartmaca/park_rnaseq
#SBATCH --mail-user=cajus.hartmann@gmx.net
#SBATCH --mail-type=BEGIN,END,FAIL

#output files
#SBATCH -o /work/%u/job_logs/%x-%j.out
#SBATCH -e /work/%u/job_logs/%x-%j.err

#loading modules
module load Nextflow/24.10.1

#path
GENOME=/data/laub-rna/hiwi/cajus/arabid_genome

#command
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
nextflow run nf-core/rnaseq \
    --input /work/hartmaca/park_rnaseq/samplesheet_rnaseq.csv\
    --outdir /work/hartmaca/park_rnaseq \
    --fasta $GENOME/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz \
    --gtf  $GENOME/Arabidopsis_thaliana.TAIR10.60.gtf.gz \
    -profile singularity