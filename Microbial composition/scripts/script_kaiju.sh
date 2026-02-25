#!/bin/bash

#SBATCH -J kaiju
#SBATCH --time=00:45:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=100G

# output files
#SBATCH -o /work/%u/job_logs/%x-%A-%a.out
#SBATCH -e /work/%u/job_logs/%x-%A-%a.err

# loading modules
module load Miniforge3/24.3.0-0
source activate /gpfs1/schlecker/global/apps/kaiju/1.10.1

# get list of input files from first argument
input_files="$1"

# get the specific input file for this task
# this is the n-th (where n is current task ID) line of the file
input_file=$(awk "NR==$SLURM_ARRAY_TASK_ID" "$input_files")

# file directory
DATABASE=/path/to/db_kaiju/
WD=/path/to/work_directory/
READS=/path/to/unmmaped_reads/

# command
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

kaiju -t $DATABASE/nodes.dmp -e 5 -s 50 -z $SLURM_CPUS_PER_TASK \
    -f $DATABASE/kaiju_db_nr_euk.fmi \
    -i $READS/"$input_file".fastq \
    -o $WD/"$input_file"_kaiju.out

kaiju2table -t $DATABASE/nodes.dmp -n $DATABASE/names.dmp -r genus \
    -l superkingdom,phylum,class,order,family,genus \
    -o $WD/"$input_file"_kaiju.tsv $WD/"$input_file"_kaiju.out
