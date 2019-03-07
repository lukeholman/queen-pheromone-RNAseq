#!/bin/bash
#SBATCH --job-name=trinity
#SBATCH --partition=compute
#SBATCH --mem=80G
#SBATCH --time=3-00:00:00
#SBATCH --cpus-per-task=12
##SBATCH --mail-user=%u@oist.jp
##SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --input=none
#SBATCH --output=%j.out
##SBATCH --error=job_%j.err
. $HOME/.bashrc

species=lf
refdir=../data/raw_reads/merged
module load Trinity/2.1.1  bowtie/1.1.0
Trinity --seqType fq --max_memory 75G --left $refdir/"$species"_1.fastq.gz \
 --right  $refdir/"$species"_2.fastq.gz  \
    --CPU 12 --output ../data/assembly/trinity_"$species"
