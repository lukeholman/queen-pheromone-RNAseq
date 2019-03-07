#!/bin/bash
#SBATCH --job-name=rsem
#SBATCH --partition=compute
#SBATCH --mem=4G
#SBATCH --time=2-00:00:00
#SBATCH --cpus-per-task=10
##SBATCH --mail-user=%u@oist.jp
##SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --input=none
#SBATCH --output=%j.out
##SBATCH --error=job_%j.err
. $HOME/.bashrc


module load rsem bowtie2

cd data/rsem
refdir=../raw_reads/trimmed
a=($refdir/*_1.fastq.gz) #43
b=($refdir/*_2.fastq.gz)
f=${a["SLURM_ARRAY_TASK_ID"]}
r=${b["SLURM_ARRAY_TASK_ID"]}
species=`basename $f | cut -c-2`
base=`basename $f _1.fastq.gz`

if [ "$species" == "bt" ]; then
    ref=../../ref/bt
elif [ "$species" == "am" ]; then
    ref=../../ref/am
elif [ "$species" == "lf" ]; then
    ref=../../ref/lf
elif [ "$species" == "ln" ]; then
    ref=../../ref/ln
fi

echo rsem-calculate-expression  -p 8 --paired-end --bowtie2 $f $r $ref $base
/apps/unit/MikheyevU/sasha/RSEM/rsem-calculate-expression -p 10 --bowtie2  --paired-end $f $r $ref $base
