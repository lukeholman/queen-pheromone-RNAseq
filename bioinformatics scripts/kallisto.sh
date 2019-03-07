#!/bin/bash
#SBATCH --job-name=kalisto
#SBATCH --partition=compute
#SBATCH --mem=3G
#SBATCH --time=5:00:00
#SBATCH --cpus-per-task=1
##SBATCH --mail-user=%u@oist.jp
##SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --input=none
#SBATCH --output=%j.out
##SBATCH --error=job_%j.err
. $HOME/.bashrc


module load kallisto

##species=lf #am, bt, ln
##SLURM_ARRAY_TASK_ID=0
#for sbatch: --export=species=am
#10 for bt, 6 for am, 12 for ln, 15 lf
refdir=./data/raw_reads/trimmed
a=($refdir/$species*_1.fastq.gz) 
b=($refdir/$species*_2.fastq.gz)
f=${a["SLURM_ARRAY_TASK_ID"]}
r=${b["SLURM_ARRAY_TASK_ID"]}
base=`basename $f _1.fastq.gz`

if [ "$species" == "bt" ]; then
    ref=./ref/GCF_000214255.1_Bter_1.0_rna
elif [ "$species" == "am" ]; then
    ref=./ref/GCF_000002195.4_Amel_4.5_rna
elif [ "$species" == "lf" ]; then
	ref=./ref/lf
elif [ "$species" == "ln" ]; then
    ref=./ref/ln
fi

kallisto quant -i $ref -o data/kallisto/$base -b 100 $f $r 