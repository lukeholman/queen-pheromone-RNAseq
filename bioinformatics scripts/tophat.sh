#!/bin/bash
#SBATCH --job-name=tophat
#SBATCH --partition=compute
##SBATCH --mem=1G
#SBATCH --time=1-00:00:00
##SBATCH --cpus-per-task=1
##SBATCH --mail-user=%u@oist.jp
##SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --input=none
#SBATCH --output=%j.out
##SBATCH --error=job_%j.err
. $HOME/.bashrc

#species=ln  #am, bt, ln
#for sbatch: --export=species=am
refdir=./data/raw_reads/trimmed
a=($refdir/$species*_1.fastq.gz) #10 for bt, 6 for am, 12 for ln
b=($refdir/$species*_2.fastq.gz)
f=${a["SLURM_ARRAY_TASK_ID"]}
r=${b["SLURM_ARRAY_TASK_ID"]}
base=`basename $f _1.fastq.gz`

module load bowtie2/2.2.6 tophat/2.1.1 cufflinks/2.2.1

if [ "$species" == "bt" ]; then
    gff=./ref/GCF_000214255.1_Bter_1.0_genomic.gff
    ref=./ref/bt
elif [ "$species" == "am" ]; then
    gff=./ref/GCF_000002195.4_Amel_4.5_genomic.gff
    ref=./ref/am
else
    gff=./ref/GCA_001045655.1_ASM104565v1_genomic.gff 
    ref=./ref/ln
fi

echo $SLURM_ARRAY_TASK_ID $species $gff $ref

echo  $species $ref  $base $gff $f $r 
#tophat2 -p 1 -G $gff -o ./data/assembly/tophat/$base $ref $f $r
cufflinks -p 1 -g $gff ./data/assembly/tophat/$base/accepted_hits.bam -o ./data/assembly/tophat/$base