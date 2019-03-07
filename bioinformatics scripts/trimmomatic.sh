#!/bin/bash
#SBATCH --job-name=trim
#SBATCH --partition=compute
#SBATCH --mem=20G
#SBATCH --cpus-per-task=8
##SBATCH --mail-user=%u@oist.jp
##SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --input=none
#SBATCH --output=%j.out
##SBATCH --error=job_%j.err

trimmomatic=/apps/unit/MikheyevU/sasha/trimmomatic
##SLURM_ARRAY_TASK_ID=1
a=(*/*_1.fastq.gz) # 43
b=(*/*_2.fastq.gz) # 43
f=${a["SLURM_ARRAY_TASK_ID"]}   
r=${b["SLURM_ARRAY_TASK_ID"]}   
base=$(basename $f _1.fastq.gz)
java -jar $trimmomatic/trimmomatic-0.32.jar PE -threads 8 -phred33 $f $r trimmed/$base"_1.fastq.gz" trimmed/$base"_unpaired1.fastq.gz" trimmed/$base"_2.fastq.gz" trimmed/$base"_unpaired2.fastq.gz" ILLUMINACLIP:$trimmomatic/adapters/TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:4:15 MINLEN:25