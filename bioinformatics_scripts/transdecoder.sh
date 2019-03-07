#!/bin/bash
#SBATCH --job-name=transdecoder
#SBATCH --partition=compute
#SBATCH --mem=10G
#SBATCH --cpus-per-task=10
#SBATCH --time=1-00:00:00
#SBATCH --ntasks=1
##SBATCH --mail-user=%u@oist.jp
##SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --input=none
#SBATCH --output=%j.out
##SBATCH --error=job_%j.err

##SLURM_ARRAY_TASK_ID=2

module load Trinity/2.1.1 

infile=./data/assembly/trinity_lf/Trinity.fasta

base=lf

TransDecoder -t $infile --reuse --workdir ./output_"$base" --search_pfam /apps/unit/MikheyevU/sasha/TransDecoder_r20140704/pfam/Pfam-AB.hmm.bin --CPU 10 -v