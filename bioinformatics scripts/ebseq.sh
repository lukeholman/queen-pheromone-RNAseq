#!/bin/bash
#SBATCH --job-name=ebseq
#SBATCH --partition=compute
#SBATCH --mem=500M
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=1
##SBATCH --mail-user=%u@oist.jp
##SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --input=none
#SBATCH --output=%j.out
##SBATCH --error=job_%j.err
. $HOME/.bashrc


rsem=/apps/unit/MikheyevU/sasha/RSEM
refdir=./data/rsem
a=($(ls -1  ./data/raw_reads/trimmed/*_1.fastq.gz |awk -F/ '{print $NF}' |cut -c-2 |uniq)) #4
species=${a["SLURM_ARRAY_TASK_ID"]}
echo $species

# contrasts are relative to control treatments, i.e., control goes first

if [ "$species" == "bt" ]; then
	samples=`for i in 2 4 5 8 10 1 3 6 7 9; do echo -ne $refdir/$species$i".isoforms.results ";done`
	$rsem/rsem-generate-data-matrix $samples > $refdir/$species.matrix
    conditions="5,5"
elif [ "$species" == "am" ]; then    
	samples=`for i in 2 4 6 1 3; do echo -ne $refdir/$species$i".isoforms.results ";done`
	$rsem/rsem-generate-data-matrix $samples > $refdir/$species.matrix
    conditions="3,2"
elif [ "$species" == "lf" ]; then
	samples=`for i in 4 6 8 10 12 14 16 1 3 5 7 11 13; do echo -ne $refdir/$species$i".isoforms.results " ;done`
	echo $samples
	$rsem/rsem-generate-data-matrix $samples > $refdir/$species.matrix
    conditions="7,6"
elif [ "$species" == "ln" ]; then
    samples=`for i in 1 3 5 7 11 2 4 6 8  12; do echo -ne $refdir/$species$i".isoforms.results ";done`
	$rsem/rsem-generate-data-matrix $samples > $refdir/$species.matrix
    conditions="5,5"
fi

$rsem/rsem-generate-data-matrix `echo $samples | sed 's/isoforms/genes/g'` > $refdir/$species.genes.matrix
$rsem/rsem-run-ebseq $refdir/$species.genes.matrix $conditions $refdir/$species".genes.ebseq"
$rsem/rsem-control-fdr $refdir/$species".genes.ebseq" .05 $refdir/$species".genes.padj.ebseq"

$rsem/rsem-run-ebseq --ngvector ./ref/$species.ngvec $refdir/$species.matrix $conditions $refdir/$species".ebseq"
$rsem/rsem-control-fdr $refdir/$species".ebseq" .05 $refdir/$species".padj.ebseq"

