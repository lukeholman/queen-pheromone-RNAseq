# Bioinformatic analysis starting from raw reads.

## Workflow

There are already sequenced genomes of _A. mellifera_ and _B. terrestris_ and _L. niger_ in NCBI. We'll have to make something for _L. flavus_ from scratch.

We also need to determine possible isoforms for the genes in the data set, which we will do using tophat for the genomes that have references. 

### _L. niger_ workflow

This specis has no isoform data in NCBI, so we create our own. First we create bowtie2 indexed for the reference genomes using `bowtie2-build`, and then tophat and cufflinks are executed using `tophat.sh` using reference-based mapping to identify alternative splicing.


### _L flavus_ workflow
`trinity.sh` assemble lf from raw reads. We are not using a genome guided assembly, since the ln genome is fragmented, as per [recommendation](https://github.com/trinityrnaseq/trinityrnaseq/wiki/Genome-Guided-Trinity-Transcriptome-Assembly) of the trinity authors
`transdecoder.sh` predict proteins from the transcripts. 


### Gene expression analysis

**kallisto** We use kallisto for gene expression analysis of transcripts from the NCBI data bases for _A. mellifera_ and _B. terrestis_, for the TopHat assembly of _L. niger_ and for the Trinity assembly of _L. flavus_. 

**rsem and ebseq** This is the more traditional approach using the same data sources. We prepare references as ngvector files from the predicted transcripts, as per instructions.

## Orthology

We handle this by reciprocal blastp hit on the proteins vs the honey bee genome (the best annotated of the bunch) _e.g.,_ `blastp -num_threads 12 -query lf.fa -db amel -outfmt 6 -evalue 1e-4 -max_target_seqs 1`

First, we have to generate proteins with the same names as the genes to simplify matching of orthologs later

	gffread GCF_000002195.4_Amel_4.5_genomic.gff -g GCF_000002195.4_Amel_4.5_genomic.fna -F -y -  |grep ">" | perl -ne 'print "$1\t$2\n" if /transcript_id=(.._\w+\..).*protein_id=(.._\w+\..)/'  > transcript2proteinApis.txt # 21764 transcipts
	join -1 1 -2 2 <(sort -k1,1 transcript2proteinApis.txt) <(sort -k 2,2 knownIsoformApis.txt) > transcript2protein3geneApis.txt
	awk 'NR==FNR {q[$2]=$3; next} $0!~/^>/ {print; next} {sub(/>/,"",$1); if ($1 in q) print ">"q[$1]; else print ">"$1}  ' transcript2protein3geneApis.txt GCF_000002195.4_Amel_4.5_protein.faa | awk 'BEGIN{RS=">"}NR>1{sub("\n","\t"); gsub("\n",""); print RS$0}' | awk '!seen[$1]++' | sed 's/.$//'  | tr "\t" "\n" > am_prot.fa	

	gffread GCF_000214255.1_Bter_1.0_genomic.gff -g GCF_000214255.1_Bter_1.0_genomic.fna -F -y - | perl -ne '(/GeneID:(\w+)/ && print ">$1\n") || print' | awk 'BEGIN{RS=">"}NR>1{sub("\n","\t"); gsub("\n",""); print RS$0}' | awk '!seen[$1]++' | sed 's/.$//' | grep -v "\." | tr "\t" "\n" > bt_prot.fa
	gffread GCA_001045655.1_ASM104565v1_genomic.gff -g GCA_001045655.1_ASM104565v1_genomic.fna -F -y - | perl -ne '(/gene=(\w+)/ && print ">$1\n") || print' | awk 'BEGIN{RS=">"}NR>1{sub("\n","\t"); gsub("\n",""); print RS$0}' | awk '!seen[$1]++' | sed 's/.$//' | grep -v "\." | tr "\t" "\n" > ln_prot.fa
	cat ../data/assembly/trinity_lf/Trinity.fasta.transdecoder.pep |sed 's/_i.*//' | tr -d '*' > lf_prot.fa


	for i in  *.txt; do cat $i | awk -v OFS=, '{print $1"|"$2,$(NF-1)}' | sort -t , -k1,2 -g | tr "|" "," |  awk -F, '!x[$1$2]++' > `basename $i txt`csv; done

