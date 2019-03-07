1 (A): BeeBase 'Release1 of the OGS' GB identifier

2 (B): UCSC Jan 2005 assembly 'NCBI Genes' identifier (mapped to GB IDs using
NCBI geneinfo and gene2refseq databases (Note that if multiple rows exist for
ID, only first was taken. Some may have NM and XM IDs, possibly resulting in
some very small loss of data.)

3-6 (C-F): Sequence lengths (total coding exon, total intron, total untranslated
exon)for UCSC Jan 2005 assembly NCBI Genes

7 (G): Coding exon total length for OGS Release1

8-12 (H-L): BAGEL normalized expression values from Grozinger et al. 2007 cDNA
microarray. (12) log2Q/W expression ratio.

13-18 (M-R): P-values for queen (q), sterile-worker (N), and reproductive worker
(O) gene expression comparisons from Grozinger et al. 2007 cDNA microarray.  

19 (S): Significance designation of expression difference between queen and sterile-worker
indicating upregulated caste or nonsignificance, from Grozinger et al. 2007
supplementary materials.

20-24 (T-X): Ortholog and branch length data for Amel and Nvit from alignments of 3-5
insects (muscle and phyml used; from Stefan - revision1).  FlyID is flybase ID for Dmel.

26 (Z): Codon adaptation index. CAI was calculated by first running Amel OGS
release1 through the Emboss tool cusp to calculate a codon frequency table. The codon
frequency table and Amel OGS release1 were analyzed by Emboss cai to calculate
the Codon Adaptation Index.

28 (AB): Effective number of codons for Amel OGS release1. ENC (Nc; Wright 1990)
was calculated using codonW and confirmed with the EMBOSS tool chips.

30-49 (AD-AW): codonW base composition analysis of Amel OGS release1.  Overall
GC content was confirmed with the EMBOSS tool geecee and codon position GC
content was confirmed with a small sample file. GC3s (33) is the G+C content of 3rd
position of synonymous codons. This is the fraction of codons, that are
synonymous at the third codon position, which have either a guanine of cytosine
at that third codon position.  GC3 (37) is the raw G+C content at the third codon
position.

51 (AY): G+C content of introns (UCSC Jan 2005 assembly NCBI Genes) calculated by
the EMBOSS tool geecee.

52-68 (AZ-BP): CpG analysis of UCSC Jan 2005 assembly NCBI gene bodies (introns
+ exons) used for Elango et al. 2009 methylation analysis.

69 (BQ): ABS value of log2Q/W expression ratio

