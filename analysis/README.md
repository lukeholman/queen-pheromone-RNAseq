## Click [here](https://mikheyev.github.io/queen-pheromone) to view the HTML report


## Main .Rmd file
All of our R statistical analyses are described in the `queen_pheromone.Rmd` file, which can be used to recreate all the figures, tables and analyses presented in the paper and supplementary material. To reproduce our analysis, clone this Github repository, open the .Rmd file, install all the R packages you don't already have, and then Run All or press the "knit" button in R studio.


## Database containing most of the necessary data
We created a sqlite3 database using `import_data.sh` prior to analysis in R. The component spreadsheets (tables) from this database are also available as .csv files in the directory `component spreadsheets of queen_pheromone.db` if you prefer to access them that way.

The spreadsheet has the following tables:

- **am2bt, am2lf, am2ln, bt2am, bt2lf, bt2ln, lf2am, lf2bt, lf2ln, ln2am, ln2bt, ln2lf**: These tables give the BLAST hits of species1 genes in species2, with the corresponding evalue. For example, am2bt BLASTs Apis genes against Bombus. Can be used to identify genes that are each other's reciprocal best BLAST hit.

- **bee_go, bee_kegg** Lists the GO nad KEGG terms associated with each Apis gene. See the file `Script to set up for GO analyses.R`, which created the bee_kegg file from data on Entrez.  

- **ebseq_gene_am, ebseq_gene_bt, ebseq_gene_lf, ebseq_gene_ln** Output of EB-seq analysis of gene-level expression data for each species. The PostFC is probably what you want. This list includes all the genes, significant or non-significant.

- **ebseq_padj_gene_am, ebseq_padj_gene_bt, ebseq_padj_gene_lf, ebseq_padj_gene_ln** Same, except that the list of genes has been culled to only include those showing significant response to treatment at p < 0.05, after Benjamini-Hochberg FDR correction.

- **ebseq_padj_isoform_am, ebseq_padj_isoform_bt, ebseq_padj_isoform_lf, ebseq_padj_isoform_ln** Output of EB-seq analysis of isoform-level expression data for each species. The PostFC is probably what you want. This list only includes isoforms showing a significant response to treatment at p < 0.05, after Benjamini-Hochberg FDR correction.

- **isoforms_am, isoforms_bt, isoforms_lf, isoforms_ln** List of mappings of isoforms to genes for each species.

- **rsem_am, rsem_bt, rsem_lf, rsem_ln** Gene expression values (RSEM) for each gene in each species.

- **treatments** Gives the treament, colony and species for each RNA-seq library. 

## Other files

### Files in `apis_gene_comparisons` directory

This directory contains files kindly shared with us by other research groups, all of which were used to compare our _Apis_ gene-level variables with other, previously-measured variables such as DNA methylation. 

- `apis_gene_methyl_CG_OE.csv` CpG and gene body methylation data provided by Soojin Yi and Xin Wu. The data are from Galbraith et al. 2016 _PNAS_.

- `harpur_etal_gamma.txt` Gamma values (i.e. strength of positive selection) from Harpur et al. 2014 _PNAS_. Provided by Brock Harpur.

- `Amel_AllData_012709.txt` Various gene-level metrics provided by Brendan Hunt. This is the source of the queen vs sterile/reproductive worker data (measured in Grozinger et al. 2007 _Mol. Ecol._), as well as the Codon Adaptation Index.  

- `am.gene_info.txt` Contains mappings between Entrez IDs, old Beebase IDs, and new Beebase IDs. Needed to relate different gene lists to each other.

### Custom gene set collection files

These were made by Luke following instructions at http://bioconductor.org/packages/2.11/bioc/vignettes/GOstats/inst/doc/GOstatsForUnsupportedOrganisms.pdf. One file is for GO, one for KEGG. Needed for GSEA in the GOstats package.

- `gene_set_collection.RData`
- `gene_set_collection_kegg.RData`


### Files for comparing our results to Morandin et al.

- `Morandin to Holman orthology.csv` Orthology for placing our genes into the OGGs in Morandin et al. _Genome Biology_ and vice versa. 

- `Morandin module membership.csv` Gives the module membership for each OGG in Morandin et al. _Genome Biology_.

- `Morandin module membership.csv` Gives the caste bias of each module in Morandin et al. _Genome Biology_.


### Additional scripts

- `Script to set up for GO analyses.R` Used to get the KEGG terms associated with each Apis gene.

- `orthodb.py` A Python script to access orthodb and get gene orthology information. 
