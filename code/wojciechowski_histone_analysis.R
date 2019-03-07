library(ape)
library(readr)
library(stringr)
library(dplyr)
library(ggplot2)
library(EBSeq)
library(reshape2)

########### Analysis of the ChIPseq data (reference number is GSE110640, https://www.ncbi.nlm.nih.gov/bioproject/?term=GSE110640)

# Open the ChIPseq data from wojciechowski et al. 2018 Genome Biology
files <- list.files("data/apis_gene_comparisons/wojciechowski_histone_data/GSE110640_RAW/")

# Make a list of all the ChIPseq samples
chipseq_samples <- data.frame(sample = 1:length(files), strsplit(files, split = "_") %>% 
  do.call("rbind", .) %>% as.data.frame() %>% 
    mutate(V3 = gsub("[.]txt[.]gz", "", V3),
           V2 = substr(V2,1,4))) %>%
  rename(sample=V1, caste=V2, histone=V3)

# Open all the ChIP files into a list
x <- lapply(list.files("data/apis_gene_comparisons/wojciechowski_histone_data/GSE110640_RAW/", full.names = TRUE), 
            function(x) read_tsv(x, col_names = c("seq_id", "start", "end", "value")))
# Make sure there is one observation per sample per site, by averaging across any duplicates
x <- lapply(x, function(x) group_by(x, seq_id, start, end) %>% summarise(value = mean(value)))

# Concatenate all the ChIP-seq data into one data frame using full_join
histone_data <- full_join(x[[1]], x[[2]], by = c("seq_id", "start", "end")) %>%
  dplyr::rename(sample1 = value.x, sample2 = value.y) 
for(i in 3:length(files)){
  histone_data <- full_join(histone_data, x[[i]], by = c("seq_id", "start", "end")) 
  names(histone_data)[names(histone_data) == "value"] <- paste("sample", i, sep="")
}


# Open the A. mellifera 4.5 genome annotation (same one as used by wojciechowski et al)
apis_gff <- read.gff("data/apis_gene_comparisons/wojciechowski_histone_data/GCF_000002195.4_Amel_4.5_genomic.gff.gz") %>% 
  filter(type == "gene") %>% select(seqid, start, end, attributes) %>%
  mutate(attributes = lapply(
    lapply(str_extract_all(attributes, "GB[:digit:]+"), unique), 
    function(x) ifelse(length(x) > 0, x, NA))) %>%
  rename(seq_id = seqid, gene = attributes) %>%
  filter(!is.na(gene), end - start < 50000) # only keep measurements where we know a gene nearby

# Remove any contigs that are in the ChIP-seq data but not the genome annotation
histone_data <- histone_data[histone_data$seq_id %in% apis_gff$seq_id, ] %>% ungroup()

# Annotate the ChIPseq data with the overlapping gene ID or IDs
histone_data <- parallel::mclapply(1:nrow(histone_data), function(i, apis_gff, histone_data){
  foc_histone <- histone_data[i,]
  foc_gff <- apis_gff %>% filter(seq_id == foc_histone$seq_id)
  hits <- which(
    (foc_histone$start > foc_gff$start & foc_histone$start < foc_gff$end) |
      (foc_histone$end > foc_gff$start & foc_histone$end < foc_gff$end) 
  )
  num_hits <- length(hits)
  if(num_hits == 1) return(data.frame(gene = as.character(apis_gff$gene[hits]), histone_data[i, ] %>% select(starts_with("sample"))))
  if(num_hits > 1)  return(data.frame(gene = as.character(apis_gff$gene[hits]), histone_data[rep(i, num_hits), ] %>% select(starts_with("sample"))))
  else return(NULL)
}, mc.cores = 8, apis_gff = apis_gff, histone_data = histone_data) %>% 
  bind_rows() %>%
  filter(!is.na(gene)) # get rid of epimarks with no known gene nearby

# Get the data into "tidy" format, for all three types of histone modification
H3K4me3_data <- histone_data %>% 
  select(gene, sample1, sample2, sample7, sample8) %>% 
  filter(!(is.na(sample1) | is.na(sample2) | is.na(sample7) | is.na(sample8))) %>% 
  gather(sample, normalised_hm, sample1, sample2,  sample7,  sample8) %>%
  mutate(caste = "W",
         caste = replace(caste, sample %in% c("sample1", "sample2"), "Q")) 

H3K27ac_data <- histone_data %>% 
  select(gene, sample3, sample4, sample9, sample10) %>% 
  filter(!(is.na(sample3) | is.na(sample4) | is.na(sample9) | is.na(sample10))) %>% 
  gather(sample, normalised_hm, sample3, sample4, sample9, sample10) %>%
  mutate(caste = "W",
         caste = replace(caste, sample %in% c("sample3", "sample4"), "Q")) 

H3K36me3_data <- histone_data %>% 
  select(gene, sample5, sample6, sample11, sample12) %>% 
  filter(!(is.na(sample5) | is.na(sample6) | is.na(sample11) | is.na(sample12))) %>% 
  gather(sample, normalised_hm, sample5, sample6, sample11, sample12) %>%
  mutate(caste = "W",
         caste = replace(caste, sample %in% c("sample5", "sample6"), "Q")) 

parse_histone_data <- function(df){
  group_by(df, gene, caste) %>% 
    summarise(mean_histone_score = mean(normalised_hm)) %>%
    group_by(gene) %>% 
    summarise(mean_Q = mean_histone_score[1],
              mean_W = mean_histone_score[2],
              caste_difference = mean_histone_score[1] - mean_histone_score[2]) %>%
    arrange(-caste_difference)
}

# Find the average level of each modification per gene, for both castes (and the caste difference)
H3K4me3 <- parse_histone_data(H3K4me3_data) 
H3K27ac <- parse_histone_data(H3K27ac_data)
H3K36me3 <- parse_histone_data(H3K36me3_data)

my_db$con %>% db_drop_table(table = "H3K4me3")
copy_to(my_db, H3K4me3, "H3K4me3", temporary = FALSE)

my_db$con %>% db_drop_table(table = "H3K27ac")
copy_to(my_db, H3K27ac, "H3K27ac", temporary = FALSE)

my_db$con %>% db_drop_table(table = "H3K36me3")
copy_to(my_db, H3K36me3, "H3K36me3", temporary = FALSE)



# I took this section out because the RNAseq expression data from wojciechowski et al are archived in TPM format, which is a bit dubious to analyse, and I didn't feel like re-mapping all the reads just for this one test. The below code runs ebseq on RNAseq TPM data, which I later found out is a bad idea (one should instead analyse *count* data, e.g. as produced by RSEM). I did contact the corresponding author to ask for their gene expression results (e.g. fold change and p-values), and did not hear back, so decided to give up.

# ########### Analysis of the RNA-seq data (reference number is GSE110641, https://www.ncbi.nlm.nih.gov/bioproject/?term=GSE110641)
# 
# my_db <- src_sqlite("data/queen_pheromone.db")
# 
# files <- list.files("data/apis_gene_comparisons/wojciechowski_histone_data/GSE110641_RAW/")
# sampleID <- data.frame(sample = 1:length(files), strsplit(files, split = "_") %>% 
#                          do.call("rbind", .) %>% as.data.frame()) %>% 
#   mutate(V2 = substr(V2,1,7)) %>% rename(id = V1, caste = V2)
# 
# expression_data_isoforms <- lapply(list.files("data/apis_gene_comparisons/wojciechowski_histone_data/GSE110641_RAW/", full.names = TRUE), 
#                           function(x) read_tsv(x) %>% select(target_id, tpm)) %>% 
#   purrr::reduce(left_join, by = "target_id") 
# 
# expression_data_isoforms <- expression_data_isoforms[rowSums(expression_data_isoforms[,-1]) > 10^-5, ]
# 
# expression_data_genes <- left_join(expression_data_isoforms,
#                                    tbl(my_db, "isoforms_am") %>% 
#                                      collect(), 
#                                    by = c("target_id" = "isoform")) %>%
#   select(-target_id) %>%
#   filter(!is.na(gene)) %>%
#   group_by(gene) %>%
#   summarise_all(sum)
# 
# gene <- expression_data_genes$gene
# expression_data_genes <- as.matrix(expression_data_genes[,-1])
# rownames(expression_data_genes) <- gene
# colnames(expression_data_genes) <- sampleID$id
# 
# isoforms <- expression_data_isoforms$target_id
# expression_data_isoforms <- as.matrix(expression_data_isoforms[,-1])
# rownames(expression_data_isoforms) <- isoforms
# colnames(expression_data_isoforms) <- sampleID$id
# 
# 
# # Gene-level
# Sizes <- MedianNorm(expression_data_genes)
# EBOut <- EBTest(Data=expression_data_genes, 
#                 Conditions = as.factor(rep(c("QLH","WLH"), times = c(5,4))),
#                 sizeFactors=Sizes, maxround=50)
# wojciechowski_ebseq_gene <- GetDEResults(EBOut, FDR = 0.05)
# GeneFC <- PostFC(EBOut)
# wojciechowski_ebseq_gene <- left_join(melt(GeneFC$PostFC) %>% tibble::rownames_to_column() %>% rename(PostFC = value),
#                              melt(GeneFC$RealFC) %>% tibble::rownames_to_column() %>% rename(RealFC = value),  by = "rowname") %>%
#   left_join(as.data.frame(wojciechowski_ebseq_gene[[2]]) %>% tibble::rownames_to_column(),  by = "rowname") %>%
#   rename(gene = rowname) %>% 
#   filter(!is.na(PPEE)) %>%
#   arrange(-PPDE)
# 
# # Isoform level
# isoforms_am <- tbl(my_db, "isoforms_am") %>% collect()
# NgList <- GetNg(isoforms_am$isoform, isoforms_am$gene)
# IsoNgTrun <- NgList$IsoformNgTrun
# IsoSizes <- MedianNorm(expression_data_isoforms)
# 
# IsoEBOut <- EBTest(Data = expression_data_isoforms, 
#                    NgVector = IsoNgTrun,
#                    Conditions = as.factor(rep(c("QLH","WLH"), times = c(5,4))),
#                    sizeFactors = IsoSizes, maxround = 50)
# wojciechowski_ebseq_isoform <- GetDEResults(IsoEBOut, FDR = 0.05)
# IsoformFC <- PostFC(IsoEBOut)
# wojciechowski_ebseq_isoform  <- left_join(melt(IsoformFC$PostFC) %>% tibble::rownames_to_column() %>% rename(PostFC = value),
#                                  melt(IsoformFC$RealFC) %>% tibble::rownames_to_column() %>% rename(RealFC = value),  by = "rowname") %>%
#   left_join(as.data.frame(wojciechowski_ebseq_isoform[[2]]) %>% tibble::rownames_to_column(),  by = "rowname") %>%
#   rename(gene = rowname) %>% arrange(-PPDE)





# my_db$con %>% db_drop_table(table = "wojciechowski_ebseq_gene")
# copy_to(my_db, wojciechowski_ebseq_gene, "wojciechowski_ebseq_gene", temporary = FALSE)

# my_db$con %>% db_drop_table(table = "wojciechowski_ebseq_isoform")
# copy_to(my_db, wojciechowski_ebseq_isoform, "wojciechowski_ebseq_isoform", temporary = FALSE)