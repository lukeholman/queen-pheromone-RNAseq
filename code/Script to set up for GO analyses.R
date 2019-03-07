# script to make the tables 'go_meanings', 'kegg_meanings' and 'bee_kegg' in the database
# Prepared separately because the GO.db package causes many annoying conflicts with other packages, e.g. dplyr, that we use a lot.
library(clusterProfiler)
library(GO.db)
library(dplyr)
my_db <- dplyr::src_sqlite("data/queen_pheromone.db", create = FALSE)

library(AnnotationHub)
hub <- AnnotationHub()
# query(hub, c("mellifera","sqlite")) # Get the number for mellifera
# Useful database of Apis orthologs xx <- hub[["AH10452"]]
zz <- hub[["AH62534"]] # Apis mellifera
saveRDS(zz, file = "data/apis.org.db")
test <- dbconn(zz)
db_list_tables(test)
bee_GO <- tbl(test, "go") %>% left_join(tbl(test, "genes")) %>% dplyr::select(GID, GO, ONTOLOGY) %>% collect()

names(geneList) <- entrez.tbl$entrez.id[match(names(geneList), gsub("BEEBASE:", "", entrez.tbl$beebase2))]
geneList <- geneList[!is.na(names(geneList))]

x <- gseGO(geneList, ont = "CC", OrgDb = apis.db, pvalueCutoff = 1, minGSSize = 3, maxGSSize = Inf)@result[,c(1, 2,5,6,7)] %>% arrange(-NES)
x <- x[x$pvalue < 0.05, ]
split(paste(x$ID, x$Description),  x$NES) %>% purrr::keep(function(x) length(x) > 2) 
split(x$Description, x$NES)

annotationTable <- data.frame(termID = bee_GO %>% filter(ONTOLOGY == "BP") %>% .$GO, 
                             geneID = bee_GO %>% filter(ONTOLOGY == "BP") %>% .$GID,
                             termName = "filler",
                             dbName = "apis_GO_BP",
                             description = "filler",
                             stringsAsFactors = FALSE)
annotationTable$termName <- annotationTable %>% left_join(tbl(dbconn(GO.db), "go_term") %>% filter(go_id %in% annotationTable$termID) %>% collect(), by = c("termID" = "go_id")) %>% .$term
annotationTable$description <- annotationTable %>% left_join(tbl(dbconn(GO.db), "go_term") %>% filter(go_id %in% annotationTable$termID) %>% collect(), by = c("termID" = "go_id")) %>% .$definition

gene_universe <- entrez.tbl[,1]
ranked_genes <- names(geneList)[names(geneList) %in% gene_universe]
collection <- buildSetCollection(annotationTable, referenceSet = gene_universe)
options(mc.cores=1)
out <- setRankAnalysis(ranked_genes[1:100], collection, setPCutoff = .99, fdrCutoff = 1, use.ranks = TRUE) 
out %>% igraph::get.data.frame()


geneID2GO <- split(bee_GO$GO, bee_GO$GID)
topDiffGenes <- function(allScore) {
  return(allScore > 0.1)
}
GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList, 
              annot = annFUN.gene2GO, geneSel = topDiffGenes, gene2GO = geneID2GO[names(geneID2GO) %in% names(geneList)], nodeSize = 3)
test.stat <- new("weight01Score", testStatistic = GOKSTest, name = "KS tests")
resultKS <- getSigGroups(GOdata, test.stat)
GenTable(GOdata, ks = runTest(GOdata, algorithm = "elim", statistic = "ks"), orderBy = "ks", topNodes = 200) %>% filter(Significant > Expected)


#Get the meanings for each GO term ID from the GO.db database, and save as a dataframe 
go.meanings <- suppressMessages(
  AnnotationDbi::select(GO.db, 
                        unique((tbl(my_db, "bee_go") %>% 
                                  dplyr::select(GO) %>% 
                                  collect(n=Inf) %>% 
                                  as.data.frame())[,1]), c("GOID", "ONTOLOGY", "TERM")))
names(go.meanings) <- c("GO", "ontology", "term")

# Make a table mapping Entrez gene names to Beebase gene names
entrez.tbl <- read.delim("./data/apis_gene_comparisons/am.gene_info.txt", stringsAsFactors = FALSE)[,c(2,5,6)]
names(entrez.tbl) <- c("entrez.id", "beebase1", "beebase2")
entrez.tbl$beebase1 <- str_extract(entrez.tbl$beebase1, "GB[:digit:]*")
entrez.tbl$beebase2 <- gsub("BEEBASE:", "", entrez.tbl$beebase2)
entrez.tbl$beebase2[entrez.tbl$beebase2 == "-"] <- NA
entrez.tbl <- entrez.tbl %>%
  filter(!(is.na(beebase1) & is.na(beebase2))) %>%
  filter(!is.na(entrez.id))
entrez.tbl <- tbl(my_db, "rsem_am") %>%
  select(gene) %>% collect(n=Inf) %>%
  left_join(entrez.tbl, by = c("gene" = "beebase2")) %>%
  select(gene, entrez.id)
entrez.tbl$entrez.id[is.na(entrez.tbl$entrez.id)] <- entrez.tbl$gene[is.na(entrez.tbl$entrez.id)]

# Get the KEGG annotations for bee genes (mapped to Entrez gene names)
bee_kegg_download <- clusterProfiler::download_KEGG("ame", keggType = "KEGG", keyType = "kegg")
kegg_meanings <- bee_kegg_download[[2]] %>% rename(kegg = from, name = to)
bee_kegg <- bee_kegg_download[[1]] %>% rename(kegg = from, entrez.id = to)

# Replace the Entrez names with Beebase names
bee_kegg <- bee_kegg %>% left_join(entrez.tbl, by = "entrez.id") %>% select(gene, kegg)

# Add/overwrite the new tables to the database
my_db$con %>% db_drop_table(table = "go_meanings")
my_db$con %>% db_drop_table(table = "kegg_meanings")
my_db$con %>% db_drop_table(table = "bee_kegg")
copy_to(my_db, go.meanings, "go_meanings", temporary = FALSE)
copy_to(my_db, kegg_meanings, "kegg_meanings", temporary = FALSE)
copy_to(my_db, bee_kegg, "bee_kegg", temporary = FALSE)