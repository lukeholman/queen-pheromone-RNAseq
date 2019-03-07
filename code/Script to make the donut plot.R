# R script used to make the component parts of the "goldeneye" or "donut ring" in Figure 1

library(ggdendro)
library(ggplot2)
library(gridExtra)
library(ape)
options(expressions=10000)

# Define the dendrogram
hc <- network[[1]]$dendrograms[[1]]

# Make a table mapping the OGGs genes to the 9 modules, and sort by module
module_membership <- data.frame(module = network[[1]]$colors, 
                                gene = OGGs[[2]]$am, 
                                row.names = NULL,
                                stringsAsFactors = FALSE) %>% 
  arrange(module) %>%
  mutate(module = factor(module, unique(module)))
hc$labels <- NULL

# Make the dataframe holding information on each gene's module membership, for the outer-most donut
module_donut <- melt(table(module_membership$module)) %>%
  mutate(percent_share = 100 * value / sum(value),
         ymax = cumsum(percent_share),
         ymin = c(0, ymax[1:(n()-1)])) %>%
  rename(Module = Var1, nGenes = value) %>%
  mutate(Module = factor(Module, Module))

# Make a dataframe holding all the data for the 4 inner-most donuts
gene_donut <- rbind
(
  tbl(my_db, "ebseq_gene_am") %>% # Get the Apis log fold change EBSeq results
    select(gene, PostFC) %>% collect(n=Inf) %>% 
    filter(gene %in% OGGs[[2]]$am) %>% mutate(phero_sensitivity = log2(PostFC)) %>% select(-PostFC) %>% 
    left_join(module_membership, by = "gene") %>% 
    arrange(module) %>%
    rename(am = gene) %>%
    mutate(percent_share = 1 / nrow(OGGs[[2]]), 
           ymax = cumsum(percent_share), 
           ymin = c(0, ymax[1:(n()-1)]),
           xmin = 5, xmax = 6),
  
  tbl(my_db, "ebseq_gene_bt") %>% # Get the Bombus and Lasius log fold change EBSeq results, and match the gene names to their Apis orthologs
    select(gene, PostFC) %>% collect(n=Inf) %>% 
    filter(gene %in% OGGs[[2]]$bt) %>% 
    left_join(OGGs[[2]] %>% select(am, bt), by = c("gene" = "bt")) %>% select(-gene) %>%
    mutate(phero_sensitivity = log2(PostFC)) %>% select(-PostFC) %>% 
    left_join(module_membership, by = c("am" = "gene")) %>%
    arrange(module) %>%
    mutate(percent_share = 1 / nrow(OGGs[[2]]), 
           ymax = cumsum(percent_share), 
           ymin = c(0, ymax[1:(n()-1)]),
           xmin = 4, xmax = 5),
  
  tbl(my_db, "ebseq_gene_lf") %>% 
    select(gene, PostFC) %>% collect(n=Inf) %>% 
    filter(gene %in% OGGs[[2]]$lf) %>% 
    left_join(OGGs[[2]] %>% select(am, lf), by = c("gene" = "lf")) %>% select(-gene) %>%
    mutate(phero_sensitivity = log2(PostFC)) %>% select(-PostFC) %>% 
    left_join(module_membership, by = c("am" = "gene")) %>%
    arrange(module) %>%
    mutate(percent_share = 1 / nrow(OGGs[[2]]), 
           ymax = cumsum(percent_share), 
           ymin = c(0, ymax[1:(n()-1)]),
           xmin = 3, xmax = 4),
  
  tbl(my_db, "ebseq_gene_ln") %>% 
    select(gene, PostFC) %>% collect(n=Inf) %>% 
    filter(gene %in% OGGs[[2]]$ln) %>% 
    left_join(OGGs[[2]] %>% select(am, ln), by = c("gene" = "ln")) %>% select(-gene) %>%
    mutate(phero_sensitivity = log2(PostFC)) %>% select(-PostFC) %>% 
    left_join(module_membership, by = c("am" = "gene")) %>%
    arrange(module) %>%
    mutate(percent_share = 1 / nrow(OGGs[[2]]), 
           ymax = cumsum(percent_share), 
           ymin = c(0, ymax[1:(n()-1)]),
           xmin = 2, xmax = 3)  
) %>% 
  mutate(module = as.numeric(gsub("Module ", "", module)),
         species = unlist(mapply(rep, c("am", "bt", "lf", "ln"), 
                                 each = c(3465, 3465, 3465, 3462)))) %>% 
  as.data.frame() 

# Plot the inner donut rings showing pheromone sensitivity
log_fold_rings <- gene_donut %>%
  mutate(phero_sensitivity_category = 
           cut(abs(phero_sensitivity), c(0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, Inf))) %>% # Round off data to categories for faster plotting (still very slow)
  ggplot() + 
  geom_rect(aes(fill=phero_sensitivity_category, ymin=ymin, ymax=ymax, xmin=xmin, xmax=xmax)) + # draw the 4 inner donuts
  geom_vline(xintercept = 2) + geom_vline(xintercept = 3) + geom_vline(xintercept = 4) + # add black concentric lines
  geom_vline(xintercept = 5) + geom_vline(xintercept = 6) +
  coord_polar(theta="y") + xlim(c(0, 7)) + theme_dendro() + 
  scale_fill_brewer(palette = "YlOrRd", type = "div") +
  theme(legend.position = "none")

# Plot the outer donut ring showing modules
module_ring <- module_donut %>%  
  mutate(Module = factor(Module, paste("Module", c(1:9,0)))) %>%
  ggplot() + geom_rect(aes(fill = Module, ymin=ymin, ymax=ymax, xmin=6.2, xmax=7)) +
  geom_segment(x = 6.2, xend = 7, aes(y = ymax, yend = ymax)) + # Add radial black lines between the modules
  coord_polar(theta="y") + xlim(c(0, 7)) + theme_dendro() +
  theme(legend.position = "none")

# Plot the dendrogram
dendro <- OGGs[[1]] %>% remove.effects.combat() %>% .[[1]] %>% t() %>% # Get the gene expression data
  dist() %>% hclust()
dendro$labels <- as.numeric(gsub("Module ", "", module_membership$module[match(dendro$labels, module_membership$gene)]))
start <- 1
for(i in 0:9) { # Sort the tips of the phylogeny by module, as in the donut ring data
  length_foc <- sum(dendro$labels == i)
  dendro$order[dendro$labels == i] <- start:(start+length_foc-1)
  start <- start + length_foc
}

# Save results to disk for editing using Inkscape or Illustrator (I couldn't find a good way to combine all of these in R)
donuts <- grid.arrange(log_fold_rings, module_ring) 
ggsave(donuts, file = "figures/donuts.svg")
svg(file = "figures/phylo.svg")
plot(as.phylo(dendro), type = "fan", edge.width = 1.5, show.tip.label = FALSE)
dev.off()