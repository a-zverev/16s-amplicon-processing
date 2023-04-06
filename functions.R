library(dada2)
library(Biostrings)
library(DECIPHER)
library(phyloseq)
library(ggplot2)
library(data.table)
library(ape)
library(dplyr)
library(DESeq2)

# Dada2 pipeline follows to sequence table
reads_to_seqtable_by_dada2 <- function(raw_files_path, trimLeft, truncLen,
                                   pool=TRUE, cores=TRUE){
  require(dada2)
  
  # following dada2 pipeline
  fnFs <- sort(list.files(raw_files_path, pattern="_R1_001.fastq", full.names = TRUE))
  fnRs <- sort(list.files(raw_files_path, pattern="_R2_001.fastq", full.names = TRUE))
  sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
  
  # show some info about data
  print(sample.names)
  print(plotQualityProfile(fnRs, aggregate = T))
  print(plotQualityProfile(fnRs, aggregate = T))
  
  filtFs <- file.path(raw_files_path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
  filtRs <- file.path(raw_files_path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
  
  # trim data
  out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=truncLen, 
                       trimLeft=trimLeft, maxN=0, maxEE=c(2,5), 
                       rm.phix=TRUE, compress=TRUE, multithread=cores)
  
  # error rates
  errF <- learnErrors(filtFs, multithread=cores)
  errR <- learnErrors(filtRs, multithread=cores)
  
  # dereplication
  derepFs <- derepFastq(filtFs, verbose=F)
  derepRs <- derepFastq(filtRs, verbose=F)
  
  names(derepFs) <- sample.names
  names(derepRs) <- sample.names
  
  # apply error rate and merge
  dadaFs <- dada(derepFs, err=errF, multithread=cores, pool=pool)
  dadaRs <- dada(derepRs, err=errR, multithread=cores, pool=pool)
  mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
  seqtab <- makeSequenceTable(mergers)
  
  # get some stats
  dim(seqtab)
  table(nchar(getSequences(seqtab)))
  getN <- function(x) sum(getUniques(x))
  seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=cores, verbose=TRUE)
  
  # logs
  
  getN <- function(x) sum(getUniques(x))
  track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
  colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
  rownames(track) <- sample.names
  print(track)
  write.csv(track, "processing.log")
  
  return(seqtab.nochim)
}

# Read metadata from mapfile
read_metadata <- function(filename, sample.names, ...){
  
  # read metadata and pass names of samples to rownames
  metadata <- read.csv(filename, ...)
  rownames(metadata) <- metadata[,sample.names]
  metadata
}

# Assign taxonomy
assign_taxonomy <- function(seqtab, set_train_path, train_set_species_path, cores = TRUE){
  require(dada2)
  
  taxa.dada2 <- assignTaxonomy(seqtab, set_train_path , multithread=cores)
  taxa.dada2 <- addSpecies(taxa.dada2, train_set_species_path)
  return(taxa.dada2)
}

# Rename sequence table by our names
rename_seqtab_by_metadata <- function(seqtab, metadata, old.names){
  require(dplyr)
  
  # sort metadata by %seqtab% rownames order
  metadata <- metadata %>% arrange(factor(metadata[,old.names], levels = rownames(seqtab)))
  
  # if names of %seqtab% and %old.names% equal and in correct order, rename
  if (all(metadata[,old.names] == rownames(seqtab))) {
    rownames(seqtab) <- rownames(metadata)
  } else {print('Wrong column names')}
  return(seqtab)
}

# Combine a plyloseq object (without a tree)
assemble_phyloseq <- function(seqtab, metadata, taxonomy, filter.organells = T, write_fasta = TRUE){
  require(Biostrings)
  require(phyloseq)
  require(dplyr)
  
  # create phyloseq
  ps_object <- phyloseq(otu_table(seqtab, taxa_are_rows=FALSE), 
                 sample_data(metadata), 
                 tax_table(taxa))
  
  # move references from %ps_object% taxa names to dedicated data
  dna <- Biostrings::DNAStringSet(taxa_names(ps_object))
  names(dna) <- taxa_names(ps_object)
  ps_object <- merge_phyloseq(ps_object, dna)
  taxa_names(ps_object) <- paste0("ASV", seq(ntaxa(ps_object)))
  
  # write fastas, if required
  if (write_fasta == TRUE){
    writeXStringSet(ps_object@refseq, format = 'fasta', filepath = 'refseqs.fasta')
  }
  
  # filter mitochondria and chloroplasts, if requied. NA protected by masking
  if (filter.organells == TRUE){
    
    # filter non-attributed
    ps_object <- subset_taxa(ps_object, Phylum != "NA")
    
    ps_object@tax_table[is.na(ps_object@tax_table)] <- TRUE
    ps_object <- subset_taxa(ps_object,
                      !(Family  == "Mitochondria" |
                          Class   == "Chloroplast" |
                          Order   == "Chloroplast"))
    ps_object@tax_table <- dplyr::na_if(ps_object@tax_table, TRUE)
  }
  return(ps_object)
}

# # Root the tree and add to ps object
# add_tree_to_ps <- function(ps, tree){
#   require(phyloseq)
#   require(ape)
#   require(data.table)
#   
#   # find longest branch and return it's label
#   pick_new_outgroup <- function(tree.unrooted){
#     # tablify parts of tree that we need.
#     treeDT <- 
#       cbind(
#         data.table(tree.unrooted$edge),
#         data.table(length = tree.unrooted$edge.length)
#       )[1:Ntip(tree.unrooted)] %>% 
#       cbind(data.table(id = tree.unrooted$tip.label))
#     # Take the longest terminal branch as outgroup
#     new.outgroup <- treeDT[which.max(length)]$id
#     return(new.outgroup)
#   }
#   
#   m.fasttree <- read_tree(tree)
#   m.fasttree <- ape::root(m.fasttree, outgroup=pick_new_outgroup(m.fasttree), resolve.root=TRUE)
#   
#   ps <- merge_phyloseq(ps, phy_tree(m.fasttree))
# }

# Write ASVs table from ps-object
write_ASVs_table <- function(ps_object, filename){
  write.csv(cbind(ps_object@otu_table %>% t() %>%  data.frame(),
                  ps_object@tax_table %>% data.frame()),
            filename)
}

# Draw barplot of relative abundance by taxa level
bargraph <- function(ps, rank, threshold=0.05, percents=FALSE){
  require(dplyr)
  require(ggplot2)
  require(phyloseq)
  
  ps <- prune_taxa(taxa_sums(ps) > 0, ps)
  ps2 <- tax_glom(ps, taxrank = rank)
  ps3 = transform_sample_counts(ps2, function(x) x / sum(x) )
  data <- psmelt(ps3) # create dataframe from phyloseq object
  data$Plot <- as.character(data[,rank]) # convert to character
  data$Plot[data$Abundance < threshold] <- paste0("<", threshold, " abund.")
  medians <- data %>% group_by(Plot) %>% mutate(median=median(data$Abundance))
  remainder <- medians[medians$median <= threshold,]$Plot
  data$Percentage = ifelse(data$Plot != paste0("<", threshold, " abund."),
                           round(data$Abundance, 3)*100, NA)
  
  # create palette long enough for our data
  base.palette <- c("darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "lightskyblue", 
                    "darkgreen", "deeppink", "khaki2", "firebrick", "brown1", "darkorange1", "cyan1", "royalblue4", 
                    "darksalmon", "dodgerblue3", "steelblue1", "darkgoldenrod1", "brown1", "cyan1", "darkgrey")
  required.colors <- nlevels(factor(data$Plot))
  repeats = required.colors %/% length(base.palette) + 1
  palette <- rep(base.palette, length.out = repeats * length(base.palette))
  
  p <- ggplot(data=data, aes(x=Sample, y=Abundance, fill=Plot))
  p + geom_bar(aes(), stat="identity", position="stack") + theme_light() +
    scale_fill_manual(values = palette) +
    theme(legend.position="bottom") + guides() +
    theme(axis.text.x = element_text(angle = 90)) +
    if (percents) {
      geom_text(aes(label = Percentage),
                position = position_stack(vjust = 0.5), size = 1.5)
    }
  
}

# Calculate several alpha-diversity indexes, return one dataframe
alpha_div_table <- function(ps, metric, group){
  require(phyloseq)
  
  ps <- prune_taxa(taxa_sums(ps) > 0, ps)
  
  obs_sim <- estimate_richness(ps, split = TRUE, measures = metric)
  alpha <- cbind(sample_data(ps)[,group], obs_sim)
  return(alpha)
}

# Plot alpha-diversity by selected metric
plot_alpha <- function(ps, group, metric) {
  require(phyloseq)
  require(ggplot2)
  
  ps_a <- prune_taxa(taxa_sums(ps) > 0, ps)
  plot_richness(ps_a, x=group, measures=metric) + 
    geom_boxplot() +
    geom_point(size=1.2, alpha=0.3) + 
    theme_light() + scale_color_brewer(palette="Dark2") +
    theme(axis.text.x = element_text(angle = 90),
          axis.title.x=element_blank()) +
    labs(y=paste(metric, "index")) +
    ggtitle(metric)
}

# Plot beta-diversity
beta_plot <- function(ps, method, distance, ...){
  require(phyloseq)
  require(ggplot2)
  
  ps.prop <- transform_sample_counts(ps, function(x) x/sum(x))
  ords <- ordinate(ps.prop, method=method, distance=distance)
  plot_ordination(ps.prop, ords, title=deparse(substitute(ps)), ...) +
    geom_point(size=3, alpha=0.7) + 
    theme_light()
}

# Gets table of significant differs from ps and formula
sig_table <- function(ps, formula, threshold){
  require(DESeq2)
  require(dplyr)
  
  ds <- phyloseq_to_deseq2(ps, formula)
  ds = estimateSizeFactors(ds, type="poscounts")
  ds = estimateDispersions(ds, fitType = "local")
  ds = DESeq(ds)
  #mcols(ds, use.names=TRUE)
  res = results(ds, cooksCutoff = FALSE)
  alpha = 0.05
  sigtab = res[which(res$padj < alpha), ]
  sigtab <- sigtab %>% data.frame() %>% dplyr::filter(baseMean >= threshold[1] &
                                                        abs(log2FoldChange) >= threshold[2])
  if (nrow(sigtab) == 0) {
    return(NA)
  }
  sigtab = cbind(as(sigtab, "data.frame"),
                 as(tax_table(ps)[rownames(sigtab), ], "matrix")
  )
  return(sigtab)
}

# Draw a plot by %sig_table% data
draw_sig_table <- function(sig_table, rank){
  require(ggplot2)
  
  sig_table$Class <-  ifelse(is.na(sig_table$Class), 
                             paste(sig_table$Phylum, "// NA"), 
                             paste(sig_table$Class))
  sig_table$Order <-  ifelse(is.na(sig_table$Order), 
                             paste(sig_table$Class, "// NA"), 
                             paste(sig_table$Order))
  sig_table$Family <- ifelse(is.na(sig_table$Family), 
                             paste(sig_table$Order, "// NA"), 
                             paste(sig_table$Family))
  sig_table$Genus <- ifelse(is.na(sig_table$Genus), 
                            paste(sig_table$Class, "// NA"), 
                            paste(sig_table$Genus))
  
  sig_table[sig_table == "Burkholderia-Caballeronia-Paraburkholderia"] <- "Burkholderia"
  sig_table[sig_table == "Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium"] <- "Pararhizobium"
  
  ggplot(sig_table, aes(y = sig_table[,rank], x=log2FoldChange, size = log(baseMean))) + 
    geom_point(aes(color=sig_table[,rank])) + guides(colour=FALSE) +
    theme_light() + facet_grid(rows = vars(Phylum), space = 'free_y', scales = 'free') +
    theme(strip.text.y = element_text(angle=0)) +
    theme(legend.position = "top") + ylab(rank)
}

# Plot heatmap on selected
plot_heatmap <- function(ps, group = "SampleID", log.transform = TRUE){
  require(dplyr)
  require(phyloseq)
  require(ggplot2)
  
  ps <- prune_taxa(taxa_sums(ps) > 0, ps)
  ps <- tax_glom(ps, "Genus")
  sig.taxa.long <- psmelt(ps) %>%
    arrange(Phylum) %>% 
    mutate(row = row_number()) %>% 
    mutate_at(c('Genus', 'Family', 'Order', 'Class', 'Phylum'), as.character)  # fuckin magic for ifelse tier
  
  sig.taxa.long$Abundance <- as.numeric(sig.taxa.long$Abundance)
  sig.taxa.long$Taxa <- ifelse(is.na(sig.taxa.long$Genus),
                               ifelse(is.na(sig.taxa.long$Family), 
                                      ifelse(is.na(sig.taxa.long$Order), 
                                             ifelse(is.na(sig.taxa.long$Class), 
                                                    paste(sig.taxa.long$Phylum, "// NA"), 
                                                    paste(sig.taxa.long$Class, "// NA")),
                                             paste(sig.taxa.long$Order, "// NA")),
                                      paste(sig.taxa.long$Family, "// NA")), 
                               sig.taxa.long$Genus)
  sig.taxa.long[sig.taxa.long == "Burkholderia-Caballeronia-Paraburkholderia"] <- "Burkholderia"
  sig.taxa.long[sig.taxa.long == "Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium"] <- "Pararhizobium"
  
  
  ggplot(sig.taxa.long, aes(x = sig.taxa.long[,group], y = reorder(Taxa, row))) + 
    {if(log.transform)geom_tile(aes(fill=log(Abundance)))} +
    {if(!log.transform)geom_tile(aes(fill=Abundance))} +
    scale_fill_distiller(palette = "OrRd", trans = "reverse") +
    facet_grid(rows = vars(Phylum), scales = "free_y", space = "free") +
    theme(strip.text.y = element_text(angle = 0),
          panel.spacing = unit(0.05,'lines')) +
    xlab("") + ylab("") +
    theme(axis.text.x = element_text(angle = 90))
}