library(dada2)
library(Biostrings)
library(DECIPHER)
library(phyloseq)
library(ggplot2)
library(ape)
library(plyr)
library(dplyr)
library(DESeq2)

dada2_processing_reads <- function(raw_files_path, cores=TRUE, trimLeft = c(19, 20)){
  
  #following dada2 pipeline
  fnFs <- sort(list.files(raw_files_path, pattern="_R1_001.fastq", full.names = TRUE))
  fnRs <- sort(list.files(raw_files_path, pattern="_R2_001.fastq", full.names = TRUE))
  sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
  
  print(sample.names)
  
  filtFs <- file.path(raw_files_path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
  filtRs <- file.path(raw_files_path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
  
  out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(220,180), 
                       trimLeft=trimLeft, maxN=0, maxEE=c(2,5), 
                       rm.phix=TRUE, compress=TRUE, multithread=cores)
  
  errF <- learnErrors(filtFs, multithread=cores)
  errR <- learnErrors(filtRs, multithread=cores)
  
  derepFs <- derepFastq(filtFs, verbose=TRUE)
  derepRs <- derepFastq(filtRs, verbose=TRUE)
  
  names(derepFs) <- sample.names
  names(derepRs) <- sample.names
  
  dadaFs <- dada(derepFs, err=errF, multithread=cores, pool=TRUE)
  dadaRs <- dada(derepRs, err=errR, multithread=cores, pool=TRUE)
  
  mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
  seqtab <- makeSequenceTable(mergers)
  
  dim(seqtab)
  table(nchar(getSequences(seqtab)))
  getN <- function(x) sum(getUniques(x))
  seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=cores, verbose=TRUE)
  
  return(seqtab.nochim)
}


dada2_assign_taxonomy <- function(seqtab, set_train_path, train_set_species_path, cores = TRUE){
  taxa.dada2 <- assignTaxonomy(seqtab, set_train_path , multithread=cores)
  taxa.dada2 <- addSpecies(taxa.dada2, train_set_species_path)
  
  return(taxa.dada2)
}


rename_by_metadata <- function(seqtab, metadata){
  rownames(metadata) <-metadata$SampleID
  metadata <- metadata %>% arrange(factor(Filename, levels = rownames(seqtab)))
  if (all(metadata$Filename == rownames(seqtab))) {
    rownames(seqtab) <- rownames(metadata)
  } else {print('Wrong column names')}
  return(seqtab)
}


create_ASV_references <- function(ps_object, write = TRUE){
  dna <- Biostrings::DNAStringSet(taxa_names(ps_object))
  names(dna) <- taxa_names(ps_object)
  ps_object <- merge_phyloseq(ps_object, dna)
  taxa_names(ps_object) <- paste0("ASV", seq(ntaxa(ps_object)))
  if (write == TRUE){
    writeXStringSet(ps_object@refseq, format = 'fasta', filepath = 'refseqs.fasta')
  }
  return(ps_object)
}


# Normalise data
normalize_phyloseq <- function(ps_object, method='rarefaction'){
  if (method == 'rarefaction'){ return(rarefy_even_depth(ps_object, verbose = FALSE)) }
  if (method == 'relative'){return(transform_sample_counts(ps_object, function(x) x/sum(x)))}
} 


# Draw barplot of relative abundance by taxa level
bargraph <- function(ps_object, rank, threshold){
  ph_ps <- tax_glom(ps_object, taxrank = rank)
  physeq3 = transform_sample_counts(ph_ps, function(x) x / sum(x) )
  
  data <- psmelt(physeq3) # create dataframe from phyloseq object
  data$Plot <- as.character(data[,rank]) #convert to character
  data$Plot[data$Abundance < threshold] <- paste0("<", threshold, " abund.")
  medians <- ddply(data, ~Plot, function(x) c(median=median(x$Abundance)))
  remainder <- medians[medians$median <= threshold,]$Plot
  p <- ggplot(data=data, aes(x=Sample, y=Abundance, fill=Plot))
  p + geom_bar(aes(), stat="identity", position="stack") + theme_light() +
    scale_fill_manual(values = c("darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "lightskyblue", "darkgreen", "deeppink", "khaki2", "firebrick", "brown1", "darkorange1", "cyan1", "royalblue4", "darksalmon", "darkblue",
                                 "royalblue4", "dodgerblue3", "steelblue1", "lightskyblue", "darkseagreen", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "brown1", "darkorange1", "cyan1", "darkgrey")) +
    theme(legend.position="bottom") + guides() +
    theme(axis.text.x = element_text(angle = 90))
}


# Calculate several alpha-diversity indexes, return one dataframe
alpha_div <- function(ps, metric, group){
  ps <- prune_taxa(taxa_sums(ps) > 0, ps)
  
  obs_sim <- estimate_richness(ps, split = TRUE, measures = metric)
  alpha <- cbind(sample_data(ps)[,group], obs_sim)
  return(alpha)
}


# Plot alpha-diversity by selected metric
plot_alpha <- function(ps, metric, group) {
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


# Gets table of significant differs from ps and formula
sig_table <- function(ps_object, formula){
  ds <- phyloseq_to_deseq2(ps_object, formula)
  ds = estimateSizeFactors(ds, type="poscounts")
  ds = estimateDispersions(ds, fitType = "local")
  ds = DESeq(ds)
  #mcols(ds, use.names=TRUE)
  res = results(ds, cooksCutoff = FALSE)
  alpha = 0.05
  sigtab = res[which(res$padj < alpha), ]
  if (nrow(sigtab) == 0) {
    return(NA)
  }
  sigtab = cbind(as(sigtab, "data.frame"),
                 as(tax_table(ps_object)[rownames(sigtab), ], "matrix")
  )
  return(sigtab)
}


draw_sig_table <- function(sig_table, taxa){
  ggplot(sig_table, aes(x = sig_table[,taxa], y=log2FoldChange, color=sig_table[,taxa])) + 
    geom_point(size = log(sig_table$baseMean)) + 
    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) +
    theme(legend.position = "none") + labs(x = taxa)
}


