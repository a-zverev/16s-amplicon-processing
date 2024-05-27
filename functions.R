library(dada2)
library(Biostrings)
library(DECIPHER)
library(phyloseq)
library(ggplot2)
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
  if (exists("write_fasta")){
    writeXStringSet(ps_object@refseq, format = 'fasta', filepath = write_fasta)
  }
  
  # filter mitochondria and chloroplasts, if requied. NA protected by masking
  if (filter.organells == TRUE){
    
    # filter non-attributed
    ps_object <- subset_taxa(ps_object, Phylum != "NA")
    
    asvs.keep <-ps_object@tax_table %>% 
      data.frame() %>%  
      filter((Family != "Mitochondria" & Order != "Chloroplast") %>% 
      tidyr::replace_na(TRUE))  %>% 
      rownames()
    ps_object <- prune_taxa(asvs.keep, ps_object)
  }
  return(ps_object)
}


# Write ASVs table from ps-object
write_ASVs_table <- function(ps, filename){
  write.csv(cbind(ps@otu_table %>% t() %>%  data.frame(),
                  ps@tax_table %>% data.frame()),
            filename)
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

