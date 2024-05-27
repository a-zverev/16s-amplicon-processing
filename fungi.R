require(dada2)
require(Biostrings)
require(DECIPHER)
require(phyloseq)
require(ggplot2)
library(stringr)
require(seqinr)
require(plyr)
library(ape)
library(phangorn)

setwd('/home/alexey/Analysis/RW/Field_RW/fungi/R')

path <- '../raw'
list.files

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)


#Primers and orientations (btw, this is ITS2)
FWD <- "GCATCGATGAAGAACGCAGC" #its3
REV <- "TCCTCCGCTTATTGATATGC"  #its4

allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients

fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE)

#Cut the primers
cutadapt <- '/home/alexey/.conda/envs/py_work/bin/cutadapt'
path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC) 
# Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i])) # input files
}

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))

# Forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(path.cut, pattern = "_R1_001.fastq", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "_R2_001.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format:
get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))
sample.names


# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

seqtab <- makeSequenceTable(mergers)
dim(seqtab)

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

saveRDS(seqtab.nochim, 'seqtab.nochim.RData')
seqtab.nochim <- readRDS('seqtab.nochim.RData')

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, 
                                                                       getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", 
                     "nonchim")
rownames(track) <- sample.names


#########################################################################

## Read map, make a df to rename asv matrix and rename seqtab.nochim via cbind()
mdat <- read.csv2('../new_mapITS.csv', sep = '\t', header = T, row.names = 1)
pairs <- DataFrame(row.names = mdat$Filename, Names = rownames(mdat))
write.csv(merge(pairs, track, by=0, all=T), file = "processing.log")

rownames(seqtab.nochim) <- sapply(strsplit(basename(rownames(seqtab.nochim)), "_"), `[`, 1)
tab <- merge(pairs, seqtab.nochim, by=0, all=T)
rownames(tab) <- tab$Names
tab <- tab[, -(1:2)]
matr.tab <- data.matrix(tab)

#Make first phyloseq object,
ps <- phyloseq(otu_table((matr.tab), taxa_are_rows=F),
               sample_data(mdat))


#Make reference reads for taxa assignment
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

# Taxonomy? 
taxa <- assignTaxonomy(ps@refseq, "/home/alexey/tax_n_refs/unite/sh_general_release_dynamic_s_04.02.2020.fasta", multithread=TRUE)
rownames(taxa) <- colnames(ps@otu_table)
ps <- merge_phyloseq(ps, tax_table(taxa))


#Phylogeny? By fasttree from QIIME2-plugin
writeXStringSet(ps@refseq, format = 'fasta', filepath = 'refseqs.fasta')

#In QIIME2 run this:

#  qiime tools import --input-path refseqs.fasta --output-path sequences.qza --type 'FeatureData[Sequence]'
#  qiime alignment mafft --i-sequences sequences.qza --o-alignment aligned-rep-seqs.qza
#  qiime alignment mask --i-alignment aligned-rep-seqs.qza --o-masked-alignment masked-aligned-rep-seqs.qza

#  qiime phylogeny fasttree --i-alignment masked-aligned-rep-seqs.qza --o-tree masked-fasttree.qza
#  extract them from archives

m.fasttree <- read_tree('tree.nwk')
ps <- merge_phyloseq(ps, phy_tree(m.fasttree))

#################################################################################

ps@tax_table[,'Phylum'] <- sapply(strsplit(basename(matrix(ps@tax_table[,'Phylum'])), "_"), `[`, 3)
ps@tax_table[,'Class'] <- sapply(strsplit(basename(matrix(ps@tax_table[,'Class'])), "_"), `[`, 3)
ps@tax_table[,'Order'] <- sapply(strsplit(basename(matrix(ps@tax_table[,'Order'])), "_"), `[`, 3)
ps@tax_table[,'Family'] <- sapply(strsplit(basename(matrix(ps@tax_table[,'Family'])), "_"), `[`, 3)
ps@tax_table[,'Genus'] <- sapply(strsplit(basename(matrix(ps@tax_table[,'Genus'])), "_"), `[`, 3)

ps = subset_taxa(ps, Phylum != 'NA')

#################################################################################

#save ps-object
saveRDS(ps, "fungi_ps.RData")

####################################################################

fungi.PS <- subset_samples(ps, Soil %in% "PS")
fungi.PS <- filter_taxa(fungi.PS, function(x) mean(x) > 0, TRUE)

fungi.BS <- subset_samples(ps, Soil %in% "BS")
fungi.BS <- filter_taxa(fungi.BS, function(x) mean(x) > 0, TRUE)



#Save otu-table and taxa from ps-object
save_otu_table <- function(ps_object){
  otu_df = as.data.frame(t(as.matrix(ps_object@otu_table)))
  taxa = as.data.frame(as(tax_table(ps_object), "matrix"))
  res <- cbind(otu_df, taxa)
  write.table(res, file = paste(deparse(substitute(ps_object)), "_otu_table.txt", sep=""), 
              sep = '\t', quote = FALSE, col.names=NA)
}

save_otu_table(fungi.BS)
save_otu_table(fungi.PS)

saveRDS(fungi.BS, "fungi.BS.RData")
saveRDS(fungi.PS, "fungi.PS.RData")

#######################################################################