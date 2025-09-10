library(dada2)
library(Biostrings)
library(DECIPHER)
library(phyloseq)
library(ggplot2)
library(dplyr)
library(DESeq2)

# Dada2 pipeline follows to sequence table
reads_to_seqtable_16s <- function(raw_files_path,
                                  trimLeft,
                                  truncLen,
                                  log_file,
                                  pool = TRUE,
                                  cores = TRUE,
                                  maxEE = c(2, 5),
                                  pattern = c("_R1_001.fastq", "_R2_001.fastq"),
                                  plot_profile = TRUE) {
  require(dada2)
  
  # List forward and reverse read files using provided pattern
  fnFs <- sort(list.files(raw_files_path, pattern = pattern[1], full.names = TRUE))
  fnRs <- sort(list.files(raw_files_path, pattern = pattern[2], full.names = TRUE))
  
  # Extract sample names from filenames
  sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
  
  # Print sample names
  print(sample.names)
  
  # Plot quality profiles if enabled
  if (plot_profile) {
    plotQualityProfile(fnFs, aggregate = TRUE)
    plotQualityProfile(fnRs, aggregate = TRUE)
  }
  
  # Define filtered file paths
  filtFs <- file.path(raw_files_path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
  filtRs <- file.path(raw_files_path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
  
  # Filter and trim reads
  out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen = truncLen,
                       trimLeft = trimLeft, maxN = 0, maxEE = maxEE,
                       rm.phix = TRUE, compress = TRUE, multithread = cores)
  
  # Learn error rates from filtered reads
  errF <- learnErrors(filtFs, multithread = cores)
  errR <- learnErrors(filtRs, multithread = cores)
  
  # Dereplicate filtered reads
  derepFs <- derepFastq(filtFs, verbose = FALSE)
  derepRs <- derepFastq(filtRs, verbose = FALSE)
  
  names(derepFs) <- sample.names
  names(derepRs) <- sample.names
  
  # Apply dada algorithm and merge paired reads
  dadaFs <- dada(derepFs, err = errF, multithread = cores, pool = pool)
  dadaRs <- dada(derepRs, err = errR, multithread = cores, pool = pool)
  mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose = TRUE)
  seqtab <- makeSequenceTable(mergers)
  
  # Remove chimeras
  seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus",
                                      multithread = cores, verbose = TRUE)
  
  # Track reads through the pipeline
  getN <- function(x) sum(getUniques(x))
  track <- cbind(out,
                 sapply(dadaFs, getN),
                 sapply(dadaRs, getN),
                 sapply(mergers, getN),
                 rowSums(seqtab.nochim))
  colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
  rownames(track) <- sample.names
  
  print(track)
  
  # Write log to file (log_file is mandatory)
  write.csv(track, log_file)
  
  return(seqtab.nochim)
}

trim_ITS_primers <- function(raw_files_path, FWD, REV, cutadapt) {
  require(dada2)
  require(ShortRead)
  require(Biostrings)
  
  # Locate all FASTQ pairs
  fnFs <- sort(list.files(raw_files_path, pattern = "_1\\.f(ast)?q(.gz)?$", full.names = TRUE))
  fnRs <- sort(list.files(raw_files_path, pattern = "_2\\.f(ast)?q(.gz)?$", full.names = TRUE))
  
  if (length(fnFs) == 0 || length(fnRs) == 0) stop("No matching FASTQ files found in path.")
  if (length(fnFs) != length(fnRs)) stop("Mismatched number of forward and reverse reads.")
  
  # Filter N-containing reads for all samples
  filt_path <- file.path(raw_files_path, "filtN")
  dir.create(filt_path, showWarnings = FALSE)
  fnF_filt <- file.path(filt_path, basename(fnFs))
  fnR_filt <- file.path(filt_path, basename(fnRs))
  
  message(">> Filtering all reads with Ns...")
  filterAndTrim(fnFs, fnF_filt, fnRs, fnR_filt, maxN = 0, multithread = TRUE, verbose = FALSE)
  
  # Use only first sample to assess primer orientation
  testF <- fnF_filt[1]
  testR <- fnR_filt[1]
  
  allOrients <- function(primer) {
    dna <- DNAString(primer)
    orients <- c(Forward = dna,
                 Complement = complement(dna),
                 Reverse = reverse(dna),
                 RevComp = reverseComplement(dna))
    return(sapply(orients, toString))
  }
  FWD.orients <- allOrients(FWD)
  REV.orients <- allOrients(REV)
  
  primerHits <- function(primer, fn) {
    nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
    return(sum(nhits > 0))
  }
  
  message(">> Primer hits BEFORE trimming (1st sample):")
  print(rbind(
    FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = testF),
    FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = testR),
    REV.ForwardReads = sapply(REV.orients, primerHits, fn = testF),
    REV.ReverseReads = sapply(REV.orients, primerHits, fn = testR)
  ))
  
  # Prepare cutadapt output dir
  outdir <- file.path(raw_files_path, "cutadapt")
  dir.create(outdir, showWarnings = FALSE)
  
  rc_FWD <- dada2::rc(FWD)
  rc_REV <- dada2::rc(REV)
  
  # Run cutadapt on all files
  message(">> Running cutadapt on all samples...")
  for (i in seq_along(fnFs)) {
    outF <- file.path(outdir, basename(fnFs[i]))
    outR <- file.path(outdir, basename(fnRs[i]))
    
    cmd <- paste(
      shQuote(cutadapt),
      "-g", shQuote(FWD),
      "-G", shQuote(REV),
      "-a", shQuote(rc_REV),
      "-A", shQuote(rc_FWD),
      "-o", shQuote(outF),
      "-p", shQuote(outR),
      shQuote(fnF_filt[i]),
      shQuote(fnR_filt[i])
    )
    
    system(paste(cmd, "> /dev/null 2>&1"))
  }
  
  # Check primer hits after trimming (again only for first sample)
  trimmedF <- file.path(outdir, basename(fnFs[1]))
  trimmedR <- file.path(outdir, basename(fnRs[1]))
  
  message(">> Primer hits AFTER trimming (1st sample):")
  print(rbind(
    FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = trimmedF),
    FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = trimmedR),
    REV.ForwardReads = sapply(REV.orients, primerHits, fn = trimmedF),
    REV.ReverseReads = sapply(REV.orients, primerHits, fn = trimmedR)
  ))
  
  message(">> Done. Trimmed files are in: ", outdir)
  return(outdir)
}

reads_to_seqtable_ITS <- function(raw_files_path,
                                  trimLeft = c(0, 0),
                                  truncLen = c(0, 0),
                                  log_file,
                                  pool = TRUE,
                                  cores = TRUE,
                                  maxEE = c(2, 5),
                                  pattern = c("_R1_001.fastq", "_R2_001.fastq"),
                                  plot_profile = TRUE,
                                  justConcatenate = TRUE,
                                  remove_chimeras = FALSE) {
  require(dada2)
  
  # 1. List reads
  fnFs <- sort(list.files(raw_files_path, pattern = pattern[1], full.names = TRUE))
  fnRs <- sort(list.files(raw_files_path, pattern = pattern[2], full.names = TRUE))
  sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
  
  if (plot_profile) {
    print(plotQualityProfile(fnFs, aggregate = TRUE))
    print(plotQualityProfile(fnRs, aggregate = TRUE))
  }
  
  # 2. Filtering
  filtFs <- file.path(raw_files_path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
  filtRs <- file.path(raw_files_path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
  out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                       truncLen = truncLen,
                       trimLeft = trimLeft,
                       maxN = 0, maxEE = maxEE,
                       rm.phix = TRUE,
                       compress = TRUE, multithread = cores)
  
  # 3. Learn errors
  errF <- learnErrors(filtFs, multithread = cores)
  errR <- learnErrors(filtRs, multithread = cores)
  
  # 4. Dereplication
  derepFs <- derepFastq(filtFs, verbose = FALSE)
  derepRs <- derepFastq(filtRs, verbose = FALSE)
  names(derepFs) <- sample.names
  names(derepRs) <- sample.names
  
  # 5. Denoising
  dadaFs <- dada(derepFs, err = errF, multithread = cores, pool = pool)
  dadaRs <- dada(derepRs, err = errR, multithread = cores, pool = pool)
  
  # 6. Merging with optional concatenation
  mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs,
                        verbose = TRUE, justConcatenate = justConcatenate)
  
  seqtab <- makeSequenceTable(mergers)
  
  # 7. Chimera removal (optional)
  if (remove_chimeras) {
    seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = cores)
  } else {
    seqtab.nochim <- seqtab
  }
  
  # 8. Track pipeline
  getN <- function(x) sum(getUniques(x))
  track <- cbind(out,
                 sapply(dadaFs, getN),
                 sapply(dadaRs, getN),
                 sapply(mergers, getN),
                 rowSums(seqtab.nochim))
  colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
  rownames(track) <- sample.names
  
  print(track)
  write.csv(track, log_file)
  
  return(seqtab.nochim)
}

# Read metadata from mapfile
read_metadata <- function(filename, sample.names, sep = ",", check_files = NULL, ...) {
  # Read the metadata file with specified separator and additional arguments
  metadata <- read.csv(filename, sep = sep, ...)
  
  # Check if the specified sample.names column exists
  if (!sample.names %in% colnames(metadata)) {
    stop(paste("Column", sample.names, "not found in the file."))
  }
  
  # Extract sample IDs and check for uniqueness
  sample_ids <- metadata[[sample.names]]
  if (anyDuplicated(sample_ids)) {
    stop("Values in the 'sample.names' column are not unique.")
  }
  
  # Set row names to sample IDs
  rownames(metadata) <- sample_ids
  
  # If a directory is provided, check for file presence per sample
  if (!is.null(check_files)) {
    if (!dir.exists(check_files)) {
      stop(paste("Directory", check_files, "does not exist."))
    }
    
    # List all files in the directory
    all_files <- list.files(check_files)
    
    # Check that each sample has 2 associated files
    missing_or_incomplete <- sample_ids[sapply(sample_ids, function(id) {
      sum(startsWith(all_files, id)) != 2
    })]
    
    if (length(missing_or_incomplete) > 0) {
      warning(paste(
        "The following samples do not have 2 corresponding files in the directory:",
        paste(missing_or_incomplete, collapse = ", ")
      ))
    }
  }
  
  return(metadata)
}

# Assign taxonomy
assign_taxonomy <- function(seqtab, set_train_path, cores = TRUE){
  require(dada2)
  
  taxa.dada2 <- assignTaxonomy(seqtab, set_train_path , multithread=cores)
  return(taxa.dada2)
}

# Combine a plyloseq object (without a tree)
assemble_phyloseq_16s <- function(seqtab, metadata, taxonomy, 
                                  filter.organelles = TRUE, write_fasta = NA) {
  library(Biostrings)
  library(phyloseq)
  library(dplyr)
  library(tidyr)
  
  # Create a phyloseq object from OTU table, metadata, and taxonomy
  ps_object <- phyloseq(
    otu_table(seqtab, taxa_are_rows = FALSE),
    sample_data(metadata),
    tax_table(taxonomy)
  )
  
  # Create a DNAStringSet from taxa names and attach it to phyloseq object
  dna <- DNAStringSet(taxa_names(ps_object))
  names(dna) <- taxa_names(ps_object)
  ps_object <- merge_phyloseq(ps_object, dna)
  
  # Rename ASVs to standard names: ASV1, ASV2, ...
  taxa_names(ps_object) <- paste0("ASV", seq(ntaxa(ps_object)))
  
  # Log the initial number of ASVs
  n_initial <- ntaxa(ps_object)
  message("Initial number of ASVs: ", n_initial)
  
  # Filter out ASVs with NA in Phylum
  taxonomy_df <- as.data.frame(tax_table(ps_object))
  asvs_no_na_phylum <- taxonomy_df %>%
    filter(!is.na(Phylum)) %>%
    rownames()
  ps_object <- prune_taxa(asvs_no_na_phylum, ps_object)
  n_after_phylum <- ntaxa(ps_object)
  message("ASVs removed due to NA in Phylum: ", n_initial - n_after_phylum)
  
  # Filter to keep only Bacteria and Archaea
  taxonomy_df <- as.data.frame(tax_table(ps_object))  # update taxonomy
  asvs_domain_filtered <- taxonomy_df %>%
    filter(Kingdom %in% c("Bacteria", "Archaea")) %>%
    rownames()
  ps_object <- prune_taxa(asvs_domain_filtered, ps_object)
  n_after_kingdom <- ntaxa(ps_object)
  message("ASVs removed due to non-Bacteria/Archaea Kingdom: ", 
          n_after_phylum - n_after_kingdom)
  
  # Optionally filter out mitochondrial and chloroplast sequences
  if (isTRUE(filter.organelles)) {
    taxonomy_df <- as.data.frame(tax_table(ps_object))  # update taxonomy again
    asvs_keep <- taxonomy_df %>%
      filter(
        is.na(Family) | Family != "Mitochondria",
        is.na(Order)  | Order  != "Chloroplast"
      ) %>%
      rownames()
    
    ps_object <- prune_taxa(asvs_keep, ps_object)
    n_final <- ntaxa(ps_object)
    message("ASVs removed as organelles (Mitochondria/Chloroplast): ", 
            n_after_kingdom - n_final)
  } else {
    n_final <- ntaxa(ps_object)
  }
  
  message("Final number of ASVs: ", n_final)
  
  # Write FASTA file if a valid file path is provided (not NA)
  if (!is.na(write_fasta) && is.character(write_fasta)) {
    writeXStringSet(ps_object@refseq, format = "fasta", filepath = write_fasta)
    message("FASTA file written to: ", write_fasta)
  }
  
  return(ps_object)
}

assemble_phyloseq_ITS <- function(seqtab, metadata, taxonomy) {
  library(Biostrings)
  library(phyloseq)
  library(dplyr)
  library(tidyr)
  
  # Remove UNITE-style prefixes (like k__, p__, etc.) if present
  taxonomy_df <- as.data.frame(taxonomy)
  has_unite_prefixes <- any(grepl("^[a-z]__", unlist(taxonomy_df)))
  
  if (has_unite_prefixes) {
    message("UNITE-style prefixes detected. Removing prefixes (e.g., 'k__', 'p__', etc.)...")
    taxonomy_df <- taxonomy_df %>%
      mutate(across(everything(), ~ sub("^[a-z]__+", "", .)))
  }
  
  # Reconstruct taxonomy as matrix
  taxonomy_clean <- as.matrix(taxonomy_df)
  
  # Create phyloseq object
  ps_object <- phyloseq(
    otu_table(seqtab, taxa_are_rows = FALSE),
    sample_data(metadata),
    tax_table(taxonomy_clean)
  )
  
  # Add DNA sequences as refseq
  dna <- Biostrings::DNAStringSet(taxa_names(ps_object))
  names(dna) <- taxa_names(ps_object)
  ps_object <- merge_phyloseq(ps_object, dna)
  
  # Rename ASVs to ASV1, ASV2, ...
  taxa_names(ps_object) <- paste0("ASV", seq_len(ntaxa(ps_object)))
  
  # Initial count
  n_initial <- ntaxa(ps_object)
  message("Initial number of ASVs: ", n_initial)
  
  # Filter ASVs: keep only those assigned to Kingdom Fungi
  taxonomy_df <- as.data.frame(tax_table(ps_object))
  asvs_fungi <- taxonomy_df %>%
    filter(!is.na(Kingdom) & Kingdom == "Fungi") %>%
    rownames()
  ps_object <- prune_taxa(asvs_fungi, ps_object)
  
  n_final <- ntaxa(ps_object)
  message("ASVs removed due to non-Fungi Kingdom or NA: ", n_initial - n_final)
  message("Final number of ASVs: ", n_final)
  
  return(ps_object)
}

# Write ASVs table from ps-object
# write_ASVs_table <- function(ps, filename, tax_level = NA) {
#   require(phyloseq)
#   require(dplyr)
#   
#   # If tax_level is specified, aggregate using tax_glom
#   if (!is.na(tax_level)) {
#     if (!(tax_level %in% colnames(tax_table(ps)))) {
#       stop("Specified tax_level not found in tax_table.")
#     }
#     ps <- tax_glom(ps, taxrank = tax_level, NArm = FALSE)
#   }
#   
#   # Extract abundance data
#   otu_df <- as.data.frame(t(otu_table(ps)))
#   
#   # Extract taxonomy
#   tax_df <- as.data.frame(tax_table(ps))
#   
#   # Define taxon names
#   tax_names <- rownames(otu_df)
#   
#   # Add taxon IDs
#   otu_df$Taxon <- rownames(otu_df)
#   tax_df$Taxon <- rownames(tax_df)
#   
#   # Join abundance and taxonomy
#   merged_df <- left_join(otu_df, tax_df, by = "Taxon")
#   
#   # Move Taxon column to front
#   merged_df <- merged_df[, c("Taxon", setdiff(colnames(merged_df), "Taxon"))]
#   
#   # Write to CSV
#   write.csv(merged_df, file = filename, row.names = FALSE)
#   message("Table written to: ", filename)
# }

# Write ASVs table from ps-object
write_ASVs_table <- function(ps, filename, tax_level = NA) {
  require(phyloseq)
  require(dplyr)
  
  # If tax_level is specified, aggregate using tax_glom
  if (!is.na(tax_level)) {
    if (!(tax_level %in% colnames(tax_table(ps)))) {
      stop("Specified tax_level not found in tax_table.")
    }
    ps <- tax_glom(ps, taxrank = tax_level, NArm = FALSE)
  }
  
  # Extract abundance data
  otu_df <- as.data.frame(t(otu_table(ps)))
  
  # Extract taxonomy
  tax_df <- as.data.frame(tax_table(ps))
  
  # Add taxon IDs
  otu_df$Taxon <-  rownames(otu_df)
  tax_df$Taxon <-  rownames(tax_df)
  
  # Join abundance and taxonomy
  merged_df <- left_join(otu_df, tax_df, by = "Taxon")
  
  # Move Taxon column to front
  merged_df <- merged_df[, c("Taxon", setdiff(colnames(merged_df), "Taxon"))]
  
  # Write to CSV
  write.csv(merged_df, file = filename, row.names = FALSE)
  message("Table written to: ", filename)
}

alpha_div_table <- function(ps, cols_to_keep, sample.size = NULL, metric = c("Observed", "Shannon", "Simpson")) {
  require(phyloseq)
  require(picante)
  require(dplyr)
  require(tidyr)
  
  # Helper: root tree by longest branch
  root_by_longest_branch <- function(tree) {
    dist_root <- node.depth.edgelength(tree)
    tip_index <- which.max(dist_root[1:Ntip(tree)])
    tip_label <- tree$tip.label[tip_index]
    tree_rooted <- root(tree, outgroup = tip_label, resolve.root = TRUE)
    return(tree_rooted)
  }
  
  # Prune empty taxa
  ps <- prune_taxa(taxa_sums(ps) > 0, ps)
  
  # Rarefaction
  if (is.null(sample.size)) {
    sample.size <- min(sample_sums(ps))
    message("No rarefaction depth specified. Using minimum depth: ", sample.size)
  } else {
    message("Rarefying to depth: ", sample.size)
  }
  
  ps.r <- rarefy_even_depth(ps, sample.size = sample.size, verbose = FALSE, replace = FALSE)
  
  n_total <- nsamples(ps)
  n_kept <- nsamples(ps.r)
  message("Samples retained after rarefaction: ", n_kept, "/", n_total)
  
  # Calculate basic alpha diversity metrics
  obs_sim <- estimate_richness(ps.r, split = TRUE, measures = metric)
  
  # Initialize alpha table
  alpha <- cbind(as(sample_data(ps.r), "data.frame")[, cols_to_keep, drop = FALSE], obs_sim) %>%
    tibble::rownames_to_column("SampleID")
  
  # Initialize additional metrics list
  extra_metrics <- list()
  
  # Approximate MPD and PD if tree is available
  if (!is.null(ps.r@phy_tree)) {
    otu <- ps.r@otu_table
    if (taxa_are_rows(otu)) otu <- t(otu)
    otu_df <- as.data.frame(otu)
    
    # Calculate cumulative abundance and filter top ASVs (80% of total reads)
    asv_totals <- colSums(otu_df)
    sorted_asvs <- sort(asv_totals, decreasing = TRUE)
    cumulative_sum <- cumsum(sorted_asvs)
    total_sum <- sum(sorted_asvs)
    keep_asvs <- names(sorted_asvs)[cumulative_sum <= 0.8 * total_sum]
    
    if (length(keep_asvs) >= 2) {
      otu_sub <- otu_df[, keep_asvs, drop = FALSE]
      shared_asvs <- intersect(keep_asvs, ps.r@phy_tree$tip.label)
      
      if (length(shared_asvs) >= 2) {
        otu_sub <- otu_sub[, shared_asvs, drop = FALSE]
        tree_sub <- ape::keep.tip(ps.r@phy_tree, shared_asvs)
        
        # Root if needed
        if (is.null(attr(tree_sub, "rooted")) || !isTRUE(attr(tree_sub, "rooted"))) {
          tree_sub <- root_by_longest_branch(tree_sub)
          attr(tree_sub, "rooted") <- TRUE
        }
        
        # MPD
        tree_dist <- cophenetic(tree_sub)
        mpd_vals <- mpd(otu_sub, tree_dist, abundance.weighted = TRUE)
        extra_metrics$MPD <- mpd_vals
      }
    }
    
    # PD (on full tree)
    pd_vals <- pd(as(ps.r@otu_table, "matrix"), ps.r@phy_tree, include.root = TRUE)
    extra_metrics$PD <- pd_vals$PD
  }
  
  # Combine results
  if (length(extra_metrics) > 0) {
    alpha <- cbind(alpha, as.data.frame(extra_metrics))
  }
  
  # Reshape to long format
  all_metrics <- c(metric, names(extra_metrics))
  alpha_long <- alpha %>%
    select(any_of(c("SampleID", cols_to_keep, all_metrics))) %>%
    pivot_longer(cols = all_of(all_metrics), names_to = "Metric", values_to = "Value")
  
  return(alpha_long)
}

#Plot alpha-diversity by all metrics
plot_alpha <- function(alpha_table, X, colour) {
  ggplot(alpha_table, aes_string(x = X, y = "Value", color = colour)) +
    geom_boxplot() +
    geom_point(position = position_dodge(0.75)) +
    facet_wrap(~Metric, scales = "free_y", nrow = 1) +
    theme_light() +
    xlab('') + ylab('') +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
}

plot_heatmap <- function(ps, X, taxa = "Genus", log.transform = TRUE){
  require(dplyr)
  require(phyloseq)
  require(ggplot2)
  
  # Remove taxa with zero total abundance
  ps <- prune_taxa(taxa_sums(ps) > 0, ps)
  
  tax_ranks <- colnames(ps@tax_table)
  
  # Prepare taxonomy as data.frame of characters
  tax_df <- as.data.frame(ps@tax_table, stringsAsFactors = FALSE)
  tax_df[] <- lapply(tax_df, function(x) if (is.factor(x)) as.character(x) else as.character(x))
  
  # Fill chosen rank with fallback "NA // <nearest higher rank>" BEFORE tax_glom
  fill_one <- function(row_vec) {
    rr <- as.list(row_vec)
    val <- rr[[taxa]]
    if (!is.null(val) && !is.na(val) && nzchar(val)) return(val)
    
    idx <- match(taxa, tax_ranks)
    if (!is.na(idx) && idx > 1) {
      for (hr in rev(tax_ranks[1:(idx - 1)])) {
        vh <- rr[[hr]]
        if (!is.null(vh) && !is.na(vh) && nzchar(vh)) {
          return(paste("NA //", vh))
        }
      }
    }
    return("NA")
  }
  tax_df[[taxa]] <- apply(tax_df, 1, fill_one)
  
  # Replace taxonomy in phyloseq with modified version
  tax_table(ps) <- as.matrix(tax_df)
  
  # Collapse to the selected taxonomic rank (no ASV lost now)
  ps <- tax_glom(ps, taxa)
  
  # Melt phyloseq object
  sig.taxa.long <- psmelt(ps) %>%
    arrange(Phylum) %>% 
    mutate(row = row_number())
  
  sig.taxa.long$Abundance <- as.numeric(sig.taxa.long$Abundance)
  sig.taxa.long$Taxa <- sig.taxa.long[, taxa]
  
  # Plot heatmap
  ggplot(sig.taxa.long, aes(x = !!sym(X), y = reorder(Taxa, row))) +
    {if(log.transform) geom_tile(aes(fill = log(Abundance)))} +
    {if(!log.transform) geom_tile(aes(fill = Abundance))} +
    scale_fill_distiller(palette = "OrRd", trans = "reverse") +
    facet_grid(rows = vars(Phylum), scales = "free", space = "free") +
    theme_light() +
    theme(
      strip.text.y = element_text(angle = 0),
      panel.spacing = unit(0.02,'lines'),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
    ) +
    xlab("") + ylab("")
}

