---
title: "Base R template for 16S amplicon processing"
output: 
  html_document: 
    keep_md: yes
    toc: yes
    toc_float: yes
    toc_collapsed: yes
---




``` r
# Load required libraries
library(dada2)
library(Biostrings)
library(DECIPHER)
library(phyloseq)
library(tidyverse)
library(telegramNotify) # custom notification package

# Load helper functions
source('functions.R')
set.seed(5678)
```

# Read metadata

Read metadata, assign sample names as rownames (important for downstream analysis).


``` r
mdat <- read_metadata(filename = "16s_map.csv",
                      sample.names = "Filename",
                      sep = ",",
                      check_files = "16s_raw/")
mdat %>% head()
```

```
##      Filename Type      Group Treatment  Dosage
## S246     S246  16S Experiment   Control control
## S247     S247  16S Experiment   Control control
## S248     S248  16S Experiment   Control control
## S249     S249  16S Experiment   Control control
## S250     S250  16S Experiment   Control control
## S251     S251  16S Experiment   Control control
```

# ASV picking

Run DADA2 pipeline for 16S sequences.


``` r
notified_reads_to_seqtable_16s <- run_with_notification(fun = reads_to_seqtable_16s, 
                                          task_name = "16S ASV picking",
                                          notify_success = TRUE)

seqtab <- notified_reads_to_seqtable_16s(raw_files_path = "16s_raw",
                                trimLeft = c(20, 19),
                                truncLen = c(240, 230), 
                                log_file = "16S_processing.log",
                                pool = "pseudo",
                                cores = TRUE,
                                maxEE = c(2,5),
                                pattern = c("_L2_1.fq.gz", "_L2_2.fq.gz"),
                                plot_profile = TRUE)
```

```
##  [1] "S246" "S247" "S248" "S249" "S250" "S251" "S255" "S256" "S257" "S258"
## [11] "S259" "S260" "S264" "S265" "S266" "S267" "S268" "S269" "S273" "S274"
## [21] "S275" "S276" "S277" "S278" "S282" "S283" "S284" "S285" "S286" "S287"
## [31] "S291" "S292" "S293" "S294" "S295" "S296" "S300" "S301" "S302" "S303"
## [41] "S304" "S305"
## 844140 total bases in 3837 reads from 42 samples will be used for learning the error rates.
## 809607 total bases in 3837 reads from 42 samples will be used for learning the error rates.
## Sample 1 - 96 reads in 75 unique sequences.
## Sample 2 - 89 reads in 65 unique sequences.
## Sample 3 - 94 reads in 74 unique sequences.
## Sample 4 - 90 reads in 68 unique sequences.
## Sample 5 - 85 reads in 68 unique sequences.
## Sample 6 - 90 reads in 75 unique sequences.
## Sample 7 - 90 reads in 71 unique sequences.
## Sample 8 - 95 reads in 78 unique sequences.
## Sample 9 - 88 reads in 76 unique sequences.
## Sample 10 - 96 reads in 75 unique sequences.
## Sample 11 - 91 reads in 74 unique sequences.
## Sample 12 - 95 reads in 77 unique sequences.
## Sample 13 - 92 reads in 75 unique sequences.
## Sample 14 - 85 reads in 62 unique sequences.
## Sample 15 - 95 reads in 71 unique sequences.
## Sample 16 - 91 reads in 61 unique sequences.
## Sample 17 - 94 reads in 71 unique sequences.
## Sample 18 - 93 reads in 61 unique sequences.
## Sample 19 - 92 reads in 81 unique sequences.
## Sample 20 - 86 reads in 71 unique sequences.
## Sample 21 - 96 reads in 79 unique sequences.
## Sample 22 - 89 reads in 73 unique sequences.
## Sample 23 - 91 reads in 72 unique sequences.
## Sample 24 - 94 reads in 74 unique sequences.
## Sample 25 - 89 reads in 61 unique sequences.
## Sample 26 - 86 reads in 56 unique sequences.
## Sample 27 - 86 reads in 58 unique sequences.
## Sample 28 - 90 reads in 64 unique sequences.
## Sample 29 - 91 reads in 71 unique sequences.
## Sample 30 - 94 reads in 72 unique sequences.
## Sample 31 - 95 reads in 75 unique sequences.
## Sample 32 - 88 reads in 65 unique sequences.
## Sample 33 - 96 reads in 79 unique sequences.
## Sample 34 - 94 reads in 70 unique sequences.
## Sample 35 - 94 reads in 74 unique sequences.
## Sample 36 - 91 reads in 74 unique sequences.
## Sample 37 - 95 reads in 78 unique sequences.
## Sample 38 - 92 reads in 67 unique sequences.
## Sample 39 - 92 reads in 74 unique sequences.
## Sample 40 - 89 reads in 61 unique sequences.
## Sample 41 - 93 reads in 74 unique sequences.
## Sample 42 - 85 reads in 69 unique sequences.
## 
##    selfConsist step 2Sample 1 - 96 reads in 76 unique sequences.
## Sample 2 - 89 reads in 66 unique sequences.
## Sample 3 - 94 reads in 69 unique sequences.
## Sample 4 - 90 reads in 62 unique sequences.
## Sample 5 - 85 reads in 63 unique sequences.
## Sample 6 - 90 reads in 74 unique sequences.
## Sample 7 - 90 reads in 70 unique sequences.
## Sample 8 - 95 reads in 77 unique sequences.
## Sample 9 - 88 reads in 73 unique sequences.
## Sample 10 - 96 reads in 73 unique sequences.
## Sample 11 - 91 reads in 72 unique sequences.
## Sample 12 - 95 reads in 76 unique sequences.
## Sample 13 - 92 reads in 74 unique sequences.
## Sample 14 - 85 reads in 63 unique sequences.
## Sample 15 - 95 reads in 69 unique sequences.
## Sample 16 - 91 reads in 62 unique sequences.
## Sample 17 - 94 reads in 69 unique sequences.
## Sample 18 - 93 reads in 62 unique sequences.
## Sample 19 - 92 reads in 77 unique sequences.
## Sample 20 - 86 reads in 68 unique sequences.
## Sample 21 - 96 reads in 76 unique sequences.
## Sample 22 - 89 reads in 75 unique sequences.
## Sample 23 - 91 reads in 72 unique sequences.
## Sample 24 - 94 reads in 74 unique sequences.
## Sample 25 - 89 reads in 64 unique sequences.
## Sample 26 - 86 reads in 51 unique sequences.
## Sample 27 - 86 reads in 61 unique sequences.
## Sample 28 - 90 reads in 67 unique sequences.
## Sample 29 - 91 reads in 70 unique sequences.
## Sample 30 - 94 reads in 81 unique sequences.
## Sample 31 - 95 reads in 74 unique sequences.
## Sample 32 - 88 reads in 63 unique sequences.
## Sample 33 - 96 reads in 82 unique sequences.
## Sample 34 - 94 reads in 67 unique sequences.
## Sample 35 - 94 reads in 77 unique sequences.
## Sample 36 - 91 reads in 73 unique sequences.
## Sample 37 - 95 reads in 80 unique sequences.
## Sample 38 - 92 reads in 61 unique sequences.
## Sample 39 - 92 reads in 75 unique sequences.
## Sample 40 - 89 reads in 62 unique sequences.
## Sample 41 - 93 reads in 74 unique sequences.
## Sample 42 - 85 reads in 59 unique sequences.
## 
##    selfConsist step 2
```

```
##      input filtered denoisedF denoisedR merged nonchim
## S246   100       96        49        43     37      37
## S247   100       89        47        48     38      38
## S248   100       94        37        47     32      32
## S249   100       90        52        63     43      43
## S250   100       85        41        49     38      38
## S251   100       90        29        31     26      26
## S255   100       90        46        48     42      42
## S256   100       95        43        41     26      26
## S257   100       88        38        40     20      20
## S258   100       96        50        54     36      36
## S259   100       91        37        44     32      32
## S260   100       95        40        37     30      30
## S264   100       92        44        44     28      28
## S265   100       85        46        48     39      39
## S266   100       95        45        45     37      37
## S267   100       91        58        59     36      36
## S268   100       94        48        51     35      35
## S269   100       93        68        67     57      57
## S273   100       92        34        44     29      29
## S274   100       86        37        36     24      24
## S275   100       96        42        45     34      34
## S276   100       89        40        38     28      28
## S277   100       91        50        51     36      36
## S278   100       94        54        54     41      41
## S282   100       89        61        60     50      50
## S283   100       86        56        63     43      43
## S284   100       86        55        53     41      41
## S285   100       90        45        39     33      33
## S286   100       91        44        49     39      39
## S287   100       94        50        31     29      29
## S291   100       95        41        40     26      26
## S292   100       88        59        60     52      52
## S293   100       96        39        34     31      31
## S294   100       94        47        51     37      37
## S295   100       94        40        36     30      30
## S296   100       91        43        47     34      34
## S300   100       95        40        34     31      31
## S301   100       92        51        57     44      44
## S302   100       92        49        42     37      37
## S303   100       89        66        57     55      55
## S304   100       93        44        44     34      34
## S305   100       85        38        45     31      31
```

``` r
send_telegram_file("16S_processing.log")
```

```
## Response [https://api.telegram.org/bot8070906504:AAHLUS70wr6hRCs00ARZUfiNaHiwNxpq0Sc/sendDocument]
##   Date: 2025-09-25 13:05
##   Status: 200
##   Content-Type: application/json
##   Size: 483 B
```

# Taxonomic assignment

Assign taxonomy using SILVA database.


``` r
notified_assign_taxonomy <- run_with_notification(fun = assign_taxonomy, 
                                          task_name = "16S Taxonomy",
                                          notify_success = TRUE)

taxa <- notified_assign_taxonomy(seqtab = seqtab,
                        set_train_path = "/home/alexey/tax_n_refs/silva_nr_v138_train_set.fa.gz",
                        cores = TRUE)
taxa[1:10]
```

```
##  [1] "Bacteria" "Bacteria" "Bacteria" "Bacteria" "Bacteria" "Archaea" 
##  [7] "Bacteria" "Bacteria" "Bacteria" "Bacteria"
```

# Phyloseq object assembly

Combine ASV table, taxonomy, and metadata into a `phyloseq` object. Optionally, filter out organelle sequences and export reference FASTA.


``` r
ps <- assemble_phyloseq_16s(seqtab = seqtab,
                            metadata = mdat, 
                            taxonomy = taxa,
                            filter.organelles = TRUE, 
                            write_fasta = 'refseqs.fasta')
```

```
## Initial number of ASVs: 392
```

```
## ASVs removed due to NA in Phylum: 46
```

```
## ASVs removed due to non-Bacteria/Archaea Kingdom: 0
```

```
## ASVs removed as organelles (Mitochondria/Chloroplast): 3
```

```
## Final number of ASVs: 343
```

```
## FASTA file written to: refseqs.fasta
```

``` r
ps
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 343 taxa and 42 samples ]
## sample_data() Sample Data:       [ 42 samples by 5 sample variables ]
## tax_table()   Taxonomy Table:    [ 343 taxa by 6 taxonomic ranks ]
## refseq()      DNAStringSet:      [ 343 reference sequences ]
```

# Phylogenetic tree construction

Build tree from representative sequences and merge into `phyloseq` object.


``` r
system2(c('./build_tree.sh', 'refseqs.fasta'))
tree <- ape::read.tree('tree.nwk')
ps <- merge_phyloseq(ps, tree)

saveRDS(ps, "bact.RData")
```

# Basic stats

Explore the dataset: number of taxa, reads per sample, taxonomy and ASV tables. Save `phyloseq` object to `.RData` file.



``` r
sample_names(ps) # Names of samples
```

```
##  [1] "S246" "S247" "S248" "S249" "S250" "S251" "S255" "S256" "S257" "S258"
## [11] "S259" "S260" "S264" "S265" "S266" "S267" "S268" "S269" "S273" "S274"
## [21] "S275" "S276" "S277" "S278" "S282" "S283" "S284" "S285" "S286" "S287"
## [31] "S291" "S292" "S293" "S294" "S295" "S296" "S300" "S301" "S302" "S303"
## [41] "S304" "S305"
```

``` r
sample_sums(ps)  # Number of reads per sample
```

```
## S246 S247 S248 S249 S250 S251 S255 S256 S257 S258 S259 S260 S264 S265 S266 S267 
##   27   28   26   39   33   24   42   20   16   33   29   26   26   39   28   36 
## S268 S269 S273 S274 S275 S276 S277 S278 S282 S283 S284 S285 S286 S287 S291 S292 
##   31   51   29   24   30   28   36   41   50   40   33   33   36   29   14   52 
## S293 S294 S295 S296 S300 S301 S302 S303 S304 S305 
##   29   35   27   31   26   42   37   52   34   31
```

``` r
tax_table(ps)[1:5, 1:4] # First rows of taxonomy table
```

```
## Taxonomy Table:     [5 taxa by 4 taxonomic ranks]:
##        Kingdom   Phylum          Class             Order              
## ASV98  "Archaea" "Crenarchaeota" "Nitrososphaeria" "Nitrososphaerales"
## ASV17  "Archaea" "Crenarchaeota" "Nitrososphaeria" "Nitrososphaerales"
## ASV6   "Archaea" "Crenarchaeota" "Nitrososphaeria" "Nitrososphaerales"
## ASV60  "Archaea" "Crenarchaeota" "Nitrososphaeria" "Nitrososphaerales"
## ASV171 "Archaea" "Crenarchaeota" "Nitrososphaeria" "Nitrososphaerales"
```

``` r
otu_table(ps)[1:4, 1:5] # First rows of ASV table
```

```
## OTU Table:          [5 taxa and 4 samples]
##                      taxa are columns
##      ASV98 ASV17 ASV6 ASV60 ASV171
## S246     0     0    0     0      0
## S247     0     0    0     0      0
## S248     0     0    0     0      0
## S249     0    10    0     0      0
```

``` r
ps@sam_data # Metadata
```

```
##      Filename Type      Group          Treatment  Dosage
## S246     S246  16S Experiment            Control control
## S247     S247  16S Experiment            Control control
## S248     S248  16S Experiment            Control control
## S249     S249  16S Experiment            Control control
## S250     S250  16S Experiment            Control control
## S251     S251  16S Experiment            Control control
## S255     S255  16S Experiment     Poultry manure  single
## S256     S256  16S Experiment     Poultry manure  single
## S257     S257  16S Experiment     Poultry manure  single
## S258     S258  16S Experiment     Poultry manure  single
## S259     S259  16S Experiment     Poultry manure  single
## S260     S260  16S Experiment     Poultry manure  single
## S264     S264  16S Experiment     Poultry manure  double
## S265     S265  16S Experiment     Poultry manure  double
## S266     S266  16S Experiment     Poultry manure  double
## S267     S267  16S Experiment     Poultry manure  double
## S268     S268  16S Experiment     Poultry manure  double
## S269     S269  16S Experiment     Poultry manure  double
## S273     S273  16S Experiment Rapid fermentation  single
## S274     S274  16S Experiment Rapid fermentation  single
## S275     S275  16S Experiment Rapid fermentation  single
## S276     S276  16S Experiment Rapid fermentation  single
## S277     S277  16S Experiment Rapid fermentation  single
## S278     S278  16S Experiment Rapid fermentation  single
## S282     S282  16S Experiment Rapid fermentation  double
## S283     S283  16S Experiment Rapid fermentation  double
## S284     S284  16S Experiment Rapid fermentation  double
## S285     S285  16S Experiment Rapid fermentation  double
## S286     S286  16S Experiment Rapid fermentation  double
## S287     S287  16S Experiment Rapid fermentation  double
## S291     S291  16S Experiment  Long fermentation  single
## S292     S292  16S Experiment  Long fermentation  single
## S293     S293  16S Experiment  Long fermentation  single
## S294     S294  16S Experiment  Long fermentation  single
## S295     S295  16S Experiment  Long fermentation  single
## S296     S296  16S Experiment  Long fermentation  single
## S300     S300  16S Experiment  Long fermentation  double
## S301     S301  16S Experiment  Long fermentation  double
## S302     S302  16S Experiment  Long fermentation  double
## S303     S303  16S Experiment  Long fermentation  double
## S304     S304  16S Experiment  Long fermentation  double
## S305     S305  16S Experiment  Long fermentation  double
```
