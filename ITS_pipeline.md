---
title: "Base R template for ITS amplicon processing"
output: 
  html_document: 
    keep_md: yes
    toc: yes
    toc_float: yes
    toc_collapsed: yes
---





# Read metadata

`read_metadata(filename, sample.names, ...)`

Read metadata, assign sample names as rownames (important for downstream analysis).


```
## Warning in read_metadata(filename = "its_map.csv", sample.names = "Filename", :
## The following samples do not have 2 corresponding files in the directory: PE6
```

```
##      Filename Type      Group          Treatment  Dosage
## PE6       PE6  ITS Experiment            Control control
## PE7       PE7  ITS Experiment            Control control
## PE8       PE8  ITS Experiment            Control control
## PE9       PE9  ITS Experiment            Control control
## PE10     PE10  ITS Experiment            Control control
## PE11     PE11  ITS Experiment            Control control
## PE15     PE15  ITS Experiment     Poultry manure  single
## PE16     PE16  ITS Experiment     Poultry manure  single
## PE17     PE17  ITS Experiment     Poultry manure  single
## PE18     PE18  ITS Experiment     Poultry manure  single
## PE19     PE19  ITS Experiment     Poultry manure  single
## PE20     PE20  ITS Experiment     Poultry manure  single
## PE24     PE24  ITS Experiment     Poultry manure  double
## PE25     PE25  ITS Experiment     Poultry manure  double
## PE26     PE26  ITS Experiment     Poultry manure  double
## PE27     PE27  ITS Experiment     Poultry manure  double
## PE28     PE28  ITS Experiment     Poultry manure  double
## PE29     PE29  ITS Experiment     Poultry manure  double
## PE33     PE33  ITS Experiment Rapid fermentation  single
## PE34     PE34  ITS Experiment Rapid fermentation  single
## PE35     PE35  ITS Experiment Rapid fermentation  single
## PE36     PE36  ITS Experiment Rapid fermentation  single
## PE37     PE37  ITS Experiment Rapid fermentation  single
## PE38     PE38  ITS Experiment Rapid fermentation  single
## PE42     PE42  ITS Experiment Rapid fermentation  double
## PE43     PE43  ITS Experiment Rapid fermentation  double
## PE44     PE44  ITS Experiment Rapid fermentation  double
## PE45     PE45  ITS Experiment Rapid fermentation  double
## PE46     PE46  ITS Experiment Rapid fermentation  double
## PE47     PE47  ITS Experiment Rapid fermentation  double
## PE51     PE51  ITS Experiment  Long fermentation  single
## PE52     PE52  ITS Experiment  Long fermentation  single
## PE53     PE53  ITS Experiment  Long fermentation  single
## PE54     PE54  ITS Experiment  Long fermentation  single
## PE55     PE55  ITS Experiment  Long fermentation  single
## PE56     PE56  ITS Experiment  Long fermentation  single
## PE60     PE60  ITS Experiment  Long fermentation  double
## PE61     PE61  ITS Experiment  Long fermentation  double
## PE62     PE62  ITS Experiment  Long fermentation  double
## PE63     PE63  ITS Experiment  Long fermentation  double
## PE64     PE64  ITS Experiment  Long fermentation  double
## PE65     PE65  ITS Experiment  Long fermentation  double
```

# Primer trimming

Use `cutadapt` to remove ITS primers.


```
## Loading required package: ShortRead
```

```
## Loading required package: BiocParallel
```

```
## Loading required package: Rsamtools
```

```
## Loading required package: GenomicAlignments
```

```
## 
## Attaching package: 'GenomicAlignments'
```

```
## The following object is masked from 'package:dplyr':
## 
##     last
```

```
## 
## Attaching package: 'ShortRead'
```

```
## The following object is masked from 'package:dplyr':
## 
##     id
```

```
## The following object is masked from 'package:purrr':
## 
##     compose
```

```
## The following object is masked from 'package:tibble':
## 
##     view
```

```
## >> Filtering all reads with Ns...
```

```
## >> Primer hits BEFORE trimming (1st sample):
```

```
##                  Forward Complement Reverse RevComp
## FWD.ForwardReads      93          0       0       0
## FWD.ReverseReads       0          0       0       0
## REV.ForwardReads       0          0       0       0
## REV.ReverseReads      92          0       0       0
```

```
## >> Running cutadapt on all samples...
```

```
## >> Primer hits AFTER trimming (1st sample):
```

```
##                  Forward Complement Reverse RevComp
## FWD.ForwardReads       0          0       0       0
## FWD.ReverseReads       0          0       0       0
## REV.ForwardReads       0          0       0       0
## REV.ReverseReads       0          0       0       0
```

```
## >> Done. Trimmed files are in: its_raw//cutadapt
```

```
## [1] "its_raw//cutadapt"
```

# ASV picking

Run DADA2 pipeline for ITS sequences.

![](ITS_pipeline_files/figure-html/unnamed-chunk-4-1.png)<!-- -->![](ITS_pipeline_files/figure-html/unnamed-chunk-4-2.png)<!-- -->

```
## 791061 total bases in 3424 reads from 42 samples will be used for learning the error rates.
## 788080 total bases in 3424 reads from 42 samples will be used for learning the error rates.
## Sample 1 - 73 reads in 67 unique sequences.
## Sample 2 - 89 reads in 55 unique sequences.
## Sample 3 - 78 reads in 66 unique sequences.
## Sample 4 - 77 reads in 67 unique sequences.
## Sample 5 - 86 reads in 61 unique sequences.
## Sample 6 - 84 reads in 54 unique sequences.
## Sample 7 - 85 reads in 36 unique sequences.
## Sample 8 - 83 reads in 65 unique sequences.
## Sample 9 - 83 reads in 54 unique sequences.
## Sample 10 - 79 reads in 60 unique sequences.
## Sample 11 - 83 reads in 67 unique sequences.
## Sample 12 - 75 reads in 55 unique sequences.
## Sample 13 - 81 reads in 57 unique sequences.
## Sample 14 - 85 reads in 57 unique sequences.
## Sample 15 - 47 reads in 47 unique sequences.
## Sample 16 - 87 reads in 48 unique sequences.
## Sample 17 - 86 reads in 55 unique sequences.
## Sample 18 - 79 reads in 60 unique sequences.
## Sample 19 - 88 reads in 53 unique sequences.
## Sample 20 - 86 reads in 45 unique sequences.
## Sample 21 - 87 reads in 43 unique sequences.
## Sample 22 - 78 reads in 38 unique sequences.
## Sample 23 - 84 reads in 55 unique sequences.
## Sample 24 - 90 reads in 52 unique sequences.
## Sample 25 - 89 reads in 39 unique sequences.
## Sample 26 - 82 reads in 57 unique sequences.
## Sample 27 - 87 reads in 52 unique sequences.
## Sample 28 - 84 reads in 56 unique sequences.
## Sample 29 - 87 reads in 39 unique sequences.
## Sample 30 - 86 reads in 34 unique sequences.
## Sample 31 - 85 reads in 55 unique sequences.
## Sample 32 - 81 reads in 54 unique sequences.
## Sample 33 - 81 reads in 63 unique sequences.
## Sample 34 - 83 reads in 59 unique sequences.
## Sample 35 - 81 reads in 45 unique sequences.
## Sample 36 - 83 reads in 49 unique sequences.
## Sample 37 - 88 reads in 47 unique sequences.
## Sample 38 - 80 reads in 48 unique sequences.
## Sample 39 - 52 reads in 48 unique sequences.
## Sample 40 - 86 reads in 55 unique sequences.
## Sample 41 - 83 reads in 68 unique sequences.
## Sample 42 - 73 reads in 62 unique sequences.
## 
##    selfConsist step 2Sample 1 - 73 reads in 58 unique sequences.
## Sample 2 - 89 reads in 60 unique sequences.
## Sample 3 - 78 reads in 65 unique sequences.
## Sample 4 - 77 reads in 63 unique sequences.
## Sample 5 - 86 reads in 59 unique sequences.
## Sample 6 - 84 reads in 55 unique sequences.
## Sample 7 - 85 reads in 40 unique sequences.
## Sample 8 - 83 reads in 64 unique sequences.
## Sample 9 - 83 reads in 51 unique sequences.
## Sample 10 - 79 reads in 49 unique sequences.
## Sample 11 - 83 reads in 51 unique sequences.
## Sample 12 - 75 reads in 33 unique sequences.
## Sample 13 - 81 reads in 58 unique sequences.
## Sample 14 - 85 reads in 57 unique sequences.
## Sample 15 - 47 reads in 35 unique sequences.
## Sample 16 - 87 reads in 47 unique sequences.
## Sample 17 - 86 reads in 50 unique sequences.
## Sample 18 - 79 reads in 48 unique sequences.
## Sample 19 - 88 reads in 58 unique sequences.
## Sample 20 - 86 reads in 44 unique sequences.
## Sample 21 - 87 reads in 45 unique sequences.
## Sample 22 - 78 reads in 38 unique sequences.
## Sample 23 - 84 reads in 55 unique sequences.
## Sample 24 - 90 reads in 44 unique sequences.
## Sample 25 - 89 reads in 35 unique sequences.
## Sample 26 - 82 reads in 48 unique sequences.
## Sample 27 - 87 reads in 48 unique sequences.
## Sample 28 - 84 reads in 55 unique sequences.
## Sample 29 - 87 reads in 32 unique sequences.
## Sample 30 - 86 reads in 28 unique sequences.
## Sample 31 - 85 reads in 54 unique sequences.
## Sample 32 - 81 reads in 53 unique sequences.
## Sample 33 - 81 reads in 58 unique sequences.
## Sample 34 - 83 reads in 52 unique sequences.
## Sample 35 - 81 reads in 35 unique sequences.
## Sample 36 - 83 reads in 38 unique sequences.
## Sample 37 - 88 reads in 47 unique sequences.
## Sample 38 - 80 reads in 44 unique sequences.
## Sample 39 - 52 reads in 27 unique sequences.
## Sample 40 - 86 reads in 44 unique sequences.
## Sample 41 - 83 reads in 66 unique sequences.
## Sample 42 - 73 reads in 51 unique sequences.
## 
##    selfConsist step 2
```

```
## 24 paired-reads (in 9 unique pairings) successfully merged out of 24 (in 9 pairings) input.
```

```
## 68 paired-reads (in 11 unique pairings) successfully merged out of 68 (in 11 pairings) input.
```

```
## 33 paired-reads (in 14 unique pairings) successfully merged out of 33 (in 14 pairings) input.
```

```
## 26 paired-reads (in 9 unique pairings) successfully merged out of 26 (in 9 pairings) input.
```

```
## 51 paired-reads (in 13 unique pairings) successfully merged out of 51 (in 13 pairings) input.
```

```
## 56 paired-reads (in 7 unique pairings) successfully merged out of 56 (in 7 pairings) input.
```

```
## 73 paired-reads (in 10 unique pairings) successfully merged out of 73 (in 10 pairings) input.
```

```
## 50 paired-reads (in 11 unique pairings) successfully merged out of 50 (in 11 pairings) input.
```

```
## 77 paired-reads (in 8 unique pairings) successfully merged out of 77 (in 8 pairings) input.
```

```
## 60 paired-reads (in 11 unique pairings) successfully merged out of 60 (in 11 pairings) input.
```

```
## 57 paired-reads (in 8 unique pairings) successfully merged out of 57 (in 8 pairings) input.
```

```
## 69 paired-reads (in 7 unique pairings) successfully merged out of 69 (in 7 pairings) input.
```

```
## 65 paired-reads (in 16 unique pairings) successfully merged out of 65 (in 16 pairings) input.
```

```
## 58 paired-reads (in 14 unique pairings) successfully merged out of 58 (in 14 pairings) input.
```

```
## 24 paired-reads (in 7 unique pairings) successfully merged out of 24 (in 7 pairings) input.
```

```
## 72 paired-reads (in 10 unique pairings) successfully merged out of 72 (in 10 pairings) input.
```

```
## 68 paired-reads (in 9 unique pairings) successfully merged out of 68 (in 9 pairings) input.
```

```
## 60 paired-reads (in 15 unique pairings) successfully merged out of 60 (in 15 pairings) input.
```

```
## 69 paired-reads (in 10 unique pairings) successfully merged out of 69 (in 10 pairings) input.
```

```
## 73 paired-reads (in 9 unique pairings) successfully merged out of 73 (in 9 pairings) input.
```

```
## 73 paired-reads (in 10 unique pairings) successfully merged out of 73 (in 10 pairings) input.
```

```
## 67 paired-reads (in 11 unique pairings) successfully merged out of 67 (in 11 pairings) input.
```

```
## 67 paired-reads (in 10 unique pairings) successfully merged out of 67 (in 10 pairings) input.
```

```
## 76 paired-reads (in 11 unique pairings) successfully merged out of 76 (in 11 pairings) input.
```

```
## 87 paired-reads (in 5 unique pairings) successfully merged out of 87 (in 5 pairings) input.
```

```
## 59 paired-reads (in 11 unique pairings) successfully merged out of 59 (in 11 pairings) input.
```

```
## 67 paired-reads (in 12 unique pairings) successfully merged out of 67 (in 12 pairings) input.
```

```
## 61 paired-reads (in 14 unique pairings) successfully merged out of 61 (in 14 pairings) input.
```

```
## 80 paired-reads (in 10 unique pairings) successfully merged out of 80 (in 10 pairings) input.
```

```
## 81 paired-reads (in 7 unique pairings) successfully merged out of 81 (in 7 pairings) input.
```

```
## 69 paired-reads (in 15 unique pairings) successfully merged out of 69 (in 15 pairings) input.
```

```
## 58 paired-reads (in 12 unique pairings) successfully merged out of 58 (in 12 pairings) input.
```

```
## 45 paired-reads (in 16 unique pairings) successfully merged out of 45 (in 16 pairings) input.
```

```
## 61 paired-reads (in 13 unique pairings) successfully merged out of 61 (in 13 pairings) input.
```

```
## 75 paired-reads (in 9 unique pairings) successfully merged out of 75 (in 9 pairings) input.
```

```
## 72 paired-reads (in 9 unique pairings) successfully merged out of 72 (in 9 pairings) input.
```

```
## 79 paired-reads (in 16 unique pairings) successfully merged out of 79 (in 16 pairings) input.
```

```
## 75 paired-reads (in 14 unique pairings) successfully merged out of 75 (in 14 pairings) input.
```

```
## 32 paired-reads (in 8 unique pairings) successfully merged out of 32 (in 8 pairings) input.
```

```
## 68 paired-reads (in 4 unique pairings) successfully merged out of 68 (in 4 pairings) input.
```

```
## 43 paired-reads (in 12 unique pairings) successfully merged out of 43 (in 12 pairings) input.
```

```
## 43 paired-reads (in 8 unique pairings) successfully merged out of 43 (in 8 pairings) input.
```

```
##      input filtered denoisedF denoisedR merged nonchim
## PE10    99       73        27        44     24      24
## PE11   100       89        70        70     68      68
## PE15   100       78        43        43     33      33
## PE16   100       77        35        40     26      26
## PE17   100       86        53        58     51      51
## PE18   100       84        59        64     56      56
## PE19   100       85        76        73     73      73
## PE20    98       83        59        57     50      50
## PE24   100       83        79        77     77      77
## PE25    99       79        63        67     60      60
## PE26    98       83        60        69     57      57
## PE27   100       75        69        72     69      69
## PE28    98       81        67        67     65      65
## PE29   100       85        61        60     58      58
## PE33   100       47        24        41     24      24
## PE34   100       87        79        74     72      72
## PE35    99       86        75        72     68      68
## PE36   100       79        62        64     60      60
## PE37   100       88        73        72     69      69
## PE38   100       86        75        77     73      73
## PE42   100       87        77        73     73      73
## PE43    99       78        71        68     67      67
## PE44   100       84        70        69     67      67
## PE45   100       90        78        83     76      76
## PE46   100       89        87        87     87      87
## PE47   100       82        59        62     59      59
## PE51   100       87        71        67     67      67
## PE52   100       84        65        68     61      61
## PE53   100       87        81        80     80      80
## PE54   100       86        81        81     81      81
## PE55   100       85        73        70     69      69
## PE56   100       81        63        60     58      58
## PE6     99       81        49        47     45      45
## PE60   100       83        67        63     61      61
## PE61   100       81        75        76     75      75
## PE62   100       83        72        74     72      72
## PE63   100       88        83        80     79      79
## PE64   100       80        77        75     75      75
## PE65   100       52        32        44     32      32
## PE7    100       86        73        72     68      68
## PE8    100       83        48        52     43      43
## PE9     98       73        46        49     43      43
```

# Taxonomic assignment

Assign taxonomy using UNITE database.


```
## UNITE fungal taxonomic reference detected.
```

# Phyloseq object assembly

Combine ASV table, taxonomy, and metadata into a `phyloseq` object.


```
## UNITE-style prefixes detected. Removing prefixes (e.g., 'k__', 'p__', etc.)...
```

```
## Initial number of ASVs: 167
```

```
## ASVs removed due to non-Fungi Kingdom or NA: 27
```

```
## Final number of ASVs: 140
```

# Basic stats

Explore the dataset: number of taxa, reads per sample, taxonomy and ASV tables. Save `phyloseq` object to `.RData` file.


```
##  [1] "PE10" "PE11" "PE15" "PE16" "PE17" "PE18" "PE19" "PE20" "PE24" "PE25"
## [11] "PE26" "PE27" "PE28" "PE29" "PE33" "PE34" "PE35" "PE36" "PE37" "PE38"
## [21] "PE42" "PE43" "PE44" "PE45" "PE46" "PE47" "PE51" "PE52" "PE53" "PE54"
## [31] "PE55" "PE56" "PE6"  "PE60" "PE61" "PE62" "PE63" "PE64" "PE65" "PE7" 
## [41] "PE8"  "PE9"
```

```
## PE10 PE11 PE15 PE16 PE17 PE18 PE19 PE20 PE24 PE25 PE26 PE27 PE28 PE29 PE33 PE34 
##   15   61   20   21   37   56   73   39   76   57   23   63   58   56   21   70 
## PE35 PE36 PE37 PE38 PE42 PE43 PE44 PE45 PE46 PE47 PE51 PE52 PE53 PE54 PE55 PE56 
##   64   59   69   70   71   66   67   76   87   55   21   61   80   10   67   56 
##  PE6 PE60 PE61 PE62 PE63 PE64 PE65  PE7  PE8  PE9 
##   38   59   72   67   77   73   32   66   29   42
```

```
## Taxonomy Table:     [5 taxa by 4 taxonomic ranks]:
##      Kingdom Phylum          Class             Order             
## ASV1 "Fungi" "Ascomycota"    "Sordariomycetes" "Sordariales"     
## ASV2 "Fungi" "Ascomycota"    "Pezizomycetes"   "Pezizales"       
## ASV3 "Fungi" "Basidiomycota" "Tremellomycetes" "Trichosporonales"
## ASV4 "Fungi" "Basidiomycota" "Tremellomycetes" "Trichosporonales"
## ASV5 "Fungi" "Ascomycota"    "Sordariomycetes" "Sordariales"
```

```
## OTU Table:          [5 taxa and 4 samples]
##                      taxa are columns
##      ASV1 ASV2 ASV3 ASV4 ASV5
## PE10    0    0    0    0    0
## PE11    0   43    0    0    3
## PE15    3    0    0    0    0
## PE16    0    0    0    1    0
```

```
##      Filename Type      Group          Treatment  Dosage
## PE10     PE10  ITS Experiment            Control control
## PE11     PE11  ITS Experiment            Control control
## PE15     PE15  ITS Experiment     Poultry manure  single
## PE16     PE16  ITS Experiment     Poultry manure  single
## PE17     PE17  ITS Experiment     Poultry manure  single
## PE18     PE18  ITS Experiment     Poultry manure  single
## PE19     PE19  ITS Experiment     Poultry manure  single
## PE20     PE20  ITS Experiment     Poultry manure  single
## PE24     PE24  ITS Experiment     Poultry manure  double
## PE25     PE25  ITS Experiment     Poultry manure  double
## PE26     PE26  ITS Experiment     Poultry manure  double
## PE27     PE27  ITS Experiment     Poultry manure  double
## PE28     PE28  ITS Experiment     Poultry manure  double
## PE29     PE29  ITS Experiment     Poultry manure  double
## PE33     PE33  ITS Experiment Rapid fermentation  single
## PE34     PE34  ITS Experiment Rapid fermentation  single
## PE35     PE35  ITS Experiment Rapid fermentation  single
## PE36     PE36  ITS Experiment Rapid fermentation  single
## PE37     PE37  ITS Experiment Rapid fermentation  single
## PE38     PE38  ITS Experiment Rapid fermentation  single
## PE42     PE42  ITS Experiment Rapid fermentation  double
## PE43     PE43  ITS Experiment Rapid fermentation  double
## PE44     PE44  ITS Experiment Rapid fermentation  double
## PE45     PE45  ITS Experiment Rapid fermentation  double
## PE46     PE46  ITS Experiment Rapid fermentation  double
## PE47     PE47  ITS Experiment Rapid fermentation  double
## PE51     PE51  ITS Experiment  Long fermentation  single
## PE52     PE52  ITS Experiment  Long fermentation  single
## PE53     PE53  ITS Experiment  Long fermentation  single
## PE54     PE54  ITS Experiment  Long fermentation  single
## PE55     PE55  ITS Experiment  Long fermentation  single
## PE56     PE56  ITS Experiment  Long fermentation  single
## PE6       PE6  ITS Experiment            Control control
## PE60     PE60  ITS Experiment  Long fermentation  double
## PE61     PE61  ITS Experiment  Long fermentation  double
## PE62     PE62  ITS Experiment  Long fermentation  double
## PE63     PE63  ITS Experiment  Long fermentation  double
## PE64     PE64  ITS Experiment  Long fermentation  double
## PE65     PE65  ITS Experiment  Long fermentation  double
## PE7       PE7  ITS Experiment            Control control
## PE8       PE8  ITS Experiment            Control control
## PE9       PE9  ITS Experiment            Control control
```
