# Network drawer
# For a threshold, we use significance of correlation by permutations

#create cor graph by qgraph
require(viridis)
require(qgraph)
library(scales)

some_ps <- ps.f
otus.pr.Cd <- some_ps@otu_table@.Data
tax.Cd.d <- as.data.frame(some_ps@tax_table@.Data)
lev <- levels(tax.Cd.d$Order)
len_lev <- length(lev)
pal <- hue_pal(h = c(0, 360), c = 100, l = 30, direction = -1)(len_lev)
vertex.col <- pal[tax.Cd.d$Order]
qgraph(cor(otus.pr.Cd), layout = "spring", minimum = "sig", alpha=0.05, groups = tax.Cd.d$Order,
       sampleSize = 16, graph = "cor", threshold = "BH", vsize = 3, label.cex = 1.8)
```

# Simple test for network drawer ability

```{r}

library(vegan)

some_ps <- ps.f.i.Cdt.wet
otus.dt <- as.data.frame(some_ps@otu_table@.Data)
mean(vegdist(otus.dt, method = "jaccard"))

seq.list <- seq(0, 500, by=10)
jac.steppo <- function(some_number){
  some_ps <- prune_taxa(taxa_sums(some_ps) > some_number , some_ps) 
  otus.dt <- as.data.frame(some_ps@otu_table@.Data)
  jacco.mean <- mean(vegdist(otus.dt, method = "jaccard"))
  return(jacco.mean)
}

some_ps <- ps.f.i
seq.list <- seq(0, 600, by=1)
jacco.list <- sapply(seq.list, jac.steppo)
plot(jacco.list)
```
