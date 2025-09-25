plot_heatmap <- function(ps, taxa = "Genus", log.transform = TRUE){
  require(dplyr)
  require(phyloseq)
  require(ggplot2)
  
  ps <- prune_taxa(taxa_sums(ps) > 0, ps)
  ps <- tax_glom(ps, taxa)
  sig.taxa.long <- psmelt(ps) %>%
    arrange(Phylum) %>% 
    mutate(row = row_number())
  
  sig.taxa.long$Abundance <- as.numeric(sig.taxa.long$Abundance)
  sig.taxa.long$Taxa <- sig.taxa.long[,taxa]
  
  sig.taxa.long[sig.taxa.long == "Burkholderia-Caballeronia-Paraburkholderia"] <- "Burkholderia"
  sig.taxa.long[sig.taxa.long == "Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium"] <- "Pararhizobium"
  
  ggplot(sig.taxa.long, aes(x = SampleID, y = reorder(Taxa, row))) +
    {if(log.transform) geom_tile(aes(fill=log(Abundance)))} +
    {if(!log.transform) geom_tile(aes(fill=Abundance))} +
    scale_fill_distiller(palette = "OrRd", trans = "reverse") +
    facet_grid(rows = vars(Phylum), scales = "free", space = "free") +
    theme_light() +
    theme(strip.text.y = element_text(angle = 0),
          panel.spacing = unit(0.02,'lines'),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    xlab("") + ylab("")
}


ps.major <- prune_taxa(taxa_sums(ps) > 100, ps)

plot_heatmap(ps.major, taxa = "Family") + 
  facet_grid(rows = vars(Phylum), 
             cols = vars(Group, BioRepeat), scales = "free", space = "free")