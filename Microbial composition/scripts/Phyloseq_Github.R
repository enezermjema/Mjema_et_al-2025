library("phyloseq")
library("ggplot2")
library("dplyr")
library("Polychrome")

# Loading tax and asv tables, and metadata
seqtab.nochim <- read.csv2("otu_table.csv", row.names = 1)
seqtab.nochim[is.na(seqtab.nochim)] <- 0

taxa <- read.csv2("tax_table.csv", row.names = 1)

taxa <- as.matrix(taxa)

metadata_new <- read.csv2("metadata_final.csv", stringsAsFactors = TRUE, row.names = 1)

# Making phyloseq object
ps = phyloseq(
  otu_table(seqtab.nochim, taxa_are_rows=FALSE),
  sample_data(metadata_new),
  tax_table(taxa))
ps

# Analysis with ONLY classified reads
# Filtering out unclassified and not-non-virus, respectively
unique(tax_table(ps)[, "Kingdom"])

ps_filt <- subset_taxa(ps, Kingdom != "unclassified" & Kingdom != "cannot be assigned to a (non-viral) genus")
unique(tax_table(ps_filt)[, "Kingdom"])

# Taxonomy plots - top 15 Genus
# Agglomerate taxa at genus level
ps_filt_glommed_gen <- tax_glom(ps_filt, taxrank = "Genus")

# Top N taxa
N <- 15
top <- names(sort(taxa_sums(ps_filt_glommed_gen), decreasing = TRUE))[1:N]

# Subset object to top N taxa
ps_filt_glommed_gen.top <- prune_taxa(top, ps_filt_glommed_gen)

# Calculating the percentage of reads in each sample considering only the top 15 ASVs
sort(sample_sums(ps_filt_glommed_gen))
sort(sample_sums(prune_taxa(top, ps_filt_glommed_gen)))
boxplot(sort(sample_sums(prune_taxa(top, ps_filt_glommed_gen))/sample_sums(ps_filt_glommed_gen)))

# Building color palette
pallet1 <- createPalette(15,  c("#ff0000", "#00ff00", "#0000ff"))
swatch(pallet1)
names(pallet1) <- NULL

# Plotting
pdf("Figure_3b.pdf", width = 28, height = 8)
plot_bar(ps_filt_glommed_gen.top, fill = "Genus") +
  geom_bar(aes(), stat="identity", position="stack", colour=NA) + 
  scale_fill_manual(values = pallet1) + 
  facet_wrap(~Season, scales = "free_x") +
  scale_x_discrete(name = "Samples") +
  theme(axis.text.x = element_blank(), axis.ticks.x=element_blank(), 
        text = element_text(size = 22), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 
dev.off()
