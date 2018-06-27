# Session configuration ####

setwd("~/Documents/PHP/fungi")

library(ggplot2)
library(data.table)
library(phyloseq)
library(gplots)
library(VennDiagram)
library(vegan)
library(DESeq2)
library(doMC)
library(foreach)
library("pheatmap")
library(dplyr)
library(RColorBrewer)
options(warnings=-1)

# post-processing and error correction of metadata ####
# correct version is imported below
# this can be disregarded, but included for historical considerations

#phpmeta2<-read.table("Master_Label_Table_PHPv1.csv", header = T, na.strings = "", sep = ",")
#phpmeta2[2]<-NULL
#phpmeta3<-phpmeta2[,-1]
#rownames(phpmeta3)<-phpmeta2[,1]
#saveRDS(phpmeta2, file = "fungi_metadata.rds")

#phpmeta$System.loc<-as.character(phpmeta$System.loc)
#phpmeta$System.loc[phpmeta$System.loc=="CT-FF"]<-"CT" #other transformations done as well
#saveRDS(phpmeta, file = "phpmeta_corrected.rds")

# Construct initial phyloseq object ####
# create phyloseq object and change taxa names to something shorter

seqtab.new<-readRDS("fungi_seqs2.rds")
phpfuntax.new<-readRDS(file = "fungi_taxa2.rds")
phpfuntax.new <- cbind( phpfuntax.new, "seq" = rownames(phpfuntax.new), "fasta_name" = paste0("phpfung", seq(1:nrow(phpfuntax.new))))
phpmeta<-readRDS("phpmeta_corrected.rds")

# Fix Stoneville System.loc names
phpmeta$System.loc <- as.character(phpmeta$System.loc)
phpmeta$Herbicide_History <- as.character(phpmeta$Herbicide_History)
phpmeta$System.loc[phpmeta$Location == "Stoneville"] <- paste(phpmeta$System.loc[phpmeta$Location == "Stoneville"], phpmeta$Herbicide_History[phpmeta$Location == "Stoneville"], sep = "_")
phpmeta$System.loc <- as.factor(phpmeta$System.loc)
phpmeta$Herbicide_History <- as.factor(phpmeta$Herbicide_History)
phpmeta$year<-as.factor(phpmeta$year)

# Make alternate sampleID and variables
phpmeta$group <- factor(paste(phpmeta$System.loc, phpmeta$Glyphosphate_Treatment, sep = "_"))
phpmeta$q_id <- factor(paste(phpmeta$System.loc, phpmeta$Glyphosphate_Treatment, phpmeta$Soil_Zone, phpmeta$year, sep = "_"))

ps.new <- phyloseq(otu_table(seqtab.new, taxa_are_rows=FALSE), 
                   sample_data(phpmeta), 
                   tax_table(phpfuntax.new))
ntaxa(ps.new)

taxa_names(ps.new)[1:5]

if (identical(taxa_names(ps.new), as.character(phpfuntax.new[,"seq"]))) {
  taxa_names(ps.new) <- as.character(phpfuntax.new[,"fasta_name"])
}

taxa_names(ps.new)[1:5]

ps.new <- subset_taxa(ps.new, Kingdom == "k__Fungi")
ntaxa(ps.new)

# Ordination for all samples, all sites ####
# At this point in the workflow there is no filtering of low abundance taxa
# You might want to do that later, but you will need to build a filtered
# phyloseq object 
# remove Urbana sites and non-RR genotype samples first
ps.newREL <- subset_samples(ps.new, Location != "Urbana" & genotype == "RR")
ps.newREL <- prune_taxa(taxa_sums(ps.newREL) > 10, ps.newREL)

ps.newREL  = transform_sample_counts(ps.newREL, function(x) x / sum(x) )
set.seed(1978)
all.ord.nmds.bray <- ordinate(ps.newREL, method="PCoA", distance="bray") #, maxit = 30000, sratmax = 0.99999999, sfgrmin = 1e-10)#, sratmax = 0.9999999)#, k = 2, sfgrmin = 1e-8, sratmax = 0.999999, maxit = 30000)
all.plot<-plot_ordination(ps.newREL, all.ord.nmds.bray, type="samples", color="System.loc", shape = "crop") + 
  geom_hline(yintercept = 0, color = "grey90") + 
  geom_vline(xintercept = 0, color = "grey90") +
  geom_point(size = 2.5) + 
  theme_classic() + 
  scale_color_manual(values = c("#FF7F00", "#E31A1C", "#6A3D9A", "#B15928", "#33A02C","#1F78B4")) + 
  ggtitle("All sites")
all.plot
ggsave("all_sites_pcoa_reltrans_taxa.pdf", plot = all.plot, device = pdf, height = 6, width = 8)


# Run standard workflow in parallel ####
# create subsets for combinations of sample data values
# this version is for Location, crop 
# Genotype is hardcoded for this version. Search "genotype.row".

unique.meta <- as.data.frame(unique(phpmeta[c("Location", "crop")]))
unique.meta <- unique.meta[(unique.meta$Location != "Urbana"),]

# Set custom color palettes for each site when plotting PCoA
# These colors are used throughout all graphics. 
# Darker shades for System.loc, lighter for those that have been sprayed 
unique.meta$palette[unique.meta$Location == "Beltsville"] <- list(c("#FF7F00", "#FDBF6F", "#E31A1C", "#FB9A99", "#33A02C", "#B2DF8A", "#1F78B4", "#A6CEE3"))
unique.meta$palette[unique.meta$Location == "Stoneville"] <- list(c("#6A3D9A", "#CAB2D6", "#B15928", "#EEAD0E"))

registerDoMC(cores = 8)

ptm <- proc.time()
foreach(a=1:nrow(unique.meta), .packages = c("phyloseq")) %dopar% { 
  # get variable values
  loc <- as.character(unique.meta["Location"][a,])
  crop <- as.character(unique.meta["crop"][a,])
  pal <- unique.meta$palette[[a]]
  
  # make a name to recycle
  thing <- paste(loc, crop, sep = "_")
  
  # make folder tree to save output
  out.dir <- paste(getwd(),"vst",loc,thing, sep = "/")
  dir.create(out.dir, recursive = T)
  
  # determine the samples to work on
  #loc.list <- as.character(get_variable(ps.new, "Location")) == loc
  loc.row <- row.names(sample_data(ps.new)[sample_data(ps.new)$Location %in% loc, ])
  crop.row <- row.names(sample_data(ps.new)[sample_data(ps.new)$crop %in% crop, ])
  genotype.row <- row.names(phpmeta[phpmeta$genotype %in% "RR", ])
  common.row <- Reduce(intersect, list(loc.row, crop.row, genotype.row))
  
  # prune phyloseq object to desired samples
  ps.0 <- prune_samples(common.row, ps.new)
  ps.0 <- prune_taxa(taxa_sums(ps.0) > 10, ps.0)
  saveRDS(ps.0, file.path(out.dir, paste(thing,".rds", sep = "")))
  
  # estimate richness, selected metrics
  rich <- data.frame(estimate_richness(ps.0), measures = c("Observed", "Chao1", "Shannon", "Simpson"))
  row.names(rich) <- sub("X", "", row.names(rich))
  rich$samples <- row.names(rich)
  
  # DESeq2 variance stablization
  ps.ds <- phyloseq_to_deseq2(ps.0, ~ group + year + Sampling_date + group:Sampling_date)
  ps.ds <- estimateSizeFactors(ps.ds)
  ps.ds <- estimateDispersions(ps.ds)
  psVST <- getVarianceStabilizedData(ps.ds)
  saveRDS(psVST, file.path(out.dir, paste(thing, "_vst", ".rds", sep = "")))
  
  # perform ordinations
  psVST[psVST < 0.0] <- 0.0
  otu_table(ps.0) <- otu_table(psVST, taxa_are_rows = TRUE)
  bray <- ordinate(ps.0, method="PCoA", distance="bray")
  ps.plot<-plot_ordination(ps.0, bray, type="samples", color="group", shape="year") +
    ggtitle(paste(thing, "samples:", nsamples(ps.0), "taxa:", ntaxa(ps.0), sep = " ")) + 
    theme_classic() +
    geom_hline(yintercept = 0, color = "grey90") + 
    geom_vline(xintercept = 0, color = "grey90") +
    geom_point(size = 3) + 
    scale_color_manual(values = pal)
  ggsave(filename = file.path(out.dir, paste(thing,"pdf", sep = ".")), 
         plot = ps.plot, width = 20, height = 14, units = "cm", 
         device = "pdf", dpi = 300)
  
  # save vectors from PCoA and richness estimates for use in qiime2 longitudinal test
  axisvals <- data.frame(bray$vectors)[,1:3]
  axisvals$samples <- row.names(axisvals)
  
  meta <- data.frame(sample_data(ps.0))
  meta$samples <- row.names(meta)
  qiime <- left_join(meta, rich, by = "samples") %>% 
    left_join(., axisvals, by = "samples")
  write.table(qiime, sep = "\t", 
              file = file.path(out.dir, paste(thing, "qiime", "txt", sep = ".")), 
              row.names = F, quote = F)
  
  # perform PERMANOVA 
  out <- adonis(t(otu_table(ps.0)) ~ System.loc * year * Soil_Zone * Glyphosphate_Treatment * Sampling_date, 
                strata = sample_data(ps.0)$Loc_plot_ID, as(sample_data(ps.0), "data.frame"))
  
  write.table(as.matrix(cbind(rownames(out$aov.tab), out$aov.tab)), 
              file = file.path(out.dir, paste(thing,"perm", "txt", sep = ".")), 
              row.names = F, 
              col.names = c("Variable", colnames(out$aov.tab)),
              quote = F, sep = "\t")
}
proc.time() - ptm

# Unix and Perl one-liner to combine all PERMANOVA output. 
# Execute from top level directory for each site. Edit as LOCATION necessary.
# tail -n +1 LOCATION_*/*.perm.txt | perl -p -e 's/^.*\/(\w+)\.perm\.txt.*$/$1/g' > all_permanova.txt

# Fusarium Figures ####

#Stoneville heatmap palette

sv_colors <- list(
  Glyphosphate_Treatment = c(no_spray = "#9BCD9B", spray = "#CDAA7D"),
  System.loc = c(NT_none = "#B15928", NT_15yrs = "#6A3D9A"),
  Sampling_date = c(post = "grey80", pre = "grey50")
)

# Stoneville Soy Fusarium 
sv.soy <- readRDS("vst/Stoneville/Stoneville_soy/Stoneville_soy.rds")
sv.soy.VST <- readRDS("vst/Stoneville/Stoneville_soy/Stoneville_soy_vst.rds")

sv.soy.VST[sv.soy.VST < 0.0] <- 0.0
otu_table(sv.soy) <- otu_table(sv.soy.VST, taxa_are_rows = TRUE)
sv.soy <- subset_samples(sv.soy, Sampling_date == "post")
sv.soy <- prune_taxa(taxa_sums(sv.soy) > 0, sv.soy)

df <- data.frame(sample_data(sv.soy)[,c("System.loc", "Glyphosphate_Treatment")])
#df$Sampling_date <- ordered(df$Sampling_date, levels = c("pre", "post"))
df <- df[with(df, order(Glyphosphate_Treatment,System.loc)), ]
select <- taxa_names(subset_taxa(sv.soy, Genus == "g__Fusarium"))
vsttest <- as.data.frame(sv.soy.VST[select,])
vsttest <-vsttest[,rownames(df)]
heat.tax <- as.data.frame(tax_table(sv.soy)[row.names(vsttest)])
heat.tax <- lapply(heat.tax, function(x) gsub("^[a-z]_\\w", "", x, perl = T))
row.names(vsttest) <- paste(heat.tax$Genus, heat.tax$Species, heat.tax$fasta_name, sep = "_")
vsttest <- vsttest[order(row.names(vsttest)),]

pheatmap(vsttest, cluster_rows=FALSE, show_rownames=T, 
         show_colnames = F,
         breaks = NA,
         cluster_cols=F, 
         annotation_col=df,
         main = "Stoneville Soy Fusarium",
         gaps_col = 32,
         annotation_colors = sv_colors,
         border_color = NA,
         filename = "vst/Stoneville/Stoneville_soy/sv_soy_vst_fusarium_heat2.pdf")

# Stoneville Corn Fusarium
sv.corn <- readRDS("vst/Stoneville/Stoneville_corn/Stoneville_corn.rds")
sv.corn.VST <- readRDS("vst/Stoneville/Stoneville_corn/Stoneville_corn_vst.rds")

sv.corn.VST[sv.corn.VST < 0.0] <- 0.0
otu_table(sv.corn) <- otu_table(sv.corn.VST, taxa_are_rows = TRUE)
sv.corn <- subset_samples(sv.corn, Sampling_date == "post")
sv.corn <- prune_taxa(taxa_sums(sv.corn) > 0, sv.corn)

df <- data.frame(sample_data(sv.corn)[,c("System.loc", "Glyphosphate_Treatment")])
#df$Sampling_date <- ordered(df$Sampling_date, levels = c("pre", "post"))
df <- df[with(df, order(Glyphosphate_Treatment,System.loc)), ]
select <- taxa_names(subset_taxa(sv.corn, Genus == "g__Fusarium"))
vsttest <- as.data.frame(sv.corn.VST[select,])
vsttest <-vsttest[,rownames(df)]
heat.tax <- as.data.frame(tax_table(sv.corn)[row.names(vsttest)])
heat.tax <- lapply(heat.tax, function(x) gsub("^[a-z]_\\w", "", x, perl = T))
row.names(vsttest) <- paste(heat.tax$Genus, heat.tax$Species, heat.tax$fasta_name, sep = "_")
vsttest <- vsttest[order(row.names(vsttest)),]

pheatmap(vsttest, cluster_rows=FALSE, show_rownames=T, 
         show_colnames = F,
         breaks = NA,
         cluster_cols=F, 
         annotation_col=df,
         main = "Stoneville Corn Fusarium",
         annotation_colors = sv_colors,
         gaps_col = 32,
         border_color = NA,
         filename = "vst/Stoneville/Stoneville_corn/sv_corn_vst_fusarium_heat2.pdf")

#Beltsville heatmap palette

fsp_colors <- list(
  Glyphosphate_Treatment = c(no_spray = "#9BCD9B", spray = "#CDAA7D"),
  System.loc = c(CT = "#FF7F00", NT = "#E31A1C", Org_3 = "#33A02C", Org_6 = "#1F78B4"),
  Sampling_date = c(post = "grey80", pre = "grey50")
)

# Beltsville Soy Fusarium
fsp.soy <- readRDS("vst/Beltsville/Beltsville_soy/Beltsville_soy.rds")
fsp.soy.VST <- readRDS("vst/Beltsville/Beltsville_soy/Beltsville_soy_vst.rds")
fsp.soy <- subset_samples(fsp.soy, Sampling_date == "post")
fsp.soy <- prune_taxa(taxa_sums(fsp.soy) > 0, fsp.soy)
fsp.soy.VST[fsp.soy.VST < 0.0] <- 0.0
otu_table(fsp.soy) <- otu_table(fsp.soy.VST, taxa_are_rows = TRUE)

df <- data.frame(sample_data(fsp.soy)[,c("System.loc", "Glyphosphate_Treatment")])
df <- df[with(df, order(Glyphosphate_Treatment,System.loc)), ]
select <- taxa_names(subset_taxa(fsp.soy, Genus == "g__Fusarium"))
vsttest <- as.data.frame(fsp.soy.VST[select,])
vsttest <-vsttest[,rownames(df)]
heat.tax <- as.data.frame(tax_table(fsp.soy)[row.names(vsttest)])
heat.tax <- lapply(heat.tax, function(x) gsub("^[a-z]_\\w", "", x, perl = T))
row.names(vsttest) <- paste(heat.tax$Genus, heat.tax$Species, heat.tax$fasta_name, sep = "_")
vsttest <- vsttest[order(row.names(vsttest)),]

pheatmap(vsttest, cluster_rows=FALSE, show_rownames=T, 
         show_colnames = F,
         breaks = NA,
         cluster_cols=F, 
         annotation_col=df,
         main = "Beltsville Soy Fusarium",
         annotation_colors = fsp_colors,
         gaps_col = 64,
         border_color = NA,
         filename = "vst/Beltsville/Beltsville_soy/fsp_soy_vst_fusarium_heat2.pdf")

# Beltsville Corn Fusarium
fsp.corn <- readRDS("vst/Beltsville/Beltsville_corn/Beltsville_corn.rds")
fsp.corn.VST <- readRDS("vst/Beltsville/Beltsville_corn/Beltsville_corn_vst.rds")
fsp.corn <- subset_samples(fsp.corn, Sampling_date == "post")
fsp.corn <- prune_taxa(taxa_sums(fsp.corn) > 0, fsp.corn)

fsp.corn.VST[fsp.corn.VST < 0.0] <- 0.0
otu_table(fsp.corn) <- otu_table(fsp.corn.VST, taxa_are_rows = TRUE)

df <- data.frame(sample_data(fsp.corn)[,c("System.loc", "Glyphosphate_Treatment")])
df <- df[with(df, order(Glyphosphate_Treatment,System.loc)), ]
select <- taxa_names(subset_taxa(fsp.corn, Genus == "g__Fusarium"))
vsttest <- as.data.frame(fsp.corn.VST[select,])
vsttest <-vsttest[,rownames(df)]
heat.tax <- as.data.frame(tax_table(fsp.corn)[row.names(vsttest)])
heat.tax <- lapply(heat.tax, function(x) gsub("^[a-z]_\\w", "", x, perl = T))
row.names(vsttest) <- paste(heat.tax$Genus, heat.tax$Species, heat.tax$fasta_name, sep = "_")
vsttest <- vsttest[order(row.names(vsttest)),]

pheatmap(vsttest, cluster_rows=FALSE, show_rownames=T, 
         show_colnames = F,
         breaks = NA,
         cluster_cols=F, 
         annotation_col=df,
         main = "Beltsville Corn Fusarium",
         annotation_colors = fsp_colors,
         gaps_col = 64,
         border_color = NA,
         filename = "vst/Beltsville/Beltsville_corn/fsp_corn_vst_fusarium_heat2.pdf")

# Differentially abundant taxa ####
res <- results(diagdds, contrast = c("group", "CT_spraypost", "CT_no_spraypost"))
alpha = 0.01
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(fsp.corn)[rownames(sigtab), ], "matrix"))
#head(sigtab)

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Order, function(x) max(x))
x = sort(x, TRUE)
sigtab$Order = factor(as.character(sigtab$Order), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
ggplot(sigtab, aes(x=log2FoldChange, y=Genus, color=Order)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) +
  ggtitle("CT_spraypost CT_no_spraypost")

#### Rarefied dataset ####
# Import the rarefied phyloseq object ####

ps.rare <- readRDS("fungi_RAREcomDat.rds")
stone <- data.frame(sample_data(ps.rare))
stone$System.loc <- as.character(stone$System.loc)
stone$Herbicide_History <- as.character(stone$Herbicide_History)
stone$System.loc[stone$Location == "Stoneville"] <- paste(stone$System.loc[stone$Location == "Stoneville"], stone$Herbicide_History[stone$Location == "Stoneville"], sep = "_")
stone$System.loc <- as.factor(stone$System.loc)
stone$Herbicide_History <- as.factor(stone$Herbicide_History)
stone$year<-as.factor(stone$year)
stone$group <- factor(paste(stone$System.loc, stone$Glyphosphate_Treatment, sep = "_"))
stone$q_id <- factor(paste(stone$System.loc, stone$Glyphosphate_Treatment, stone$Soil_Zone, stone$year, sep = "_"))

sample_data(ps.rare) <- stone
#sample_data(ps.rare)$location_id <- paste(sample_data(ps.rare)$Location, sample_data(ps.rare)$Loc_plot_ID, sep = "_")

# NMDS for all samples, all sites ####
 
# remove Urbana sites and non-RR genotype samples first
ps.rareREL <- subset_samples(ps.rare, Location != "Urbana" & genotype == "RR")
ps.rareREL <- prune_taxa(taxa_sums(ps.rareREL) > 10, ps.rareREL)

ps.rareREL  = transform_sample_counts(ps.rareREL, function(x) x / sum(x) )
set.seed(1978)
all.ord.nmds.bray <- ordinate(ps.rareREL, method="PCoA", distance="bray") #, maxit = 30000, sratmax = 0.99999999, sfgrmin = 1e-10)#, sratmax = 0.9999999)#, k = 2, sfgrmin = 1e-8, sratmax = 0.999999, maxit = 30000)
all.plot<-plot_ordination(ps.rareREL, all.ord.nmds.bray, type="samples", color="System.loc", shape = "crop") + 
  geom_hline(yintercept = 0, color = "grey90") + 
  geom_vline(xintercept = 0, color = "grey90") +
  geom_point(size = 2.5) + 
  theme_classic() + 
  scale_color_manual(values = c("#FF7F00", "#E31A1C", "#6A3D9A", "#B15928", "#33A02C","#1F78B4")) + 
  ggtitle("All sites")
all.plot
ggsave("all_sites_pcoa_reltrans_rare_taxa.pdf", plot = all.plot, device = pdf, height = 6, width = 8)


# Run standard workflow in parallel ####
# create subsets for combinations of sample data values
# this version is for Location, crop 
# Genotype is hardcoded for this version. Search "genotype.row".
meta <-data.frame(sample_data(ps.rare))
unique.meta <- unique(meta[c("Location", "crop")])
unique.meta <- unique.meta[(unique.meta$Location != "Urbana"),]

# Set custom color palettes for each site when plotting PCoA
# These colors are used throughout all graphics. 
# Darker shades for System.loc, lighter for those that have been sprayed 
unique.meta$palette[unique.meta$Location == "Beltsville"] <- list(c("#FF7F00", "#FDBF6F", "#E31A1C", "#FB9A99", "#33A02C", "#B2DF8A", "#1F78B4", "#A6CEE3"))
unique.meta$palette[unique.meta$Location == "Stoneville"] <- list(c("#6A3D9A", "#CAB2D6", "#B15928", "#EEAD0E"))

registerDoMC(cores = 8)

ptm <- proc.time()
foreach(a=1:nrow(unique.meta), .packages = c("phyloseq")) %dopar% { 
  
  # get variable values
  loc <- as.character(unique.meta["Location"][a,])
  crop <- as.character(unique.meta["crop"][a,])
  pal <- unique.meta$palette[[a]]
  
  # make a name to recycle
  thing <- paste(loc, crop, sep = "_")
  
  # make folder tree to save output
  out.dir <- paste(getwd(),"grid",loc,thing, sep = "/")
  dir.create(out.dir, recursive = T)
  
  # determine the samples to work on
  #loc.list <- as.character(get_variable(ps.rare, "Location")) == loc
  loc.row <- row.names(sample_data(ps.rare)[sample_data(ps.rare)$Location %in% loc, ])
  crop.row <- row.names(sample_data(ps.rare)[sample_data(ps.rare)$crop %in% crop, ])
  genotype.row <- row.names(meta[meta$genotype %in% "RR", ])
  common.row <- Reduce(intersect, list(loc.row, crop.row, genotype.row))
  
  # prune phyloseq object to desired samples
  ps.0 <- prune_samples(common.row, ps.rare)
  ps.0 <- prune_taxa(taxa_sums(ps.0) > 10, ps.0)
  saveRDS(ps.0, file.path(out.dir, paste(thing,".rds", sep = "")))
  
  # estimate richness, selected metrics
  rich <- data.frame(estimate_richness(ps.0, measures = c("Observed", "Chao1", "Shannon", "Simpson")))
  row.names(rich) <- sub("X", "", row.names(rich))
  rich$samples <- row.names(rich)
  sample_data(ps.0)$Shannon <- rich$Shannon
  
  rich.grid <- ggplot(data = sample_data(ps.0), aes(x = Sampling_date, y = Shannon, group = Glyphosphate_Treatment)) + 
    #geom_point(aes(color = year, shape = Soil_Zone)) + 
    geom_point(position = position_jitter(width = 0.2, height = 0.2), alpha = 0.2) + 
    #stat_summary(fun.data = 'mean_sdl', geom = 'errorbar', width = 0.2, size = 1) +
    stat_summary(fun.y = mean, geom = 'point', size = 3, color = 'red') +
    stat_summary(fun.y = mean, geom = 'line', size = 1, color = 'red') +
    facet_grid(Glyphosphate_Treatment ~ System.loc) +
    theme_bw() +
    ggtitle(thing) +
    xlab("Sampling Date") +
    ylab("Richness") + 
    scale_x_discrete(limits=c("pre", "post")) +
    theme(axis.text.x=element_text(angle=90,hjust=1))
  
  ggsave(filename = file.path(out.dir, paste(thing,"sys_loc_rich","pdf", sep = ".")), 
         plot = rich.grid, width = 20, height = 10, units = "cm", 
         device = "pdf", dpi = 300)
  
  # plot richness after merging all samples per system
  ps.0.st <- merge_samples(ps.0, "System.loc")
  rich.p <- plot_richness(ps.0.st, measures = c("Observed", "Chao1", "Shannon", "Simpson"))
  ggsave(filename = file.path(out.dir, paste(thing,"rare_rich","pdf", sep = ".")), 
         plot = rich.p, width = 20, height = 14, units = "cm", 
         device = "pdf", dpi = 300)
  rm(ps.0.st, rich.p)
  
  # perform ordinations
  bray <- ordinate(ps.0, method="PCoA", distance="bray")
  ps.plot<-plot_ordination(ps.0, bray, type="samples", color="group", shape="year") +
    ggtitle(paste(thing, "samples:", nsamples(ps.0), "taxa:", ntaxa(ps.0), sep = " ")) + 
    theme_classic() +
    geom_hline(yintercept = 0, color = "grey90") + 
    geom_vline(xintercept = 0, color = "grey90") +
    geom_point(size = 3) + 
    scale_color_manual(values = pal)
  ggsave(filename = file.path(out.dir, paste(thing,"pdf", sep = ".")), 
         plot = ps.plot, width = 20, height = 14, units = "cm", 
         device = "pdf", dpi = 300)
  
  # save vectors from PCoA and richness estimates for use in qiime2 longitudinal test
  axisvals <- data.frame(bray$vectors)[,1:3]
  axisvals$samples <- row.names(axisvals)
  
  meta <- data.frame(sample_data(ps.0))
  meta$samples <- row.names(meta)
  qiime <- left_join(meta, rich, by = "samples") %>% 
    left_join(., axisvals, by = "samples")
  write.table(qiime, sep = "\t", 
              file = file.path(out.dir, paste(thing, "qiime", "txt", sep = ".")), 
              row.names = F, quote = F)
  
  # perform PERMANOVA 
  out <- adonis(otu_table(ps.0) ~ System.loc * year * Soil_Zone * Glyphosphate_Treatment * Sampling_date, 
                strata = sample_data(ps.0)$Loc_plot_ID, 
                data = as(sample_data(ps.0), "data.frame"),
                method = "bray")
  
  write.table(as.matrix(cbind(rownames(out$aov.tab), out$aov.tab)), 
              file = file.path(out.dir, paste(thing,"perm", "txt", sep = ".")), 
              row.names = F, 
              col.names = c("Variable", colnames(out$aov.tab)),
              quote = F, sep = "\t")
}
proc.time() - ptm

# Unix and Perl one-liner to combine all PERMANOVA output. 
# Execute from top level directory for each site. Edit as necessary.
# tail -n +1 Stoneville_*/*.perm.txt | perl -p -e 's/^.*\/(\w+)\.perm\.txt.*$/$1/g' > all_permanova.txt

# END ####