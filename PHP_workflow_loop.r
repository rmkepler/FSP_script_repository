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
options(warnings=-1)

# creation of phyloseq object ####

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

# Import output from dada2 and metadata ####
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

# Run standard workflow in parallel ####
# create subsets for combinations of sample data values
# this version is for Location, Soil_Zone, crop 
# Genotype is hardcoded for this version. See line BLAH

unique.meta <- as.data.frame(unique(phpmeta[c("Location", "crop")]))
unique.meta <- unique.meta[(unique.meta$Location == "Stoneville"),]

registerDoMC(cores = 8)

ptm <- proc.time()
foreach(a=1:nrow(unique.meta), .packages = c("phyloseq")) %dopar% { 
  # get variable values
  loc <- as.character(unique.meta["Location"][a,])
  crop <- as.character(unique.meta["crop"][a,])
  
  # make a name to recycle
  thing <- paste(loc, crop, sep = "_")
  
  # make folder tree to save output
  out.dir <- paste(getwd(),loc,thing, sep = "/")
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
  
  # DESeq2 variance stablization
  ps.ds <- phyloseq_to_deseq2(ps.0, ~ group + year + Sampling_date + group:Sampling_date)
  ps.ds = estimateSizeFactors(ps.ds)
  ps.ds = estimateDispersions(ps.ds)
  psVST = getVarianceStabilizedData(ps.ds)
  saveRDS(psVST, file.path(out.dir, paste(thing, "_vst", ".rds", sep = "")))
  
  # perform ordinations
  psVST[psVST < 0.0] <- 0.0
  otu_table(ps.0) <- otu_table(psVST, taxa_are_rows = TRUE)
  bray <- ordinate(ps.0, method="PCoA", distance="bray")
  ps.plot<-plot_ordination(ps.0, bray, type="samples", color="group", shape="year") +
    ggtitle(paste(thing, "samples:", nsamples(ps.0), "taxa:", ntaxa(ps.0), sep = " ")) + 
    geom_point(size = 3)
  ggsave(filename = file.path(out.dir, paste(thing,"pdf", sep = ".")), plot = ps.plot, width = 20, height = 14, units = "cm", device = "pdf", dpi = 300)
  
  # save vectors from PCoA for use in qiime2 longitudinal test
  axisvals <- bray$vectors
  qiime <- merge(as.data.frame(sample_data(ps.0)), axisvals, by = "row.names")
  write.table(qiime, sep = "\t", file = file.path(out.dir, paste(thing, "qiime", "txt", sep = ".")), row.names = F, quote = F)
  
  # perform PERMANOVA 
  out <- adonis(t(otu_table(ps.0)) ~ System.loc * year * Soil_Zone * Glyphosphate_Treatment * Sampling_date, strata = sample_data(ps.0)$Loc_plot_ID, as(sample_data(ps.0), "data.frame"))
  
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

# Test some of the saved subsets ####

sv.soy <- readRDS("Stoneville/Stoneville_soy/Stoneville_soy.rds")
sv.soy.VST <- readRDS("Stoneville/Stoneville_soy/Stoneville_soy_vst.rds")

sv.soy.VST[sv.soy.VST < 0.0] <- 0.0
otu_table(sv.soy) <- otu_table(sv.soy.VST, taxa_are_rows = TRUE)

# taxon abundance across samples ####

sv.soy.DS <-phyloseq_to_deseq2(sv.soy, ~ group + year + Sampling_date + group:Sampling_date)

#diagdds = DESeq(fsp.soy.rhizDS, test="Wald", fitType="parametric")

#vsd <- vst(diagdds)
select <- order(rowMeans(counts(sv.soy.DS,normalized=F)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(sv.soy.DS)[,c("group","Sampling_date")])

df <- df[with(df, order(group,Sampling_date)), ]
vsdASSAY <- sv.soy.VST[select,]

vsdASSAY <- vsdASSAY[,rownames(df)]

heat.tax <- as.data.frame(tax_table(sv.soy)[row.names(sv.soy.VST)])
heat.tax <- lapply(heat.tax, function(x) gsub("^[a-z]_\\w", "", x, perl = T))

row.names(sv.soy.VST) <- paste(heat.tax$fasta_name, heat.tax$Genus, heat.tax$Species, sep = "_")

#assay(vsd)[select,]
pheatmap(vsdASSAY, cluster_rows=FALSE, show_rownames=T,
         cluster_cols=F, annotation_col=df) #,
         #filename = "sv_soy_heat.pdf")

# fsp heat maps ####
fsp.corn <- readRDS("Beltsville/Beltsville_corn/Beltsville_corn.rds")

#fsp.corn.VST <- readRDS("Beltsville/Beltsville_corn/Beltsville_corn_vst.rds")

#fsp.corn.VST[fsp.corn.VST < 0.0] <- 0.0
#otu_table(fsp.corn) <- otu_table(fsp.corn.VST, taxa_are_rows = TRUE)

# taxon abundance across samples ####
sample_data(fsp.corn)$group <- factor(paste0(sample_data(fsp.corn)$group, sample_data(fsp.corn)$Sampling_date))
fsp.corn.DS <-phyloseq_to_deseq2(fsp.corn, ~ group)
diagdds = DESeq(fsp.corn.DS, test="Wald", fitType="parametric")

vsd <- vst(diagdds)
select <- order(rowMeans(counts(diagdds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
select <- taxa_names(subset_taxa(fsp.corn, Order == "o__Hypocreales"))
df <- as.data.frame(colData(diagdds)[,c("group","Sampling_date")])

df <- df[with(df, order(group,Sampling_date)), ]
vsdASSAY <-assay(vsd)[select,]

vsdASSAY <- vsdASSAY[,rownames(df)]

heat.tax <- as.data.frame(tax_table(fsp.corn)[row.names(vsdASSAY)])
heat.tax <- lapply(heat.tax, function(x) gsub("^[a-z]_\\w", "", x, perl = T))

row.names(vsdASSAY) <- paste(heat.tax$Genus, heat.tax$Species, heat.tax$fasta_name, sep = "_")
vsdASSAY <- vsdASSAY[order(row.names(vsdASSAY)), ]
#assay(vsd)[select,]
pheatmap(vsdASSAY, cluster_rows=FALSE, show_rownames=T, 
         show_colnames = F,
         cluster_cols=F, annotation_col=df) #,
         #filename = "fsp_corn_vst_nectriaceae_heat.pdf")
# Differentially abundant taxa 
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
# END ####
