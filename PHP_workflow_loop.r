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

phpmeta<-readRDS("phpmeta_corrected.rds")
phpmeta$year<-as.factor(phpmeta$year)
seqtab.new<-readRDS("fungi_seqs2.rds")
phpfuntax.new<-readRDS(file = "fungi_taxa2.rds")
phpfuntax.new <- cbind( phpfuntax.new, "seq" = rownames(phpfuntax.new), "fasta_name" = paste0("phpfung", seq(1:nrow(phpfuntax.new))))
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

ntaxa(ps.new)
ps.new <- subset_taxa(ps.new, Kingdom == "k__Fungi")
ntaxa(ps.new)

# Run standard workflow in parallel ####
# create subsets for combinations of sample data values
# this version is for Location, Soil_Zone, crop 
# Genotype is hardcoded for this version. See line BLAH

unique.meta <- as.data.frame(unique(phpmeta[c("Location", "Soil_Zone", "crop")]))

unique.meta <- unique.meta[(unique.meta$Location == "Beltsville"),]

registerDoMC(cores = 8)

#ptm <- proc.time()
foreach(a=1:nrow(unique.meta), .packages = c("phyloseq")) %dopar% { 
  # get variable values
  loc <- as.character(unique.meta["Location"][a,])
  soil <- as.character(unique.meta["Soil_Zone"][a,])
  crop <- as.character(unique.meta["crop"][a,])
  
  # make a name to recycle
  thing <- paste(loc, soil, crop, sep = "_")
  
  # make folder tree to save output
  out.dir <- paste(getwd(),loc,thing, sep = "/")
  dir.create(out.dir, recursive = T)
  
  # determine the samples to work on
  #loc.list <- as.character(get_variable(ps.new, "Location")) == loc
  loc.row <- row.names(sample_data(ps.new)[sample_data(ps.new)$Location %in% loc, ])
  crop.row <- row.names(sample_data(ps.new)[sample_data(ps.new)$crop %in% crop, ])
  soil.row <- row.names(sample_data(ps.new)[sample_data(ps.new)$Soil_Zone %in% soil, ])
  genotype.row <- row.names(phpmeta[phpmeta$genotype %in% "RR", ])
  common.row <- Reduce(intersect, list(loc.row, crop.row, soil.row, genotype.row))
  
  # prune phyloseq object to desired samples
  ps.0 <- prune_samples(common.row, ps.new)
  ps.0 <- prune_taxa(taxa_sums(ps.0) > 10, ps.0)
  saveRDS(ps.0, file.path(out.dir, paste(thing,".rds", sep = "")))
  
  # DESeq2 variance stablization
  ps.ds <- phyloseq_to_deseq2(ps.0, ~ group + Sampling_date + group:Sampling_date)
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
  write.table(qiime, sep = "\t", file = file.path(out.dir, paste(thing,"txt", sep = ".")), row.names = F, quote = F)
  
  # perform PERMANOVA 
  out <- adonis(t(otu_table(ps.0)) ~ System.loc * Glyphosphate_Treatment * Sampling_date, strata = sample_data(ps.0)$Loc_plot_ID, as(sample_data(ps.0), "data.frame"))
  write.table(as.matrix(out$aov.tab), file = file.path(out.dir, paste(thing,"perm", "txt", sep = "."))  , quote = F, sep = ",")
}
#proc.time() - ptm

# heatmap for some of the saved subsets ####

sv.soy.rhiz <- readRDS("Stoneville/Stoneville_rhizosphere_soy/Stoneville_rhizosphere_soy.rds")
sv.soy.rhiz.VST <- readRDS("Stoneville/Stoneville_rhizosphere_soy/Stoneville_rhizosphere_soy_vst.rds")

sv.soy.rhiz.VST[sv.soy.rhiz.VST < 0.0] <- 0.0
otu_table(sv.soy.rhiz) <- otu_table(sv.soy.rhiz.VST, taxa_are_rows = TRUE)
out <- adonis(t(otu_table(sv.soy.rhiz)) ~ System.loc * Glyphosphate_Treatment * Sampling_date, strata = sample_data(sv.soy.rhiz)$Loc_plot_ID, as(sample_data(sv.soy.rhiz), "data.frame"))


# homogeneity of variances ####
fsr.dist <- dist(otu_table(fsp.soy.rhiz), method = "euclidean")
fsr.beta <- betadisper(fsr.dist, sample_data(fsp.soy.rhiz)$Sampling_date)
anova(fsr.beta)
fsr.beta
 

gpt <- prune_taxa(names(sort(taxa_sums(fsp.soy.rhiz),TRUE)[1:300]), fsp.soy.rhiz)
plot_heatmap(gpt, sample.label="group")

# differential expression ####

fsp.soy.rhizDS <-phyloseq_to_deseq2(fsp.soy.rhiz, ~ group + Sampling_date + group:Sampling_date)

diagdds = DESeq(fsp.soy.rhizDS, test="Wald", fitType="parametric")

vsd <- vst(diagdds)
select <- order(rowMeans(counts(diagdds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(diagdds)[,c("group","Sampling_date")])

df <- df[with(df, order(group,Sampling_date)), ]
vsdASSAY <- assay(vsd)[select,]

vsdASSAY <- vsdASSAY[,rownames(df)]

heat.tax <- as.data.frame(tax_table(fsp.soy.rhiz)[row.names(vsdASSAY)])
heat.tax <- lapply(heat.tax, function(x) gsub("^[a-z]_\\w", "", x, perl = T))

row.names(vsdASSAY) <- paste(heat.tax$fasta_name, heat.tax$Genus, heat.tax$Species, sep = "_")

#assay(vsd)[select,]
pheatmap(vsdASSAY, cluster_rows=FALSE, show_rownames=T,
         cluster_cols=F, annotation_col=df)

# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))
#####
fsp.corn.rhiz <- readRDS("Beltsville/Beltsville_rhizosphere_corn/Beltsville_rhizosphere_corn.rds")
fsp.corn.rhizDS <-phyloseq_to_deseq2(fsp.corn.rhiz, ~ group + Sampling_date + group:Sampling_date)

diagdds = DESeq(fsp.corn.rhizDS, test="Wald", fitType="parametric")

vsd <- vst(diagdds)
select <- order(rowMeans(counts(diagdds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(diagdds)[,c("group","Sampling_date")])

df <- df[with(df, order(group,Sampling_date)), ]
vsdASSAY <- assay(vsd)[select,]

vsdASSAY <- vsdASSAY[,rownames(df)]

heat.tax <- as.data.frame(tax_table(fsp.corn.rhiz)[row.names(vsdASSAY)])
heat.tax <- lapply(heat.tax, function(x) gsub("^[a-z]_\\w", "", x, perl = T))

row.names(vsdASSAY) <- paste(heat.tax$fasta_name, heat.tax$Genus, heat.tax$Species, sep = "_")

#assay(vsd)[select,]
pheatmap(vsdASSAY, cluster_rows=FALSE, show_rownames=T,
         cluster_cols=F, annotation_col=df)


# END ####

fsp.RRsoy <- subset_samples(ps.new, Location == "Beltsville" & genotype == "RR" & crop == "corn" & Soil_Zone =="bulk") #  & Sampling_date == "pre" & crop == "soy"  & Glyphosphate_Treatment == "no_spray")
fsp.RRsoy <- prune_taxa(taxa_sums(fsp.RRsoy) > 10, fsp.RRsoy)
ntaxa(fsp.RRsoy)
nsamples(fsp.RRsoy)

# DESeq routines
# make DESeq object
fspDS <- phyloseq_to_deseq2(fsp.RRsoy, ~ Glyphosphate_Treatment)

# Set the model to something meaningful
design(fspDS) <- formula(~ group + Sampling_date + group:Sampling_date)

# Perform variance stabilization
# Start the clock
ptm <- proc.time()
fspDS = estimateSizeFactors(fspDS)
fspDS = estimateDispersions(fspDS)
fspVST = getVarianceStabilizedData(fspDS)
# Stop the clock
proc.time() - ptm

# use VST counts
# backup original FSP-Soy-RR phyloseq object
fsp.RRsoy0 <- fsp.RRsoy
fspVST0 <- fspVST
fspVST[fspVST < 0.0] <- 0.0

otu_table(fsp.RRsoy) <- otu_table(fspVST, taxa_are_rows = TRUE)
ntaxa(fsp.RRsoy)
fsp.ord.pca.bray <- ordinate(fsp.RRsoy, method="PCoA", distance="bray")
fsp.RRsoy.eigen <- eigen(fsp.ord.pca.bray)

fsp.plot<-plot_ordination(fsp.RRsoy, fsp.ord.pca.bray, type="samples", color="System.loc", shape="Glyphosphate_Treatment") 
fsp.plot + ggtitle("FSP-Soy") +
  geom_point(size = 3)

fsp.RRsoy.axisvals <- fsp.ord.pca.bray$vectors

fsp.soy.qiime <- merge(as.data.frame(sample_data(fsp.RRsoy)), fsp.RRsoy.axisvals, by = "row.names")
write.table(fsp.soy.qiime, sep = "\t", file = "fsp.soy.qiime.txt", row.names = F, quote = F)
