# session configuration ####

setwd("~/Documents/PHP/fungi")

library(ggplot2)
library(data.table)
library(phyloseq)
library(gplots)
library(VennDiagram)
library(DESeq2)

# creation of phyloseq object ####

# post-processing and error correction of metadata.
# correct version is imported below
##
#phpmeta2<-read.table("Master_Label_Table_PHPv1.csv", header = T, na.strings = "", sep = ",")
#phpmeta2[2]<-NULL
#phpmeta3<-phpmeta2[,-1]
#rownames(phpmeta3)<-phpmeta2[,1]
#saveRDS(phpmeta2, file = "fungi_metadata.rds")

#phpmeta$System.loc<-as.character(phpmeta$System.loc)
#phpmeta$System.loc[phpmeta$System.loc=="CT-FF"]<-"CT" #other transformations done as well
#saveRDS(phpmeta, file = "phpmeta_corrected.rds")

phpmeta<-readRDS("phpmeta_corrected.rds")
phpmeta$year<-as.factor(phpmeta$year)
seqtab.new<-readRDS("fungi_seqs2.rds")
phpfuntax.new<-readRDS(file = "fungi_taxa2.rds")

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

#rm(seqtab.nochim, phpmeta, phpfuntax)


# Raw data summaries ####

# determine shared taxa between sites 
bv.names <- taxa_names(prune_taxa(taxa_sums(subset_samples(ps.new, Location=="Beltsville")) > 0, ps.new))
sv.names <- taxa_names(prune_taxa(taxa_sums(subset_samples(ps.new, Location=="Stoneville")) > 0, ps.new))
ur.names <- taxa_names(prune_taxa(taxa_sums(subset_samples(ps.new, Location=="Urbana")) > 0, ps.new))

sharetax <- venn.diagram(list(Beltsville = bv.names, Stoneville = sv.names, Urbana = ur.names), fill = 2:4, alpha = 0.3, filename = NULL)
grid.draw(sharetax)

# NMDS for all samples, all sites 
# At this point in the workflow there is no filtering of low abundance taxa
# You might want to do that later, but you will need to build a filtered
# phyloseq object 

ps.newREL  = transform_sample_counts(ps.new, function(x) x / sum(x) )
set.seed(1977)
all.ord.nmds.bray <- ordinate(ps.newREL, k = 4, method="NMDS", distance="bray", maxit = 30000, sratmax = 0.9999999, sfgrmin = 1e-10)#, sratmax = 0.9999999)#, k = 2, sfgrmin = 1e-8, sratmax = 0.999999, maxit = 30000)
all.plot<-plot_ordination(ps.newREL, all.ord.nmds.bray, type="samples", color="Location")#, shape="Location") 
all.plot + ggtitle("All sites Relative Abundance Transformation NMDS (k=4, sratmax, sfgrmin) 101602 taxa") +
  geom_point(size = 2.5)#+ geom_polygon(geom_point(size=5) + ggtitle("samples"))
ggsave("all_sites_nmds_reltrans_101602_taxa.pdf", plot = last_plot(), device = pdf, height = 6, width = 8)

# plot sequening depth per plate
sdt = data.table(as(sample_data(ps.new), "data.frame"),
                 TotalReads = sample_sums(ps.new), keep.rownames = TRUE)
setnames(sdt, "rn", "SampleID")
pSeqDepth = ggplot(sdt, aes(TotalReads)) + geom_histogram(binwidth = 1000) + ggtitle("Sequencing Depth")
pSeqDepth + facet_wrap(Location~ITS_seq_plate) + geom_vline(xintercept =2000)

#plot richness accumulation with sequencing
log(estimate_richness(ps.new, measures = "Observed")) -> sdt$richness
pRichAccum<-ggplot(sdt, aes(x = TotalReads, y = richness)) + 
  geom_point(aes(colour = factor(ITS_seq_plate)))
pRichAccum

# Undo the commenet below to clean up things saved in memory
#rm(sdt, pSeqDepth, pRichAccum)

# plot OTUs lost after filtering at different count totals
tdt = data.table(tax_table(ps.new),
                 TotalCounts = taxa_sums(ps.new),
                 OTU = taxa_names(ps.new))
ggplot(tdt, aes(TotalCounts)) + 
  geom_histogram(binwidth = 1000) + 
  ggtitle("Histogram of Total Counts")

# taxa cumulative sum
taxcumsum = tdt[, .N, by = TotalCounts]
setkey(taxcumsum, TotalCounts)
taxcumsum[, CumSum := cumsum(N)]
# Define the plot
pCumSum = ggplot(taxcumsum, aes(TotalCounts, CumSum)) + 
  geom_point() +
  xlab("Filtering Threshold, Minimum Total Counts") +
  ylab("OTUs Filtered") +
  ggtitle("OTUs that would be filtered vs. the minimum count threshold")
pCumSum

pCumSum + xlim(0, 3000)

# Undo the commenet below to clean up things saved in memory
#rm(tdt, taxcumsum, pCumSum)

# filter low abundance sequence variants #### 
# This filters the entire dataset. Ignore unless you want to run analyses for all sites
psF.new <- filter_taxa(ps.new, function(x) sum(x) > 10, TRUE)
ntaxa(ps.new)

# Make phyloseq object for taxa shared with defined sites ####
universe <- unique(c(bv.names, sv.names, ur.names))

bvtaxa.l <- universe %in% bv.names
svtaxa.l <- universe %in% sv.names
urtaxa.l <- universe %in% ur.names

venntest <- function(x) (x %in% bv.names) & (x %in% sv.names) & (x %in% ur.names)
php_overlap <- universe[venntest(universe)]

# make sure you have the right phyloseq object here 
psF.shared <- subset_taxa(psF, rownames(tax_table(psF)) %in% php_overlap)
psF.sharedR <- transform_sample_counts(psF.shared, function(x) x / sum(x))
ntaxa(psF.sharedR)

# FSP specific ####
# Should be run from unfiltered data
# Make FSP-RR-Soy only phyloseq object from unfiltered original 
fsp.RRsoy <- subset_samples(ps.new, Location == "Beltsville" & genotype == "RR" & crop == "soy") #  & Sampling_date == "pre" & crop == "soy"  & Glyphosphate_Treatment == "no_spray")
fsp.RRsoy <- prune_taxa(taxa_sums(fsp.RRsoy) > 10, fsp.RRsoy)
ntaxa(fsp.RRsoy)
nsamples(fsp.RRsoy)

# DESeq routines

fspDS <- phyloseq_to_deseq2(fsp.RRsoy, ~ Glyphosphate_Treatment)

fspDS$group <- factor(paste(fspDS$System.loc, fspDS$Glyphosphate_Treatment))

levels(fspDS$group)<-sub(" ", "", levels(fspDS$group))

design(fspDS) <- formula(~ group + Sampling_date + Soil_Zone + group:Sampling_date)

fspDS = estimateSizeFactors(fspDS)
fspDS = estimateDispersions(fspDS)
fspVST = getVarianceStabilizedData(fspDS)

# This is a different type of transformation. I don't think it is appropriate for the fungi.
#fspRLOG <- rlog(fspDS)
#plotPCA(fspRLOG, intgroup = "group")
#fspRLOGmatrix <- assay(fspRLOG)

fspDS <- DESeq(fspDS, test="Wald", fitType="parametric", parallel = T)
resFSP<-results(fspDS)

resFSPgroup <- results(fspDS, contrast = c("group", "Org_6spray", "Org_6no_spray"))
plotDispEsts(fspDS)
#fspRLD <- rlog(fspDS)

fspVST2 <- vst(fspDS, blind = T)

plotPCA(fspVST2, intgroup = "group")
relibrary("pheatmap")
select <- order(rowMeans(counts(fspDS,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(fspDS)[,c("group","Sampling_date")])

#dd[with(dd, order(-z, b)), ]
select <- order(rowMeans(counts(fspDS,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(fspDS)[,c("genotype", "System.loc")])
newdf<-df
row.names(newdf)->newdfr
population[order(population$age),]
newdf<-df[order(df$System.loc, df$genotype),]
fspheat<-assay(fspVST2)
fspheat<-fspheat[,newdfr]

pheatmap(fspheat[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=F, annotation_col=newdf)

# END of DESeq routines

# FSP Soy Ordination. Interesting, but the DESeq PCoA seem to be prefereable. 
#You could always run a PCoA instead of NMDS here...  You do you.

fsp.RRsoy0 <- fsp.RRsoy

otu_table(fsp.RRsoy) <- otu_table(fspVST, taxa_are_rows = TRUE)
fspsoyREL  = transform_sample_counts(fsp.RRsoy0, function(x) x / sum(x) )
rm(fspsoyRELfilt)
sample_data(fspsoyREL)<-phpmeta
set.seed(1977)
fsp.ord.nmds.bray <- ordinate(fspsoyREL, k = 2, method="NMDS", distance="bray", maxit = 30000)#, sratmax = 0.9999999)#, sfgrmin = 1e-8, sratmax = 0.9999999)#, k = 2, sfgrmin = 1e-8, sratmax = 0.999999, maxit = 30000)
fsp.plot<-plot_ordination(fspsoyREL, fsp.ord.nmds.bray, type="samples", color="System.loc", shape="year") 
fsp.plot + ggtitle("FSP-Soy Relative transform NMDS (k=2) 16092 taxa") +
  geom_point(size = 3)#+ geom_polygon(geom_point(size=5) + ggtitle("samples"))
ggsave("fsp_nmds_reltrans_16092_taxa.pdf", plot = last_plot(), device = pdf, height = 6, width = 8)

#rlog transformed (with model)

fsp.RRsoy0 <- fsp.RRsoy
fspRLOGmatrix0 <- fspRLOGmatrix
fspRLOGmatrix[fspRLOGmatrix < 0.0] <- 0.0
otu_table(fsp.RRsoy) <- otu_table(fspRLOGmatrix, taxa_are_rows = TRUE)
library(vegan)
set.seed(1977)
metadata <- as(sample_data(fsp.RRsoy), "data.frame")
fsp.dist <- phyloseq::distance(fsp.RRsoy, "bray")
adonis(fsp.dist ~ Glyphosphate_Treatment + 
         year + 
         System.loc + 
         System.loc:year + 
         Sampling_date + 
         Soil_Zone + 
         year:Sampling_date + 
         Glyphosphate_Treatment:Sampling_date +
         Glyphosphate_Treatment:Sampling_date:System.loc +
         Glyphosphate_Treatment:Soil_Zone + 
         Glyphosphate_Treatment:Sampling_date:System.loc:year, 
       as(sample_data(fsp.RRsoy), "data.frame"))

fspsoyREL  = transform_sample_counts(fsp.RRsoy0, function(x) x / sum(x) )
rm(fspsoyRELfilt)
sample_data(fspsoyREL)<-phpmeta
set.seed(1977)
fsp.ord.nmds.bray <- ordinate(fspsoyREL, k = 2, method="NMDS", distance="bray", maxit = 30000)#, sratmax = 0.9999999)#, sfgrmin = 1e-8, sratmax = 0.9999999)#, k = 2, sfgrmin = 1e-8, sratmax = 0.999999, maxit = 30000)
fsp.plot<-plot_ordination(fspsoyREL, fsp.ord.nmds.bray, type="samples", color="System.loc", shape="year") 
fsp.plot + ggtitle("FSP-Soy Relative transform NMDS (k=2) 16092 taxa") +
  geom_point(size = 3)#+ geom_polygon(geom_point(size=5) + ggtitle("samples"))
ggsave("fsp_nmds_reltrans_16092_taxa.pdf", plot = last_plot(), device = pdf, height = 6, width = 8)

# FSP Soy PERMANOVA
fspVST[fspVST < 0.0] <- 0.0
otu_table(fsp.RRsoy) <- otu_table(fspVST, taxa_are_rows = TRUE)
library(vegan)
set.seed(1977)
metadata <- as(sample_data(fsp.RRsoy), "data.frame")
fsp.dist <- phyloseq::distance(fsp.RRsoy, "bray")
adonis(fsp.dist ~ System.loc + 
         year + 
         Glyphosphate_Treatment + 
         Sampling_date + 
         Soil_Zone + 
         Glyphosphate_Treatment + 
         Glyphosphate_Treatment:Soil_Zone + 
         Glyphosphate_Treatment:Sampling_date +
         Glyphosphate_Treatment:Sampling_date:System.loc:year, 
       as(sample_data(fsp.RRsoy), "data.frame"))


# FSP-RR-Corn ####
# phyloseq object from unfiltered original
fsp.RRcorn <- subset_samples(ps.new, Location == "Beltsville" & genotype == "RR" & crop == "corn") #  & Sampling_date == "pre" & crop == "corn"  & Glyphosphate_Treatment == "no_spray")
fsp.RRcorn <- prune_taxa(taxa_sums(fsp.RRcorn) > 10, fsp.RRcorn)
ntaxa(fsp.RRcorn)
nsamples(fsp.RRcorn)

# DESeq routines
# inherits phyloseq object from FSP-only 

fspDS <- phyloseq_to_deseq2(fsp.RRcorn, ~ Glyphosphate_Treatment)

fspDS$group <- factor(paste(fspDS$System.loc, fspDS$Glyphosphate_Treatment))

levels(fspDS$group)<-sub(" ", "", levels(fspDS$group))

design(fspDS) <- formula(~ group + Sampling_date + Soil_Zone + group:Sampling_date)

fspDS = estimateSizeFactors(fspDS)
fspDS = estimateDispersions(fspDS)
fspVST = getVarianceStabilizedData(fspDS)

#dists <- dist(t(assay(fspVST2)))

fspDS <- DESeq(fspDS, test="Wald", fitType="parametric", parallel = T)
resFSP<-results(fspDS)

resFSPgroup <- results(fspDS, contrast = c("group", "Org_6spray", "Org_6no_spray"))
plotDispEsts(fspDS)
#fspRLD <- rlog(fspDS)

fspVST2 <- vst(fspDS, blind = T)
fspVST3 <- vst(fspDS)
plotPCA(fspVST3, intgroup = "group")
relibrary("pheatmap")
select <- order(rowMeans(counts(fspDS,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(fspDS)[,c("group","Sampling_date")])

#dd[with(dd, order(-z, b)), ]
df1 <- df[with(df, order(group,Sampling_date)), ]
fspASSAY <- assay(fspVST)[select,]

pheatmap(assay(fspVST)[select,], cluster_rows=FALSE, show_rownames=F,
         cluster_cols=F, annotation_col=df)

# END of DESeq routines

# FSP Ordination
fspcornREL  = transform_sample_counts(fsp.RRcorn0, function(x) x / sum(x) )
#rm(fspcornRELfilt)
#sample_data(fspcornREL)<-phpmeta
set.seed(1977)
fsp.ord.nmds.bray <- ordinate(fspcornREL, k = 2, method="NMDS", distance="bray", maxit = 30000)#, sratmax = 0.9999999)#, sfgrmin = 1e-8, sratmax = 0.9999999)#, k = 2, sfgrmin = 1e-8, sratmax = 0.999999, maxit = 30000)
fsp.plot<-plot_ordination(fspcornREL, fsp.ord.nmds.bray, type="samples", color="System.loc", shape="year") 
fsp.plot + ggtitle("FSP-corn Relative transform NMDS (k=2) 18094 taxa") +
  geom_point(size = 3)#+ geom_polygon(geom_point(size=5) + ggtitle("samples"))
ggsave("fsp_corn_nmds_reltrans_18094_taxa.pdf", plot = last_plot(), device = pdf, height = 6, width = 8)

# FSP Corn PERMANOVA
fspVST[fspVST < 0.0] <- 0.0
fsp.RRcorn0 <- fsp.RRcorn
otu_table(fsp.RRcorn) <- otu_table(fspVST, taxa_are_rows = TRUE)
library(vegan)
set.seed(1977)
metadata <- as(sample_data(fsp.RRcorn), "data.frame")
fsp.dist <- phyloseq::distance(fsp.RRcorn, "bray")
adonis(fsp.dist ~ System.loc + 
         year + 
         Glyphosphate_Treatment + 
         Sampling_date + 
         Soil_Zone + 
         Glyphosphate_Treatment + 
         Glyphosphate_Treatment:Soil_Zone + 
         Glyphosphate_Treatment:Sampling_date +
         Glyphosphate_Treatment:Sampling_date:System.loc:year, 
       as(sample_data(fsp.RRcorn), "data.frame"))


#stressplot(fsp.ord.nmds.bray)
ggsave("fsp_nmds_converged.pdf", plot = last_plot(), device = pdf)
 
#Summary graphs ####
# Transform with median sample count
total = median(sample_sums(fsp.RRsoy))
standf = function(x, t=total) round(t * (x / sum(x)))
fsp.trans = transform_sample_counts(fsp.RRsoy, standf)

title = "plot_bar; Beltsville vst Classes"
#plot_bar(psF.subset, x = "System.loc", y = "Abundance", "Family", title=title)
bvfilt <- plot_bar(fsp.trans, "System.loc", "Abundance", "Class", title=title, facet_grid = "year~")
bvfilt + geom_bar(aes(color=Class, fill=Class), stat="identity", position="stack")

ggsave("FSP_taxa_barchart_vst_class.pdf", plot = last_plot(), device = pdf, height = 20, width = 20)


# STUFF
fsp.tax = data.table(tax_table(fsp.slice),
                     TotalCounts = taxa_sums(fsp.slice),
                     OTU = taxa_names(fsp.slice))

sdt.slice = data.table(as(sample_data(fsp.slice), "data.frame"),
                       TotalReads = sample_sums(fsp.slice), keep.rownames = TRUE)
setnames(sdt.slice, "rn", "SampleID")
log(estimate_richness(fsp.slice, measures = "Observed")) -> sdt.slice$richness
pRichAccum<-ggplot(sdt.slice, aes(x = TotalReads, y = richness)) +
  geom_point(aes(colour = factor(Soil_Zone)))
pRichAccum

fsp.filt = filter_taxa(fsp.slice, function(x) sd(x)/mean(x) > 3.0, TRUE)
ntaxa(fsp.filt)

fsp.top = names(sort(taxa_sums(fsp.slice), decreasing = TRUE)[1:100])
ex1 = prune_taxa(fsp.top, fsp.slice)

# Urbana ####

#Make a special phyloseq object due to labeling
# This could definitely stand to be cleaned up or modified based on what your metadata looks like
# Previously I had removed things like nt_bahn etc.
phpmeta2<-read.table("Master_Label_Table_PHPv1.csv", header = T, na.strings = "", sep = ",")
phpmeta2[2]<-NULL
phpmeta3<-phpmeta2[,-1]
rownames(phpmeta3)<-phpmeta2[,1]
#phpmeta3$Loc_plot_ID<-as.character(phpmeta3$Loc_plot_ID)
#phpmeta3$ITS_seq_plate<-as.character(phpmeta3$ITS_seq_plate)
phpmeta3$year<-as.factor(phpmeta3$year)
seqtab.new<-readRDS("fungi_seqs2.rds")
phpfuntax.new<-readRDS(file = "fungi_taxa2.rds")
psU <- phyloseq(otu_table(seqtab.new, taxa_are_rows=FALSE), 
                sample_data(phpmeta3), 
                tax_table(phpfuntax.new))
ur.RRcorn <- subset_samples(psU, Location == "Urbana" & genotype == "RR" & crop == "corn") #  & Sampling_date == "pre" & crop == "corn"  & Glyphosphate_Treatment == "no_spray")
ur.RRcorn <- prune_taxa(taxa_sums(ur.RRcorn) > 10, ur.RRcorn)
ntaxa(ur.RRcorn)
nsamples(ur.RRcorn)

# DESeq routines
# inherits phyloseq object from ur-only 

urDS <- phyloseq_to_deseq2(ur.RRcorn, ~ Glyphosphate_Treatment)
levels(urDS$System.loc) <- sub("-","_", levels(urDS$System.loc))

urDS$group <- factor(paste(urDS$System.loc, urDS$Glyphosphate_Treatment))

levels(urDS$group)<-sub(" ", "", levels(urDS$group))

design(urDS) <- formula(~ group + Sampling_date + Soil_Zone + group:Sampling_date)

urDS = estimateSizeFactors(urDS)
urDS = estimateDispersions(urDS)
urVST = getVarianceStabilizedData(urDS)

#dists <- dist(t(assay(urVST2)))

urDS <- DESeq(urDS, test="Wald", fitType="parametric", parallel = T)
resur<-results(urDS)

resurgroup <- results(urDS, contrast = c("group", "Org_6spray", "Org_6no_spray"))
plotDispEsts(urDS)
#urRLD <- rlog(urDS)

urVST2 <- vst(urDS, blind = T)
urVST3 <- vst(urDS)
plotPCA(urVST2, intgroup = "group")
relibrary("pheatmap")
select <- order(rowMeans(counts(urDS,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(urDS)[,c("group","Sampling_date")])

#dd[with(dd, order(-z, b)), ]
df1 <- df[with(df, order(group,Sampling_date)), ]
urASSAY <- assay(urVST)[select,]

pheatmap(assay(urVST)[select,], cluster_rows=FALSE, show_rownames=F,
         cluster_cols=F, annotation_col=df)

# END of DESeq stuff

# ur Ordination
ur.RRcorn0 <- ur.RRcorn
urcornREL  = transform_sample_counts(ur.RRcorn0, function(x) x / sum(x) )
#rm(urcornRELfilt)
#sample_data(urcornREL)<-phpmeta
set.seed(1977)
ur.ord.nmds.bray <- ordinate(urcornREL, k = 2, method="NMDS", distance="bray", maxit = 30000)#, sratmax = 0.9999999)#, sfgrmin = 1e-8, sratmax = 0.9999999)#, k = 2, sfgrmin = 1e-8, sratmax = 0.999999, maxit = 30000)
ur.plot<-plot_ordination(urcornREL, ur.ord.nmds.bray, type="samples", color="System.loc", shape="year") 
ur.plot + ggtitle("ur-corn Relative transform NMDS (k=2) 9150 taxa") +
  geom_point(size = 3)#+ geom_polygon(geom_point(size=5) + ggtitle("samples"))
ggsave("ur_corn_nmds_reltrans_9150_taxa.pdf", plot = last_plot(), device = pdf, height = 6, width = 8)

# ur Corn PERMANOVA
urVST[urVST < 0.0] <- 0.0
otu_table(ur.RRcorn) <- otu_table(urVST, taxa_are_rows = TRUE)
library(vegan)
set.seed(1977)
metadata <- as(sample_data(ur.RRcorn), "data.frame")
ur.dist <- phyloseq::distance(ur.RRcorn, "bray")
adonis(ur.dist ~ System.loc + 
         year + 
         Glyphosphate_Treatment + 
         Sampling_date + 
         Soil_Zone + 
         Glyphosphate_Treatment + 
         Glyphosphate_Treatment:Soil_Zone + 
         Glyphosphate_Treatment:year + 
         Glyphosphate_Treatment:Sampling_date +
         Glyphosphate_Treatment:Sampling_date:System.loc:year, 
       as(sample_data(ur.RRcorn), "data.frame"))


ur.slice <- subset_samples(psU, Location == "Urbana") # & Sampling_date == "pre" & crop == "soy" & genotype == "Non_RR" & Glyphosphate_Treatment == "no_spray")
ur.slice <- prune_taxa(taxa_sums(ur.slice) > 0, ur.slice)
ntaxa(ur.slice)
ur.slice <- filter_taxa(ur.slice, function(x) sum(x > 3) > (0.2*length(x)), TRUE)
ntaxa(ur.slice)

# Ordination
ur.ord.nmds.bray <- ordinate(ur.slice, method="NMDS", distance="bray", maxit = 10000)
ur.plot<-plot_ordination(ur.slice, ur.ord.nmds.bray, type="samples", color="year", shape="ITS_seq_plate") 
ur.plot + ggtitle("Urbana") #+ geom_polygon(geom_point(size=5) + ggtitle("samples"))

total = median(sample_sums(ur.slice))
standf = function(x, t=total) round(t * (x / sum(x)))
ur.trans = transform_sample_counts(ur.slice, standf)

# Diversity barchart
title = "plot_bar; Urbana taxa filtered"
#plot_bar(psF.subset, x = "System.loc", y = "Abundance", "Family", title=title)
ur.barplot <- plot_bar(ur.trans, "System.loc", "Abundance", "Order", title=title, facet_grid = "~crop")
ur.barplot + geom_bar(aes(color=Order, fill=Order), stat="identity", position="stack")

ggsave("Urbana_taxa_diversity.pdf", plot = last_plot(), device = pdf, height = 10, width = 20)

ur.meta <- data.table(as(sample_data(ur.slice), "data.frame"),
                      TotalReads = sample_sums(ps), keep.rownames = TRUE)
setnames(ur.meta, "rn", "SampleID")

ur.tax = data.table(tax_table(ur.slice),
                    TotalCounts = taxa_sums(ur.slice),
                    OTU = taxa_names(ur.slice))
unique(ur.tax[, list(Order)])->ur.orders
length(ur.orders[[1]])

intersect(list[[1]], ur.orders[[1]])

# Stoneville ####

phpmeta2<-read.table("Master_Label_Table_PHPv1.csv", header = T, na.strings = "", sep = ",")
phpmeta2[2]<-NULL
phpmeta3<-phpmeta2[,-1]
rownames(phpmeta3)<-phpmeta2[,1]
phpmeta3$Loc_plot_ID<-as.character(phpmeta3$Loc_plot_ID)
phpmeta3$ITS_seq_plate<-as.character(phpmeta3$ITS_seq_plate)

psU <- phyloseq(otu_table(seqtab.new, taxa_are_rows=FALSE), 
                sample_data(phpmeta3), 
                tax_table(phpfuntax.new))

stone.RRcorn <- subset_samples(psU, Location == "Stoneville" & genotype == "RR" & crop == "corn") #  & Sampling_date == "pre" & crop == "corn"  & Glyphosphate_Treatment == "no_spray")
stone.RRcorn <- prune_taxa(taxa_sums(stone.RRcorn) > 10, stone.RRcorn)
ntaxa(stone.RRcorn)
nsamples(stone.RRcorn)

# DESeq routines
# inherits phyloseq object from stone-only 

stoneDS <- phyloseq_to_deseq2(stone.RRcorn, ~ Glyphosphate_Treatment)
levels(stoneDS$System.loc) <- sub("-","_", levels(stoneDS$System.loc))

stoneDS$group <- factor(paste(stoneDS$System.loc, stoneDS$Glyphosphate_Treatment))

levels(stoneDS$group)<-sub(" ", "", levels(stoneDS$group))

design(stoneDS) <- formula(~ group + Sampling_date + Soil_Zone + group:Sampling_date)

stoneDS = estimateSizeFactors(stoneDS)
stoneDS = estimateDispersions(stoneDS)
stoneVST = getVarianceStabilizedData(stoneDS)

#dists <- dist(t(assay(stoneVST2)))

stoneDS <- DESeq(stoneDS, test="Wald", fitType="parametric", parallel = T)
resstone<-results(stoneDS)

resstonegroup <- results(stoneDS, contrast = c("group", "Org_6spray", "Org_6no_spray"))
plotDispEsts(stoneDS)
#stoneRLD <- rlog(stoneDS)

stoneVST2 <- vst(stoneDS, blind = T)
stoneVST3 <- vst(stoneDS)
plotPCA(stoneVST2, intgroup = "group")
relibrary("pheatmap")
select <- order(rowMeans(counts(stoneDS,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(stoneDS)[,c("group","Sampling_date")])

#dd[with(dd, order(-z, b)), ]
df1 <- df[with(df, order(group,Sampling_date)), ]
stoneASSAY <- assay(stoneVST)[select,]

pheatmap(assay(stoneVST)[select,], cluster_rows=FALSE, show_rownames=F,
         cluster_cols=F, annotation_col=df)

# END of DESeq stuff

# stone Ordination
stone.RRcorn0 <- stone.RRcorn
stonecornREL  = transform_sample_counts(stone.RRcorn0, function(x) x / sum(x) )
#rm(stonecornRELfilt)
#sample_data(stonecornREL)<-phpmeta
set.seed(1977)
stone.ord.nmds.bray <- ordinate(stonecornREL, k = 2, method="NMDS", distance="bray", maxit = 30000)#, sratmax = 0.9999999)#, sfgrmin = 1e-8, sratmax = 0.9999999)#, k = 2, sfgrmin = 1e-8, sratmax = 0.999999, maxit = 30000)
stone.plot<-plot_ordination(stonecornREL, stone.ord.nmds.bray, type="samples", color="System.loc", shape="year") 
stone.plot + ggtitle("stone-corn Relative transform NMDS (k=2) 11355 taxa") +
  geom_point(size = 3)#+ geom_polygon(geom_point(size=5) + ggtitle("samples"))
ggsave("stone_corn_nmds_reltrans_11355_taxa.pdf", plot = last_plot(), device = pdf, height = 6, width = 8)

# stone Corn PERMANOVA
stoneVST[stoneVST < 0.0] <- 0.0
otu_table(stone.RRcorn) <- otu_table(stoneVST, taxa_are_rows = TRUE)
library(vegan)
set.seed(1977)
metadata <- as(sample_data(stone.RRcorn), "data.frame")
stone.dist <- phyloseq::distance(stone.RRcorn, "bray")
adonis(stone.dist ~ System.loc + 
         year + 
         Glyphosphate_Treatment + 
         Sampling_date + 
         Soil_Zone + 
         Glyphosphate_Treatment + 
         Glyphosphate_Treatment:Soil_Zone + 
         Glyphosphate_Treatment:year + 
         Glyphosphate_Treatment:Sampling_date +
         Glyphosphate_Treatment:Sampling_date:System.loc:year, 
       as(sample_data(stone.RRcorn), "data.frame"))

# Stoneville Soy

phpmeta2<-read.table("Master_Label_Table_PHPv1.csv", header = T, na.strings = "", sep = ",")
phpmeta2[2]<-NULL
phpmeta3<-phpmeta2[,-1]
rownames(phpmeta3)<-phpmeta2[,1]
phpmeta3$Loc_plot_ID<-as.character(phpmeta3$Loc_plot_ID)
phpmeta3$ITS_seq_plate<-as.character(phpmeta3$ITS_seq_plate)

psU <- phyloseq(otu_table(seqtab.new, taxa_are_rows=FALSE), 
                sample_data(phpmeta3), 
                tax_table(phpfuntax.new))

stone.RRsoy <- subset_samples(psU, Location == "Stoneville" & genotype == "RR" & crop == "soy") #  & Sampling_date == "pre" & crop == "soy"  & Glyphosphate_Treatment == "no_spray")
stone.RRsoy <- prune_taxa(taxa_sums(stone.RRsoy) > 10, stone.RRsoy)
ntaxa(stone.RRsoy)
nsamples(stone.RRsoy)

# DESeq routines
# inherits phyloseq object from stone-only 

stoneDS <- phyloseq_to_deseq2(stone.RRsoy, ~ Glyphosphate_Treatment)
levels(stoneDS$System.loc) <- sub("-","_", levels(stoneDS$System.loc))

stoneDS$group <- factor(paste(stoneDS$System.loc, stoneDS$Glyphosphate_Treatment))

levels(stoneDS$group)<-sub(" ", "", levels(stoneDS$group))

design(stoneDS) <- formula(~ group + Sampling_date + Soil_Zone + group:Sampling_date)

stoneDS = estimateSizeFactors(stoneDS)
stoneDS = estimateDispersions(stoneDS)
stoneVST = getVarianceStabilizedData(stoneDS)

#dists <- dist(t(assay(stoneVST2)))

stoneDS <- DESeq(stoneDS, test="Wald", fitType="parametric", parallel = T)
resstone<-results(stoneDS)

resstonegroup <- results(stoneDS, contrast = c("group", "Org_6spray", "Org_6no_spray"))
plotDispEsts(stoneDS)
#stoneRLD <- rlog(stoneDS)

stoneVST2 <- vst(stoneDS, blind = T)
stoneVST3 <- vst(stoneDS)
plotPCA(stoneVST2, intgroup = "group")
relibrary("pheatmap")
select <- order(rowMeans(counts(stoneDS,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(stoneDS)[,c("group","Sampling_date")])

#dd[with(dd, order(-z, b)), ]
df1 <- df[with(df, order(group,Sampling_date)), ]
stoneASSAY <- assay(stoneVST)[select,]

pheatmap(assay(stoneVST)[select,], cluster_rows=FALSE, show_rownames=F,
         cluster_cols=F, annotation_col=df)

# END of DESeq stuff

# stone Ordination
stone.RRsoy0 <- stone.RRsoy
stonesoyREL  = transform_sample_counts(stone.RRsoy0, function(x) x / sum(x) )
#rm(stonesoyRELfilt)
#sample_data(stonesoyREL)<-phpmeta
set.seed(1977)
stone.ord.nmds.bray <- ordinate(stonesoyREL, k = 4, method="NMDS", distance="bray", maxit = 30000, sratmax = 0.9999999)#, sfgrmin = 1e-8, sratmax = 0.9999999)#, k = 2, sfgrmin = 1e-8, sratmax = 0.999999, maxit = 30000)
stone.plot<-plot_ordination(stonesoyREL, stone.ord.nmds.bray, type="samples", color="System.loc", shape="year") 
stone.plot + ggtitle("stone-soy Relative transform NMDS (k=4) 11590 taxa") +
  geom_point(size = 3)#+ geom_polygon(geom_point(size=5) + ggtitle("samples"))
ggsave("stone_soy_nmds_reltrans_11590_taxa.pdf", plot = last_plot(), device = pdf, height = 6, width = 8)

# stone soy PERMANOVA
stoneVST[stoneVST < 0.0] <- 0.0
otu_table(stone.RRsoy) <- otu_table(stoneVST, taxa_are_rows = TRUE)
library(vegan)
set.seed(1977)
metadata <- as(sample_data(stone.RRsoy), "data.frame")
stone.dist <- phyloseq::distance(stone.RRsoy, "bray")
adonis(stone.dist ~ System.loc + 
         year + 
         Glyphosphate_Treatment + 
         Sampling_date + 
         Soil_Zone + 
         Glyphosphate_Treatment + 
         Glyphosphate_Treatment:Soil_Zone + 
         Glyphosphate_Treatment:year + 
         Glyphosphate_Treatment:Sampling_date +
         Glyphosphate_Treatment:Sampling_date:System.loc:year, 
       as(sample_data(stone.RRsoy), "data.frame"))

stone.slice <- subset_samples(psU, Location == "Stoneville") # & Sampling_date == "pre" & crop == "soy" & genotype == "Non_RR" & Glyphosphate_Treatment == "no_spray")
stone.slice <- prune_taxa(taxa_sums(stone.slice) > 0, stone.slice)
ntaxa(stone.slice)
stone.slice <- filter_taxa(stone.slice, function(x) sum(x > 3) > (0.2*length(x)), TRUE)
ntaxa(stone.slice)

stone.ord.nmds.bray <- ordinate(stone.slice, method="NMDS", distance="bray", maxit = 10000)
stone.plot<-plot_ordination(stone.slice, stone.ord.nmds.bray, type="samples", label = "Sample_ID", color="System.loc", shape="crop") 
stone.plot + ggtitle("Stoneville: notice those errors!") +
  xlim(0.25, 0.5) +
  ylim(0.25, 0.5)#+ geom_polygon(geom_point(size=5) + ggtitle("samples"))

ggsave("Stoneville_ordination.pdf", plot = last_plot(), device = pdf, height = 10, width = 20)
# Diversity barchart
total = median(sample_sums(stone.slice))
standf = function(x, t=total) round(t * (x / sum(x)))
stone.trans = transform_sample_counts(stone.slice, standf)

title = "plot_bar; Stoneville taxa filtered"
#plot_bar(psF.subset, x = "System.loc", y = "Abundance", "Family", title=title)
stone.barplot <- plot_bar(stone.trans, "System.loc", "Abundance", "Order", title=title, facet_grid = "~crop")
stone.barplot + geom_bar(aes(color=Order, fill=Order), stat="identity", position="stack")

ggsave("Stoneville_taxa_diversity.pdf", plot = last_plot(), device = pdf, height = 10, width = 20)

ur.meta <- data.table(as(sample_data(ur.slice), "data.frame"),
                      TotalReads = sample_sums(ps), keep.rownames = TRUE)
setnames(ur.meta, "rn", "SampleID")

ur.tax = data.table(tax_table(ur.slice),
                    TotalCounts = taxa_sums(ur.slice),
                    OTU = taxa_names(ur.slice))
unique(ur.tax[, list(Order)])->ur.orders
length(ur.orders[[1]])

intersect(list[[1]], ur.orders[[1]])
#### ####
total = median(sample_sums(psf))
standf = function(x, t=total) round(t * (x / sum(x)))
psfs = transform_sample_counts(psf, standf)
phpF = filter_taxa(psfs, function(x) sd(x)/mean(x) > 3.0, TRUE)
bvsfh <- subset_taxa(bvsf, Order == "o__Hypocreales")

shyp <- subset_taxa(sville, Order == "o__Hypocreales")
shyp1 <-  prune_samples(sample_sums(shyp)>=20, shyp)



####END####
plot_richness(ps, x="System.loc", measures=c("Shannon", "Simpson"), color="Soil_Zone") + theme_bw()

ord.nmds.bray <- ordinate(ps, method="NMDS", distance="bray")

top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
plot_bar(ps, x="crop", fill="Family") + facet_wrap(~Location, scales="free_x")
