#!/software/apps/r/gcc/64/3.4.1/bin/Rscript

library(dada2); packageVersion("dada2")

args = commandArgs(trailingOnly=TRUE)
outdir <- args[1]
pathF <-args[2]
pathR <- args[3]
run <- args[4]

filtpathF <- file.path(outdir, "filtered_for") # make path for forward filtered output
filtpathR <- file.path(outdir, "filtered_rev") # make path for reverse filtered output 

fastqFs <- sort(list.files(pathF, pattern="fastq.gz"))
fastqRs <- sort(list.files(pathR, pattern="fastq.gz"))

if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")
# Filtering: THESE PARAMETERS ARENT OPTIMAL FOR ALL DATASETS
filterAndTrim(fwd=file.path(pathF, fastqFs), filt=file.path(filtpathF, fastqFs),
	      rev=file.path(pathR, fastqRs), filt.rev=file.path(filtpathR, fastqRs),
	      minLen=50, maxEE=2, truncQ=11, maxN=0, rm.phix=TRUE,
	      compress=TRUE, verbose=TRUE, multithread=TRUE)

#unlink(pathF, recursive = TRUE)
#unlink(pathR, recursive = TRUE)

filtFs <- list.files(filtpathF, pattern="fastq.gz", full.names = TRUE)
filtRs <- list.files(filtpathR, pattern="fastq.gz", full.names = TRUE)

sample.names <- sapply(strsplit(basename(filtFs), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
sample.namesR <- sapply(strsplit(basename(filtRs), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz

if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")

names(filtFs) <- sample.names
names(filtRs) <- sample.names

set.seed(100)
# Learn forward error rates
errF <- learnErrors(filtFs, nread=1e6, multithread=TRUE)
pdf(file = paste(outdir,"/",run,"_for_error.pdf", sep = ""))
plotErrors(errF, nominalQ=TRUE)
dev.off()
# Learn reverse error rates
errR <- learnErrors(filtRs, nread=1e6, multithread=TRUE)
pdf(file = paste(outdir,"/",run,"_rev_error.pdf", sep = ""))
plotErrors(errR, nominalQ=TRUE)
dev.off()
# Sample inference and merger of paired-end reads
mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names
for(sam in sample.names) {
	cat("Processing:", sam, "\n")
	derepF <- derepFastq(filtFs[[sam]])
	ddF <- dada(derepF, err=errF, multithread=TRUE)
	derepR <- derepFastq(filtRs[[sam]])
	ddR <- dada(derepR, err=errR, multithread=TRUE)
	merger <- mergePairs(ddF, derepF, ddR, derepR)
	mergers[[sam]] <- merger
}

rm(derepF); rm(derepR)
# Construct sequence table and remove chimeras
seqtab <- makeSequenceTable(mergers)

saveRDS(seqtab, file = paste(outdir,"/",run,"_seqtab.RDS", sep = "")) # CHANGE ME to where you want sequence table saved
