### Compare LRS and SRS SVs
library(sveval)
library(GenomicRanges)

# Check if the command-line argument is provided
if (length(commandArgs(trailingOnly = TRUE)) < 2) {
  stop("Usage: Rscript svevalOl.R lrs.vcf srs.vcf outfile.tsv")
}

chromosomes <- c(paste0("chr", 1:22), "chrX", "chrY")

# Get the files from the command-line argument
lr.vcf <- commandArgs(trailingOnly = TRUE)[1]
sr.vcf <- commandArgs(trailingOnly = TRUE)[2]
outfile.tsv <- commandArgs(trailingOnly = TRUE)[3]

## load SR SVs
sr = readSVvcf(sr.vcf)
## keep only autosomes and sex-chromosomes
sr <- subset(sr, seqnames(sr) %in% chromosomes)
sr = subset(sr, type!='BND' & type!='CPX')
## convert DUP into INS
sr$type = ifelse(sr$type == 'DUP', 'INS', sr$type)
## keep only PASS filters
sr = subset(sr, filter=='PASS')

## load LR SVs
lr = readSVvcf(lr.vcf)
## keep only autosomes and sex-chromosomes
lr <- subset(lr, seqnames(lr) %in% chromosomes)
lr = subset(lr, type!='BND')
## convert DUP:TANDEM into INS
lr$type = ifelse(lr$type == 'DUP:TANDEM', 'INS', lr$type)

## load simple repeats
simprep = read.table('simpleRepeat_GRCh38.bed.gz', as.is=TRUE)

## make into a genomic ranges object
simprep = GRanges(simprep[,1], IRanges(simprep[,2], simprep[,3]))

## "evaluation" using LR as truth set
sveval.o = svevalOl(sr, lr, min.ol=.1, method='coverage', max.ins.dist=100, simprep=simprep)

## summary table
print(sveval.o$eval)

write.table(sveval.o$eval, file = outfile.tsv, sep = "\t", row.names = FALSE)

## opens an app to visualize the FP/FN/TP and check if they look real or a matching problem from sveval
#explore_eval_svs(sveval.o)


