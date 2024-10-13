######## LRS and SRS only coverage regions in a karyoplot ############ 
## ---------------------------------------------------------------------------------------------------------------------------- ##
## ------------------------------------------------For T2T-CHM13v2.0 Karyotype------------------------------------------------- ##
## ---------------------------------------------------------------------------------------------------------------------------- ##

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("karyoploteR")

library(karyoploteR)
library(GenomicRanges)
library(ggplot2)
library(dplyr)
library(biomaRt)
library(regioneR)
library("data.table")

# Load regions (after removing simulated centromeric regions, and only keep regions on min 1kb length. Only plots for one sample)
inlrs_notsrs <- read.table(file = "inlrs_notsrs.callable.intersectBED.bed", sep = "\t", header = FALSE) %>% 
  filter(V5>=1000)
lrs_only <- GRanges(seqnames = inlrs_notsrs$V1, ranges = IRanges(start = inlrs_notsrs$V2, end = inlrs_notsrs$V3))
insrs_notlrs <- read.table(file = "insrs_notlrs.callable.intersectBED.bed", sep = "\t", header = FALSE) %>%
  filter(V5>=1000)
srs_only <- GRanges(seqnames = insrs_notlrs$V1, ranges = IRanges(start = insrs_notlrs$V2, end = insrs_notlrs$V3))

# Create karyotype with regions
CYTO_chm13 <- toGRanges(fread("../databases/chm13v2.0_cytobands_allchrs.bed", col.names = c("chr","start","end", "name","gieStain")))
kp <- plotKaryotype(genome="hs1", plot.type=2, cytobands = CYTO_chm13, labels.plotter = NULL);
kpAddChromosomeNames(kp, srt=0, cex=1)
kpDataBackground(kp, r0=0, r1=0.3, color = "white")
kpDataBackground(kp, r0=0.4, r1=0.7, color = "white");

kpPlotRegions(kp, data=lrs_only, col="#3182bd", border="#3182bd", layer.margin = 0.01, r0=0.1, r1=0.6, avoid.overlapping=FALSE, data.panel=1);
kpPlotRegions(kp, data=extendRegions(srs_only), col="#c51b8a", border="#c51b8a", layer.margin = 0.01, r0=.1, r1=.6, avoid.overlapping=FALSE, data.panel=2);

# Legend
legend(x = "bottomright", fill = c("#c51b8a", "#3182bd"),
       legend = c("Regions > 1kb callable by SRS,\nbut not by LRS", 
                  "Regions > 1kb callable by LRS,\nbut not by SRS"),
       cex=0.8,
       box.lty = 0)

ggsave(filename = "Figure2B_CHM13.pdf", width = 9, height = 6, units = "in", dpi=300)


## ---------------------------------------------------------------------------------------------------------------------------- ##
## ------------------------------------------------For GRCh38 Karyotype------------------------------------------------- ##
## ---------------------------------------------------------------------------------------------------------------------------- ##

# Load regions (after removing simulated centromeric regions, and only keep regions on min 1kb length)
inlrs_notsrs <- read.table(file = "inlrs_notsrs.callable.intersectBED.noCent.1kb.bed", sep = "\t", header = FALSE) %>% 
  filter(V5>=1000)
lrs_only <- GRanges(seqnames = inlrs_notsrs$V1, ranges = IRanges(start = inlrs_notsrs$V2, end = inlrs_notsrs$V3))
insrs_notlrs <- read.table(file = "insrs_notlrs.callable.intersectBED.noCent.1kb.bed", sep = "\t", header = FALSE) %>%
  filter(V5>=1000)
srs_only <- GRanges(seqnames = insrs_notlrs$V1, ranges = IRanges(start = insrs_notlrs$V2, end = insrs_notlrs$V3))

# Create karyotype with regions
kp <- plotKaryotype(genome="hg38", plot.type=2, labels.plotter = NULL);
kpAddChromosomeNames(kp, srt=0, cex=1)
kpDataBackground(kp, r0=0, r1=0.3, color = "white")
kpDataBackground(kp, r0=0.4, r1=0.7, color = "white");

kpPlotRegions(kp, data=lrs_only, col="#3182bd", border="#3182bd", layer.margin = 0.01, r0=0.1, r1=0.6, avoid.overlapping=FALSE, data.panel=1);
kpPlotRegions(kp, data=extendRegions(srs_only), col="#c51b8a", border="#c51b8a", layer.margin = 0.01, r0=.1, r1=.6, avoid.overlapping=FALSE, data.panel=2);

# Legend
legend(x = "right", fill = c("#c51b8a", "#3182bd"),
       legend = c("Regions > 1kb well-covered by SRS,\nbut not by LRS", 
                  "Regions > 1kb well-covered by LRS,\nbut not by SRS"),
       cex=0.8,
       box.lty = 0)

ggsave(filename = "Figure2B_GRCh38.pdf", width = 8, height = 6, units = "in", dpi=300)


## ---------------------------------------------------------------------------------------------------------------------------- ##
# Annotate Mendelian disease genes well-covered by LRS-only (ZOOM IN)
gene.symbols <- c("C4A","C4B","CBS","CRYAA","ICOSLG","PRODH","SMN1","SMN2")
#ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", mirror = 'www')
genes <- toGRanges(getBM(attributes=c('chromosome_name', 'start_position', 'end_position', 'hgnc_symbol'),
                         filters = 'hgnc_symbol', values =gene.symbols, mart = ensembl))
seqlevelsStyle(genes) <- "UCSC"
head(genes)
kpPlotMarkers(kp, data=genes, labels=genes$hgnc_symbol, text.orientation = "horizontal", avoid.overlapping=FALSE,
              r1=0.5, cex=0.8)

# List of chromosomes
chromosomes <- c("chr5", "chr6", "chr21", "chr22")
# Create a plotKaryotype object
kp <- plotKaryotype(genome="hg38", chromosomes=chromosomes, plot.type=2, labels.plotter = NULL)
kpAddChromosomeNames(kp, srt=0, cex=1)
kpDataBackground(kp, r0=0.1, r1=0.4, color = "white")
kpDataBackground(kp, r0=0.1, r1=0.4, color = "white")

# Create a plot with multiple panels for each chromosome
kpPlotRegions(kp, data=lrs_only[lrs_only@seqnames %in% chromosomes], col="#3182bd", border="#3182bd", layer.margin = 0.01, r0=0, r1=0.4, avoid.overlapping=FALSE, data.panel = 1)
kpPlotRegions(kp, data=extendRegions(srs_only[srs_only@seqnames %in% chromosomes]), col="#c51b8a", border="#c51b8a", layer.margin = 0.01, r0=0, r1=.4, avoid.overlapping=FALSE, data.panel = 2)
kpPlotMarkers(kp, data=genes, labels=genes$hgnc_symbol, text.orientation = "horizontal",
              r1=0.9, cex=1, ignore.chromosome.ends=TRUE)
# Legend
legend(x = "bottomright", fill = c("#c51b8a", "#3182bd"), 
       legend = c("Regions > 1kb well-covered by SRS,\nbut not by LRS\n", 
                  "Regions > 1kb well-covered by LRS,\nbut not by SRS"),
       cex=0.9,
       box.lty = 0)
