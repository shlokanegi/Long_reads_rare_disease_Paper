suppressPackageStartupMessages(library(dplyr))
library(stringr)
library(ggplot2)
library(knitr)
library(tidyr)
suppressPackageStartupMessages(library(VariantAnnotation))
library(GenomicRanges)

### Prepare a GRanges object with general variant information
prepareGRanges <- function(vcf.o){  # vcf.o is a VCF object generated somewhere from readVCF()
  vars = rowRanges(vcf.o)   ## from VariantAnnotation package.
  vars$gt = geno(vcf.o)$GT[,1]
  vars$gq = geno(vcf.o)$GQ[,1]
  vars$AF = unlist(info(vcf.o)$AF)
  vars$AF = ifelse(is.na(vars$AF), 0, vars$AF)
  vars$CADD = unlist(info(vcf.o)$dbNSFP_CADD_phred)
  return(vars)
}

######## Make a ANNotation data.frame
## to shorten long allelic sequence (for SVs it can be thousands of bases long)
shorten <- function(x, max.length=30){
  x.l = nchar(x)
  long.x = which(x.l>max.length)
  if(length(long.x)>0){
    x[long.x] = paste0(substr(x[long.x], 1, max.length), '_', as.integer(x.l[long.x]-max.length))
  }
  return(x)
}
## convert the GRange object to a data.frame more efficiently by shortening the sequence of the alleles 
## and flattening the lists
grToDF <- function(gr){
  df = tibble(ref=shorten(as.character(gr$REF)),
              alt=shorten(as.character(unlist(gr$ALT))))
  df = cbind(df, mcols(gr)[,setdiff(colnames(mcols(gr)), c('REF','ALT', 'paramRangeID'))])
  mcols(gr) = NULL
  cbind(as.data.frame(gr)[,c('seqnames', 'start', 'end')], df)
}
## parse the ANN field and extract the different information (gene name, effect, impact, etc)
formatAnn <- function(df){
  ann.l = strsplit(df$ann, split='\\|')
  df$ann = NULL
  cbind(df, tibble(
    allele=unlist(lapply(ann.l, '[', 1)),
    effect=unlist(lapply(ann.l, '[', 2)),
    impact=unlist(lapply(ann.l, '[', 3)),
    gene=unlist(lapply(ann.l, '[', 4)),
    gene_type=unlist(lapply(ann.l, '[', 8))))
}
## use the functions above to create the data.frame from the GRanges and VCF objects
makeANNdataframe <-function(vars, vcf.o){
  ann = info(vcf.o)$ANN
  vars.ann = vars[rep(1:length(ann), unlist(lapply(ann, length)))]
  vars.ann$ann = unlist(ann)
  names(vars.ann) = NULL
  ann.df = grToDF(vars.ann)
  ann.df$impact = NULL
  unique(formatAnn(ann.df))
}

### Low mappability and low complexity BED files from GIAB GRCh38 stratification set
lowmap <- read.table("../databases/GRCh38_lowmappabilityall.bed", sep = "\t", header = FALSE)
lowmap <- GRanges(seqnames = lowmap$V1,
                  ranges = IRanges(start = lowmap$V2, end = lowmap$V3))
lowcomp <- read.table("../databases/GRCh38_AllTandemRepeats_Homopols_satellites.bed", sep = "\t", header = FALSE)
lowcomp <- GRanges(seqnames = lowcomp$V1,
                   ranges = IRanges(start = lowcomp$V2, end = lowcomp$V3))

annotateRepeats <- function(vcf.o, sd, gaps, homopol, map, lowmap, lowcomp){
  ## segmental duplications
  rowRanges(vcf.o)$sd = suppressWarnings(overlapsAny(vcf.o, sd))
  rowRanges(vcf.o)$sd99 = suppressWarnings(overlapsAny(vcf.o, subset(sd, fracMatch>=.99)))
  ## assembly gaps, centromere, telomeres
  dist.o = suppressWarnings(distanceToNearest(vcf.o, gaps)) %>% as.data.frame
  rowRanges(vcf.o)$gapd = Inf
  rowRanges(vcf.o)$gapd[dist.o$queryHits] = dist.o$distance
  ## homopolymers
  rowRanges(vcf.o)$homopolymer = suppressWarnings(overlapsAny(vcf.o, homopol))
  ## mappable regions (UMAP and GIAB_lowmap)
  rowRanges(vcf.o)$mappable = suppressWarnings(overlapsAny(vcf.o, map))
  rowRanges(vcf.o)$giab_lowmap = suppressWarnings(overlapsAny(vcf.o, lowmap))
  ## GIAB low-complexity (All TRs, homs, satellites)
  rowRanges(vcf.o)$giab_lowcomp = suppressWarnings(overlapsAny(vcf.o, lowcomp))
  ## return annotated vcf
  vcf.o
}

loadVCF <- function(path, label){
  ## download and read VCF
  system(str_glue('rm -rf lrs.variants.vcf.gz 2>&1'), intern = TRUE)
  system(str_glue('cp {path} lrs.variants.vcf.gz 2>&1'), intern = TRUE)
  vcf.o = readVcf("lrs.variants.vcf.gz")
  ## annotate repeats
  vcf.o = annotateRepeats(vcf.o, sd, gaps, homopol, map, lowmap, lowcomp)
  ## prepare data.frames
  vars = prepareGRanges(vcf.o)
  ann.df = makeANNdataframe(vars, vcf.o)
  # condense al.df to match vep's/srs TSV
  ann.df <- ann.df %>% mutate(label=label) %>%
    group_by(label, seqnames, start, end, ref, alt) %>%
    mutate(
      alltranscript.gene = toString(gene),
      alltranscript.impact = toString(impact),
      alltranscript.effect = toString(effect)
    ) %>% ungroup() %>% distinct(across(c("label", "seqnames", "start", "end", "ref", "alt")), .keep_all = TRUE)
  ## remove rownames
  ann.df = ann.df %>% mutate(locus=paste0(seqnames,":",start))
  rownames(ann.df) <- NULL
  
  return(ann.df)
}

############################################################################################################################
###########################################################################################################################################################
############################### LRS FAV analysis ############################

############################################################################################################################
############################################################################################################################

############### READ all VCFs
probands <- c('RGP_696_3', 'RGP_1608_3', 'RGP_1219_3', 'RGP_25_3', 'RGP_1219_4', 'RGP_1081_3', 'RGP_686_3', 'RGP_607_3', 'RGP_123_3', 'RGP_1360_3', 'RGP_1238_3', 'RGP_1395_3', 'RGP_1630_3', 'RGP_858_3', 'RGP_731_3', 'RGP_876_3', 'RGP_12_3', 'RGP_2040_3', 'RGP_558_3', 'RGP_1770_3', 'RGP_1316_3')

if(LOAD_CACHED_DATAFRAME){
  load('annotated.lrsfav.gq20.df.rdata')
} else {
  ann.df = lapply(probands, function(x){
    loadVCF(paste0('lrs_only_FAVs_gq20/', x, '.fp.vcf.gz'), x)
  }) %>% bind_rows
  save(ann.df, file='annotated.lrsfav.gq20.df.rdata')
}

## Without GQ filtering
LOAD_CACHED_DATAFRAME = TRUE
if(LOAD_CACHED_DATAFRAME){
  load('annotated.lrsfav.df.rdata')
} else {
  ann.df = lapply(probands, function(x){
    loadVCF(paste0('lrs_only_FAVs/', x, '.fp.vcf.gz'), x)
  }) %>% bind_rows
  save(ann.df, file='annotated.lrsfav.df.rdata')
}

## How many LRS-only FAVs are rare? How many of these overlap disease genes? What are number of hets and homs.
ann.df$label <- sample_to_family[ann.df$label]
ann.df <- ann.df %>% group_by(seqnames,start,end) %>% mutate(AC=n())
median((ann.df %>% group_by(label) %>% summarise(count=n()))$count)

ann.rare.df <- ann.df %>% filter(AF<0.001 & AC<=1)
median((ann.rare.df %>% group_by(label) %>% summarise(count=n()))$count)

## How many in segdups and low mappability regions
print(paste(median((ann.rare.df %>% filter(sd==TRUE) %>% group_by(label) %>% summarise(sd_favs=n()))$sd_favs), median((ann.rare.df %>% filter(sd==TRUE) %>% group_by(label) %>% summarise(count=n()))$count)*100/median((ann.rare.df %>% group_by(label) %>% summarise(count=n()))$count)))
print(paste(median((ann.rare.df %>% filter(giab_lowmap==TRUE) %>% group_by(label) %>% summarise(count=n()))$count), median((ann.rare.df %>% filter(giab_lowmap==TRUE) %>% group_by(label) %>% summarise(count=n()))$count)*100/median((ann.rare.df %>% group_by(label) %>% summarise(count=n()))$count)))
print(paste(median((ann.rare.df %>% filter(giab_lowcomp==TRUE) %>% group_by(label) %>% summarise(count=n()))$count), median((ann.rare.df %>% filter(giab_lowcomp==TRUE) %>% group_by(label) %>% summarise(count=n()))$count)*100/median((ann.rare.df %>% group_by(label) %>% summarise(count=n()))$count)))
print(paste(median((ann.rare.df %>% filter(homopolymer==TRUE) %>% group_by(label) %>% summarise(count=n()))$count), median((ann.rare.df %>% filter(homopolymer==TRUE) %>% group_by(label) %>% summarise(count=n()))$count)*100/median((ann.rare.df %>% group_by(label) %>% summarise(count=n()))$count)))

## How many of these were INDEL contributions
ann.rare.df <- ann.rare.df %>% mutate(type=ifelse((str_length(ref)==1 & str_length(alt)==1), "snv", "indel"))
median((ann.rare.df %>% filter(type=="indel") %>% group_by(label) %>% summarise(indels=n()))$indels)
median((ann.rare.df %>% filter(type=="snv") %>% group_by(label) %>% summarise(snvs=n()))$snvs)

x <- ann.rare.df %>% filter(gene %in% dis.genes) %>% mutate(type=ifelse((str_length(ref)==1 & str_length(alt)==1), "snv", "indel"))
median((x %>% filter(type=="indel") %>% group_by(label) %>% summarise(indels=n()))$indels)
median((x %>% filter(type=="snv") %>% group_by(label) %>% summarise(snvs=n()))$snvs)


## How many LRS-only genes are disease genes?
## disease genes
dis.genes = unique(c(clingen.genes$gene, gencc.genes$gene_symbol))
omim.genes <- read.csv(file="databases/omim_disease_genes.txt", sep="\t", header = TRUE)
omim.genes <- omim.genes %>% filter(Approved.Gene.Symbol!="")
dis.genes <- unique(c(dis.genes,omim.genes$Approved.Gene.Symbol))

## dominant disease genes
dominant.genes = unique(c(
  subset(omim.genes, grepl('AD', moi) | grepl('XLD', moi))$Approved.Gene.Symbol,
  subset(clingen.genes, inheritance == 'AD')$gene,
  subset(gencc.genes, grepl('dominant', moi_title))$gene_symbol))
## recessive disease genes
recessive.genes = unique(c(
  subset(omim.genes, grepl('AR', moi))$Approved.Gene.Symbol,
  subset(clingen.genes, inheritance == 'AR')$gene,
  subset(gencc.genes, grepl('recessive', moi_title))$gene_symbol))

lrs_only_genes <- unique(ann.rare.df$gene)
print(paste("No. of LRS-only genes in disease genes:", sum(lrs_only_genes %in% dis.genes)))
print(paste("No. of LRS-only genes in recessive genes:", sum(lrs_only_genes %in% recessive.genes)))

lrs_only_genes[which(lrs_only_genes %in% dis.genes)]

############### Barplot representing LRS-only FAVs, LRS-only-SRS-nocov, rare, and variants in disease genes ##############
# Make a single dataframe with counts of a 4 categories, per proband
lro_favs_ALL <- ann.df %>% group_by(label) %>% summarise(ALL=n())
lro_favs_RARE <- ann.df %>% group_by(label) %>% summarise(RARE=count(AC<=1 & AF<0.001))
lro_favs_OMIM <- ann.df %>% group_by(label) %>% filter(AF<0.001 & AC<=1) %>% summarise(OMIM=count(gene%in%dis.genes))

comb.df <- merge(merge(lro_favs_ALL, lro_favs_RARE, by="label", all=TRUE), lro_favs_OMIM, by="label", all=TRUE)
# Reshape the data for plotting
comb.df_long <- tidyr::pivot_longer(comb.df, cols = c(ALL, RARE, OMIM), names_to = "category", values_to = "value")
# Change the order of levels in the "category" variable
comb.df_long$category <- factor(comb.df_long$category, levels = c("ALL", "RARE", "OMIM"))
comb.df_long <- comb.df_long %>% mutate(plot=ifelse(category %in% c("ALL"), "p1", "p2"))

########## MEDIANS
comb.df_long %>% group_by(category) %>% summarise(n=median(value))

## Recaculate intersection
comb.df <- comb.df %>%
  mutate(ALL=ALL-RARE,
         RARE=RARE-OMIM)

# Reshape the data for plotting
comb.df_long <- tidyr::pivot_longer(comb.df, cols = c(ALL, RARE, OMIM), names_to = "category", values_to = "value")
# Change the order of levels in the "category" variable
comb.df_long$category <- factor(comb.df_long$category, levels = c("ALL", "RARE", "OMIM"))
comb.df_long <- comb.df_long %>% mutate(plot=ifelse(category %in% c("ALL"), "p1", "p2"))

### PLOT EVERYTHING
library(paletteer)
paletteer_d("MapPalettes::tealberry_pie")
colors <- c("#82A2BFFF", "#D4D2E5FF", "#A30248FF")
legend_labels <- c("ALL" = "Common in gnomAD or cohort",
                   "RARE" = "Rare",
                   "OMIM" = "Rare affecting disease genes")
# Reorder the labels for the y-axis (to match order with other figure in paper)
comb.df_long$label <- factor(comb.df_long$label, levels = rev(c("RGP_558_3", "RGP_2040_3", "RGP_731_3", "RGP_858_3",
                                                            "RGP_1630_3", "RGP_1219_4", "RGP_1395_3", "RGP_1238_3",
                                                            "RGP_1360_3", "RGP_123_3", "RGP_607_3", "RGP_12_3",
                                                            "RGP_686_3", "RGP_1081_3", "RGP_25_3", "RGP_1608_3",
                                                            "RGP_1219_3", "RGP_696_3", "RGP_1316_3", "RGP_1770_3",
                                                            "RGP_876_3")))
comb.df_long <- comb.df_long %>% arrange(label)

b <- ggplot(comb.df_long, aes(x = label, y = value, fill = category)) +
  coord_flip() +
  geom_bar(stat = "identity", position = "stack", width = 1, color="black", size=0.2) +
  geom_text(aes(label = value), position = position_stack(vjust = 0.5), hjust=-1.3, size = 3, color = "black") +  # Add text labels inside bars
  labs(x = "", y = "", fill = "Category") +
  scale_fill_manual(values = colors, labels = legend_labels) +
  #scale_y_continuous(limits = c(0, 4001), breaks = seq(0, 5000, 500)) +
  scale_y_continuous(limits = c(0, 1801), breaks = seq(0, 1900, 300)) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.direction = "horizontal", legend.text = element_text(size = 12),
        legend.title = element_blank(),
        axis.text.x = element_text(size = 10, hjust = 1, vjust = 1), 
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14), 
        plot.title = element_text(size = 16, hjust = 0.5),
        axis.line = element_line(), axis.ticks = element_line(),
        strip.text = element_blank()) +
  ggtitle("Additional functionally-annotated variant yield with long reads") + xlab("") + ylab("Number of LRS only FAVs")

ggsave(filename = "Figure3B.pdf", plot = b, width = 10, height = 6, units = "in", dpi=300)
