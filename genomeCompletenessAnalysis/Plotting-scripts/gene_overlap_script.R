
### All files generated using the bash script coveval_grch38.sh and gene_coverage.sh

# SV annotation database
load("../databases/figures-tables_sv_annotation_database.RData")

############################### USING GENCODE ANNOTATIONS ###########################
################################# (all protein-coding) ##############################
dis.genes = unique(c(clingen.genes$gene, gencc.genes$gene_symbol))
omim.genes <- read.csv(file="../databases/omim_disease_genes.txt", sep="\t", header = TRUE)
omim.genes <- new.omim.genes %>% filter(Approved.Gene.Symbol!="")
dis.genes <- unique(c(dis.genes,new.omim.genes$Approved.Gene.Symbol))

######## Plot overlap % v/s number of genes whose atleast x% coding sequence is covered by LRS or SRS-only.
# Number of overlap fractions for calculating cumulative features
n_ovlfracs <- 10
probands <- c("M11AO","M11BG","M11AJ","M11AU","M11AV","M11AR","M11BN","M11B1","M11B4","M11BA","M11BD","M11BJ","M11BM","M11BO","M11BR","M11BU","M11AH","M11B7","M11C1","M11CV","M11AI")

# Create an empty list to store the plots
plot_list <- list()

###### FOR CDS overlaps
calculate_counts <- function(df_ovl, output_name) {
  vec <- numeric(0)  # Initialize empty vector
  for (x in seq(0, 1, length.out = n_ovlfracs)) {
    count <- df_ovl %>% filter(tot.ovl.perc >= x) 
    count <- length(unique(count$gene_name))
    vec <- c(vec, count)  # Append count to x vector
  }
  assign(output_name, vec, envir = .GlobalEnv)  # Assign lr_vec to global environment with specified output_name
}

###### FOR GENE overlaps
calculate_gene_counts <- function(df_ovl, output_name) {
  vec <- numeric(0)  # Initialize empty vector
  for (x in seq(0, 1, length.out = n_ovlfracs)) {
    count <- df_ovl %>% filter(ovl.perc >= x) %>% nrow()
    vec <- c(vec, count)  # Append count to x vector
  }
  assign(output_name, vec, envir = .GlobalEnv)  # Assign lr_vec to global environment with specified output_name
}

################## Overlaps with disease genes ###################
calculate_disgenes <- function(df_ovl, output_name) {
  vec <- numeric(0)  # Initialize empty vector
  for (x in seq(0, 1, length.out = n_ovlfracs)) {
    d1 <- df_ovl %>% filter(ovl.perc >= x)
    count <- length(d1$gene_name[which(d1$gene_name %in% dis.genes)])
    vec <- c(vec, count)  # Append count to x vector
  }
  assign(output_name, vec, envir = .GlobalEnv)  # Assign lr_vec to global environment with specified output_name
}

library(openxlsx)
# Create a new workbook
wb <- createWorkbook()
# Add a new sheet
addWorksheet(wb, sheetName = "overlap_analysis_output")
hs1 <- createStyle(
  fgFill = "#DCE6F1", halign = "CENTER", textDecoration = "italic",
  border = "Bottom"
)

master.df <- data.frame()

#### CDS
## Location of CDS
meta.cds <- read.csv("subset_hg38.gencode.CDS.merged.sort.bed", sep="\t", header = FALSE, col.names = c('chrom','start','end','gene_name'))
# label CDS number
meta.cds <- meta.cds %>% group_by(gene_name) %>% mutate(cds_number=row_number())

for (sample in probands){
  
  print(sample)
  
  ############### LRS CDS
  lr_cds <- read.csv(paste0("proband_files/",sample,".inlrs_notsrs.intersectBED.cds.int.bed"), sep="\t", header = FALSE, 
                     col.names=c('chrom', 'start', 'end', 'gene_name', 'chrom_cov', 'start_cov', 'end_cov', 'label', 'length', 'ovl'))
  lr_cds$chrom <- as.character(lr_cds$chrom)
  lr_cds$gene_name <- as.character(lr_cds$gene_name)
  lr_cds <- lr_cds %>% left_join(meta.cds, by=c('chrom', 'start', 'end', 'gene_name'))
  # only keep rows which overlapped genes
  lr_cds <- lr_cds %>% distinct() %>% filter(label=="NO_COVERAGE") %>% mutate(cds.length=(end-start))
  # find total overlapped bases per gene, and percentage of overlap
  lr_cds_ovl <- lr_cds %>% mutate(ovl.perc=(ovl/cds.length))
  # If CDS is the same with many LRS-only callable segments, then add ovl.perc
  lr_cds_ovl <- lr_cds_ovl %>% group_by(chrom,start,end,gene_name) %>%
    mutate(tot.ovl.perc=sum(ovl.perc))
  # Only keep needed columns
  lr_cds_ovl_fin <- lr_cds_ovl %>% dplyr::select(chrom,start,end,gene_name,cds.length,tot.ovl.perc,cds_number) %>%
    distinct()
  calculate_counts(lr_cds_ovl_fin, "lr_cds_vec")
  # Only keep first coding exon
  lr_cds_ovl_fin_c1 <- lr_cds_ovl_fin %>% filter(cds_number==1)
  calculate_counts(lr_cds_ovl_fin_c1, "lr_cds_c1_vec")
  
  ################ SRS CDS
  sr_cds <- read.csv(paste0("proband_files/",sample,".insrs_notlrs.intersectBED.cds.int.bed"), sep="\t", header = FALSE, 
                     col.names=c('chrom', 'start', 'end', 'gene_name', 'chrom_cov', 'start_cov', 'end_cov', 'label', 'length', 'ovl'))
  sr_cds$chrom <- as.character(sr_cds$chrom)
  sr_cds$gene_name <- as.character(sr_cds$gene_name)
  sr_cds <- sr_cds %>% left_join(meta.cds, by=c('chrom', 'start', 'end', 'gene_name'))
  # only keep rows which overlapped genes
  sr_cds <- sr_cds %>% distinct() %>% filter(label=="NO_COVERAGE") %>% mutate(cds.length=(end-start))
  # find total overlapped bases per gene, and percentage of overlap
  sr_cds_ovl <- sr_cds %>% mutate(ovl.perc=(ovl/cds.length))
  sr_cds_ovl <- sr_cds_ovl %>% group_by(chrom,start,end,gene_name) %>%
    mutate(tot.ovl.perc=sum(ovl.perc))
  # Only keep needed rows
  sr_cds_ovl_fin <- sr_cds_ovl %>% dplyr::select(chrom,start,end,gene_name,cds.length,tot.ovl.perc,cds_number) %>%
    distinct()
  calculate_counts(sr_cds_ovl_fin, "sr_cds_vec")
  # Only keep first coding exon
  sr_cds_ovl_fin_c1 <- sr_cds_ovl_fin %>% filter(cds_number==1)
  calculate_counts(sr_cds_ovl_fin_c1, "sr_cds_c1_vec")
  
  ############### LRS GENE
  # This is LRS intersection file with overlaps
  lr <- read.csv(paste0("proband_files/",sample,".inlrs_notsrs.intersectBED.int.bed"), sep="\t", header = FALSE, 
                 col.names=c('chrom', 'start', 'end', 'gene_name', 'chrom_cov', 'start_cov', 'end_cov', 'label', 'length', 'ovl'))
  # only keep rows which overlapped genes
  lr <- lr %>% filter(label=="NO_COVERAGE") %>% mutate(gene.length=(end-start))
  # find total overlapped bases per gene, and percentage of overlap
  lr_ovl <- lr %>% group_by(gene_name, gene.length) %>% summarise(total.ovl=sum(ovl)) %>%
    mutate(ovl.perc=(total.ovl/gene.length))
  calculate_gene_counts(lr_ovl, "lr_gene")
  
  ################ SRS GENE
  # This is SRS intersection file with overlaps
  lr <- read.csv(paste0("proband_files/",sample,".insrs_notlrs.intersectBED.int.bed"), sep="\t", header = FALSE, 
                 col.names=c('chrom', 'start', 'end', 'gene_name', 'chrom_cov', 'start_cov', 'end_cov', 'label', 'length', 'ovl'))
  # only keep rows which overlapped genes
  sr <- sr %>% filter(label=="NO_COVERAGE") %>% mutate(gene.length=(end-start))
  # find total overlapped bases per gene, and percentage of overlap
  sr_ovl <- sr %>% group_by(gene_name, gene.length) %>% summarise(total.ovl=sum(ovl)) %>%
    mutate(ovl.perc=(total.ovl/gene.length))
  calculate_gene_counts(sr_ovl, "sr_gene")
  
  ########### DIS-GENES
  calculate_disgenes(lr_ovl, "lr_disgene")
  calculate_disgenes(sr_ovl, "sr_disgene")
  
  ##### if you want to plot disease genes as well ######
  tot <- c(lr_cds_vec, lr_cds_c1_vec, lr_disgene, lr_gene, sr_cds_vec, sr_cds_c1_vec, sr_disgene, sr_gene)
  ovlap_frac <- seq(0, 1, length.out = n_ovlfracs)
  df <- data.frame(
    overlap_fraction = rep(ovlap_frac,8),
    cumulative_features = tot,
    Feature=c(rep("Coding Exon",n_ovlfracs), rep("First Coding Exon",n_ovlfracs), rep("Disease Gene",n_ovlfracs), rep("Gene",n_ovlfracs), rep("Coding Exon",n_ovlfracs), rep("First Coding Exon",n_ovlfracs), rep("Disease Gene",n_ovlfracs), rep("Gene",n_ovlfracs)),
    tech = c(rep("LRS-only",n_ovlfracs*4), rep("SRS-only",n_ovlfracs*4))
  )
  
  df <- df %>% mutate(sample=sample)
  
  ## append to master.df
  master.df <- rbind(master.df, df)
  
  ###############################################################
  my_colors <- c("Coding Exon" = "#D9565CFF", "First Coding Exon" = "#F1A226", "Gene" = "#449DB3FF", "Disease Gene" = "#A559AA")
  
  p <- ggplot(df, aes(x=overlap_fraction, y=cumulative_features, group=interaction(tech, Feature), color=Feature)) +
    geom_point(stat="identity", size=0.3) +
    geom_line(linewidth=1) + 
    theme_minimal() +
    facet_grid(.~tech) +
    labs(x = "Overlap fraction", y = "Cumulative feature counts") +
    ggtitle(paste0(sample)) +
    scale_x_continuous(breaks=seq(0,1.1,0.2), limits=c(0,1)) +
    scale_y_continuous(breaks=seq(0,700,100), limits=c(0,700)) +
    scale_color_manual(values=my_colors) +
    theme_bw() +
    theme(legend.position = "none")
  
  # Store each plot in the list
  plot_list[[sample]] <- p
  
}

# Print all plots in a grid layout
library(gridExtra)
grid.arrange(grobs = plot_list, ncol = 6)

# Save the workbook
# Keep updating into an excel workbook
# Write the data to the sheet
writeData(wb, sheet = "overlap_analysis_output", x = master.df, 
          borders = "rows", headerStyle = hs1, borderStyle = "dashed")
saveWorkbook(wb, "gene-overlap-results.xlsx", overwrite = TRUE)



