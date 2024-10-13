library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(ggpubr)
library(tidyverse)
library(patchwork)
library(readxl)

# Read all files
loadTSV <- function(path, label){
  df = read.table(path, sep="\t", header = TRUE)
  df <- df %>% arrange(desc(phase_block_size)) %>% 
    mutate(phase_block_size = phase_block_size/1e6,
           cumcov = cumsum(as.numeric(phase_block_size)),
           cumcov.perc = cumcov/3088.269832,
           label=label)
}

probands <- c('M11AO', 'M11AU', 'M11AR', 'M11AV', 'M11BN', 'M11B1', 'M11B4', 'M11BA', 'M11BD', 'M11BG', 'M11BJ', 'M11BM', 'M11BO', 'M11BR', 'M11BU', 'M11AH', 'M11B7', 'M11C1', 'M11CV', 'M11AI', 'M11AJ', 'R10_11110', 'R10_68680', 'R10_1461460', 'R10_1641640', 'R10_1671670', 'CVG60720', 'DEN63290', 'R10_86860','ARB20220','CVG60130','DEN63190','DEN63270','DEN63360','DCA70650','LGA80170','LGA80420','LGA80510','LGA80520','ORD90210','ORD90570','ARB21022','HRFD01250')
parents <- c('M11AN','M11AM','M11AT','M11AS','M11AW','M11AX','M11AP','M11AQ','M11AZ','M11AY','M11B3','M11B2','M11B9','M11B8','M11BC','M11BB','M11BF','M11BE','M11BI','M11BH','M11BL','M11BK','M11BQ','M11BP','M11BT','M11BS','M11BW','M11BV','M11D1','M11CZ','M11B6','M11B5','M11C3','M11C2','M11CU','M11CT','M11CX','M11CW','M11AL','M117Z','R10_139111','R10_12112','R10_70681','R10_69682','R10_1571461','R10_1471462','R10_1661641','R10_1651642','R10_1691671','R10_1681672','CVG60721','CVG60722','DEN63291','DEN63292','R10_80862')

all_samples <- c(probands,parents)

###########################################################################################
##############################  Raw Sequencing Stats  ###############################
###########################################################################################
###########################################################################################
seq.df <- read_excel('raw_seqdata_stats.xlsx', sheet = "Sheet1")

seq.cut.df <- seq.df %>% dplyr::select(c(2,3,4)) %>% mutate(read.n50=read.n50/1e3)
# For only samples with 1 FC
seq.cut.df <- seq.df %>% filter(FCs==1) %>% dplyr::select(c(2,3,4)) %>% mutate(read.n50=read.n50/1e3)

# Reshape the dataframe to long format
seq.cut.df.long <- seq.cut.df %>%
  pivot_longer(cols = c("read.n50", "total.gbp", "mean.identity"), 
               names_to = "metric", values_to = "value")

# Split the data into three separate dataframes
seq.cut.read.n50 <- seq.cut.df.long %>% filter(metric == "read.n50")
seq.cut.total.gbp <- seq.cut.df.long %>% filter(metric == "total.gbp")
seq.cut.mean.identity <- seq.cut.df.long %>% filter(metric == "mean.identity")

# Custom x-tick labels
custom_labels <- c("read.n50" = "", "total.gbp" = "", "mean.identity" = "")

# Create individual plots
plot_read_n50 <- ggplot(seq.cut.read.n50, aes(x=metric, y=value)) +
  geom_violin(trim=FALSE, scale = "width", alpha=0.5, fill="#CBB3BFFF") +
  geom_boxplot(width = 0.2, outlier.shape = NA, fill="#CBB3BFFF") +
  stat_summary(fun=median, geom="point", fill="blue", shape=21, size=1) + 
  theme_bw() +
  labs(title = "", x = " ", y = "read length N50 (Kbp)", fill = " ") +
  scale_x_discrete(labels = custom_labels) +
  theme(legend.position = "none")

plot_total_gbp <- ggplot(seq.cut.total.gbp, aes(x=metric, y=value)) +
  geom_violin(trim=FALSE, scale = "width", alpha=0.5, fill="#FFC857FF") +
  geom_boxplot(width = 0.2, outlier.shape = NA, fill="#FFC857FF") +
  stat_summary(fun=median, geom="point", fill="blue", shape=21, size=1) + 
  theme_bw() +
  labs(title = "", x = " ", y = "total sequence (Gbp)", fill = " ") +
  scale_y_continuous(sec.axis = sec_axis(~ . / 3.1, name = "genome coverage")) +
  scale_x_discrete(labels = custom_labels) +
  theme(legend.position = "none")

plot_mean_identity <- ggplot(seq.cut.mean.identity, aes(x=metric, y=value)) +
  geom_violin(trim=FALSE, scale = "width", alpha=0.5, fill="#119DA4FF") +
  geom_boxplot(width = 0.2, outlier.shape = NA,fill="#119DA4FF") +
  stat_summary(fun=median, geom="point", fill="blue", shape=21, size=1) + 
  theme_bw() +
  labs(title = "", x = " ", y = "identity", fill = " ") +
  scale_x_discrete(labels = custom_labels) +
  theme(legend.position = "none")

# Combine the plots using patchwork
combined_plot <- plot_total_gbp + plot_spacer() + plot_read_n50 + plot_spacer() + plot_mean_identity +
  plot_layout(ncol = 5, widths = c(1, 0.1, 1, 0.1, 1))
  
# Display the combined plot
combined_plot

## Save as PDF for AI
ggsave(filename = "Figure1B_1.pdf", plot = plot_read_n50, width = 6, height = 8, units = "cm", dpi=300)
ggsave(filename = "Figure1B_2.pdf", plot = plot_total_gbp, width = 6, height = 8, units = "cm", dpi=300)
ggsave(filename = "Figure1B_3.pdf", plot = plot_mean_identity, width = 6, height = 8, units = "cm", dpi=300)


###########################################################################################
##############################  Phase block stats  ###############################
###########################################################################################
###########################################################################################
seq.df <- read_excel('raw_seqdata_stats.xlsx', sheet = "Sheet1")
seq.df <- seq.df %>% mutate(read.n50=read.n50/1e3)

al.df = lapply(all_samples, function(x){
  loadTSV(paste0('vcfstats_outputs_phasing/', x, '_hvcf.phased.phase_block_stats.tsv'), x)
}) %>% bind_rows
m.df <- merge(al.df,seq.df, by=c('label'))

## Plotting Phase block NGx
library(wesanderson)
library(scico)
a <- m.df %>% ggplot(aes(x = cumcov.perc, y = phase_block_size, group=label, color=read.n50)) +
  geom_line(alpha=0.6, size=0.8) +
  geom_vline(xintercept = 0.5, linetype=2) +
  labs(x = "Cumulative Coverage", y = "Phase Block Size (In Mbp)", title = "Phase Block NGx") +
  scale_x_continuous(limits = c(0,1), breaks = seq(0,1,0.2)) +
  scale_y_continuous(limits = c(0,40), breaks = seq(0,40,5)) +
  scale_color_scico(palette = "hawaii") + 
  theme_bw()
ggsave(filename = "SuppFig16.pdf", plot = a, width = 17, height = 12, units = "cm", dpi=300)


### Scatter plots between Ng50s and raw sequencing stats (Read N50, coverage)
# Custom function to calculate N50
calculate_n50 <- function(sizes) {
  sorted_sizes <- sort(sizes)
  half_total_size <- sum(sorted_sizes) / 2
  cumulative_sum <- 0
  n50 <- 0
  
  for (size in sorted_sizes) {
    cumulative_sum <- cumulative_sum + size
    if (cumulative_sum >= half_total_size) {
      n50 <- size
      break
    }
  }
  
  return(n50)
}

p.n50.df <- al.df %>%
  group_by(label) %>%
  summarise(
    half_tot_size = sum(phase_block_size) / 2,
    phase.n50 = calculate_n50(phase_block_size)
  )

m.df <- merge(p.n50.df,seq.df, by=c('label'))

## Correlation between pb.n50 and read.N50, colored by read coverage
m.df %>% mutate(FCs=as.factor(FCs)) %>%
  ggplot(aes(x=phase.n50, y=read.n50, colour = cov)) +
  geom_point() +
  scale_y_continuous(limits = c(0,50), breaks = seq(0,50,10)) +
  #geom_text(aes(label = label)) +
  labs(x='phase.n50 (in Mbp)', y='read.n50 (in Kbp)') +
  theme_bw()

######## Number of variants per phase-block.
al.df %>% ggplot(aes(x = phase_block_size, y = variant_count)) +
  geom_line() +
  facet_wrap(.~label) +
  theme_bw() +
  theme(legend.position = "none")

al.df %>% filter(label=="M117Z") %>% 
  ggplot(aes(x = phase_block_size, y = variant_count)) +
  geom_point() +
  theme_bw() +
  theme(legend.position = "none")


##########################################################################################
#######################       Genes in phase-blocks analysis    ##########################
##########################################################################################
# Read phased_genes files
loadBED <- function(path, label){
  df = read.table(path, sep="\t", header = FALSE, col.names = c('chr','start','end','gene','n_phase_blocks'))
  df <- df %>% mutate(label=label,
                      n_phase_blocks=as.factor(n_phase_blocks))
}
## All samples
## GENES
pb.genes.df = lapply(all_samples, function(x){
  loadBED(paste0('phased_genes/', x, '_phase_block_genes.bed'), x)
}) %>% bind_rows

pb.genes.count.df <- pb.genes.df %>%
  group_by(label,n_phase_blocks) %>%
  summarise(gene_counts=n()) %>%
  mutate(perc.pcgenes=gene_counts/20019)

dataMedian <- pb.genes.count.df %>%
  group_by(n_phase_blocks) %>% summarise(MD=median(gene_counts), q95=quantile(gene_counts,probs=0.95))
pb.genes.count.df %>%
  ggplot(aes(x = n_phase_blocks, y=gene_counts, fill=n_phase_blocks)) +
  geom_violin(trim=FALSE, scale = "width", alpha=0.5) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +  # Add boxplot for better visualization
  geom_text(data=dataMedian, aes(y=q95, label=round(MD,0)), position=position_dodge(width = 0.8), vjust=-3, hjust=1.2,size = 4, color="black") +
  stat_summary(fun=median, geom="point", fill="blue", shape=21, size=1) + 
  theme_bw() +
  labs(title = "", x = " ", y = "Number of protein-coding genes", fill = "Number of phase blocks") +
  theme(legend.position = "bottom")

# in fraction
dataMedian <- pb.genes.count.df %>%
  group_by(n_phase_blocks) %>% summarise(MD=median(perc.pcgenes), q95=quantile(perc.pcgenes,probs=0.95))
pb.genes.count.df %>%
  ggplot(aes(x = n_phase_blocks, y=perc.pcgenes, fill=n_phase_blocks)) +
  geom_violin(trim=FALSE, scale = "width", alpha=0.5) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +  # Add boxplot for better visualization
  geom_text(data=dataMedian, aes(y=q95, label=round(MD,2)), position=position_dodge(width = 0.8), vjust=-3, hjust=1.2,size = 4, color="black") +
  stat_summary(fun=median, geom="point", fill="blue", shape=21, size=1) + 
  theme_bw() +
  labs(title = "", x = " ", y = "Fraction of protein-coding genes", fill = "Number of phase blocks") +
  theme(legend.position = "bottom")

## All samples
## CDS
pb.cds.df = lapply(all_samples, function(x){
  loadBED(paste0('phased_cds/', x, '_phase_block_cds.bed'), x)
}) %>% bind_rows

pb.cds.count.df <- pb.cds.df %>%
  group_by(label,n_phase_blocks) %>%
  summarise(cds_counts=n())

dataMedian <- pb.cds.count.df %>%
  group_by(n_phase_blocks) %>% summarise(MD=median(cds_counts), q95=quantile(cds_counts,probs=0.95))
pb.cds.count.df %>%
  ggplot(aes(x = n_phase_blocks, y=cds_counts, fill=n_phase_blocks)) +
  geom_violin(trim=FALSE, scale = "width", alpha=0.8) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +  # Add boxplot for better visualization
  geom_text(data=dataMedian, aes(y=q95, label=MD), position=position_dodge(width = 0.8), vjust=-3, hjust=1.2,size = 4, color="black") +
  stat_summary(fun=median, geom="point", fill="blue", shape=21, size=1) + 
  theme_bw() +
  labs(title = "", x = " ", y = "Number of conding exons", fill = "Number of phase blocks") +
  theme(legend.position = "bottom")

pb.genes.df.single <- pb.genes.df %>% filter(n_phase_blocks==1)

##########################################################################################
##########################################################################################
# Read phased_genes_overlap files
loadBED <- function(path, label){
  df = read.table(path, sep="\t", header = FALSE, col.names = c('chr','start','end','gene','p_chrom','p_start','p_end','p_size', 'n_vars', 'bases.ovl'))
  df <- df %>% mutate(label=label,g_size=end-start) %>% dplyr::select(label,chr,start,end,g_size,gene,p_chrom,p_start,p_end,p_size,bases.ovl)
  df <- df %>% group_by(label,chr,start,end,g_size,gene) %>% mutate(n_phase_blocks=n()) %>% ungroup() %>%
    group_by(label,p_chrom,p_start,p_end,p_size) %>% mutate(n_genes=ifelse(p_start==-1,NA,n()))
  df <- df %>% mutate(gene.ovl.frac=bases.ovl/g_size)
}

## All samples
## GENES
pb.genes.df = lapply(all_samples, function(x){
  loadBED(paste0('phased_genes_with_overlap/', x, '_phase_block_genes.bed'), x)
}) %>% bind_rows

## For each gene, only choose one best phase_block, i.e. one with largest overlap fraction.
# In case of a tie, keep any one of the phase-blocks
pb.genes.filt.df <- pb.genes.df %>%
  group_by(label, gene) %>%
  arrange(desc(gene.ovl.frac)) %>%
  distinct(label, gene, .keep_all = TRUE)

## Distribution of gene.ovl.fraction per sample
pb.genes.filt.df %>%
  ggplot(aes(x = gene.ovl.frac, color = label)) +
  geom_histogram(aes(y = ..count..), position = "identity", alpha = 0.4, binwidth = 0.05, fill="lightgrey") +  # Adjust binwidth as needed
  #geom_density(alpha = 0.7) +
  labs(x = "Gene Overlap Fraction", y = "Count", title = "Histogram with Counts by Label") +
  theme_bw() +
  theme(legend.position = "none")  # Remove legend


pb.genes.count.df <- pb.genes.filt.df %>%
  group_by(label) %>% summarise(n_0=sum(gene.ovl.frac==0),
                                n_0_25=sum(gene.ovl.frac>0 & gene.ovl.frac<0.25),
                                n_25_50=sum(gene.ovl.frac>=0.25 & gene.ovl.frac<0.5),
                                n_50_75=sum(gene.ovl.frac>=0.5 & gene.ovl.frac<0.75),
                                n_75_100=sum(gene.ovl.frac>=0.75 & gene.ovl.frac<1),
                                n_100=sum(gene.ovl.frac==1))

# Reshape the data for plotting
pb.genes.count.df_long <- pb.genes.count.df %>%
  pivot_longer(cols = starts_with("n_"), names_to = "Overlap_Fraction_Bin", values_to = "Count")
pb.genes.count.df_long$Overlap_Fraction_Bin <- factor(pb.genes.count.df_long$Overlap_Fraction_Bin, levels = c("n_0", "n_0_25", "n_25_50", "n_50_75", "n_75_100", "n_100"))

dataMedian <- pb.genes.count.df_long %>%
  group_by(Overlap_Fraction_Bin) %>%
  summarise(MD = round(median(Count),0),
            q95 = round(quantile(Count, probs = 0.95),0))

# Distribution plot
a <- pb.genes.count.df_long %>%
  ggplot(aes(x = Overlap_Fraction_Bin, y = Count, fill = Overlap_Fraction_Bin)) +
  geom_violin(trim = FALSE, scale = "width", alpha = 0.5) +
  geom_boxplot(width = 0.2, outlier.shape = NA, position = position_dodge(width = 0.75)) +
  geom_point(data = dataMedian, aes(y = MD), color = "blue", size = 1, position = position_dodge(width = 0.75)) +
  geom_text(data = dataMedian, aes(label = round(MD, 1), y = q95 + 0.1), color = "black", size = 3, vjust = -4, position = position_dodge(width = 0.75)) +  # Add median text
  scale_x_discrete(labels = NULL) +
  labs(title = "Violin Plot of Gene Counts by Overlap Fraction Bin",
       y = "Gene Count",
       fill = "Overlap Fraction Bin") +
  theme_bw()
a

################## Let's make a better plot #################
calculate_cumulative_counts <- function(df, n_bins, output_name) {
  # Create bin values
  bin_values <- seq(0, 1, length.out = n_bins)
  
  # Expand the dataframe to include all bin values for each label
  expanded_df <- df %>%
    crossing(bin_value = bin_values) %>%
    group_by(label, bin_value) %>%
    summarize(cumulative_count = sum(gene.ovl.frac >= bin_value), .groups = 'drop')
  
  # Assign the resulting dataframe to the global environment with the specified name
  assign(output_name, expanded_df, envir = .GlobalEnv)
}

calculate_cumulative_counts(pb.genes.filt.df, 50, "cumulative_counts_df")

# Merge phase-block N50 column
cumulative_counts_df <- merge(cumulative_counts_df,p.n50.df, by="label")

###### Check if 0 phase-block gene numbers differ between probands and parents
probands <- c('M11AO', 'M11AU', 'M11AR', 'M11AV', 'M11BN', 'M11B1', 'M11B4', 'M11BA', 'M11BD', 'M11BG', 'M11BJ', 'M11BM', 'M11BO', 'M11BR', 'M11BU', 'M11AH', 'M11B7', 'M11C1', 'M11CV', 'M11AI', 'M11AJ', 'R10_11110', 'R10_68680', 'R10_1461460', 'R10_1641640', 'R10_1671670', 'CVG60720', 'DEN63290', 'R10_86860','ARB20220','CVG60130','DEN63190','DEN63270','DEN63360','DCA70650','LGA80170','LGA80420','LGA80510','LGA80520','ORD90210','ORD90570','ARB21022','HRFD01250')
parents <- c('M11AN','M11AM','M11AT','M11AS','M11AW','M11AX','M11AP','M11AQ','M11AZ','M11AY','M11B3','M11B2','M11B9','M11B8','M11BC','M11BB','M11BF','M11BE','M11BI','M11BH','M11BL','M11BK','M11BQ','M11BP','M11BT','M11BS','M11BW','M11BV','M11D1','M11CZ','M11B6','M11B5','M11C3','M11C2','M11CU','M11CT','M11CX','M11CW','M11AL','M117Z','R10_139111','R10_12112','R10_70681','R10_69682','R10_1571461','R10_1471462','R10_1661641','R10_1651642','R10_1691671','R10_1681672','CVG60721','CVG60722','DEN63291','DEN63292','R10_80862')

pb.genes.count.df <- pb.genes.filt.df %>%
  group_by(label) %>% summarise(n_0=sum(gene.ovl.frac==0),
                                n_100=sum(gene.ovl.frac==1))

test.df <- pb.genes.count.df %>% dplyr::select(label,n_0) %>% 
  mutate(stype=ifelse(label %in% probands, "proband", "parent"))

dataMedian <- test.df %>% group_by(stype) %>%
  summarise(MD = round(median(n_0),0),
            q95 = round(quantile(n_0, probs = 0.95),0))
colors <- c('blue','red')
test.df %>% group_by(stype) %>%
  ggplot(aes(x = stype, y = n_0, fill = stype)) +
  geom_violin(trim = FALSE, scale = "width", alpha = 0.5, width=0.5) +
  geom_boxplot(width = 0.2, outlier.shape = NA, position = position_dodge(width = 0.75)) +
  geom_point(data = dataMedian, aes(y = MD), color = "white", size = 1, position = position_dodge(width = 0.75)) +
  geom_text(data = dataMedian, aes(label = round(MD, 1), y = q95 + 0.1), color = "black", size = 4, vjust = -1, hjust = 2.5, position = position_dodge(width = 0.75)) +  # Add median text
  scale_fill_manual(values = colors) +  # Apply custom colors
  labs(y = "Protein-coding genes\noverlapped by 0 phase-blocks") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()
  )  # Remove x-axis text

# Plot the data using ggplot2
phased_genes <- ggplot(cumulative_counts_df, aes(x = bin_value, y = cumulative_count, color = phase.n50, group = label)) +
  geom_line(alpha=0.5, size=0.6) +
  #geom_point(alpha=0.5) +
  labs(x = "Overlap Fraction", y = "Cumulative protein-coding gene counts", color = "Phase-Block\nNG50") +
  theme_bw() +
  scale_color_scico(palette = "hawaii") + 
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(13000,20500), breaks = c(seq(13000,21000,1000))) +
  #coord_flip() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.title = element_text(),
        axis.text.x = element_text(vjust=0.5))

phased_genes

ggsave(filename = "Figure6A.pdf", plot = phased_genes, width = 6, height = 6, units = "in", dpi=300)

library(gridExtra)
grid.arrange(phased_genes, a, ncol = 2)


## BOXPLOTS
pb.genes.count.df <- pb.genes.filt.df %>%
  group_by(label) %>% summarise(n_0=sum(gene.ovl.frac==0),
                                n_100=sum(gene.ovl.frac==1))

# Reshape the data for plotting
pb.genes.count.df_long <- pb.genes.count.df %>%
  pivot_longer(cols = starts_with("n_"), names_to = "Overlap_Fraction_Bin", values_to = "Count")
pb.genes.count.df_long$Overlap_Fraction_Bin <- factor(pb.genes.count.df_long$Overlap_Fraction_Bin, levels = c("n_0", "n_0_25", "n_25_50", "n_50_75", "n_75_100", "n_100"))

dataMedian <- pb.genes.count.df_long %>%
  group_by(Overlap_Fraction_Bin) %>%
  summarise(MD = round(median(Count),0),
            q95 = round(quantile(Count, probs = 0.95),0))

# Distribution plot
## Combined
dataMedian <- pb.genes.count.df_long %>% group_by(Overlap_Fraction_Bin) %>%
  summarise(MD = round(median(Count),0),
            q95 = round(quantile(Count, probs = 0.95),0))

colors <- c('#007191','#f47a00')
c <- pb.genes.count.df_long %>% group_by(Overlap_Fraction_Bin) %>%
  ggplot(aes(x = Overlap_Fraction_Bin, y = Count, fill = Overlap_Fraction_Bin)) +
  geom_violin(trim = FALSE, scale = "width", alpha = 0.5, width=0.5) +
  geom_boxplot(width = 0.2, outlier.shape = NA, position = position_dodge(width = 0.75)) +
  geom_point(data = dataMedian, aes(y = MD), color = "blue", size = 1, position = position_dodge(width = 0.75)) +
  geom_text(data = dataMedian, aes(label = round(MD, 1), y = q95 + 0.1), color = "black", size = 4, vjust = -1, hjust = 1.5, position = position_dodge(width = 0.75)) +  # Add median text
  scale_fill_manual(values = colors) +  # Apply custom colors
  scale_x_discrete(labels = NULL) +
  labs(y = "Protein-coding genes") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()
        )  # Remove x-axis text
c

ggsave(filename = "Figure6A_2.pdf", plot = c, width = 6, height = 4, units = "in", dpi=300)


## Separate
### FOR n=0
dataMedian <- pb.genes.count.df_long %>% filter(Overlap_Fraction_Bin=="n_0") %>%
  summarise(MD = round(median(Count),0),
            q95 = round(quantile(Count, probs = 0.95),0))
a <- pb.genes.count.df_long %>% filter(Overlap_Fraction_Bin=="n_0") %>%
  ggplot(aes(x = 1, y = Count, fill = "red")) +
  geom_violin(trim = FALSE, scale = "width", alpha = 0.5, fill = "#eddca5") +
  geom_boxplot(width = 0.2, outlier.shape = NA, position = position_dodge(width = 0.75), fill = "#c99b38") +
  geom_point(data = dataMedian, aes(y = MD), color = "blue", size = 1, position = position_dodge(width = 0.75)) +
  #geom_text(data = dataMedian, aes(label = round(MD, 1), y = q95 + 0.1), color = "black", size = 3, vjust = -1, position = position_dodge(width = 0.75)) +  # Add median text
  scale_x_discrete(labels = NULL) +
  labs(y = "Protein-coding genes") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        legend.position = "none")  # Remove x-axis text

ggsave(filename = "Figure6A_0.pdf", plot = a, width = 6, height = 4, units = "in", dpi=300)


### FOR n=1
dataMedian <- pb.genes.count.df_long %>% filter(Overlap_Fraction_Bin=="n_100") %>%
  summarise(MD = round(median(Count),0),
            q95 = round(quantile(Count, probs = 0.95),0))
b <- pb.genes.count.df_long %>% filter(Overlap_Fraction_Bin=="n_100") %>%
  ggplot(aes(x = 1, y = Count)) +
  geom_violin(trim = FALSE, scale = "width", alpha = 0.5, fill = "#8fd7d7") +
  geom_boxplot(width = 0.2, outlier.shape = NA, position = position_dodge(width = 0.75), fill = "#00b0be") +
  geom_point(data = dataMedian, aes(y = MD), color = "blue", size = 1, position = position_dodge(width = 0.75)) +
  #geom_text(data = dataMedian, aes(label = round(MD, 1), y = q95 + 0.1), color = "black", size = 3, vjust = -1, position = position_dodge(width = 0.75)) +  # Add median text
  scale_x_discrete(labels = NULL) +
  labs(y = "Protein-coding genes") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        legend.position = "none")  # Remove x-axis text

ggsave(filename = "Figure6A_1.pdf", plot = b, width = 6, height = 4, units = "in", dpi=300)


# in fraction
dataMedian <- pb.genes.count.df %>%
  group_by(n_phase_blocks) %>% summarise(MD=median(perc.pcgenes), q95=quantile(perc.pcgenes,probs=0.95))
pb.genes.count.df %>%
  ggplot(aes(x = n_phase_blocks, y=perc.pcgenes, fill=n_phase_blocks)) +
  geom_violin(trim=FALSE, scale = "width", alpha=0.5) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +  # Add boxplot for better visualization
  geom_text(data=dataMedian, aes(y=q95, label=round(MD,2)), position=position_dodge(width = 0.8), vjust=-3, hjust=1.2,size = 4, color="black") +
  stat_summary(fun=median, geom="point", fill="blue", shape=21, size=1) + 
  theme_bw() +
  labs(title = "", x = " ", y = "Fraction of protein-coding genes", fill = "Number of phase blocks") +
  theme(legend.position = "bottom")

## All samples
## CDS
pb.cds.df = lapply(all_samples, function(x){
  loadBED(paste0('phased_cds/', x, '_phase_block_cds.bed'), x)
}) %>% bind_rows

pb.cds.count.df <- pb.cds.df %>%
  group_by(label,n_phase_blocks) %>%
  summarise(cds_counts=n())

dataMedian <- pb.cds.count.df %>%
  group_by(n_phase_blocks) %>% summarise(MD=median(cds_counts), q95=quantile(cds_counts,probs=0.95))
pb.cds.count.df %>%
  ggplot(aes(x = n_phase_blocks, y=cds_counts, fill=n_phase_blocks)) +
  geom_violin(trim=FALSE, scale = "width", alpha=0.8) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +  # Add boxplot for better visualization
  geom_text(data=dataMedian, aes(y=q95, label=MD), position=position_dodge(width = 0.8), vjust=-3, hjust=1.2,size = 4, color="black") +
  stat_summary(fun=median, geom="point", fill="blue", shape=21, size=1) + 
  theme_bw() +
  labs(title = "", x = " ", y = "Number of conding exons", fill = "Number of phase blocks") +
  theme(legend.position = "bottom")

pb.genes.df.single <- pb.genes.df %>% filter(n_phase_blocks==1)

