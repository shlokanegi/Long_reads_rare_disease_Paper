## Prepare Gencode annotation files for protein-coding genes and CDS

#GENCODE Release 45 - https://www.gencodegenes.org/human/release_45.html. Used the comprehensive gene annotation set on the reference chromosomes only

zcat gencode.v45.annotation.gtf.gz | awk 'BEGIN{FS="\t"}{split($9,a,";"); if($3~"gene") print a[1]"\t"a[3]"\t"$1":"$4"-"$5"\t"a[2]"\t"$7}' | \
    sed 's/gene_id "//' | sed 's/gene_id "//' | sed 's/gene_type "//'| sed 's/gene_name "//' | sed 's/"//g' | \
    awk 'BEGIN{FS="\t"}{split($3,a,"[:-]"); print $1"\t"$2"\t"a[1]"\t"a[2]"\t"a[3]"\t"$4"\t"$5"\t"a[3]-a[2];}' | \
    sed "1i\Geneid\tGeneSymbol\tChromosome\tStart\tEnd\tClass\tStrand\tLength" > subset_hg38.gencode.tsv
## All protein_coding genes
awk 'NR>1 && $6=="protein_coding"' subset_hg38.gencode.tsv > subset_hg38.gencode.genes.tsv
protein_coding_genes="subset_hg38.gencode.genes.tsv"

## All protein coding CDS regions
zcat gencode.v45.annotation.gtf.gz | awk 'BEGIN{FS="\t"}{split($9,a,";"); if($3~"CDS") print a[1]"\t"a[4]"\t"$1":"$4"-"$5"\t"a[3]"\t"$7}' | \
    sed 's/gene_id "//' | sed 's/gene_id "//' | sed 's/gene_type "//'| sed 's/gene_name "//' | sed 's/"//g' | \
    awk 'BEGIN{FS="\t"}{split($3,a,"[:-]"); print $1"\t"$2"\t"a[1]"\t"a[2]"\t"a[3]"\t"$4"\t"$5"\t"a[3]-a[2];}' | \
    sed "1i\Geneid\tGeneSymbol\tChromosome\tStart\tEnd\tClass\tStrand\tLength" > subset_hg38.gencode.CDS.tsv

# --------------------------------------------------------------------------------- #
cds_tsv="subset_hg38.gencode.CDS.tsv"
gene_tsv="subset_hg38.gencode.genes.tsv"

## convert to sorted BED
awk 'NR>1 {print $3"\t"$4"\t"$5"\t"$2}' $gene_tsv | sort -k1,1 -k2,2n -k3,3n > subset_hg38.gencode.genes.sort.bed
awk 'NR>1 {print $3"\t"$4"\t"$5"\t"$2}' $cds_tsv | sort -k1,1 -k2,2n -k3,3n > subset_hg38.gencode.CDS.sort.bed
# merged overlapping CDS regions per gene
bedtools merge -i subset_hg38.gencode.CDS.sort.bed -c 4 -o collapse | awk 'BEGIN{FS="\t"; OFS="\t"}{split($4,a,","); print $1, $2, $3, a[1]}' | sort -k1,1 -k2,2n -k3,3n > subset_hg38.gencode.CDS.merged.sort.bed

# --------------------------------------------------------------------------------- #
#### Commands for generating gene-coverage overlap BED files, used for technology exclusive gene coverage overlap analysis

for sample_id in $(cat probands.txt); 
do
    echo "Running $sample_id"
    
    cat $sample_id/lrs.quant.bed | awk '$4=="CALLABLE"' > lrs.cov.callable.bed
    cat $sample_id/srs.quant.bed | awk '$4=="CALLABLE"' > srs.cov.callable.bed

    ####### Generate intersected BED files to view on IGV and stratify SRS-only and LRS-only covered regions
    bedtools intersect -a $sample_id/srs.nocov.bed -b lrs.cov.callable.bed | cut -f1-4 | awk '{FS=OFS="\t"}{print $0, ($3-$2)}' > inlrs_notsrs.callable.intersectBED.bed
    bedtools intersect -a $sample_id/lrs.nocov.bed -b srs.cov.callable.bed | cut -f1-4 | awk '{FS=OFS="\t"}{print $0, ($3-$2)}' > insrs_notlrs.callable.intersectBED.bed

    gene_bed="subset_hg38.gencode.genes.sort.bed"
    # We want to find overlaps of each gene with the coverage file and amount of overlaps
    bedtools intersect -a $gene_bed -b $lrs_only -wao > $sample_id.inlrs_notsrs.intersectBED.int.bed # imported in R for data visualization and plotting
    bedtools intersect -a $gene_bed -b $srs_only -wao > $sample_id.insrs_notlrs.intersectBED.int.bed # imported in R for data visualization and plotting

    ## For CDS
    cds_bed="subset_hg38.gencode.CDS.merged.sort.bed"
    bedtools intersect -a "$cds_bed" -b "$lrs_only" -wo > $sample_id.inlrs_notsrs.intersectBED.cds.int.bed # imported in R for data visualization and plotting
    bedtools intersect -a "$cds_bed" -b "$srs_only" -wo > $sample_id.insrs_notlrs.intersectBED.cds.int.bed # imported in R for data visualization and plotting
done