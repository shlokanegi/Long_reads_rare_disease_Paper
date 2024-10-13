#!/bin/bash

'''
Author: Shloka Negi, shnegi@ucsc.edu
Purpose: This script generates callable coverage statistics for SRS and LRS BAMs aligned to GRCh38 and identifies exclusive coverage regions for each technology.
Input file requirements: Requires input LRS and SRS BAM/CRAM (aligned to GRCh38, sorted and indexed)
Usage: ./coveval_grch38.sh -l lrs.bam -s srs.bam
'''

# Using getopts to parse the two required input files
while getopts l:s: flag
do
    case "${flag}" in
        l) lrs_bam=${OPTARG};;
        s) srs_bam=${OPTARG};;
    esac
done

# Displaying the input values
echo "LRS BAM: $lrs_bam"
echo "SRS BAM: $srs_bam"

# Set the exit code of a pipeline to that of the rightmost command
# to exit with a non-zero status, or zero if all commands of the pipeline exit
set -o pipefail
# cause a bash script to exit immediately when a command fails
set -e
# cause the bash shell to treat unset variables as an error and exit immediately
set -u

pwd; hostname; date

sample="$(basename ${lrs_bam} .haplotagged.bam)"
REF="GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"

mkdir -p $sample

######## Run mosdepth
echo "run mosdepth"
cd $sample
#1. SRS
mosdepth --threads 128 -f $REF --use-median --fast-mode --mapq 10 --quantize 0:1:10:150: $sample.srs.quantized $srs_bam
zcat $sample.srs.quantized.quantized.bed.gz | awk '{FS=OFS="\t"}{print $0, ($3-$2)}' > $sample.srs.quantized.tsv
awk 'BEGIN {FS=OFS="\t"} {
    if (($1 ~ /^(chr)?[0-9]+$/) || ($1 ~ /^(chr)?[XY]$/)) { 
        if ($4 == "0:1") {
            $4 = "NO_COVERAGE"
        } else if ($4 == "1:10") {
            $4 = "LOW_COVERAGE"
        } else if ($4 == "10:150") {
            $4 = "CALLABLE"
        } else if ($4 == "150:inf") {
            $4 = "HIGH_COVERAGE"
        }
        print
    }
}' $sample.srs.quantized.tsv > srs.quant.bed

#2. LRS
mosdepth --threads 128 -f $REF --use-median --fast-mode --mapq 10 --quantize 0:1:10:150: $sample.lrs.quantized $lrs_bam
zcat $sample.lrs.quantized.quantized.bed.gz | awk '{FS=OFS="\t"}{print $0, ($3-$2)}' > $sample.lrs.quantized.tsv 
awk 'BEGIN {FS=OFS="\t"} {
    if (($1 ~ /^(chr)?[0-9]+$/) || ($1 ~ /^(chr)?[XY]$/)) { 
        if ($4 == "0:1") {
            $4 = "NO_COVERAGE"
        } else if ($4 == "1:10") {
            $4 = "LOW_COVERAGE"
        } else if ($4 == "10:150") {
            $4 = "CALLABLE"
        } else if ($4 == "150:inf") {
            $4 = "HIGH_COVERAGE"
        }
        print
    }
}' $sample.lrs.quantized.tsv > lrs.quant.bed

###### Run bedtools to generate BEDs for each coverage comparison set
cat lrs.quant.bed | awk '$4!="NO_COVERAGE"' > lrs.cov.bed
cat lrs.quant.bed | awk '$4=="NO_COVERAGE"' > lrs.nocov.bed
cat srs.quant.bed | awk '$4!="NO_COVERAGE"' > srs.cov.bed
cat srs.quant.bed | awk '$4=="NO_COVERAGE"' > srs.nocov.bed

####### Find overlaps
echo "find overlaps"
#a. Covered by LRS, but not by SRS
bedtools intersect -a srs.nocov.bed -b lrs.cov.bed -wo > inlrs_notsrs.bed
a=$(awk '{sum+=$11} END{print sum}' inlrs_notsrs.bed)
echo "Number of bases where LRS has coverage, but SRS doesn't: $a"

#b. Covered by SRS, but not by LRS
bedtools intersect -a lrs.nocov.bed -b srs.cov.bed -wo > insrs_notlrs.bed
b=$(awk '{sum+=$11} END{print sum}' insrs_notlrs.bed)
echo "Number of bases where SRS has coverage, but LRS doesn't: $b"

#c. Found in both LRS and SRS
bedtools intersect -a lrs.cov.bed -b srs.cov.bed -wo > insrs_inlrs.bed
c=$(awk '{sum+=$11} END{print sum}' insrs_inlrs.bed)
echo "Number of bases where both LRS and SRS has coverage: $c"

#d. Not found in both LRS and SRS
bedtools intersect -a lrs.nocov.bed -b srs.nocov.bed -wo > notsrs_notlrs.bed
d=$(awk '{sum+=$11} END{print sum}' notsrs_notlrs.bed)
echo "Number of bases where neither LRS not SRS has coverage: $d"


####### Get callable coverage stats and save into a TSV
awk '{FS=OFS="\t"} { category_count[$4] += $5 } END { for (category in category_count) print category, category_count[category] }' srs.quant.bed > $sample.srs.coveval.tsv
awk '{FS=OFS="\t"} { category_count[$4] += $5 } END { for (category in category_count) print category, category_count[category] }' lrs.quant.bed > $sample.lrs.coveval.tsv

######## For nosrs_nolrs coverage regions, remove assembly gaps and then get counts
gaps="/private/groups/migalab/shnegi/raredis/coverage_stats/callable_loci/gaps.bed"
# Get intervals and recalculate length
bedtools intersect -a lrs.nocov.bed -b srs.nocov.bed | cut -f1-4 | awk '{FS=OFS="\t"}{print $0, ($3-$2)}' > nocovinboth.bed
# Intersect with assembly gaps bed file
bedtools subtract -a nocovinboth.bed -b $gaps | cut -f1-4 | awk '{FS=OFS="\t"}{print $0, ($3-$2)}' > nocovinboth_outsideGaps.bed
# Get counts
e=$(awk '{sum+=$5} END{print sum}' nocovinboth_outsideGaps.bed)
echo "Number of bases where neither LRS nor SRS has coverage (OUTSIDE GAPs): $e"


####### Generate intersected BED files to view on IGV and stratify SRS-only and LRS-only covered regions
bedtools intersect -a srs.nocov.bed -b lrs.cov.bed | cut -f1-4 | awk '{FS=OFS="\t"}{print $0, ($3-$2)}' > inlrs_notsrs.intersectBED.bed
bedtools intersect -a lrs.nocov.bed -b srs.cov.bed | cut -f1-4 | awk '{FS=OFS="\t"}{print $0, ($3-$2)}' > insrs_notlrs.intersectBED.bed

echo "Finished!!"