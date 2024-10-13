#!/bin/bash

'''
Author: Shloka Negi, shnegi@ucsc.edu
Purpose: This script generates callable coverage statistics for SRS and LRS BAMs aligned to T2T-CHM13 and identifies exclusive coverage regions for each technology.
Input file requirements: Requires input LRS and SRS BAM/CRAM (aligned to T2T-CHM13, sorted and indexed)
Usage: ./coveval_t2t.sh -l lrs.bam -s srs.bam
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

sample=$(basename $sr_cram | awk -F "_" '{print $1}')
REF="chm13v2.0.fa"

mkdir -p $sample/chm13

######## Run mosdepth
echo "run mosdepth"
cd $sample/chm13
#1. SRS
mosdepth --threads 64 -f $REF --use-median --fast-mode --mapq 10 --quantize 0:1:10:80: $sample.srs.quantized $sr_cram
zcat $sample.srs.quantized.quantized.bed.gz | awk '{FS=OFS="\t"}{print $0, ($3-$2)}' > $sample.srs.quantized.tsv 
awk 'BEGIN {FS=OFS="\t"} {
    if (($1 ~ /^(chr)?[0-9]+$/) || ($1 ~ /^(chr)?[XY]$/)) { 
        if ($4 == "0:1") {
            $4 = "NO_COVERAGE"
        } else if ($4 == "1:10") {
            $4 = "LOW_COVERAGE"
        } else if ($4 == "10:80") {
            $4 = "CALLABLE"
        } else if ($4 == "80:inf") {
            $4 = "HIGH_COVERAGE"
        }
        print
    }
}' $sample.srs.quantized.tsv > srs.quant.bed
cat srs.quant.bed | awk '$4=="CALLABLE"' > srs.cov.bed
cat srs.quant.bed | awk '$4=="NO_COVERAGE"' > srs.nocov.bed

#2. LRS
echo "run mosdepth"
mosdepth --threads 64 -f $REF --use-median --fast-mode --mapq 10 --quantize 0:1:10:80: $sample.lrs.quantized $lr_bam
zcat $sample.lrs.quantized.quantized.bed.gz | awk '{FS=OFS="\t"}{print $0, ($3-$2)}' > $sample.lrs.quantized.tsv 
awk 'BEGIN {FS=OFS="\t"} {
    if (($1 ~ /^(chr)?[0-9]+$/) || ($1 ~ /^(chr)?[XY]$/)) { 
        if ($4 == "0:1") {
            $4 = "NO_COVERAGE"
        } else if ($4 == "1:10") {
            $4 = "LOW_COVERAGE"
        } else if ($4 == "10:80") {
            $4 = "CALLABLE"
        } else if ($4 == "80:inf") {
            $4 = "HIGH_COVERAGE"
        }
        print
    }
}' $sample.lrs.quantized.tsv > lrs.quant.bed

cat lrs.quant.bed | awk '$4=="CALLABLE"' > lrs.cov.bed
cat lrs.quant.bed | awk '$4=="NO_COVERAGE"' > lrs.nocov.bed


# Callable by LRS, but NO_COVERAGE in SRS
bedtools intersect -a srs.nocov.bed -b lrs.cov.bed -wo > inlrs_notsrs.bed
a=$(awk '{sum+=$11} END{print sum}' inlrs_notsrs.bed)
echo "Number of bases where LRS has callable coverage, but SRS has no coverage: $a"

# Callable by SRS, but NO_COVERAGE in LRS
bedtools intersect -a lrs.nocov.bed -b srs.cov.bed -wo > insrs_notlrs.bed
b=$(awk '{sum+=$11} END{print sum}' insrs_notlrs.bed)
echo "Number of bases where SRS has callable coverage, but LRS has no coverage: $b"

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

####### Generate intersected BED files to view on IGV and stratify SRS-only and LRS-only covered regions
bedtools intersect -a srs.nocov.bed -b lrs.cov.bed | cut -f1-4 | awk '{FS=OFS="\t"}{print $0, ($3-$2)}' > inlrs_notsrs.callable.intersectBED.bed
bedtools intersect -a lrs.nocov.bed -b srs.cov.bed | cut -f1-4 | awk '{FS=OFS="\t"}{print $0, ($3-$2)}' > insrs_notlrs.callable.intersectBED.bed
