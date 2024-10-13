#!/bin/bash
'''
Author: Shloka Negi, shnegi@ucsc.edu
Purpose: This script runs sveval (https://github.com/jmonlong/sveval) for comparing SVs called by Hapdiff (assembly-based LRS SV caller) or Sniffles (reference-based LRS SV caller) and GATK-SV (SRS SV caller)
Usage: ./run_sveval.sh
'''

# Set the exit code of a pipeline to that of the rightmost command
# to exit with a non-zero status, or zero if all commands of the pipeline exit
set -o pipefail
# cause a bash script to exit immediately when a command fails
set -e
# cause the bash shell to treat unset variables as an error and exit immediately
set -u

pwd; hostname; date

samples_file="sample_ids.txt"

##########################################################################################
################## Intersection of Hadiff LRS VCFs and SRS VCFs ################
##########################################################################################

cd hapdiff_sveval_out

pwd; hostname; date
for sample in $(cat $samples_file);
do
    lrs_vcf="hapdiff_unphased.vcf.gz"
    srs_vcf="$sample.vcf.gz"
    # Filter SRS vcf
    zcat $srs_vcf | bcftools filter -e '(GT="0/0"| GT="0|0"| GT="./.")' | bcftools filter -e 'QUAL<=1' | bcftools filter -e 'ALGORITHMS="WHAM" & GT="1/1"' -o $sample.srs.vcf.gz -Oz
    tabix -p vcf $sample.srs.vcf.gz
    # Filter LRS vcf
    zcat $lrs_vcf | bcftools norm -m -any | bcftools filter -i '(abs(SVLEN)>=50) || ((STRLEN(REF)>=50 && STRLEN(ALT)>=50) && SVTYPE=="INV") || (SVTYPE=="BND")' | bcftools filter -e '(GT="0/0"| GT="0|0"| GT="./.")' -o $sample.lrs.vcf.gz -Oz
    tabix -p vcf $sample.lrs.vcf.gz
    # Run seveval
    echo "running sveval on sample $sample"
    Rscript sveval/svevalOl.R $sample.lrs.vcf.gz $sample.srs.vcf.gz $sample'_sveval_out.tsv'
    # remove filtered files
    rm -rf $sample.srs.vcf.gz* $sample.lrs.vcf.gz*
done
echo "done running sveval on all samples"

echo "generating json with comparison stats for all samples"
cd ../
python3 gen_sveval_out_json.py hapdiff_sveval_out/ hapdiff_svevals.json Hapdiff
echo "Finished!!"

##########################################################################################
################## Intersection of Sniffles LRS VCFs and SRS VCFs ################
##########################################################################################

cd sniffles_sveval_out

for sample in $(cat $samples_file);
do
    lrs_vcf="sniffles.vcf.gz"
    srs_vcf="$sample.vcf.gz"
    # Filter SRS vcf
    zcat $srs_vcf | bcftools filter -e '(GT="0/0"| GT="0|0"| GT="./.")' | bcftools filter -e 'QUAL<=1' | bcftools filter -e 'ALGORITHMS="WHAM" & GT="1/1"' -o $sample.srs.vcf.gz -Oz
    tabix -p vcf $sample.srs.vcf.gz
    # Filter LRS vcf
    zcat $lrs_vcf | bcftools norm -m -any | bcftools filter -i 'abs(SVLEN)>=50 || SVTYPE=="BND"' | bcftools filter -e '(GT="0/0"| GT="0|0"| GT="./.")' -o $sample.lrs.vcf.gz -Oz
    tabix -p vcf $sample.lrs.vcf.gz
    # Run seveval
    echo "running sveval on sample $sample"
    Rscript ../svevalOl.R $sample.lrs.vcf.gz $sample.srs.vcf.gz $sample'_sveval_out.tsv'
    # remove filtered files
    rm -rf $sample.srs.vcf.gz*
    rm -rf $sample.lrs.vcf.gz*
done
echo "done running sveval on all samples"

echo "generating json with comparison stats for all samples"
cd ../
python3 gen_sveval_out_json.py sniffles_sveval_out/ sniffles_svevals.json Sniffles
echo "Finished!!"
