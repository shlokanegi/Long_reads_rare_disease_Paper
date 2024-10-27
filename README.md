# Long reads rare disease Paper

Collection of scripts used for genomic analyses described in the "Advancing long-read nanopore genome assembly and accurate variant calling for rare disease detection" paper.

## 1. Genome Completeness Analysis

This folder includes scripts to evaluate genome completeness based on long-read ONT sequencing and SRS Illumina sequencing data.

- `coveval_grch38.sh` and `coveval_t2t.sh`: These scripts generate callable coverage statistics for short-read sequencing (SRS) and long-read sequencing (LRS) BAMs aligned to the GRCh38 and T2T-CHM13v2.0 reference genomes, respectively. They also identify exclusive coverage regions for each technology.

USAGE:
```
## For GRCh38
./coveval_grch38.sh -l lrs.bam -s srs.bam

## For T2T-CHM13
./coveval_grch38.sh -l lrs.bam -s srs.bam
```

- `gene_coverage.sh`: Commands for generating gene-coverage overlap BED files, which are used for analysis of gene regions covered exclusively by each sequencing technology.

- `databases/`: Contains relevant files and databases required for genome completeness analysis.

- `Plotting-scripts/`: Includes R scripts for generating Figure 2 for the paper:
    * `gene_overlap_script.R`: Analyzes and plots gene overlap data.
    * `karyotype_plotting.R`: Creates karyotype visualizations of coverage data.


## 2. Comparison of LRS and SRS variants

This folder contains scripts to perform variant comparisons between Long-Read Sequencing (LRS) and Short-Read Sequencing (SRS).

- `smallvars_analyses_scripts/`: Folder includes scripts used to perform small variant comparisons.

    * `exclusive_FAVs_LRS_stratified.R`: An R script for preprocessing and analyzing LRS-only VCF files generated using `rtgtools vcfeval` focusing on small variants exclusively identified by LRS (Figure 3B).

    * `exclusive_FAVs_LRS_stratified.R`: An R script for preprocessing and analyzing SRS-only VCF files produced via `rtgtools vcfeval`, targeting small variants exclusively identified by SRS (Supplementary Figure 3A).

    * `FAV_intersection_plot.ipynb`: A Jupyter notebook for generating intersection statistics of functionally-annotated variants (FAVs). It first creates an input JSON from `rtgtools vcfeval` outputs using the `vcfeval-intersection-json.py` script, then produces the intersection plot (Figure 3A).

- `run_sveval.sh`: This script runs [sveval](https://github.com/jmonlong/sveval) to compare Structural Variants (SVs) identified by Hapdiff (assembly-based LRS SV caller), Sniffles (reference-based LRS SV caller), and GATK-SV (SRS SV caller).

- `SV_analyses_scripts/`: Includes Python scripts to generate SV intersection statistics from `sveval` outputs, such as `gen_sveval_out_json.py` and `SV_intersection_plot.py` (Figure 5B).

- `databases/`: Contains files and databases essential for variant comparison analysis.


## 3. Phase Block Stats

This folder includes scripts for generating phase-block stats using LRS harmonized VCFs.

- `pbStats.R`: Includes functions to plot raw sequencing statistics post-basecalling, calculate the phase-block NG50, and perform gene-phase block overlap analysis.


## 4. CSS1 Episignature Analysis

This folder includes scripts and data used for performing CSS1 episignature analysis

- `data/`: Contains CpG coordinate BED files for CSS1 episignatures, both before and after liftover from hg19 to hg38. These CpG regions were established as known episignatures for CSS1, as provided by [Erfan Aref-Eshghi et al., 2018](https://www.nature.com/articles/s41467-018-07193-y)

- `scripts/`: Brief description about `extract_sites.sh` and `sb_heatmap.py` scripts.
