---
layout: page
title: Tutorial
description: ~
---
This page provides a tutorial for TWAS fine-mapping using FABIO. Before runing the example code, make sure that the FABIO software is installed successfully. For instructions on installation, please see the [Installation section](https://superggbond.github.io/FABIO/documentation/02_Installation.html).

## FABIO
The example data for FABIO tutorial can be downloaded in this [page](https://superggbond.github.io/FABIO/documentation/03_Data.html). Here are the details about the input data formats and FABIO commands. 
### 1. Formats of input data for FABIO
* Predicted GReX matrix: We require the predicted GReX matrix of the TWAS cohort built up using standard softwares like [PredXican](https://github.com/hakyimlab/MetaXcan) or [BSLMM](https://github.com/genetics-statistics/GEMMA). The input GReX matrix will have gene names as the first column, with the following columns of preicted GReX at individual-level. Each column represents GReX of genes for an individual.
* Binary phenotypes: We also require the observed binary phenotypes of the TWAS cohort. The input phenotypes should be formatted as a plain text file, coded '1' for case and '0' for control under a column named 'pheno'. (see [example](https://github.com/superggbond/FABIO/blob/main/data/pheno.txt))

### 2. Running FABIO
When four sets of GWAS summary statistics are available (one for overlapped individuals and one for non-overlapped individuals for each trait), and we designate trait 1 as the target trait, the PGS construction for the target trait can be performed using the following command
```
workdir=/your/mtPGS/directory #specify the mtPGS directory
${workdir}/src/mtPGS --summstat_int ${workdir}/data/summstat/summstat_trait_1_int.assoc.txt ${workdir}/data/summstat/summstat_trait_2_int.assoc.txt \
--summstat_ext ${workdir}/data/summstat/summstat_trait_1_ext.assoc.txt ${workdir}/data/summstat/summstat_trait_2_ext.assoc.txt \
--n_s 7000 --n_ext 3000 3000 --block ${workdir}/data/EUR_LD_Block.txt --target 0 --ref ${workdir}/data/ref --mafMax 0.8 \
--vg ${workdir}/data/v_g.txt --ve ${workdir}/data/v_e.txt --output trait_1_target_beta trait_2_relevant_beta
```
The essential inputs are:
- summstat_int: specify the GWAS summary statistics computed based on overlapped individuals.
- summstat_ext: specify the GWAS summary statistics computed based on non-overlapped individuals.
- n_s: the sample size of overlapped individuals
- n_ext: the sample size of non-overlapped individuals
- block: specify the LD block information.
- target: index for the target trait (0 if the first trait is the target trait, and 1 if the second trait is the target trait)
- ref: specify the prefix of reference panel.
- mafMax: specify the maximium of the allele frequency difference between reference panel and summary data.
- vg: specify the directory of genetic variance components file.
- ve: specify the directory of environmental variance components file.
- output: specify the prefix of output files.

### 3. FABIO output
When only two sets of summary statistics are available (one for each trait), we perform mtPGS analysis using the following command
```
workdir=/your/mtPGS/directory #specify the mtPGS directory
${workdir}/mtPGS --summstat ${workdir}/data/summstat/summstat_trait_1.assoc.txt ${workdir}/data/summstat/summstat_trait_2.txt \
--n 10000 10000 --block ${workdir}/data/block.txt --target 0 --ref ${workdir}/data/ref --mafMax 0.8 --vg v_g.txt --ve v_e.txt \
--output trait_1_target_beta trait_2_relevant_beta
```
Here, instead of summstat_int, summstat_ext and n_s, n_ext, we use
- summstat: specify the GWAS summary statistics for the two traits
- n: specify sample sizes for the two traits
