---
layout: full
homepage: true
disable_anchors: true
description: Fine-mApping of causal genes for BInary Outcomes
---
## FABIO
![FABIO\_pipeline](FABIO_scheme.png)
FABIO is a TWAS fine-mapping method that relies on a probit model to directly relate multiple genetically regulated gene expression (GReX) to binary outcome in TWAS fine-mapping. Additionally, it jointly models all genes located on the same chromosome to account for the correlation among GReX arising from cis-SNP LD and expression correlation across genomic regions. Through a Markov chain Monte Carlo (MCMC) algorithm, it obtains the posterior inclusion probability (PIP) of having a non-zero effect for each gene, as a measure of evidence for the gene’s association with the outcome trait.

### Example Analysis with FABIO: [here](https://superggbond.github.io/FABIO/documentation/04_FABIO_Example.html).
