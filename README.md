# FABIO

## Fine-mApping of causal genes for BInary Outcomes (FABIO)
![FABIO](FABIO_scheme.png)
FABIO is a TWAS fine-mapping method that relies on a probit model to directly relate multiple genetically regulated gene expression (GReX) to binary outcome in TWAS fine-mapping. Additionally, it jointly models all genes located on the same chromosome to account for the correlation among GReX arising from cis-SNP LD and expression correlation across genomic regions. Through a Markov chain Monte Carlo (MCMC) algorithm, it obtains the posterior probability of having a non-zero effect for each gene, also known as the posterior inclusion probability (PIP). PIP serves as an important measure of evidence for the gene’s association with the outcome trait. 

## How to use FABIO

See [Tutorial](https://superggbond.github.io/FABIO/) for detailed documentation and examples.

## Issues

For any questions or feedbacks on FABIO software, please email to Haihan Zhang (hhzhang@umich.edu).

## Citation

Haihan Zhang, Kevin He, Lam C. Tsoi, and Xiang Zhou#. FABIO: TWAS Fine-mapping to Prioritize Causal Genes for Binary Traits.


