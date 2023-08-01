# Re-analysis of data provided by Gihawi et al. 2023 bioRxiv ([link](https://www.biorxiv.org/content/10.1101/2023.07.28.550993v1))

This repository uses the raw data on 3 TCGA cancer types provided by Gihawi et al. (bladder [BLCA], breast [BRCA], head and neck [HNSC]) to evaluate if microbial cancer type-specific differences persist in their data. For reference, their data were downloaded from the Github link ([https://github.com/yge15/Cancer_Microbiome_Reanalyzed](https://github.com/yge15/Cancer_Microbiome_Reanalyzed)) listed in the Methods section.

Briefly, the following was done:
- Their data were loaded alongside metadata from the original [Poore et al. 2020 Nature](https://www.nature.com/articles/s41586-020-2095-1) paper
- The 3 cancer type tables containing raw counts were merged. Conservatively, we retained only overlapping genera between the 3 tables. However, we removed all counts associated with the _Homo_ genus (i.e., human data) to focus just on the microbial differences. As noted by Gihawi et al., their database for KrakenUniq only contained complete genomes of microbes.
- We separated the count and metadata into subsets with a single sequencing center (e.g., Harvard Medical School [HMS]), sequencing platform (Illumina HiSeq), and experimental strategy (WGS) in order to avoid needing to apply batch correction. All subsequent analyses were based on these raw data subsets without batch correction.
- We applied beta and alpha diversity analyses to these subsets using [Qiime 2](https://qiime2.org/) (see code within the `Qiime` folder). To avoid rarefaction for the beta diversity analyses, we employed a robust Aitchison approach called RPCA ([Martino et al. 2019](https://journals.asm.org/doi/10.1128/msystems.00016-19)). We also applied more traditional beta-diversity measures (e.g., Bray-Curtis, Jaccard) and calculated alpha diversity using data rarefied to the 1st quartile of sample read counts for each subset.
- After finding that the cancer type separation was significant, we performed mutli-class and two-class machine learning within HMS primary tumor (PT) samples and HMS blood derived normal (BDN) samples.

The figure summarizing the findings and its associated caption is below. All analyses necessary to generate these plots have code contained within this repository or can be done without code using [https://view.qiime2.org/](https://view.qiime2.org/).

![alt text](https://github.com/gregpoore/tcga_rebuttal/tree/master/Figures/assembled-figure.jpg?raw=true)
