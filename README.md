# Code used in "The immune contexture of hepatocellular carcinoma predicts clinical outcome" - Foerster et al. Scientific Reports 2018
This GitHub repository contains the code to conduct the immune cell marker gene enrichment analysis and the inference of gene signatures shown in "The immune contexture of hepatocellular carcinoma predicts clinical outcome".
"script.R" contains examples how to use the code included in "functions.R". Expression data shown in Supplementary Table 5, which is employed here, can be found in the supplementary electronic materials of the publication.

## Immune cell marker gene enrichment
Immune cell marker gene enrichment is computed based on immune cell type-specific marker genes ([Bindea et al. 2013](http://www.cell.com/immunity/fulltext/S1074-7613(13)00437-8)) and gene set enrichment analyses using rotation tests ([ROAST](https://academic.oup.com/bioinformatics/article/26/17/2176/200022)), implemented in [limma](https://www.bioconductor.org/packages/release/bioc/html/limma.html).


## Development of signatures to predict survival based on gene expression
Here a score is computed for each sample (individual) based on the expression of selected genes and their association with survival. The higher the score, the longer is the expected survival. Gene-wise association with survival is detected using univariate proportional hazards models.


