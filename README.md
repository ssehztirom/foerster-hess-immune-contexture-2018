# Code used in "The immune contexture of hepatocellular carcinoma predicts clinical outcome"
This GitHub repository contains the code to conduct the immune cell marker gene enrichment analysis and the inference of gene signatures shown in "The immune contexture of hepatocellular carcinoma predicts clinical outcome".
"script.R" contains examples how to use the code included in "functions.R"

## Immune cell marker gene enrichment

[limma](https://www.bioconductor.org/packages/release/bioc/html/limma.html)

[ROAST](https://academic.oup.com/bioinformatics/article/26/17/2176/200022)

## Development of signatures to predict survival based on gene expression
Here a score is computed for each sample (individual) based on the expression of selected genes and their association with survival. The higher the score, the longer is the expected survival. Gene-wise association with survival is detected using univariate proportional hazards models.
