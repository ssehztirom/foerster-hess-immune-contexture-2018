####################################################################
## Demonstrates the immune cell marker gene enrichment analysis and the development of gene signatures shown in "The immune contexture of hepatocellular carcinoma predicts clinical outcome"
## V 1.0
## (c) Moritz Hess 2018
#####################################################################
### Load Libraries #####
library(openxlsx)
source("functions.R")
library(gplots)
library(limma)
library(edgeR)
library(survival)
cyanmagenta <- function(x){colorpanel(x,'cyan','white','magenta')}
nrot = 999 ## Increase for improved accuracy of the Enrichment p-values 


#####################################################################
#### TCGA Data ######################################################
#####################################################################
supporting.document <- "41598_2018_21937_MOESM2_ESM.xlsx" # The supplementary xlsx File MOESM2_ESM.xlsx
data <- readWorkbook(supporting.document,sheet=5,startRow = 4,rowNames=TRUE)
samples <- readWorkbook(supporting.document,sheet=3,startRow=3,rowNames=TRUE)
immunegenes <- readWorkbook(supporting.document,sheet=6,startRow=3,rowNames=FALSE)
pheno <- read.csv("pheno.csv",stringsAsFactors=FALSE)
#####################################################################
## Variables ########################################################
#####################################################################
Age = as.numeric(samples$"Age.(in.decades)")
Gender = samples$Gender
Factor = samples$Tissue.type

Replicate=pheno$Replicate
Risk_Factor <- samples$Risk.factor
Risk_Factor[Risk_Factor=='Not applicable'] <- 'AAAA'
Survival <- samples$Survival
Survival[Survival=='Not applicable'] <- 'AAAA'
Stage <- samples$Stage








###########################################################################################
############################ immune cell marker gene enrichment analysis (Figure 2A) ######
###########################################################################################
immunelist <- sapply(unique(immunegenes[,1]),function(x) immunegenes[,2][immunegenes[,1]==x])
immunelistnumbers <- sapply(immunelist,function(x) match(x,rownames(data)))
Heatmap.Order <- c('B cell',
                   'T cell',
                   'T helper cell',
                   'CD8',
                   'Cytotoxic cells',
                   'T effector memory',
                   'T central memory',
                   'TFH',
                   'Th1',
                   'Th2',
                   'Th17',
                   'Tgamma delta',
                   'Treg',
                   'NK',
                   'NK CD56bright',
                   'NK CD56dim',
                   'DC',
                   'iDC',
                   'aDC',
                   'Macrophages',
                   'Mast cells',
                   'Neutrophils',
                   'Eosinophils',
                    'MHC Class I',
                   'Cytolytic Activity',
                   'Type I IFN Reponse',
                   'Type II IFN Reponse',
                   'Co-stimulation, T cell',
                   'Co-stimulation, APC',
                   'Co-inhibition, T cell',
                   'Co-inhibition, APC'
                   )
immunelistnumbers <- immunelistnumbers[Heatmap.Order]



NewLabels <- c(
"B cells",
"T cells",
"T helper cells",
"CD8+ T cells",
"Cytotoxic T cells",
"T effector memory cells",
"T central memory cells",
"T follicular helper cells",
"Th1 cells",
"Th2 cells",
"Th17 cells",
"Tgd cells",
"CD4+ regulatory T cells",
"NK cells",
"NK CD56bright cells",
"NK CD56dim cells",
"Dendritic cells",
"Immature dendritic cells",
"Activated dendritic cells",
"Macrophages",
"Mast cells",
"Neutrophils",
"Eosinophils",
"MHC Class I",
"Cytolytic activity",
"Type I Interferon response",
"Type II Interferon response",
"Co-stimulation, T cell",
"Co-stimulation, APC",
"Co-inhibition, T cell",
"Co-inhibition, APC"
)

## Transform count data using the voom transformation and estimate intra-individual correlation (Random=TRUE). For sake of computational speed, correlation from a previous run is provided here (0.1598).
AdjustNA <- !is.na(Age) & !is.na(Gender) & !is.na(Factor) # Remove samples with missing information
Model <- model.matrix(~Factor[AdjustNA]+Age[AdjustNA] + Gender[AdjustNA]) # specify the model
colnames(Model) <- gsub('[AdjustNA]','',colnames(Model),fixed = TRUE)
Limma <- Limma.Pipe(Counts=data[,AdjustNA],Model=Model, Voom=TRUE, Random=FALSE, Block=Replicate[AdjustNA],Correlation=0.1598) # perform the transformation and the computation of the intra-individual correlation (if Random=TRUE).

## Define the contrasts that are tested.
Contrast.List <- list('HCC-NT vs. HL' = c(0,0,1,0,0),
                      'HCC-T vs. HL' = c(0,1,0,0,0),
                      'HCC-T vs. HCC-NT' = c(0,1,-1,0,0)
                      )

## Perform enrichment analysis using rotation tests:
# Wu, D., Lim, E., Vaillant, F., Asselin-Labat, M. L., Visvader, J. E., & Smyth, G. K. (2010). ROAST: rotation gene set tests for complex microarray experiments. Bioinformatics, 26(17), 2176-2182.

expression <- Limma$Transformed
correlation <- 0.1598
block = Replicate
Contrast.Results <- list()
for (Contrast.Name in names(Contrast.List)) 
{
    Current.Contrast <- Contrast.List[[Contrast.Name]]
    ROAST <- mroast(y=expression,
                    index= immunelistnumbers,
                    design = Model,
                    correlation = correlation, 
                    contrast= Current.Contrast,
                    nrot=nrot, 
                    block = Replicate[AdjustNA]
                    )
    Contrast.Results[[Contrast.Name]] <- ROAST
}
Contrast.Results.ordered <- lapply(Contrast.Results, function(x) x[Heatmap.Order,])

## Format output of ROAST for plotting
ROAST. <- do.call(cbind,Contrast.Results.ordered)
ROAST.P <- ROAST.[,grepl('PValue$',colnames(ROAST.))]
ROAST.Dir <- ROAST.[,grepl('Direction$',colnames(ROAST.))]
ROAST.Dir <- ifelse(ROAST.Dir=='Down',yes = -1, no = 1)
ROAST.Plog <- abs(log10(ROAST.P)) * ROAST.Dir
rownames(ROAST.Plog) <- NewLabels
colnames(ROAST.Plog) <- gsub(".PValue","",colnames(ROAST.Plog))
ROAST.Plog <- ROAST.Plog[,c("HCC-T vs. HCC-NT","HCC-T vs. HL","HCC-NT vs. HL")]

pdf('Figure_2A.pdf',width=10,height=15)
heatmap.2(data.matrix(ROAST.Plog),
          col='cyanmagenta',
          margins = c(10,10),
          srtCol = 90,
          notecex=1,
          key=TRUE,
          Colv =FALSE,
          Rowv = FALSE,
          dendrogram = 'none',
          trace='none',
          density.info = 'none',
          notecol='black',
          symm=FALSE,
          cexCol = 0.8,
          lmat=rbind(c(4,3),c(2,1)), lhei=c(0.1,1),lwid=c(0.3,1),
          cexRow = 0.8,
          main='Cell Type Marker Enrichment')
dev.off()



###########################################################################################
############################# Cox Regression (Figure 4B) ##################################
###########################################################################################

## Format time to event data for cox regression.
.sel <- samples$Tissue.type=="Primary_Tumor"
tod <- as.numeric(pheno[,'death_days_to'])
dead <- ifelse(!is.na(tod),yes=1,no=0)
cens <- as.numeric(pheno[,'last_contact_days_to'])
cens[cens<=0] <- NA
boolvec <- .sel & (!is.na(tod) | !is.na(cens)) & !is.na(Age) & !is.na(Gender) 
survmat = cbind(rownames(samples),Gender,Age,tod,dead,cens)
survmat = survmat[boolvec,]
survmat[is.na(survmat)] <- 0
survmat = data.frame(survmat,row.names=1,stringsAsFactors=FALSE)
survmat$tod <- as.numeric(survmat$tod)
survmat$Age <- as.numeric(survmat$Age)
survmat$cens <- as.numeric(survmat$cens)
survmat$dead <- as.numeric(survmat$dead)

survmat[['time']] <- apply(survmat[,c('cens','tod')],1,sum)
survmat[['surv']] <- Surv(survmat[['time']],survmat[['dead']])

## Voom transform expression data.
nf <- calcNormFactors(data)
expressionmat <- voom(data,lib.size = colSums(data)*nf,plot=TRUE)
expressionmat = expressionmat$E[,rownames(survmat)]


## Perform univariate time-to-event analysis using cox-regression.
Formula <- as.formula('surv ~ Expression + Age + Gender') # specify the model
cox.univariate.results <- cellcox(immunelistnumbers,expressionmat,survmat,Formula) # estimate the coefficients for each gene (univariate models)
cox.univariate.z <- lapply(cox.univariate.results, function (x) sapply(x, function(y) summary(y)[['coefficients']]['Expression',4])) # extract z-score for each gene from univariate models.
names(cox.univariate.z) <- NewLabels
cox.univariate.z <- cox.univariate.z[sapply(cox.univariate.z,length)>0]

## Calculate average z-scores per gene-set and calculate the empirical null distribution for each each gene set.
set.seed(1234)
nperm = 10000
a = permcelltest(cox.univariate.z,nperm)
a[a==0] = 1/nperm # replace an empirical p-value of zero by 1 / number of simulations.
plotzscoreenrich(rev(cox.univariate.z),rev(a),"Figure_4B") # plot the distribution of z-scores and the empirical p-values.





###########################################################################################
############################# Immune cell signature (Figure 5A) ##########################
###########################################################################################

## Infer survival in patients that differ in the expression of immune cell-type specific marker genes that were significantly associated with survival (Figure 4B).


scoregenes <- unlist(immunelist[c('T cell','Cytotoxic cells','Th2','Macrophages')]) # Fetch the immune cell-type specific marker genes that were significantly associated with survival.

## For each immune cell-type specific marker gene, determine the association with survival using cox proportional hazards models.
Formula <- as.formula('surv ~ Expression + Age + Gender')
cox.univariate.all.results <- unicoxexpression(expressionmat[scoregenes,],survmat,Formula)
cox.univariate.z.all <- unlist(lapply(cox.univariate.all.results,  function(y) summary(y)[['coefficients']]['Expression',4]))
scoregenes <- cox.univariate.z.all


## Calculate a score for each patient based on the expression of immune cell-type specific marker genes. The higher the score, the longer the expected survival.
thres <- "median"
direction <- "one"
addone <- TRUE
myscore <- score.universal(scoregenes,expressionmat,thres=thres,direction=direction,addone = addone)
## Dichotomize patients based on score.
mygroups <- ifelse(myscore>median(myscore),yes='good',no='poor')
.surv <- survmat$surv


fit <- survfit(.surv~mygroups)
fit$time <- fit$time / 365

## Plot the Kaplan - Meier curves of the two inferred patient groups.
pdf('Figure_5A.pdf')
plot(fit,col=c('blue','green'),xlim=c(0,7),xlab="Time to death (years)",ylab="Survival probability")
dev.off()

