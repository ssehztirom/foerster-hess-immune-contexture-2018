####################################################################
## Helper functions and algorithms for analysis of data shown in "The immune contexture of hepatocellular carcinoma predicts clinical outcome"
## V 1.0
## (c) Moritz Hess 2018

#####################################################################
score.universal <- function(Genes,mat,thres="median",direction="both",addone=TRUE)
{
    ## Computes a score for each sample (individual) based on the expression of selected genes and their association with survival. The higher the score, the longer is the expected survival.
    ## "Genes" is a Named Vector, indicating the gene names and the ratio of the proportional hazard: 1 = increasing hazard probability with increasing expression, -1 = decreasing hazard probability with increasing expression.
    ## "mat" is an expression matrix (rows=genes,columns=samples); approximately normally distributed eg. log transformed
    ## "score" is the returned score per individual 
    
    GeneNames <- names(Genes)
    genesdir <- Genes
    mat <- t(scale(t(mat)))
    .i <- intersect(GeneNames,rownames(mat))
    mat <- mat[.i,]
 
    score <- rep(0,ncol(mat))
    for (i in .i)
        
    {
        
        
        .med <- median(mat[i,])
        .med <- rep(.med,4)
        if (thres !="median")
        {
            .med <- quantile(mat[i,])
        }
            
        if (sign(genesdir[i]) ==1)
        {
            score.log <- mat[i,]< .med[2]
       
            
        }
        else if (sign(genesdir[i]) ==-1)
        {
            score.log <- mat[i,]> .med[4]
          
        }
        increaseby <- 1
        if (!addone)
        {
            increaseby <- genesdir[i]
        }
        score[score.log] <- score[score.log] + abs(increaseby)
        if (direction=="both")
        {
            score[!score.log] <- score[score.log] - abs(increaseby)
        }
    }
    score
}





permcelltest <- function(unip,n,externaluniverse=NULL)
{
    ## Computes the probability of an average z-score of a gene set. 
    ## "unip" is a named list containing z-scores from e.g. a cox model. Each list entry is a named vector containing the z-scores of a gene set. Gene sets have to be mutually exclusive. "n" determines the number of conducted simulations to retrieve the empirical null distribution of observing a particular average z-score per gene-set. In each simulation (S_i), for each gene set, a random number of genes (z-scores) of equal size as the gene set are sampled from all available z-scores (across all gene sets) and average z-scores are computed per gene-set. The number of simulated average absolute z-scores, exceeding the absolute value of the empirical averaged z-score per gene-set is the p-value for the respective z-score.
    gene.universe <- unlist(unip)
    if (!is.null(externaluniverse))
    {
        gene.universe = externaluniverse
    }
    
    unip.length <- sapply(unip,length)
    nrep <- n
    .bg <- lapply(unip.length, function(x) sapply(1:nrep,function(y) mean(sample(gene.universe,x))))
    .obs <- sapply(unip, function(x) mean(x))
    .obs.p.r <- sapply(names(.obs), function(x) length(which(.bg[[x]] >.obs[[x]]))/nrep)
    .obs.p.l <- sapply(names(.obs), function(x) length(which(.bg[[x]] <.obs[[x]]))/nrep)
    .obs.p <- apply(cbind(.obs.p.r,.obs.p.l),1,min)
    return (.obs.p)
}

plotzscoreenrich <- function(unip,enrichp)
{
    ## Plots the distribution of z-scores per gene-set ("unip") and the p-values from "permcelltest" ("enrichp").
    xx <- unip
    xx <- lapply(xx,function(x) x * -1)
    xx.v <- enrichp[names(unip)]
    xx.vp <- p.adjust(xx.v,method="bonferroni")
    xx.vp <- sapply(xx.vp, function(x) ifelse(x<0.05,yes="*", no=""))
    xx.v <- sapply(xx.v, function(x) sprintf("%.4f",x))
    # pdf(sprintf('%s.pdf',fnamesuffix))
    .min <- min(unlist(xx))
    .max <- max(unlist(xx))           
    plot(1000,-1000,xlim=c(.min-4,.max+4),ylim=c(0.5,length(xx)+0.5),ylab='Immune Cell Type',xlab= ' - Beta/SE(Beta),  Cox proportional hazards model ',axes=FALSE)
    axis(1)
    box <- boxplot(xx,horizontal = TRUE,add=TRUE,axes=FALSE,border='black',outline=FALSE)#,boxwex=0.7)
    parc <- par()
    .min <- parc$usr[1]
    .max <- parc$usr[2]
    avg. <- sapply(xx,mean)
    .name.pos <- ifelse(avg. > 0, yes =4,no=2)
    .digit.side <- ifelse(avg. > 0, yes =2,no=4)
    .number.pos <- ifelse(avg. > 0, yes=5,no=1)
    .number.pos. <- c()
    for (i in 1:length(.number.pos))
    {
        .number.pos. <- c(.number.pos.,box$stats[.number.pos[i],i])
    }
    text(y=1:length(xx),x=.number.pos.,pos=.name.pos,labels=names(xx),col='grey')
    .number.pos <- ifelse(avg. > 0, yes=1,no=5)
    .number.pos. <- c()
    for (i in 1:length(.number.pos))
    {
        .number.pos. <- c(.number.pos.,box$stats[.number.pos[i],i])
    }
    .lab <- apply(cbind(box$n,xx.v[box$names],xx.vp[box$names]),1,paste,collapse=" | ")
    text(y=(1:length(xx))+0.25,x=.number.pos.,labels=.lab,pos=.digit.side,cex=0.7,col='grey')
    par(xpd=FALSE)
    abline(v=0)

    # dev.off()
}


###########################
## Convenience Functions ##
###########################
unicoxexpression <- function(mat,survtable,formula)
{
    require(survival)
    cox.univariate.results <- list()
    for (id in rownames(mat))
    {
        cur.df <- data.frame(Expression = as.numeric(mat[id,]),survtable)
        
        Mod <- coxph(formula,data=cur.df)
        cox.univariate.results[[id]] <- Mod
    }
    cox.univariate.results
}

cellcox <- function(immunelist,mainmat,cox.table,Formula)
{

    cox.univariate.results <- list()
    system.time(
        
        for ( cell.type in names(immunelist))
        {
            cell.type.ids <- immunelist[[cell.type]]
            mat <- mainmat[cell.type.ids,]
            celltyperesults <- unicoxexpression(mat=mat,survtable = cox.table,formula=Formula)
            cox.univariate.results[[cell.type]] <- celltyperesults
        }
    )
    return(cox.univariate.results)
}

Limma.Pipe <- function(Counts,Model,Block=Replicate,Voom=FALSE, Random=FALSE,Correlation=0){
    library(limma)
    library(edgeR)
    
    {
        if ( Voom )
            {

                nf <- calcNormFactors(Counts)
                Self.nointeract <- voom(Counts,Model,lib.size = colSums(Counts)*nf,plot=TRUE)
            }
        else
            {
                Self.nointeract <- Counts
            }
    }
    {
        if ( Random )
            {
                corfit <- duplicateCorrelation(Self.nointeract, ndups=1, block = Block)
            }
        else
            {
                corfit <- list(consensus.correlation=Correlation)
            }
    }
    Self.nointeract.fit <- lmFit(Self.nointeract,Model,correlation = corfit$consensus.correlation,block=Block)
    Self.nointeract.fit <- eBayes(Self.nointeract.fit)
    Output <- list(Transformed = Self.nointeract, Model = Self.nointeract.fit, Correlation = corfit)
    return(Output)
}

