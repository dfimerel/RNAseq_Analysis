# This script performs differential expression analysis as a continuation from the RNAseq-pipeline.R script
# The user should define a number of things, including:
# 1. the working directory
# 2. the conditions to compare
## Some of the code was obtained from: https://github.com/trinityrnaseq/Griffithlab_rnaseq_tutorial/blob/master/scripts/Tutorial_Module4_Part4_edgeR.R


# run this script by: Rscript DiffExpression-EdgeR.R 


#!/usr/bin/env Rscript


###########
#Libraries#
###########

suppressMessages(library(edgeR))
suppressMessages(library(gplots))


############################################
#Tool/data directories for USER to define  #
############################################

workDir= [add path to working directory]


#############
 #Load data #
#############


# Read in the final table with the raw counts of all samples
rawcounts=read.table(paste(workDir,"read_counts_all/raw_counts_all.txt",sep=""), header=TRUE, stringsAsFactors=FALSE)


# Read in gene mapping between ENSEMBL and gene IDs
mapping=read.table(paste(workDir,"ENSG_ID2Name.txt",sep=""), header=FALSE, stringsAsFactors=FALSE,row.names=1)


##############
# Run edgeR  #
##############

# Make condition labels * TO BE DEFINED BY USER *
sampleCondition=c(rep("tumor",10),rep("normal", 10))


# Create color list for later plots
colorCondition=c(rep("darkred",10),rep("darkgreen", 10))


# Require at least 25% of samples to have count > 25
quant <- apply(rawcounts,1,quantile,0.75)
keep <- which((quant >= 25) == 1)
rawcounts <- rawcounts[keep,]


# Make DGEList object and normalize
y <- DGEList(counts=rawcounts, group=factor(sampleCondition))
y <- calcNormFactors(y)


# Unsupervised clustering of samples and MDS plot
lcpm=cpm(y,log=TRUE)
plotMDS(lcpm,col=colorCondition)

# Boxplot of cpm 
boxplot(lcpm,col=colorCondition,las=2,ylab="counts per million (log2)")



# Estimate dispersion and Differential expression test
y <- estimateCommonDisp(y,factor(sampleCondition))
y <- estimateTagwiseDisp(y)
et <- exactTest(y)

# Final table
finalTable <- topTags(et,n=nrow(et$table))$table


# Get DE genes with FDR <= 0.05
mat <- finalTable[finalTable$FDR<=0.05,]
mat <- mat[order(mat$logFC,decreasing=T),]
mat[,"geneNames"]=mapping[rownames(mat),1]


# Save table
write.table(mat, file=paste(workDir,"DE_genes.txt",sep=""), quote=FALSE, row.names=FALSE, sep="\t")


# Heatmap of 50 top DE genes
nn=rownames(finalTable[order(finalTable$FDR),])[1:50]
deCounts=lcpm[nn,]
deCounts <- t(scale(t(deCounts)))
hc=hclust(as.dist(1-cor(deCounts,method="spearman")),method="ward.D")
hr=hclust(as.dist(1-cor(t(deCounts),method="spearman")),method="ward.D")
heatmap.2(deCounts,col=bluered,scale="none",trace="none",Rowv=as.dendrogram(hr),Colv=as.dendrogram(hc))



######################
# Finish and quit

cat ("The pipeline is over!")
q(save='no')
