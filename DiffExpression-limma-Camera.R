# This script performs differential expression analysis using limma-voom and camera
# The user should define a number of things, including:
# 1. the file where the raw counts are
# 2. the conditions to compare


# run this script by: Rscript DiffExpression-limma-Camera.R ......


#!/usr/bin/env Rscript


###########
#Libraries#
###########

suppressMessages(library(edgeR))

############################################
#Tool/data directories for USER to define  #
############################################

workDir= [add path to working directory]

#############
# Functions #
#############

# function by https://github.com/jkim208/scRNAseq_analysis
readGMT <- function(inputFile){
  con <- file(inputFile, open = "r")
  dataList <- list()
  while (length(oneLine <- readLines(con, n = 1, warn = FALSE)) > 0) {
    myVector <- unlist(strsplit(oneLine, "\t"))
    dataList <- c(dataList, list(myVector[3:length(myVector)]))
    names(dataList)[length(dataList)] <- myVector[1]
  }
  close(con)
  return(dataList)
}



#############
 #Load data #
#############

# Read in the final table with the raw counts of all samples
rawcounts=read.table(paste(workDir,"read_counts_all/raw_counts_all.txt",sep=""), header=TRUE, stringsAsFactors=FALSE)


# Read in gene mapping between ENSEMBL and gene IDs
mapping=read.table(paste(workDir,"ENSG_ID2Name.txt",sep=""), header=FALSE, stringsAsFactors=FALSE,row.names=1)



###################
# Run limma-voom  #
###################

# Make condition labels * TO BE DEFINED BY USER *
sampleCondition=c("tumor","normal","tumor","normal","tumor","normal","tumor","tumor","normal","tumor","normal","tumor","normal","tumor","normal","tumor","normal")


# Create color list for later plots
colorCondition=sampleCondition
colorCondition[colorCondition=="tumor"]="red"
colorCondition[colorCondition=="normal"]="green"


# Require at least 25% of samples to have count > 25
quant <- apply(rawcounts,1,quantile,0.75)
keep <- which((quant >= 25) == 1)
rawcounts <- rawcounts[keep,]


# Make DGEList object and normalize
y <- DGEList(counts=rawcounts, group=factor(sampleCondition))
y <- calcNormFactors(y)
design <- model.matrix(~sampleCondition)


# Unsupervised clustering of samples
lcpm=cpm(y,log=TRUE)
plotMDS(lcpm,col=colorCondition)

# Boxplot of cpm 
boxplot(lcpm,col=colorCondition,las=2,ylab="counts per million (log2)")



# Voom transformation
v <- voom(y,design)
fit <- lmFit(v,design)
fit <- eBayes(fit)
deg <- topTable(fit,coef=ncol(design), n=nrow(v))


# Get DE genes with FDR <= 0.05
mat <- deg[deg$adj.P.Val<0.05,]
mat <- mat[order(mat$logFC,decreasing=T),]
mat[,"geneNames"]=mapping[rownames(mat),1]


# Save table
write.table(mat, file=paste(workDir,"DE_genes.txt",sep=""), quote=FALSE, row.names=FALSE, sep="\t")



# Run Camera
genesnames=as.data.frame(mapping[rownames(rawcounts),1])
c2=readGMT("c2.cp.reactome.v6.0.symbols.gmt")
lengths <- function(x) vapply(x,length,1L)
index <- ids2indices(c2, genesnames) 
gs2 <- camera(y=v, index=index, design=design,inter.gene.cor=0.01)



reacto=gs2[gs2$FDR<0.05,]
reacto[,"names"]=gsub("REACTOME_(.*)","\\1",rownames(reacto))
up=reacto[reacto$Direction=="Up",]
down=reacto[reacto$Direction=="Down",]




######################
# Finish and quit

cat ("The pipeline is over!")
q(save='no')
