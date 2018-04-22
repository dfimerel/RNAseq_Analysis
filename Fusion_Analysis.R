#############################################################################################
# R code to analyze output of gene fusion results from deFuse, tophat-Fusion or STAR-Fusion
# Input needs to be the output of these three tools and user should specify 
# which of the tool he uses
# 
# run this script by: Rscript Fusion_Analysis.R --tool defuse
# TO DO: add possibility for analyzing results from star-fusion and tophat
#############################################################################################

##
## this specific script was run for Breast Cancer data and uses subtype classification
## for some of the analysis. It assumes also that sample name and subtype is added
##If any other data, some of the plots will not make any sense
##


###########
#Libraries#
###########

library(optparse)
library(biomaRt)
library(gplots)
library(ggplot2)
library(RCircos)


###################################
# Check first for input from user #
###################################


option_list <- list(
    make_option(c("-t", "--tool"), type="character", default="defuse",
                dest = "fusion_caller", help="Choose from: defuse,star,tophat",metavar="fusion")
)

opt <- parse_args(OptionParser(option_list=option_list))
fusion_caller <- opt$fusion_caller



############################################
#Tool/data directories for USER to define  #
############################################

workDir = [add your working directory]
outputDir = paste(workDir,"fusionAnalysis/",sep="")


#############
 #Load data #
#############


defuseInput = read.csv(paste(workDir,"fusion_results_all.tsv",sep=""),header=TRUE,sep="\t",stringsAsFactors=F)

#make column with the two genes
defuseInput[,"g1g2"]=paste(defuseInput[,"gene_name1"],defuseInput[,"gene_name2"],sep="->")


# access biomart to get some additional info

ensembl.db <- useMart("ENSEMBL_MART_ENSEMBL",host="www.ensembl.org", dataset="hsapiens_gene_ensembl")
gid <- unique(c(defuseInput[,"gene1"]))
gene.an <- getBM(attributes=c("ensembl_gene_id", "entrezgene", "hgnc_symbol", "description", "strand", "band", "gene_biotype"), filters="ensembl_gene_id", values=gid, mart=ensembl.db)
gene.an <- gene.an[!duplicated(gene.an[,"ensembl_gene_id"]),]
rownames(gene.an) <- gene.an[,"ensembl_gene_id"]
colnames(gene.an) <- paste(colnames(gene.an), 1, sep="")
defuseInput <- cbind(defuseInput, gene.an[defuseInput[,"gene1"], c("entrezgene1", "description1", "strand1", "band1", "gene_biotype1")])

gid <- unique(c(defuseInput[,"gene2"]))
gene.an <- getBM(attributes=c("ensembl_gene_id", "entrezgene", "hgnc_symbol", "description", "strand", "band", "gene_biotype"), filters="ensembl_gene_id", values=gid, mart=ensembl.db)
gene.an <- gene.an[!duplicated(gene.an[,"ensembl_gene_id"]),]
rownames(gene.an) <- gene.an[,"ensembl_gene_id"]
colnames(gene.an) <- paste(colnames(gene.an), 2, sep="")
defuseInput <- cbind(defuseInput, gene.an[defuseInput[,"gene2"], c("entrezgene2", "description2", "strand2", "band2", "gene_biotype2")])



#################################################
# Perform some custom filtering for the fusions #
################################################


## filter out mitochondrial genes
defuseInput = defuseInput[defuseInput[,"gene_chromosome1"]!="MT",]
defuseInput = defuseInput[defuseInput[,"gene_chromosome2"]!="MT",]

##filter out readthrough fusions
defuseInput = defuseInput[defuseInput[,"read_through"]=="N",]

##filter out genes involved in normal sample fusions
normalFusions = defuseInput[grep("^N",defuseInput$sample),]
defuseInput = defuseInput[!defuseInput$g1g2 %in% normalFusions$g1g2,]

## filter out gene chromosome Y
defuseInput = defuseInput[defuseInput[,"gene_chromosome1"]!="Y",]
defuseInput = defuseInput[defuseInput[,"gene_chromosome2"]!="Y",]

##filter out Y_RNA genes
defuseInput = defuseInput[!grepl("^Y_RNA",defuseInput[,"gene_name1"]),]
defuseInput = defuseInput[!grepl("^Y_RNA", defuseInput[,"gene_name2"]),]

##filter out rRNA genes
defuseInput = defuseInput[!grepl("rRNA", defuseInput[,"gene_name1"]),]
defuseInput = defuseInput[!grepl("rRNA", defuseInput[,"gene_name2"]),]





#################
# Create  plots #
#################

# fusions per sample
fus=melt(table(defuseInput[,"sample"]))
fus[grep("HER2",fus[,1]),3]="HER2"
fus[grep("TN",fus[,1]),3]="TN"
fus[grep("LA",fus[,1]),3]="LA"
fus[grep("LB",fus[,1]),3]="LB"
colnames(fus)=c("sample","frequency","subtype")
pdf(file=paste(outputDir,"genesSample.pdf",sep=""))
ggplot(fus,aes(x=sample,y=frequency))+geom_bar(stat="identity",aes(fill=subtype))+coord_flip()+labs(y="Number of fusion genes",x="Sample")+ggtitle("Number of fusion genes per sample")+theme_bw()
dev.off()


# fusions per subtype
pdf(file=paste(outputDir,"genesSubtype.pdf",sep=""))
ggplot(fus,aes(x=subtype,y=frequency))+geom_boxplot(aes(fill=subtype))+labs(x="Subtype",y="Number of fusion genes")+theme_bw()+ggtitle("Number of fusion genes per subtype")
dev.off()


# deletions-inversions-eversion-interchrom
fusType=data.frame(c(length(defuseInput[defuseInput$deletion=="Y",1]),length(defuseInput[defuseInput$inversion=="Y",1]),length(defuseInput[defuseInput$eversion=="Y",1]),length(defuseInput[defuseInput$interchromosomal=="Y",1])))
fusType[,2]=c("deletion","inversion","eversion","interchromosomal")
colnames(fusType)=c("frequency","type")
pdf(file=paste(outputDir,"fusionType.pdf",sep=""))
ggplot(fusType,aes(x=type,y=frequency))+geom_bar(stat="identity",aes(fill=type))+labs(x="Type",y="Number of fusions")+ guides(fill=FALSE)+theme_bw()+ggtitle("Number of fusions per fusion type")
dev.off()


# distance of breakpoints in intra
intraDist=defuseInput[defuseInput$interchromosomal=="N",c("genomic_break_pos1","genomic_break_pos2")]
intraDist[,"dist"]=abs(intraDist[,2]-intraDist[,1])
pdf(file=paste(outputDir,"BreakDist.pdf",sep=""))
ggplot(intraDist,aes(x=dist))+geom_histogram(fill="darkblue",col="red")+labs(x="distance [Kb]",y="")+ggtitle("Breakpoint distance in intrachromosomal fusions")+theme_bw()
dev.off()


# fusion per chromosome normalized by number of samples/genes
ensembl.db <- useMart("ENSEMBL_MART_ENSEMBL",host="grch37.ensembl.org", dataset="hsapiens_gene_ensembl")
geneNumMatrix=c()
for (i in c(1:22,"X","Y")) {
  gene_num=getBM(attributes=c("ensembl_gene_id","chromosome_name"),filters="chromosome_name",values=i,mart=ensembl.db) 
  geneNumMatrix=c(geneNumMatrix,length(gene_num[,1]))
}
chr=c(defuseInput[,"gene_chromosome1"],defuseInput[,"gene_chromosome2"])
chr=factor(chr,levels=c(1:22,"X"),ordered=T)
chrNorm=melt(table(chr)/geneNumMatrix[1:23])
pdf(file=paste(outputDir,"genesChrom.pdf",sep=""))
ggplot(chrNorm,aes(x=chr,y=value))+geom_bar(stat="identity",color="red",fill="darkblue")+theme_bw()+labs(x="",y="Number of fusion genes \n(normalized by # of genes per chromosome)")+ggtitle("Number of fusion genes per chromosome")
dev.off()



# circos plot with Rcircos
data(UCSC.HG19.Human.CytoBandIdeogram)
chr.exclude <- NULL
cyto.info <- UCSC.HG19.Human.CytoBandIdeogram
tracks.inside <- 1
tracks.outside <- 0
RCircos.Set.Core.Components(cyto.info, chr.exclude,tracks.inside, tracks.outside)

linksData=defuseInput[,c("gene_chromosome1","genomic_break_pos1","genomic_break_pos1","gene_chromosome2","genomic_break_pos2","genomic_break_pos2")]
linksData[linksData[,"gene_chromosome1"]==linksData[,"gene_chromosome2"],"PlotColor"]="red"
linksData[linksData[,"gene_chromosome1"]!=linksData[,"gene_chromosome2"],"PlotColor"]="blue"

pdf(file=paste(outputDir,"circos.pdf",sep=""))
RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()
RCircos.Link.Plot(linksData,1)
dev.off()


# recurrent genes and fusions
fusSample=defuseInput[,c("g1g2","sample")]
fusSample1=fusSample[duplicated(fusSample[,1],fromLast=T) |duplicated(fusSample[,1]),]
nn=unique(fusSample1[duplicated(fusSample1[,1]),1])

recurrentFus=c()
for (i in 1:length(nn)) {
  if (length(unique(fusSample1[fusSample1[,1]==nn[i],])[,1])>1) {
    recurrentFus=rbind(recurrentFus,fusSample1[fusSample1[,1]==nn[i],])
  }
}


nn=unique(recurrentFus[,1])
plotRec=data.frame(matrix(ncol=2,nrow=length(nn)))
for (i in 1:length(nn)) {
  plotRec[i,2]=length(recurrentFus[recurrentFus[,1]==nn[i],1])
  plotRec[i,1]=nn[i]
}
plotRec[,1]=factor(plotRec[,1],levels=plotRec[order(plotRec[,2]),1])
pdf(file=paste(outputDir,"recurrentFusions.pdf",sep=""))
ggplot(plotRec,aes(x=X1,y=X2))+geom_bar(stat="identity",color="red",fill="darkblue")+coord_flip()+theme_bw()+labs(x="",y="Number of samples")+ggtitle("Recurrent fusions across samples")
dev.off()


######################
# Finish and quit

q(save='no')
