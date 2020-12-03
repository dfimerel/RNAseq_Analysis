#############################################################################################
# R code to analyze output of gene fusion results from deFuse, tophat-Fusion or STAR-Fusion
# Input needs to be the output of these three tools and user should specify 
# which of the tool he uses
# 
# run this script by: Rscript fusion_Analysis.R --tool defuse
#############################################################################################

##
## this specific script was run for Breast Cancer data and uses subtype classification
## for some of the analysis. If any other data, some of the plots will not make any sense
##


###########
#Libraries#
###########

library(optparse)
library(biomaRt)
library(gplots)
library(ggplot2)
library(RCircos)
library(reshape2)


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

workDir = [add path to working directory]
outputDir = paste(workDir,"fusionAnalysis/",sep="")




###################
# Plot functions #
###################

# fusions per sample
fusionSamplePlot <- function (inputTable) {
  fus=melt(table(inputTable[,"sample"]))
  fus[grep("HER2",fus[,1]),3]="HER2"
  fus[grep("TN",fus[,1]),3]="TN"
  fus[grep("LA",fus[,1]),3]="LA"
  fus[grep("LB",fus[,1]),3]="LB"
  colnames(fus)=c("sample","frequency","subtype")
  ggplot(fus,aes(x=sample,y=frequency))+geom_bar(stat="identity",aes(fill=subtype))+coord_flip()+labs(y="Number of fusion genes",x="Sample")+ggtitle("Number of fusion genes per sample")+theme_bw()
}


# fusions per subtype
fusionSubtypePlot <- function (inputTable) {
  fus=melt(table(inputTable[,"sample"]))
  fus[grep("HER2",fus[,1]),3]="HER2"
  fus[grep("TN",fus[,1]),3]="TN"
  fus[grep("LA",fus[,1]),3]="LA"
  fus[grep("LB",fus[,1]),3]="LB"
  colnames(fus)=c("sample","frequency","subtype")
  ggplot(fus,aes(x=subtype,y=frequency))+geom_boxplot(aes(fill=subtype))+labs(x="Subtype",y="Number of fusion genes")+theme_bw()+ggtitle("Number of fusion genes per subtype")
}



# deletions-inversions-eversion-interchrom (works only with defuse results)
fusionTypePlot <- function (inputTable) {
  fusType=data.frame(c(length(inputTable[inputTable$deletion=="Y",1]),length(inputTable[inputTable$inversion=="Y",1]),length(inputTable[inputTable$eversion=="Y",1]),length(inputTable[inputTable$interchromosomal=="Y",1])))
  fusType[,2]=c("deletion","inversion","eversion","interchromosomal")
  colnames(fusType)=c("frequency","type")
  ggplot(fusType,aes(x=type,y=frequency))+geom_bar(stat="identity",aes(fill=type))+labs(x="Type",y="Number of fusions")+ guides(fill=FALSE)+theme_bw()+ggtitle("Number of fusions per fusion type")
}


# distance of breakpoints in intra
fusionIntraDistPlot <- function (inputTable) {
  intraDist=inputTable[inputTable$gene_chromosome1==inputTable$gene_chromosome2,c("genomic_break_pos1","genomic_break_pos2")]
  intraDist[,"dist"]=abs(as.numeric(intraDist[,2])-as.numeric(intraDist[,1]))
  ggplot(intraDist,aes(x=dist))+geom_histogram(fill="darkblue",col="red")+labs(x="distance [Kb]",y="")+ggtitle("Breakpoint distance in intrachromosomal fusions")+theme_bw()
}



# fusion per chromosome normalized by number of samples/genes
fusionChromNormPlot <- function (inputTable) {
  ensembl.db <- useMart("ENSEMBL_MART_ENSEMBL",host="grch37.ensembl.org", dataset="hsapiens_gene_ensembl") #grch37.ensembl.org www.ensembl.org
  geneNumMatrix=c()
  for (i in c(1:22,"X","Y")) {
    gene_num=getBM(attributes=c("ensembl_gene_id","chromosome_name"),filters="chromosome_name",values=i,mart=ensembl.db) 
    geneNumMatrix=c(geneNumMatrix,length(gene_num[,1]))
  }
  chr=c(inputTable[,"gene_chromosome1"],inputTable[,"gene_chromosome2"])
  chr=factor(chr,levels=c(1:22,"X"),ordered=T)
  chrNorm=melt(table(chr)/geneNumMatrix[1:23])
  ggplot(chrNorm,aes(x=chr,y=value))+geom_bar(stat="identity",color="red",fill="darkblue")+theme_bw()+labs(x="",y="Number of fusion genes \n(normalized by # of genes per chromosome)")+ggtitle("Number of fusion genes per chromosome")
}



# circos plot with Rcircos 
fusionCircosPlot <- function (inputTable) {
  data(UCSC.HG19.Human.CytoBandIdeogram)
  chr.exclude <- NULL
  cyto.info <- UCSC.HG19.Human.CytoBandIdeogram
  tracks.inside <- 1
  tracks.outside <- 0
  RCircos.Set.Core.Components(cyto.info, chr.exclude,tracks.inside, tracks.outside)

  linksData=inputTable[,c("gene_chromosome1","genomic_break_pos1","genomic_break_pos1","gene_chromosome2","genomic_break_pos2","genomic_break_pos2")]
  linksData[linksData[,"gene_chromosome1"]==linksData[,"gene_chromosome2"],"PlotColor"]="red"
  linksData[linksData[,"gene_chromosome1"]!=linksData[,"gene_chromosome2"],"PlotColor"]="blue"
  RCircos.Set.Plot.Area()
  RCircos.Chromosome.Ideogram.Plot()
  RCircos.Link.Plot(linksData,1)
}



# recurrent genes and fusions
fusionRecurrentPlot <- function (inputTable) {
  fusSample=inputTable[,c("g1g2","sample")]
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
  ggplot(plotRec,aes(x=X1,y=X2))+geom_bar(stat="identity",color="red",fill="darkblue")+coord_flip()+theme_bw()+labs(x="",y="Number of samples")+ggtitle("Recurrent fusions across samples")
}




###############################
 #Load data and do the plots #
###############################


inputFusions = read.csv(paste(workDir,"tab2_final_table.tsv",sep=""),header=TRUE,sep="\t",stringsAsFactors=F) #mastertab.tsv  fusion_mastertable.tsv


if (fusion_caller=="defuse") { ## FOR THE DEFUSE RESULTS ##
  #make column with the two genes
  inputFusions[,"g1g2"]=paste(inputFusions[,"gene_name1"],inputFusions[,"gene_name2"],sep="->")
  
  # access biomart to get some additional info
  ensembl.db = useMart("ENSEMBL_MART_ENSEMBL",host="grch37.ensembl.org", dataset="hsapiens_gene_ensembl")
  gid = unique(c(inputFusions[,"gene1"]))
  gene.an = getBM(attributes=c("ensembl_gene_id", "entrezgene", "hgnc_symbol", "description", "strand", "band", "gene_biotype"), filters="ensembl_gene_id", values=gid, mart=ensembl.db)
  gene.an = gene.an[!duplicated(gene.an[,"ensembl_gene_id"]),]
  rownames(gene.an) = gene.an[,"ensembl_gene_id"]
  colnames(gene.an) = paste(colnames(gene.an), 1, sep="")
  inputFusions = cbind(inputFusions, gene.an[inputFusions[,"gene1"], c("entrezgene1", "description1", "strand1", "band1", "gene_biotype1")])

  gid = unique(c(inputFusions[,"gene2"]))
  gene.an = getBM(attributes=c("ensembl_gene_id", "entrezgene", "hgnc_symbol", "description", "strand", "band", "gene_biotype"), filters="ensembl_gene_id", values=gid, mart=ensembl.db)
  gene.an = gene.an[!duplicated(gene.an[,"ensembl_gene_id"]),]
  rownames(gene.an) = gene.an[,"ensembl_gene_id"]
  colnames(gene.an) = paste(colnames(gene.an), 2, sep="")
  inputFusions = cbind(inputFusions, gene.an[inputFusions[,"gene2"], c("entrezgene2", "description2", "strand2", "band2", "gene_biotype2")])

  
  #################################################
  # Perform some custom filtering for the fusions #
  ################################################

  ## filter out mitochondrial genes
  inputFusions = inputFusions[inputFusions[,"gene_chromosome1"]!="MT",]
  inputFusions = inputFusions[inputFusions[,"gene_chromosome2"]!="MT",]

  ##filter out readthrough fusions
  inputFusions = inputFusions[inputFusions[,"read_through"]=="N",]

  ##filter out genes involved in normal sample fusions
  normalFusions = inputFusions[grep("^N",inputFusions$sample),]
  inputFusions = inputFusions[!inputFusions$g1g2 %in% normalFusions$g1g2,]

  ## filter out gene chromosome Y
  inputFusions = inputFusions[inputFusions[,"gene_chromosome1"]!="Y",]
  inputFusions = inputFusions[inputFusions[,"gene_chromosome2"]!="Y",]

  ##filter out Y_RNA genes
  inputFusions = inputFusions[!grepl("^Y_RNA",inputFusions[,"gene_name1"]),]
  inputFusions = inputFusions[!grepl("^Y_RNA", inputFusions[,"gene_name2"]),]

  ##filter out rRNA genes
  inputFusions = inputFusions[!grepl("rRNA", inputFusions[,"gene_name1"]),]
  inputFusions = inputFusions[!grepl("rRNA", inputFusions[,"gene_name2"]),]


  
  #################
  # Create  plots #
  #################

  # fusions per sample
  pdf(file=paste(outputDir,"genesSample.pdf",sep=""))
  fusionSamplePlot(inputFusions)
  dev.off()


  # fusions per subtype
  pdf(file=paste(outputDir,"genesSubtype.pdf",sep=""))
  fusionSubtypePlot(inputFusions)
  dev.off()


  # deletions-inversions-eversion-interchrom
  pdf(file=paste(outputDir,"fusionType.pdf",sep=""))
  fusionTypePlot(inputFusions)
  dev.off()


  # distance of breakpoints in intra
  pdf(file=paste(outputDir,"BreakDist.pdf",sep=""))
  fusionIntraDistPlot(inputFusions)
  dev.off()


  # fusion per chromosome normalized by number of samples/genes
  pdf(file=paste(outputDir,"genesChrom.pdf",sep=""))
  fusionChromNormPlot(inputFusions)
  dev.off()


  # circos plot with Rcircos 
  pdf(file=paste(outputDir,"circos.pdf",sep=""))
  fusionCircosPlot(inputFusions)
  dev.off()


  # recurrent genes and fusions
  pdf(file=paste(outputDir,"recurrentFusions.pdf",sep=""))
  fusionRecurrentPlot(inputFusions)
  dev.off()


  
} else if (fusion_caller=="star") { ## FOR THE STAR-FUSION RESULTS ##

  inputFusions[,"gene_name1"]=gsub("(.*)--(.*)","\\1",inputFusions[,1])
  inputFusions[,"gene_name2"]=gsub("(.*)--(.*)","\\2",inputFusions[,1])
  inputFusions[,"g1g2"]=paste(inputFusions$gene_name1,inputFusions$gene_name2,sep="->")

  # fix the annotation so that it is easily readable
  inputFusions[,"gene_chromosome1"]=gsub("chr(.*):.*:.*","\\1",inputFusions[,"LeftBreakpoint"])
  inputFusions[,"gene_chromosome2"]=gsub("chr(.*):.*:.*","\\1",inputFusions[,"RightBreakpoint"])
  
  inputFusions[,"genomic_break_pos1"]=gsub("chr.*:(.*):.*","\\1",inputFusions[,"LeftBreakpoint"])
  inputFusions[,"genomic_break_pos2"]=gsub("chr.*:(.*):.*","\\1",inputFusions[,"RightBreakpoint"])


  #################################################
  # Perform some custom filtering for the fusions #
  ################################################

  ## filter out mitochondrial genes
  inputFusions = inputFusions[inputFusions[,"gene_chromosome1"]!="MT",]
  inputFusions = inputFusions[inputFusions[,"gene_chromosome2"]!="MT",]

  ##filter out genes involved in normal sample fusions
  normalFusions = inputFusions[grep("^N",inputFusions$sample),]
  inputFusions = defuseInput[!inputFusions$g1g2 %in% normalFusions$g1g2,]

  ## filter out gene chromosome Y
  inputFusions = inputFusions[inputFusions[,"gene_chromosome1"]!="Y",]
  inputFusions = inputFusions[inputFusions[,"gene_chromosome2"]!="Y",]

  ##filter out Y_RNA genes
  inputFusions = inputFusions[!grepl("^Y_RNA",inputFusions[,"gene_name1"]),]
  inputFusions = inputFusions[!grepl("^Y_RNA", inputFusions[,"gene_name2"]),]

  ##filter out rRNA genes
  inputFusions = inputFusions[!grepl("rRNA", inputFusions[,"gene_name1"]),]
  inputFusions = inputFusions[!grepl("rRNA", inputFusions[,"gene_name2"]),]


  #################
  # Create  plots #
  #################

  # fusions per sample
  pdf(file=paste(outputDir,"genesSample.pdf",sep=""))
  fusionSamplePlot(inputFusions)
  dev.off()


  # fusions per subtype
  pdf(file=paste(outputDir,"genesSubtype.pdf",sep=""))
  fusionSubtypePlot(inputFusions)
  dev.off()


  # distance of breakpoints in intra
  pdf(file=paste(outputDir,"BreakDist.pdf",sep=""))
  fusionIntraDistPlot(inputFusions)
  dev.off()


  # fusion per chromosome normalized by number of samples/genes
  pdf(file=paste(outputDir,"genesChrom.pdf",sep=""))
  fusionChromNormPlot(inputFusions)
  dev.off()


  # circos plot with Rcircos 
  pdf(file=paste(outputDir,"circos.pdf",sep=""))
  fusionCircosPlot(inputFusions)
  dev.off()


  # recurrent genes and fusions
  pdf(file=paste(outputDir,"recurrentFusions.pdf",sep=""))
  fusionRecurrentPlot(inputFusions)
  dev.off()

  
}




######################
# Finish and quit

cat ("The pipeline is over!")
q(save='no')
