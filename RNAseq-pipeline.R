# This script calls a number of RNA-seq pipelines
# What will be done: quality control (FastQC), alignment (STAR), gene fusion detection (STAR-Fusion), read counting (featureCounts)
# run this script by: Rscript rnaSeq-pipeline.R --name YOUR_SAMPLE_NAME --fusion Y
# *** User needs to write and change to local tool/data directories in the script ***

#!/usr/bin/env Rscript



###########
#Libraries#
###########

suppressMessages(library(Rsubread))
suppressMessages(library(optparse))



###################################
# Check first for input from user #
###################################


option_list <- list(
    make_option(c("-n","--name"), type="character",
                dest="sampleName", help="The name of the sample"),
    make_option(c("-f", "--fusion"), type="character", default="Y",
                dest = "fusion_caller", help="Should call fusions? Y or N",metavar="fusion")
)

opt <- parse_args(OptionParser(option_list=option_list))


if (length(opt$sampleName)==0) {
  stop("The sample name is missing. Check -h for help \n\n")
} else {
  sampleName <- opt$sampleName
}

fusion_caller <- opt$fusion_caller



############################################
#Tool/data directories for USER to define  #
############################################


genomeInd = [add the path to the STAR genome index (should be something like "ref_genome.fa.star.idx")]
fusionInd = [add the path to the STAR-Fusion libraries]
STAR = [add path to STAR]
STARFusion = [add path to STAR-Fusion]
FastQC = [add path to FastQC]
genomefa = [add path to the genome fasta file]
genomeGTF = [add path to the genome gtf file]
workDir = [add your working directory]
fastqDir = paste(/path/to/your/fastq/dir/,sampleName,"/",sep="") # add the correct path for the fastq here



###########
#Functions#
###########


checkFile <- function (file) {
  #Check if file exists otherwise create one
  if (!file.exists(file)) {
    cat(paste(file," does not exist! I will have to create one\n"))
    dir.create(file)
  } else {
    cat(paste(file," exists! Excellent, lets move one\n"))
  }
}


checkFileEmpty <- function (file) {
  #Check if file is empty, write custom message
  if (length(list.files(file))!=0 ) {
    cat (paste(file,"not empty, lets move to the next step\n\n\n"),sep="")
  } else {
    cat (paste(file,"seems to be empty! Something went wrong during this step\n\n\n"))
  }
}



#########################################
# first we need to create the directories

cat ("Creating directories...\n\n")

sampledir=paste(workDir,sampleName,"/",sep="")
fastqc=paste(sampledir,"fastqc/",sep="")
mapping=paste(sampledir,"mapping/",sep="")
read_counts=paste(sampledir,"read_counts/",sep="")
gene_fusion=paste(sampledir,"gene_fusion/",sep="")
diff_expression=paste(sampledir,"diff_expression/",sep="")
variants=paste(sampledir,"variants/",sep="")
read_counts_all=paste(workDir,"read_counts_all/",sep="")
fusionAnalysis=paste(workDir,"fusionAnalysisResults/",sep="")
allDirs=c(sampledir,fastqc,mapping,read_counts,gene_fusion,diff_expression,variants,read_counts_all)


invisible(lapply(allDirs,checkFile))

cat("Excellent, the directories are created, we move to the next step: STAR indexes\n\n")



###################################################################################################################
# user defines inside script where the genome/STAR indexes are. If not, we will need to create them from scratch
if (!file.exists(genomeInd)) {
    stop("I cannot find the STAR libraries, please check if your path is correct!\n If you have not created the STAR libraries, please use the next command to create one")
} else {
    cat("Ok it seems that the STAR libraries exists, lets move to the next step\n")
}


#cmd=paste(STAR," --runMode genomeGenerate --genomeDir ",STARdir ,"--genomeFastaFiles ",genomefa," --runThreadN 10 --sjdbGTFfile ",genomeGTF," --sjdbOverhang 99")
#system(cmd) 


#########################################################################################
# Run FASTQC and store results (assumes paired-end and 1 file per mate, need to fix that)

cat ("FastQC running...\n\n")

fastqs=list.files(fastqDir)
fastq1=paste(fastqDir,fastqs[1],sep="")
fastq2=paste(fastqDir,fastqs[2],sep="")

cmd=paste(FastQC," ", fastq1, " --outdir ", fastqc,sep="")
system(cmd)
cmd=paste(FastQC," ", fastq2, " --outdir ", fastqc,sep="")
system(cmd)

checkFileEmpty(fastqc)



#################################
# Align reads with STAR 2-pass

cat ("STAR running...\n\n")

cmd=paste(STAR," --genomeDir ",genomeInd," --readFilesIn ",fastq1," ",fastq2," --outReadsUnmapped Fastx --runThreadN 4 --outSAMstrandField intronMotif --outSAMtype BAM SortedByCoordinate --outSAMattrRGline ID:GRPundef --genomeLoad NoSharedMemory --twopassMode Basic --readFilesCommand 'gunzip -c' --alignSJDBoverhangMin 10 --alignMatesGapMax 100000 --alignIntronMax 100000 --alignSJstitchMismatchNmax 5 -1 5 5 --chimSegmentMin 12 --chimJunctionOverhangMin 12 --chimSegmentReadGapMax 3 --outFileNamePrefix ",mapping,sep="")
system(cmd)

checkFileEmpty(mapping)



#################################
# Perform gene fusion analysis 

cat ("STAR-Fusion running...\n\n")

if (fusion_caller=="N") {
  cat ("Fusion calling is skipped by user!\n")
} else if (fusion_caller=="Y") {
  cat ("Fusion calling...\n")

cmd=paste(STARFusion," --genome_lib_dir ",fusionInd," -J ",mapping,"Chimeric.out.junction --output_dir ",gene_fusion,sep="")
system(cmd)
}

checkFileEmpty(gene_fusion)



#################################
# Count reads
# Things to be careful:
# 1) this is for non-stranded protocol AND 2) captures unique only alignments

cat ("Counting reads...\n\n")

bamFile=paste(mapping,"Aligned.sortedByCoord.out.bam",sep="")

fcHuman <- featureCounts(files=bamFile,annot.ext=genomeGTF,isGTFAnnotationFile=T,isPairedEnd=T,nthreads=4,countChimericFragments=F,countMultiMappingReads=F,primaryOnly=F)
write.table(fcHuman$counts,file=paste(read_counts,"readCounts",sep=""),quote=F,col.names=T)

geneLength=fcHuman$annotation[,c("GeneID","Length")]
write.table(geneLength,file=paste(read_counts,"geneLength",sep=""),quote=F,col.names=T)

checkFileEmpty(read_counts)



######################
# Finish and quit


cat ("The pipeline is over! What you should do now is: Differential expression analysis and analyze the fusion data! \n\n")
q(save='no')
