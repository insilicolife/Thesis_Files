

########################################################################
#		Prepare the gbm gene expression samples Matrix	       #
########################################################################
gbmGeneExpression=read.table("gbmGeneExpression/gene_expr.tsv", header=TRUE, sep="\t")
gbmGeneInfo=gbmGeneExpression[,1:5]
gbmGeneExpressionSampleInfo=colnames(gbmGeneExpression)[6:ncol(gbmGeneExpression)]
gbmGeneExpressionMatrix=gbmGeneExpression[,6:ncol(gbmGeneExpression)]
dim(gbmGeneExpressionMatrix)
rownames(gbmGeneExpressionMatrix)=gbmGeneInfo[,4]

sampleSummeryInformation=read.table("summary.tsv", sep="\t", header=TRUE)
clinicalPatientInfo=read.table("Clinical/Biotab/nationwidechildrens.org_clinical_patient_gbm.txt", sep="\t", header=TRUE)
controlSamples=read.table("Clinical/Biotab/nationwidechildrens.org_control_cntl.txt", sep="\t", header=TRUE)
clinicalDrugInfo=read.table("Clinical/Biotab/nationwidechildrens.org_clinical_drug_gbm.txt", sep="\t", header=TRUE)
clinicalOmfInfo=read.table("Clinical/Biotab/nationwidechildrens.org_clinical_omf_v4.0_gbm.txt", sep="\t", header=TRUE)

sampleSummeryInformation$filename=gsub("-",".",sampleSummeryInformation$filename)
summaryInformationIndexForGBMSamples=na.omit(match(gbmGeneExpressionSampleInfo,sampleSummeryInformation$filename))

gbmSampleSummaryMatrix=sampleSummeryInformation[summaryInformationIndexForGBMSamples,]
#factor the samples between 01 for primary solid Tumor and 02 for Recurent solid tumor
#gbmSampleSummaryMatrix$sample_type_name=factor(gbmSampleSummaryMatrix$sample_type_name, labels=c("1TP","2TR"))

#primarySolidTumorIndex=which(gbmSampleSummaryMatrix$sample_type_name=="1TP")
#recurentSolidTumorIndex=which(gbmSampleSummaryMatrix$sample_type_name=="2TR")
source("TCGAsampleIDParser.r")
gbmGeneExpressionSampleID=TCGAsampleIdpharser(gbmSampleSummaryMatrix$barcode)
colnames(gbmGeneExpressionMatrix)=names(gbmGeneExpressionSampleID)
source("DECalculaterDESeq2.r")
#lets prepare sample lables and the data matrix to be so that we can use it in DESeq2
sampleTypes=gbmSampleSummaryMatrix$sample_type
gbmGeneExpressionSampleClassLables=data.frame(condition=sampleTypes)
rownames(gbmGeneExpressionSampleClassLables)=colnames(gbmGeneExpressionMatrix)

gbmGeneExpressionMatrixDESeq2Normalized=DESeq2Normalization(gbmGeneExpressionMatrix,gbmGeneExpressionSampleClassLables, "condition")



######################################################################
#		Prepare the novel transcript expression matrix	     #
######################################################################
gbmNovelTranscriptExpression=read.table("gbmLncRNAExpression/novel_transcript_expressions_normalized.txt", header=TRUE, sep="\t")

gbmNobelTranscriptExpressionMatrix=gbmNovelTranscriptExpression[,-1]
summaryInformationIndexForGBMSamples=na.omit(match(colnames(gbmNobelTranscriptExpressionMatrix),sampleSummeryInformation$filename))
gbmNobelTranscriptSampleSummaryMatrix=sampleSummeryInformation[summaryInformationIndexForGBMSamples,]
gbmNovelTranscriptSampleID=TCGAsampleIdpharser(gbmNobelTranscriptSampleSummaryMatrix$barcode)
colnames(gbmNobelTranscriptExpressionMatrix)=names(gbmNovelTranscriptSampleID)
rownames(gbmNobelTranscriptExpressionMatrix)=gbmNovelTranscriptExpression$Transcript

transcriptSampleType=gbmNobelTranscriptSampleSummaryMatrix$sample_type
transcriptExpressionClassLabeles=data.frame(condition=transcriptSampleType)
rownames(transcriptExpressionClassLabeles)=colnames(gbmNobelTranscriptExpressionMatrix)


gbmNobelTranscriptExpressionMatrixDESeq2Normalized=DESeq2Normalization(gbmNobelTranscriptExpressionMatrix,transcriptExpressionClassLabeles, "condition")

#######################################################################
#	Prepare the miRNA expression data from MAtti	 	      #
###############################################################3#######

#----------install R.matlab Library to access the .mat file-----------
if("R.matlab" %in% rownames(installed.packages()) == FALSE){
install.packages("R.matlab")
}
library("R.matlab")
#read teh .mat file
gbmmiRNAExpression=readMat("gbmmiRNAExprssion/mirna_expression_agilent_level3.mat")

#teh file contains the following matlab workspace variable
#mirna_expr = 
#  meta: [1x1 struct]
#  mean: [534x506 double]
#  rows: [1x1 struct]
#
#gbmmiRNAExpression
#$mirna.expr
#, , 1
#     [,1]
#meta List,42
#mean Numeric,270204
#rows List,1
#attr(,"header")
#attr(,"header")$description
#[1] "MATLAB 5.0 MAT-file, Platform: GLNXA64, Created on: Thu Apr  5 10:03:40 2012      #
#attr(,"header")$version
#[1] "5"dim
#attr(,"header")$endian
#[1] "little"
#metaData
metadata=gbmmiRNAExpression$mirna.expr[1,,]
#41 type of sample information for each samples
gbmmiRNASampleInfo=names(metadata$meta[,,1][2:42])
gbmmiRNASampleInfoValues=metadata$meta[,,1][gbmmiRNASampleInfo]
names(gbmmiRNASampleInfoValues)=gbmmiRNASampleInfo
gbmmiRNASampleInfoMatrix=matrix(nrow=506, ncol=41)
gbmmiRNASampleInfoMatrix[,1:41]=unlist(gbmmiRNASampleInfoValues[1:41])
colnames(gbmmiRNASampleInfoMatrix)=gbmmiRNASampleInfo

treatedGbmmiRNASamplesIndex=which(gbmmiRNASampleInfoMatrix[,"sample.type"]=="Treated primary GBM")
untreatedGbmmiRNASamplesIndex=which(gbmmiRNASampleInfoMatrix[,"sample.type"]=="Untreated primary (de novo) GBM")
unknownGbmmiRNASamplesIndex=which(gbmmiRNASampleInfoMatrix[,"sample.type"]=="-")

treatedGbmmiRNASamples=gbmmiRNASampleInfoMatrix[treatedGbmmiRNASamplesIndex,]
untreatedGbmmiRNASamples=gbmmiRNASampleInfoMatrix[untreatedGbmmiRNASamplesIndex,]
unknownGbmmiRNASamples=gbmmiRNASampleInfoMatrix[unknownGbmmiRNASamplesIndex,]

gbmmiRNANames=unlist(gbmmiRNAExpression$mirna.expr[[3]])

gbmmiRNAExpressionData=matrix(unlist(gbmmiRNAExpression$mirna.expr[[2]]),nrow=534,ncol=506)#534by 506
rownames(gbmmiRNAExpressionData)=gbmmiRNANames
colnames(gbmmiRNAExpressionData)=gbmmiRNASampleInfoMatrix[,"tcga.barcode"]
TreatedVsUntreatedSamplesIndex=c(treatedGbmmiRNASamplesIndex,untreatedGbmmiRNASamplesIndex)
gbmmiRNASampleInfoMatrixTreatedVsUntreated=gbmmiRNASampleInfoMatrix[TreatedVsUntreatedSamplesIndex,]

gbmmiRNAExpressionDataTreatedVsUntreatedSamples=gbmmiRNAExpressionData[,TreatedVsUntreatedSamplesIndex]
colnames(gbmmiRNAExpressionDataTreatedVsUntreatedSamples)=gbmmiRNASampleInfoMatrixTreatedVsUntreated[,"sample.id"]
dim(gbmmiRNAExpressionDataTreatedVsUntreatedSamples)
#[1] 534 487

#selecting the 169 samples that we used in the transcript and gene expression data
selected169SamplesFrommiRNADataset=na.omit(match(gbmGeneExpressionSampleID,gbmmiRNASampleInfoMatrix[,"sample.id"]))
gbmmiRNAExpressionDataFor169SampleOfInterst=gbmmiRNAExpressionData[,selected169SamplesFrommiRNADataset]


png("visualize.png")
par(mfrow=c(3,2))
plot(log(gbmGeneExpressionMatrixDESeq2Normalized), main="Scatter plot of gene expression")
hist(log(gbmGeneExpressionMatrixDESeq2Normalized), main="Histogram of gene expression")
plot(log(gbmNobelTranscriptExpressionMatrixDESeq2Normalized),main="Scatter plot of transcription expression" )
hist(log(gbmNobelTranscriptExpressionMatrixDESeq2Normalized), main="Histogram of transcription expression")
plot(log(gbmmiRNAExpressionDataFor169SampleOfInterst), main="Scatter plot of miRNA expression")
hist(log(gbmmiRNAExpressionDataFor169SampleOfInterst), main="Histogram of miRNA expression")
dev.off()

png("boxplot.png")
par(mfrow=c(1,3))
boxplot(log(gbmGeneExpressionMatrixDESeq2Normalized), main="Boxplot of gene expression")
boxplot(log(gbmNobelTranscriptExpressionMatrixDESeq2Normalized), main="Boxplot of transcription expression")
boxplot(log(gbmmiRNAExpressionDataFor169SampleOfInterst), main="Boxplot of miRNA expression")
dev.off()





#lets select only the samples that we study for the gene expression and transcript expression

sampleOfInterstSampleID=colnames(gbmGeneExpressionMatrix)
gbmmiRNASampleInfoIndexTreatedVsUntreatedForSampleOfInterst=na.omit(match(sampleOfInterstSampleID,colnames(gbmmiRNAExpressionDataTreatedVsUntreatedSamples)))
gbmmiRNASampleInfoMatrixTreatedVsUntreatedForSampleOfInterst=gbmmiRNASampleInfoMatrixTreatedVsUntreated[gbmmiRNASampleInfoIndexTreatedVsUntreatedForSampleOfInterst,]
gbmmiRNAExpressionDataTreatedVsUntreatedSamplesOfInterst=gbmmiRNAExpressionDataTreatedVsUntreatedSamples[,gbmmiRNASampleInfoIndexTreatedVsUntreatedForSampleOfInterst]
##########################################################################
#		level3 data downloaded from TCGA by myself
# lets export the TCGA barcodes of our 169 samples to a file and search the miRNA expression data for this spacific samples
##########################################################################
write(as.character(sampleOfInterstSampleID), file="TCGAsampleOfInterstBarcode.txt", sep="\n")
#after downloading the miRNA expression data lets make a data matrix for all samples in a single data matrix 
filepath="gbmmiRNAExprssion/miRNAexpressionForSampleOfInterest/Expression-miRNA/UNC__H-miRNA_8x15K/Level_3/"
file_list=list.files(filepath)
 
for (file in file_list){
       
  # if the merged dataset doesn't exist, create it
  if (!exists("miRNAExpressionMydownload")){
    	miRNAExpressionMydownload= read.table(paste(filepath,file, sep=""), header=TRUE, sep="\t")
  
	columnnames=gsub(".","-",colnames(miRNAExpressionMydownload)[2],fixed = TRUE)
}
   
  # if the merged dataset does exist, append to it
  if (exists("miRNAExpressionMydownload")){
  	temp_dataset =read.table(paste(filepath,file, sep=""), header=TRUE, sep="\t")
	newColumnnames=gsub(".","-",colnames(temp_dataset)[2], fixed = TRUE)
	columnnames=c(columnnames,newColumnnames)
    	miRNAExpressionMydownload=cbind(miRNAExpressionMydownload, temp_dataset[,-1])
    	rm(temp_dataset)
	rm(newColumnnames)
  }
 
}
rownames(miRNAExpressionMydownload)=miRNAExpressionMydownload[,1]
miRNAExpressionMydownload=miRNAExpressionMydownload[-1,-1]
myDownloadsSampleID=TCGAsampleIdpharser(columnnames)
colnames(miRNAExpressionMydownload)=myDownloadsSampleID
#getting the Meatadata infromation

myDownloadedSampleInfo=read.table("gbmmiRNAExprssion/miRNAexpressionForSampleOfInterest/METADATA/UNC__H-miRNA_8x15K/unc.edu_GBM.H-miRNA_8x15K.sdrf.txt", header=TRUE, sep="\t")

#######################################################################
#			procceed with DE Analysis		      #
#######################################################################

#From the sample infromation matrixes for gene expressiona and transcript expression we see that the samples are of two categories: Primary solid tumor and Recursive solid tumor. Therefore the DE analysis is to get Genes, Transcripts and miRNAS that are expressed differently among this two sample groups	  
#-----------------------------------------------------------------------

#-------------------DE analysis for Gene Expression Data-----------------
#source("DECalculaterDESeq2.r")

#Apply DESeq2 for DE analysis of row count data

#DEResultForGenesExpression=DECalculaterDESeq2(gbmGeneExpressionMatrix,gbmGeneExpressionSampleClassLables,"condition")
#DESeq2NormalizedAndRegularizedTransformedGeneExpressionData=regularizedlogTransform(unlist(DEResultForGenesExpression)[[1]])

#SignificantlyDifferentlyExpressedGeneIndex=which(DEResultForGenesExpression[2][[1]]$padj<.05)
#SignificantlyDifferentlyExpressedGeneMatrix=gbmGeneExpressionMatrix[SignificantlyDifferentlyExpressedGeneIndex,]
#SignificatlGenes=rownames(SignificantlyDifferentlyExpressedGeneMatrix)
#length(SignificatlGenes)
#[1] 534

#gbmGeneExpressionMatrixNormalized=regularizedlogTransform(unlist(SignificantlyDifferentlyExpressedGeneMatrix)[1])

#---------------DE analysis  based on the quartile methods-------------
source("DEquartileAnalyzer.r")
logRatioOfHighOverLowExpressedGeneIn5pcOfQuartiles=DEAnalyzerBasedOnTheQuartiles(gbmGeneExpressionMatrix,FALSE)
#take top 100 differentially expressed genes among the quartils
orderedIndex=order(logRatioOfHighOverLowExpressedGeneIn5pcOfQuartiles[[1]], decreasing = TRUE)
orderedLogRatioOfHighOverLowExpressedGeneIn5pcOfQuartiles=logRatioOfHighOverLowExpressedGeneIn5pcOfQuartiles[[1]][orderedIndex]

top61SignificantGenes=head(orderedLogRatioOfHighOverLowExpressedGeneIn5pcOfQuartiles,61)
quartileSampleGroupsLowExpression=logRatioOfHighOverLowExpressedGeneIn5pcOfQuartiles[[2]][[1]][orderedIndex]

quartileSampleGroupsHighExpression=logRatioOfHighOverLowExpressedGeneIn5pcOfQuartiles[[2]][[2]][orderedIndex]
names(quartileSampleGroupsHighExpression)=names(orderedLogRatioOfHighOverLowExpressedGeneIn5pcOfQuartiles)
names(quartileSampleGroupsLowExpression)=names(orderedLogRatioOfHighOverLowExpressedGeneIn5pcOfQuartiles)

quartileSampleGroupsLowExpressionForTop61Genes=head(quartileSampleGroupsLowExpression,61)

quartileSampleGroupsHighExpressionForTop61Genes=head(quartileSampleGroupsHighExpression,61)

names(quartileSampleGroupsHighExpressionForTop61Genes)=names(top61SignificantGenes)
names(quartileSampleGroupsLowExpressionForTop61Genes)=names(top61SignificantGenes)


top61Genes=names(top61SignificantGenes)

#---------------DE analysis for the Transcrpt expression data-----------

#transcriptSampleType=gbmNobelTranscriptSampleSummaryMatrix$sample_type
#transcriptExpressionClassLabeles=data.frame(condition=transcriptSampleType)
#rownames(transcriptExpressionClassLabeles)=colnames(gbmNobelTranscriptExpressionMatrix)
#DEResultForTranscriptExpression=DECalculaterDESeq2(gbmNobelTranscriptExpressionMatrix,transcriptExpressionClassLabeles,"condition")
#DEseqNormalixzedAndREgularizedTransformedTranscriptExpressionData=regularizedlogTransform(unlist(DEResultForTranscriptExpression)[[1]])

#SignificantlyDifferentlyExpressedTranscriptIndex=which(DEResultForTranscriptExpression[2][[1]]$padj<.05)
#SignificantlyDifferentlyExpressedTranscriptMatrix=gbmNobelTranscriptExpressionMatrix[SignificantlyDifferentlyExpressedTranscriptIndex,]
#SignificatlTranscripts=rownames(SignificantlyDifferentlyExpressedTranscriptMatrix)
#length(SignificatlTranscripts)
#[1] 8

logRatioOfTranscriptHighOverLowExpressedTranscriptIn5pcOfQuartiles=DEAnalyzerBasedOnTheQuartiles(gbmNobelTranscriptExpressionMatrix,FALSE, FALSE)
orderedIndexTranscript=order(logRatioOfTranscriptHighOverLowExpressedTranscriptIn5pcOfQuartiles[[1]],decreasing = TRUE)
orderedlogRatioOfTranscriptHighOverLowExpressedTranscriptIn5pcOfQuartiles=logRatioOfTranscriptHighOverLowExpressedTranscriptIn5pcOfQuartiles[[1]][orderedIndexTranscript]
#take top 100 differentially expressed genes among the quartils
top61SignificantTranscript=head(orderedlogRatioOfTranscriptHighOverLowExpressedTranscriptIn5pcOfQuartiles,61)

quartileTranscriptSampleGroupsLowExpression=logRatioOfTranscriptHighOverLowExpressedTranscriptIn5pcOfQuartiles[[2]][[1]][orderedIndexTranscript]

quartileTranscriptSampleGroupsHighExpression=logRatioOfTranscriptHighOverLowExpressedTranscriptIn5pcOfQuartiles[[2]][[2]][orderedIndexTranscript]
quartileSampleGroupsLowExpressionForTop61Transcript=head(quartileTranscriptSampleGroupsLowExpression,61)

quartileSampleGroupsHighExpressionForTop61Transcript=head(quartileTranscriptSampleGroupsHighExpression,61)

names(quartileSampleGroupsLowExpressionForTop61Transcript)=names(top61SignificantTranscript)
names(quartileSampleGroupsHighExpressionForTop61Transcript)=names(top61SignificantTranscript)


#XLConnect is the
if("org.Hs.eg.db" %in% rownames(installed.packages()) == FALSE){
install.packages("XLConnect")
}
library("XLConnect")
gbmNobelLincRNA = loadWorkbook("gbmLncRNAExpression/Novel Transcripts in GBM characteristics.xlsx")
setMissingValue(gbmNobelLincRNA, value = "NA")
gbm53lincRNAInfo=readWorksheet(gbmNobelLincRNA, sheet = "All 53 NTs")
gbm23nobelLincRNA=readWorksheet(gbmNobelLincRNA, sheet = "Chosen ones")

#------------DE analysis for miRNA expression data from Matti------------

logRatioOfHighOverLowExpressedmiRNAIn5pcOfQuartiles=DEAnalyzerBasedOnTheQuartiles(gbmmiRNAExpressionDataFor169SampleOfInterst,TRUE)

orderedIndexmiRNA=order(logRatioOfHighOverLowExpressedmiRNAIn5pcOfQuartiles[[1]],decreasing = TRUE)
orderedLogRatioOfHighOverLowExpressedmiRNAIn5pcOfQuartiles=logRatioOfHighOverLowExpressedmiRNAIn5pcOfQuartiles[[1]][orderedIndexmiRNA]
top61SignificantmiRNA=head(orderedLogRatioOfHighOverLowExpressedmiRNAIn5pcOfQuartiles,61)
quartilemiRNASampleGroupsLowExpression=logRatioOfHighOverLowExpressedmiRNAIn5pcOfQuartiles[[2]][[1]][orderedIndexmiRNA]

qartilemiRNASampleGroupsHighExpression=logRatioOfHighOverLowExpressedmiRNAIn5pcOfQuartiles[[2]][[2]][orderedIndexmiRNA]
quartileSampleGroupsLowExpressionForTop61miRNA=head(quartilemiRNASampleGroupsLowExpression,61)
quartileSampleGroupsHighExpressionForTop61miRNA=head(qartilemiRNASampleGroupsHighExpression,61)

names(quartileSampleGroupsLowExpressionForTop61miRNA)=names(top61SignificantmiRNA)
names(quartileSampleGroupsHighExpressionForTop61miRNA)=names(top61SignificantmiRNA)

#-------------Finding the association between the significantly differently expressed genes, transcripts and miRNAs in some of the samples------------------

#-------------Calculate the corrolation between the gene vs Transcript, gene vs miRNA and transcript vs miRNA

pdf("logratioVisualization.pdf")

#par(mfrow=c(1,3))
plot(top61SignificantGenes,top61SignificantTranscript, main="corrolation betweeen log ratio of genes and transcript", xlab="logratio of top61 genes", ylab="logratio of top61 Transcripts")
lincRNAIdentifiedIndex=na.omit(match(names(top61SignificantTranscript),gbm53lincRNAInfo$Gene))
lincRNAOfInterstIndex=na.omit(match(names(top61SignificantTranscript),gbm23nobelLincRNA$Gene))
points(top61SignificantGenes[lincRNAIdentifiedIndex],top61SignificantTranscript[lincRNAIdentifiedIndex],col="green")
points(top61SignificantGenes[match(names(LincRNASeqHg19),names(top61SignificantTranscript))],top61SignificantTranscript[match(names(LincRNASeqHg19),names(top61SignificantTranscript))],col="red")
legend('bottomright', c("53 Nobel lincRNA Transcripts","17 lincRNA of interst","Transcripts in glioma samples"),pch=c(1,1),col=c("green","red","black"))

plot(top61SignificantTranscript,top61SignificantmiRNA, main="corrolation betweeen log ratio of transcripts and miRNA", xlab="logratio of top61 Transcripts", ylab="logratio of top 61 miRNAs")

points(top61SignificantTranscript[lincRNAIdentifiedIndex],top61SignificantmiRNA[lincRNAIdentifiedIndex],col="green")
#points(top61SignificantTranscript[lincRNAOfInterstIndex],top61SignificantmiRNA[lincRNAOfInterstIndex],col="red")
points(top61SignificantTranscript[match(names(LincRNASeqHg19),names(top61SignificantTranscript))],top61SignificantmiRNA[match(names(LincRNASeqHg19),names(top61SignificantTranscript))],col="red")
legend('topleft', c("53 Nobel lincRNA Transcripts","17 lincRNA of interst","Transcripts in glioma samples"),pch=c(1,1),col=c("green","red","black"))

plot(top61SignificantGenes,top61SignificantmiRNA,main="corrolation betweeen log ratio of genes and miRNA",xlab="logratio of top61 genes",ylab="logratio of top 61 miRNAs")
dev.off()
source("corCalculater.r")

rGeneVsTranscript=corAssociation(gbmGeneExpressionMatrixDESeq2Normalized[names(top61SignificantGenes),],gbmNobelTranscriptExpressionMatrixDESeq2Normalized[names(top61SignificantTranscript),])

rGeneVsTranscript[is.na(rGeneVsTranscript)]=0
rGeneVsmiRNA=corAssociation(gbmGeneExpressionMatrixDESeq2Normalized[names(top61SignificantGenes),na.omit(match(TCGAsampleIdpharser(colnames(gbmmiRNAExpressionDataFor169SampleOfInterst)),gbmGeneExpressionSampleID))],gbmmiRNAExpressionDataFor169SampleOfInterst[names(top61SignificantmiRNA),])
rTranscriptVsmiRNA=corAssociation(gbmNobelTranscriptExpressionMatrixDESeq2Normalized[names(top61SignificantTranscript),na.omit(match(TCGAsampleIdpharser(colnames(gbmmiRNAExpressionDataFor169SampleOfInterst)),gbmGeneExpressionSampleID))],gbmmiRNAExpressionDataFor169SampleOfInterst[names(top61SignificantmiRNA),])
rTranscriptVsmiRNA[is.na(rTranscriptVsmiRNA)]=0
row.names(rGeneVsTranscript)=top61GeneNames
row.names(rGeneVsmiRNA)=top61GeneNames
rGeneVsTranscriptNet=rGeneVsTranscript>.5 | rGeneVsTranscript< -.5
rGeneVsmiRNANet=rGeneVsmiRNA>.5 | rGeneVsmiRNA< -.5
rTranscriptVsmiRNANet=rTranscriptVsmiRNA>.5 | rTranscriptVsmiRNA< -.5

rGeneVsTranscriptNetGraph=igrapher(rGeneVsTranscriptNet,"green", "red")
rGeneVsmiRNANetGraph=igrapher(rGeneVsmiRNANet,"green", "yellow")
rTranscriptVsmiRNANetGraph=igrapher(rTranscriptVsmiRNANet,"red", "yellow")
rGeneVsTranscriptNetGraph=set.edge.attribute(rGeneVsTranscriptNetGraph, "weight", index=E(rGeneVsTranscriptNetGraph), rep(50,length(E(rTranscriptVsmiRNANetGraph)) ))
rGeneVsmiRNANetGraph=set.edge.attribute(rGeneVsmiRNANetGraph, "weight", index=E(rGeneVsmiRNANetGraph), rep(70,length(E(rGeneVsmiRNANetGraph)) ))
png("corMatrix.png",  width = 600, height = 800)
par(mfrow=c(3,2))
hist(rGeneVsTranscript, main="Corolation matrix between genes and nobel transcripts")
plot(rGeneVsTranscriptNetGraph,
        layout=layout.auto,
        main="Gene-nobel transcript coexpression with r>=0.5 or r<=-0.5",
        vertex.label.dist=0,
	vertex.size=setNames(rep(35,length(V(rGeneVsTranscriptNetGraph))),V(rGeneVsTranscriptNetGraph)$name),
        vertex.frame.color="blue",
        vertex.label.color="black",
        vertex.label.font=.1,
        vertex.label=V(rGeneVsTranscriptNetGraph)$name,
        vertex.label.cex=1)

hist(rGeneVsmiRNA, main="Corolation matrix between genes and miRNA")
#E(rGeneVsmiRNANetGraph)$weight =rep(70,length(E(rGeneVsmiRNANetGraph)))
plot.igraph(rGeneVsmiRNANetGraph,
        layout=layout.fruchterman.reingold,
        main="Gene-MiRNA with r>=0.5 or r<=-0.5",
        vertex.label.dist=0,
	#vertex.size=setNames(rep(30,length(V(rGeneVsmiRNANetGraph))),V(rGeneVsmiRNANetGraph)$name),
	edge.arrow.size = 30,
	vertex.size=setNames(rep(15,length(V(rGeneVsmiRNANetGraph))),V(rGeneVsmiRNANetGraph)$name),
        vertex.frame.color="blue",
        vertex.label.color="black",
        vertex.label.font=.1,
        vertex.label=V(rGeneVsmiRNANetGraph)$name,
        vertex.label.cex=1)
hist(rTranscriptVsmiRNA, main="Corolation matrix between nobel transcripts miRNA")
plot(rTranscriptVsmiRNANetGraph,
        layout=layout.auto,
        main="Nobel transcript-MiRNA with r>=0.5 or r<=-0.5",
        vertex.label.dist=0,
	vertex.size=setNames(rep(15,length(V(rTranscriptVsmiRNANetGraph))),V(rTranscriptVsmiRNANetGraph)$name),
        vertex.frame.color="blue",
        vertex.label.color="black",
        vertex.label.font=.1,
        vertex.label=V(rTranscriptVsmiRNANetGraph)$name,
        vertex.label.cex=1)

dev.off()

lincRNAGeneCoexxpresedGraph=induced.subgraph(graph=rGeneVsTranscriptNetGraph,vids=unlist(neighborhood(graph=rGeneVsTranscriptNetGraph,order=1,nodes=V(rGeneVsTranscriptNetGraph)$name[which(V(rGeneVsTranscriptNetGraph)$name %in% LincRNA17Names)])))

top10miRNAGeneNetworkGraph=induced.subgraph(graph=rGeneVsmiRNANetGraph,vids=unlist(neighborhood(graph=rGeneVsmiRNANetGraph,order=1,nodes=V(rGeneVsmiRNANetGraph)$name[which(V(rGeneVsmiRNANetGraph)$name %in% top10DEMiRNAsNames)])))

lincRNAwithmiRNAnetworkGraph=induced.subgraph(graph=rTranscriptVsmiRNANetGraph,vids=unlist(neighborhood(graph=rTranscriptVsmiRNANetGraph,order=1,nodes=V(rTranscriptVsmiRNANetGraph)$name[which(V(rTranscriptVsmiRNANetGraph)$name %in%  LincRNA17Names)])))



# Now lets use the Mutual information based approch

top61GeneExpressionMatrix=gbmGeneExpressionMatrixDESeq2Normalized[names(top61SignificantGenes),]
top61TranscriptExpressionMatrix=gbmNobelTranscriptExpressionMatrixDESeq2Normalized[names(top61SignificantTranscript),]
top61miRNAExpressionMatrix=gbmmiRNAExpressionDataFor169SampleOfInterst[names(top61SignificantmiRNA),]

source("ARACNE.r")
png("arceneNetworks.png")
geneTranscriptNetwork=aracneAlgorithum(top61GeneExpressionMatrix,top61TranscriptExpressionMatrix,.2)
geneMiRNANetwork=aracneAlgorithum(top61GeneExpressionMatrix[,na.omit(match(TCGAsampleIdpharser(colnames(gbmmiRNAExpressionDataFor169SampleOfInterst)),gbmGeneExpressionSampleID))],top61miRNAExpressionMatrix,.4)
transcriptMiRNANetwork=aracneAlgorithum(top61TranscriptExpressionMatrix[,na.omit(match(TCGAsampleIdpharser(colnames(gbmmiRNAExpressionDataFor169SampleOfInterst)),TCGAsampleIdpharser(colnames(gbmNobelTranscriptExpressionMatrixDESeq2Normalized))))],top61miRNAExpressionMatrix,.26)

par(mfrow=c(2,2))

image(geneTranscriptNetwork, main="Gene transcript network")
image(geneMiRNANetwork, main="Gene microRNA network")
image(transcriptMiRNANetwork, main="Transcript microRNA")
dev.off()
#Network visualization with igraph package
install.packages("igraph")
library(igraph)
#Gene Transcript network visualization
row.names(geneTranscriptNetwork)=top61GeneNames
#get the coexpresses gene and transcript

igrapher=function(Network,origColor,targColor){
	net=as.data.frame(which(Network==TRUE, arr.ind=TRUE))
	origen=rownames(net)
	target=colnames(Network[,net$col])
#origen=rep(rownames(geneTranscriptNetwork), ncol(geneTranscriptNetwork))
#target=c(sapply(colnames(geneTranscriptNetwork),function(x) rep(x, nrow(geneTranscriptNetwork))))
	directed=rep("TRUE",length(target))
	NetworkGraph=data.frame(origen,target,directed)
	NetworkGraph.net=graph.data.frame(NetworkGraph,directed=F)
	V(NetworkGraph.net)$color=origColor
	V(NetworkGraph.net)[match(target,V(NetworkGraph.net)$name)]$color=targColor
#V(geneTranscriptNetworkGraph.net)$size=degree(geneTranscriptNetworkGraph.net)/10
	return(NetworkGraph.net)
}


row.names(geneTranscriptNetwork)=top61GeneNames
row.names(geneMiRNANetwork)=top61GeneNames

geneTranscriptNetwork.net=igrapher(geneTranscriptNetwork,"green", "red")
geneMiRNANetworkGraph=igrapher(geneMiRNANetwork,"green", "yellow")
transcriptMiRNANetworkGraph=igrapher(transcriptMiRNANetwork,"red", "yellow")
png("geneTranscriptInteractionNetwork.png")
#par(mfrow=c(2,2))
layout(matrix(c(1,2,3,2), 2, 2, byrow = TRUE))
plot(geneTranscriptNetwork.net,
	layout=layout.fruchterman.reingold,
	main="Gene-nobel transcript with MI >0.2",
	vertex.label.dist=0,
	edge.arrow.size = 30,
	vertex.size=25,
	vertex.frame.color="blue",
	vertex.label.color="black",
	vertex.label.font=.2,
	vertex.label=V(geneTranscriptNetwork.net)$name,
	vertex.label.cex=1)
plot(geneMiRNANetworkGraph,
        layout=layout.fruchterman.reingold,
        main="Gene-MiRNA with MI threshold of > 0.4",
        vertex.label.dist=0,
	edge.arrow.size = 30,
	vertex.size=25,
        vertex.frame.color="blue",
        vertex.label.color="black",
        vertex.label.font=.2,
        vertex.label=V(geneMiRNANetworkGraph)$name,
        vertex.label.cex=1)
plot(transcriptMiRNANetworkGraph,
        layout=layout.fruchterman.reingold,
        main="MiRNA-nobel Transcript with MI >0.26",
        vertex.label.dist=0,
	edge.arrow.size = 30,
	vertex.size=25,
        vertex.frame.color="blue",
        vertex.label.color="black",
        vertex.label.font=.2,
        vertex.label=V(transcriptMiRNANetworkGraph)$name,
        vertex.label.cex=1)

dev.off()

top10Genes=head(top61GeneNames, 11)[-4]
top10TranscriptNames=head(rownames(top61TranscriptExpressionMatrix),10)
top61Transcripts=rownames(top61TranscriptExpressionMatrix)
LincRNA17Names=names(lincRNASequencesDE)
top10DEMiRNAsNames=head(rownames(top61miRNAExpressionMatrix),10)
#subseting the graphs based on our intersts
#lets get graphs of linked with the top 10 DE genes


top10genesWithmiRNAGraph=induced.subgraph(graph=geneMiRNANetworkGraph,vids=unlist(neighborhood(graph=geneMiRNANetworkGraph,order=1,nodes=V(geneMiRNANetworkGraph)$name[which(V(geneMiRNANetworkGraph)$name %in% top10Genes)])))

top10miRNAWithGeneGraph=induced.subgraph(graph=geneMiRNANetworkGraph,vids=unlist(neighborhood(graph=geneMiRNANetworkGraph,order=1,nodes=V(geneMiRNANetworkGraph)$name[which(V(geneMiRNANetworkGraph)$name %in% top10DEMiRNAsNames)])))

lincRNAsWithGeneGraph=induced.subgraph(graph=geneTranscriptNetwork.net,vids=unlist(neighborhood(graph=geneTranscriptNetwork.net,order=1,nodes=V(geneTranscriptNetwork.net)$name[which(V(geneTranscriptNetwork.net)$name %in% LincRNA17Names)])))

allTranscriptsWithGenesGraph=induced.subgraph(graph=geneTranscriptNetwork.net,vids=unlist(neighborhood(graph=geneTranscriptNetwork.net,order=1,nodes=V(geneTranscriptNetwork.net)$name[which(V(geneTranscriptNetwork.net)$name %in%  top61Transcripts)])))

lincRNAsWithmiRNAsGraph=induced.subgraph(graph=transcriptMiRNANetworkGraph,vids=unlist(neighborhood(graph=transcriptMiRNANetworkGraph,order=1,nodes=V(transcriptMiRNANetworkGraph)$name[which(V(transcriptMiRNANetworkGraph)$name %in% LincRNA17Names)])))

top10miRNAsWithTranscripts=induced.subgraph(graph=transcriptMiRNANetworkGraph,vids=unlist(neighborhood(graph=transcriptMiRNANetworkGraph,order=1,nodes=V(transcriptMiRNANetworkGraph)$name[which(V(transcriptMiRNANetworkGraph)$name %in% top10DEMiRNAsNames)])))

top10TranscriptWithmiRNAGraph=induced.subgraph(graph=transcriptMiRNANetworkGraph,vids=unlist(neighborhood(graph=transcriptMiRNANetworkGraph,order=1,nodes=V(transcriptMiRNANetworkGraph)$name[which(V(transcriptMiRNANetworkGraph)$name %in% top10TranscriptNames)])))


png("top10GeneAndmiRNANetworks.png",  width = 160, height = 200, units = "mm", res = 300)
par(mfrow=c(3,2))
#l= layout.fruchterman.reingold(top10miRNAWithGeneGraph)*s
plot(lincRNAsWithGeneGraph,
	main="nobel lincRNA coexpressed with genes",
	edge.arrow.size=.5,
        vertex.label.font=.1,
        vertex.label.cex=1,
        layout=layout.fruchterman.reingold,
        vertex.label.dist=0,
        vertex.size=25,
        vertex.frame.color="blue",
        vertex.label.color="black")
plot(lincRNAsWithmiRNAsGraph,
	edge.arrow.size=.5,
        vertex.label.font=.1,
        vertex.label.cex=1,
        layout=layout.fruchterman.reingold,
        vertex.label.dist=0,
        vertex.size=25,
        vertex.frame.color="blue",
        vertex.label.color="black", 
	main="nobel lincRNA coexpressed with miRNA")
plot(top10genesWithmiRNAGraph,
	edge.arrow.size=.5,
	vertex.label.font=.1,
	vertex.label.cex=1,
	layout=layout.fruchterman.reingold,
	vertex.label.dist=0,
	vertex.size=25,
	vertex.frame.color="blue",
 	vertex.label.color="black",
	main="top 10 DE genes coexpressed miRNA")
#ll= layout.fruchterman.reingold(top10geneGraph)*s

plot(top10miRNAWithGeneGraph,
	edge.arrow.size=.5,
        vertex.label.font=.1,
        vertex.label.cex=1,
        layout=layout.fruchterman.reingold,
        vertex.label.dist=0,
        vertex.size=25,
        vertex.frame.color="blue",
        vertex.label.color="black",
	main="top 10 DE miRNA coexpressed with gene")
plot(lincRNAsWithmiRNAsGraph,
	edge.arrow.size=.5,
        vertex.label.font=.1,
        vertex.label.cex=1,
        layout=layout.fruchterman.reingold,
        vertex.label.dist=0,
        vertex.size=25,
        vertex.frame.color="blue",
        vertex.label.color="black",
	main="nobel lincRNAs coexpressed with miRNAs")
plot(top10miRNAsWithTranscripts, 
	edge.arrow.size=.5,
        vertex.label.font=.1,
        vertex.label.cex=1,
        layout=layout.fruchterman.reingold,
        vertex.label.dist=0,
        vertex.size=25,
        vertex.frame.color="blue",
        vertex.label.color="black",
	main="top 10 DE miRNAs coexpressed with nobel transcripts")
dev.off()
#genes that are expressed with lincRNA of interest

lincRNAsWithGeneGraph=induced.subgraph(graph=geneTranscriptNetwork.net,vids=unlist(neighborhood(graph=geneTranscriptNetwork.net,order=1,nodes=V(geneTranscriptNetwork.net)$name[which(V(geneTranscriptNetwork.net)$name %in% LincRNA17Names)])))
 
allTranscriptsWithGenesGraph=induced.subgraph(graph=geneTranscriptNetwork.net,vids=unlist(neighborhood(graph=geneTranscriptNetwork.net,order=1,nodes=V(geneTranscriptNetwork.net)$name[which(V(geneTranscriptNetwork.net)$name %in%  top61Transcripts)])))
pdf("lincRNAWithGenesGraph.pdf")
par(mfrow=c(1,2))
plot(lincRNAsWithGeneGraph, main="lincRNA expressed with genes")
plot(allTranscriptsWithGenesGraph, main="genes expressed with nobel transcripts")
dev.off()

png("lincRNAWithGenesGraph.png")
par(mfrow=c(1,2))
plot(lincRNAsWithGeneGraph, main="lincRNA expressed with genes")
plot(allTranscriptsWithGenesGraph, main="genes expressed with transcripts")
dev.off()
##transcripts that are expressed with miRNAS

top10miRNAsWithTranscripts=induced.subgraph(graph=transcriptMiRNANetworkGraph,vids=unlist(neighborhood(graph=transcriptMiRNANetworkGraph,order=1,nodes=V(transcriptMiRNANetworkGraph)$name[which(V(transcriptMiRNANetworkGraph)$name %in% top10DEMiRNAsNames)])))

lincRNAsWithmiRNAsGraph=induced.subgraph(graph=transcriptMiRNANetworkGraph,vids=unlist(neighborhood(graph=transcriptMiRNANetworkGraph,order=1,nodes=V(transcriptMiRNANetworkGraph)$name[which(V(transcriptMiRNANetworkGraph)$name %in% LincRNA17Names)])))
top10TranscriptWithmiRNAGraph=induced.subgraph(graph=transcriptMiRNANetworkGraph,vids=unlist(neighborhood(graph=transcriptMiRNANetworkGraph,order=1,nodes=V(transcriptMiRNANetworkGraph)$name[which(V(transcriptMiRNANetworkGraph)$name %in% top10TranscriptNames)])))

pdf("transcriptsWithmiRNA.pdf")
layout(matrix(c(1,2,3,2), 2, 2, byrow = TRUE))
plot(top10miRNAsWithTranscripts, main="sub-network of top 10 miRNA with nobel transcripts")
plot(lincRNAsWithmiRNAsGraph, main="sub-network of the lincRNAs with miRNAs")
plot(top10TranscriptWithmiRNAGraph, main="sub-networks of the top 10 DE nobel transcripts with miRNA")
dev.off()

png("transcriptsWithmiRNA.png", heigh=400, width=600)
layout(matrix(c(1,2,3,2), 2, 2, byrow = TRUE))
plot(top10miRNAsWithTranscripts, main="sub-network of top 10 miRNA with nobel transcripts")
plot(lincRNAsWithmiRNAsGraph, main="sub-network of the nobel lincRNAs with miRNAs")
plot(top10TranscriptWithmiRNAGraph, main="sub-networks of the top 10 DE nobel transcripts with miRNA")
dev.off()


#Gene miRNA interaction visualization
#geneMiRNAnet=which(geneMiRNANetwork==TRUE, arr.ind=TRUE)
#geneMiRNAnet=as.data.frame(geneMiRNAnet)

#origen1=rownames(geneMiRNANetwork[geneMiRNAnet$row,geneMiRNAnet$col])
#target1=colnames(geneMiRNANetwork[geneMiRNAnet$row,geneMiRNAnet$col])
#directed1=rep(TRUE, length(origen1))
#geneMiRNANetworkGraph=data.frame(origen1,target1,directed1)
#geneMiRNANetworkGraph.net=graph.data.frame(geneMiRNANetworkGraph,directed=F)
#V(geneMiRNANetworkGraph.net)$color="green"
#V(geneMiRNANetworkGraph.net)[V(geneMiRNANetworkGraph.net)$name[grep("^[h]s.*", V(geneMiRNANetworkGraph.net)$name)]]$color="yellow"


#visualizing Transcript RNA interaction

transcriptmiRNANet=which(transcriptMiRNANetwork==TRUE, arr.ind=TRUE)
transcriptmiRNANet=as.data.frame(transcriptmiRNANet)
origen2=rownames(transcriptmiRNANet[transcriptmiRNANet$row,transcriptmiRNANet$col])
target2=colnames(transcriptMiRNANetwork[transcriptmiRNANet$row,transcriptmiRNANet$col])
directed2=rep(TRUE,length(target2))
transcriptMiRNANetworkGraph=data.frame(origen2,target2,directed2)
transcriptMiRNANetworkGraphNet=graph.data.frame(transcriptMiRNANetworkGraph, directed=F)
V(transcriptMiRNANetworkGraphNet)$color="yellow"
V(transcriptMiRNANetworkGraphNet)[V(transcriptMiRNANetworkGraphNet)$name[grep("^[T]C.*", V(transcriptMiRNANetworkGraphNet)$name)]]$color="red"

png("TranscriptMiRNAInteractionNetwork.png")
plot(transcriptMiRNANetworkGraphNet,
        layout=layout.lgl,
        main="Transcript-MiRNA network",
        vertex.label.dist=0,
        vertex.frame.color="blue",
        vertex.label.color="black",
        vertex.label.font=1,
        vertex.label=V(transcriptMiRNANetworkGraphNet)$name,
        vertex.label.cex=.75)
dev.off()



################################################################################
#										#
#					Enrichment analysis			#
################################################################################

source("enrichment.r")

rankedEnsembolGeneIds=names(orderedLogRatioOfHighOverLowExpressedGeneIn5pcOfQuartiles)

rankedEntrezIds=ensembolToEntrezId(rankedEnsembolGeneIds)
rankedEntrezIdsNULLreplacedWithZero=replace(rankedEntrezIds, sapply(rankedEntrezIds, is.null), 0)
names(rankedEntrezIdsNULLreplacedWithZero)=rankedEnsembolGeneIds


ensembolToGeneSymbolConverted=ensembolToGeneSymbol(rankedEnsembolGeneIds)

#import the downloaded genesets
allGeneSets=formatGeneSets((datasets("msigdb/")))
enrichmentResults=fisherEnrichmentAnalysis(rankedEntrezIdsNULLreplacedWithZero,61,allGeneSets,.05)

max.row=max(unlist(lapply(enrichmentResults,length)))


EnrichmentResult=data.frame(KEGG_enrichment=c(enrichmentResults[[1]],rep(NA, max.row-length(enrichmentResults[[1]]))),miRNA_enrichment=c(enrichmentResults[[2]],rep(NA, max.row-length(enrichmentResults[[2]]))),tf_enrichment=c(enrichmentResults[[3]],rep(NA, max.row-length(enrichmentResults[[3]]))),GO_bpEnrichment=c(enrichmentResults[[4]],rep(NA, max.row-length(enrichmentResults[[4]]))),GOCC_enrichment=c(enrichmentResults[[5]],rep(NA, max.row-length(enrichmentResults[[5]]))),GOMF_enrichment=c(enrichmentResults[[6]],rep(NA, max.row-length(enrichmentResults[[6]]))))








top61GeneNames=head(unname(noquote(ensembolToGeneSymbolConverted)),61)





