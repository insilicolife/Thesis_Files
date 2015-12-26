if("DESeq2" %in% rownames(installed.packages()) == FALSE){
	source("http://bioconductor.org/biocLite.R")
	biocLite("DESeq2")
}
library("DESeq2")

DECalculaterDESeq2=function(countData, colData, condition){
	DEDatasetFormat= DESeqDataSetFromMatrix(countData=countData, colData=colData,design=~ condition)
	DESeqObject=DESeq(DEDatasetFormat)
	DEresults=results(DESeqObject)
	return(list(DESeqObject,DEresults,DEDatasetFormat))
}

DESignificanceDESeq2= function(DEresults,pvalThreshold){

        DEresultsAdjPvalOrdered=DEresults[order(DEresults$padj),]
        DEresultsAdjPvalOrderedSignificantlyDifferentGeneIndex=which(DEresultsAdjPvalOrdered$padj<pvalThreshold)
        DEresultsAdjPvalOrderedSignificantlyDifferentGene=DEresultsAdjPvalOrdered[DEresultsAdjPvalOrderedSignificantlyDifferentGeneIndex,]
        #dim(DEresultsAdjPvalOrderedSignificantlyDifferentGene)
        gbmGeneExpressionMatrixSignificant=gbmGeneExpressionMatrix[DEresultsAdjPvalOrderedSignificantlyDifferentGeneIndex,]

        return(list(DEresultsAdjPvalOrderedSignificantlyDifferentGene,gbmGeneExpressionMatrixSignificant))


}

DESeq2Normalization=function(countData, colData,condition){

	myDEDatasetFormat=DEDatasetFormat=DESeqDataSetFromMatrix(countData=countData, colData=colData,design=~ condition)
	myDEDatasetFormat=estimateSizeFactors(myDEDatasetFormat)
        sizeFactors(myDEDatasetFormat)
        DESeqNormalizedDataMatrix=counts(myDEDatasetFormat, normalized=TRUE)
	
	return(DESeqNormalizedDataMatrix)
}
regularizedlogTransform=function(DEResult){

        DEResultCoutDataTransformationRlog=rlog(DEResult,fast=TRUE)
        DEResultCoutDataTransformationRlogMatrix=assay(DEResultCoutDataTransformationRlog)
        return(DEResultCoutDataTransformationRlogMatrix)

}

diagnosticPlot=function(DEResult,DEresults){
	
	#mfrow(1,3)
	plotMA(DEresults,ylim=c(-3,3), main="MA-plot")
	plotDispEsts(DEResult, ylim=c(1e-6,1e1), maim="Desperssion plot")
	hist(DEresults$pvalue,col="grey", breaks=20, main="histogram of p-values")
	
	
}



