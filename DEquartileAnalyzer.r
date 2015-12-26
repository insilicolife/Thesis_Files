DEAnalyzerBasedOnTheQuartiles=function(dataMatrix,Normalization=FALSE,quartileContFilter=TRUE){
#DEAnalyzerBasedOnTheQuartiles=function(dataMatrix, geneLength){
	#topLowExpressionInFivePercentSamplesInrowCount=apply(sortedByColDataMatrix[,1:lowestExpressionWindowLenght],1,mean)

	#filtering noise
	
	#filteredGeneIndex=which(topLowExpressionInFivePercentSamplesInrowCount>100)
	#dataMatrixNoisRemoved=dataMatrix[filteredGeneIndex,]
	if(Normalization==FALSE){
		#filteredGeneIndexByRowMean=which(apply(dataMatrix,1,mean)>100)
		#dataMatrixNoisRemovedByRowMean=dataMatrix[filteredGeneIndexByRowMean,]
	#Now we can normalize the filtered dataset
	#source("playWithGbmData.r")
	#RPKMNormalizedDataMatrix=RPKMNormalization(dataMatrixNoisRemoved,geneLength[filteredGeneIndex])
		source("DECalculaterDESeq2.r")
		cond=c(rep("unknown Sample1",69),rep("unknown Sample2",100))
		colData=data.frame(condition=cond)
	#rownames(colData)=colnames(dataMatrixNoisRemoved)
		rownames(colData)=colnames(dataMatrix)
		DESeq2NormalizedDataMatrix=DESeq2Normalization(dataMatrix,colData,"condition")
	#sortedByColDataMatrixNormalizedRPKM=t(apply(RPKMNormalizedDataMatrix, 1, sort))
		
		sortedByColDataMatrixNormalizedDESeq2=t(apply(DESeq2NormalizedDataMatrix, 1, sort))
		orderedSampleIndex=t(apply(DESeq2NormalizedDataMatrix, 1, order))
	}else if(Normalization==TRUE){
	
		sortedByColDataMatrixNormalizedDESeq2=t(apply(dataMatrix, 1, sort))
		orderedSampleIndex=t(apply(dataMatrix, 1, order))
				
	}else{
		print("error")
	}
	if(quartileContFilter==FALSE){
		filteredGeneIndexByRowMean=which(apply(sortedByColDataMatrixNormalizedDESeq2,1,mean)>100)
	}
	else{
		
		filteredGeneIndexByRowMean=which(apply(sortedByColDataMatrixNormalizedDESeq2,1,mean)>0)	
	}
        dataMatrixNoisRemovedByRowMean=sortedByColDataMatrixNormalizedDESeq2[filteredGeneIndexByRowMean,]
	orderedSampleIndex=orderedSampleIndex[filteredGeneIndexByRowMean,]
	#print(orderedSampleIndex[1,])
	lowestExpressionWindowLenght=round(ncol(dataMatrixNoisRemovedByRowMean)*.05)
	
	#topLowExpressionInFivePercentSamples=apply(sortedByColDataMatrixNormalizedRPKM[,1:lowestExpressionWindowLenght],1,mean)
	
	topLowExpressionInFivePercentSamplesDESeq2=apply(dataMatrixNoisRemovedByRowMean[,1:lowestExpressionWindowLenght],1,mean)
	#noiseRemovedIndex=which(topLowExpressionInFivePercentSamplesDESeq2>100)
	#topLowExpressionInFivePercentSamplesDESeq2NoiseREmoved=topLowExpressionInFivePercentSamplesDESeq2[noiseRemovedIndex]

	highestExpressionWindowLengthStart=(ncol(dataMatrixNoisRemovedByRowMean)-lowestExpressionWindowLenght)+1
	
 	#topHighExpressionInFivePercentSamples=apply(sortedByColDataMatrixNormalizedRPKM[,highestExpressionWindowLengthStart:ncol(sortedByColDataMatrixNormalizedRPKM)],1,mean)
	topHighExpressionInFivePercentSamplesDESeq2=apply(dataMatrixNoisRemovedByRowMean[,highestExpressionWindowLengthStart:ncol(dataMatrixNoisRemovedByRowMean)],1,mean)
	#topHighExpressionInFivePercentSamplesDESeq2NoiseREmoved=topHighExpressionInFivePercentSamplesDESeq2[noiseRemovedIndex]

	#logRatioTophighOverTopLowPRKM=log((topHighExpressionInFivePercentSamples)/(topLowExpressionInFivePercentSamples),2)
	if(quartileContFilter==TRUE){
		
		quartilefilter=which(topLowExpressionInFivePercentSamplesDESeq2>100)
		orderedSampleIndex=orderedSampleIndex[quartilefilter,]
		logRatioTophighOverTopLowDESeq2=log((topHighExpressionInFivePercentSamplesDESeq2[quartilefilter])/(topLowExpressionInFivePercentSamplesDESeq2[quartilefilter]),2)
	}else{
		quartilefilter=which(topLowExpressionInFivePercentSamplesDESeq2>=0)	
		orderedSampleIndex=orderedSampleIndex[quartilefilter,]
		logRatioTophighOverTopLowDESeq2=log((topHighExpressionInFivePercentSamplesDESeq2[quartilefilter]+1)/(topLowExpressionInFivePercentSamplesDESeq2[quartilefilter]+1),2)
		
	}
		#print(dim(orderedSampleIndex))
		#print(colnames(dataMatrix)[orderedSampleIndex[1,1:lowestExpressionWindowLenght]])
		#quartileClassLables=orderedSampleIndex
		#lowestExpressedSampleGroup=
	lowExpressedSampleGroupForGene=list()
	highExpressedSampleGroupForGene=list()
	for(gene in 1:nrow(orderedSampleIndex)){

	lowExpressedSampleGroupForGene[gene]=list(colnames(dataMatrix)[orderedSampleIndex[gene,1:lowestExpressionWindowLenght]])
	highExpressedSampleGroupForGene[gene]=list(colnames(dataMatrix)[orderedSampleIndex[gene,highestExpressionWindowLengthStart:ncol(dataMatrixNoisRemovedByRowMean)]])

}

		quartileClassLables=list(lowExpressedSampleGroupForGene,highExpressedSampleGroupForGene)

		
	#logRatioTophighOverTopLowDESeq2=log((topHighExpressionInFivePercentSamplesDESeq2)/(topLowExpressionInFivePercentSamplesDESeq2),2)
	return(list(logRatioTophighOverTopLowDESeq2,quartileClassLables))
}









RPKMNormalization=function(dataMatrix,genesLength){

        RPKM=(dataMatrix/(genesLength*sum(dataMatrix)))*10e9
        return(RPKM)
}



transcriptToGeneMapping=function(geneInfoMatrix,transcriptMatrix){
	
	GeneTranscriptChromosomalLevelMatchingIndexing=match(geneInfoMatrix$CHROM, sapply(strsplit(rownames(transcriptMatrix),"-"), function(a)a[2]))
	GeneTranscriptChromosomalLevelMatching=rownames(transcriptMatrix[GeneTranscriptChromosomalLevelMatchingIndexing,])
       names(GeneTranscriptChromosomalLevelMatching)=geneInfoMatrix$NAME


	transcriptStartPosition=sapply(strsplit(GeneTranscriptChromosomalLevelMatching, "-"),function(a)a[3])
	

	geneStartPosition=geneInfoMatrix[which(geneInfoMatrix$NAME==names(GeneTranscriptChromosomalLevelMatching)),]$START

	geneEndPosition=geneInfoMatrix[which(geneInfoMatrix$NAME==names(GeneTranscriptChromosomalLevelMatching)),]$END

	geneLevelTranscriptMapping=((geneStartPosition<=transcriptStartPosition) & (transcriptStartPosition <=geneEndPosition))
	geneTranscriptOverlaped=geneLevelTranscriptMapping[which(geneLevelTranscriptMapping==TRUE)]
	geneTranscriptChromosomalMapping=GeneTranscriptChromosomalLevelMatching[which(geneLevelTranscriptMapping==TRUE)]
	return(c(geneTranscriptOverlaped, geneTranscriptChromosomalMapping))

}



