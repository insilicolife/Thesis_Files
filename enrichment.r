#download gene sets
datasets=function(directory){
#filepath="msigdb/"
	file_list=list.files(directory)
#print(file_list)
	geneSetMatrix=list()
i=1
	for(file in file_list){
	
		
		tempmaxCol=max(count.fields(paste(directory,file, sep=""), sep="\t"))
		#print(tempmaxCol)
		tempGeneSet=read.table(paste(directory,file, sep=""), header=FALSE, sep="\t", fill=TRUE, col.names=1:tempmaxCol)
		geneSetMatrix[[i]]=tempGeneSet
		rm(tempmaxCol)
		rm(tempGeneSet)
		i=i+1
	}
	return(geneSetMatrix)
}
#format gene sets
#geneSets=datasets("msigdb/")
formatGeneSets=function(geneSets){
	formatedGeneSets=list()
	i=1
	for(geneSet in geneSets){

		formatedGeneSet=geneSet[,-c(1:2)]
		rownames(formatedGeneSet)=geneSet[,1]
	
		formatedGeneSets[[i]]=formatedGeneSet	
		i=i+1
	}	
	return(formatedGeneSets)
}
#prepare rank ordered corrolations 
#source("GSEAcalculater.r")
#source("DEAnalyzer.r")

#geneLogratio=orderedLogRatioOfHighOverLowExpressedGeneIn5pcOfQuartiles
#lowExpressionSamplesForGenes=quartileSampleGroupsLowExpression
#highExpressionSamplesForGenes=quartileSampleGroupsHighExpression
#dataSet=gbmGeneExpressionMatrixDESeq2Normalized

#rankedCorrolation=phynotipicCorolation(dataSet,geneLogratio,lowExpressionSamplesForGenes,highExpressionSamplesForGenes)

#calculate enrichment score

#genes=as.list(names(rankedCorrolat
#entresGeneSymbol=function(gene){
#library('biomaRt')
#mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
#G_list=getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","entrezgene","hgnc_symbol"),values=gene,mart= mart)
#return(G_list)
#}
#ensembolGeneIds=names(rankedCorrolation)
#entresGeneSymbolList=lapply(ensembolGeneIds,entresGeneSymbol)
#entresGeneSymbols=do.call(rbind, lapply(entresGeneSymbolList, rbind))
#entresGeneSymbolDataframe=lapply(entresGeneSymbolList,rbind.fill)




ensembolToEntrezId=function(ensembolIds){

if("org.Hs.eg.db" %in% rownames(installed.packages()) == FALSE){
	source("http://bioconductor.org/biocLite.R")
	biocLite("org.Hs.eg.db")
}
	library("org.Hs.eg.db")

	entreazConversion=as.list(org.Hs.egENSEMBL2EG)
	
	#head(entreazConversion)
	return(entreazConversion[ensembolIds])
}

# get the unmapped ensembole ids and chnage it to gene symbole/id then try to find entrreze id for it
#rankedEntrezIds=ensembolToEntrezId(ensembolGeneIds)

#rankedEntrezIdsNULLreplacedWithZero=replace(rankedEntrezIds, sapply(rankedEntrezIds, is.null), 0)
#names(rankedEntrezIdsNULLreplacedWithZero)=ensembolGeneIds

ensembolToGeneSymbol=function(ensembolIds){

if("biomaRt" %in% rownames(installed.packages()) == FALSE){
source("http://bioconductor.org/biocLite.R")
biocLite("biomaRt")
}
library("biomaRt")
	ensembl=useMart("ensembl")
	ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
	filters = listFilters(ensembl)
	ensemboIdToGeneSymbol=getBM(attributes=c("ensembl_gene_id","hgnc_symbol"), filters = "ensembl_gene_id", values=ensembolIds, mart=ensembl)		
	#print(head(ensemboIdToGeneSymbol))
	ensemboIdToGeneSymbolOnly=ensemboIdToGeneSymbol[,-1]
	names(ensemboIdToGeneSymbolOnly)=ensemboIdToGeneSymbol[,-2]
	ensembolIdorderedInsex=match(ensembolIds,names(ensemboIdToGeneSymbolOnly))
	ensemboIdToGeneSymbolList=ensemboIdToGeneSymbolOnly[ensembolIdorderedInsex]
	ensemboIdToGeneSymbolList1=ensemboIdToGeneSymbolList
	names(ensemboIdToGeneSymbolList1)=ensembolIds
	return(ensemboIdToGeneSymbolList1)
}

#ensembolToGeneSymbolConverted=ensembolToGeneSymbol(ensembolGeneIds)
source("GSEAcalculater.r")
#get gene sets
#allGeneSets=formatGeneSets((datasets("msigdb/")))
#singleGeneSets=allGeneSets[[1]]
#mainES=c()
#mainPermituationScore=list()
#for(i in 1:length(singleGeneSets)){
#for(i in 1:nrow(singleGeneSets)){
	#mainES[i]=ES(rankedCorrolation,rankedEntrezIdsNULLreplacedWithZero,as.numeric(singleGeneSets[i,]))[[2]]
#	mainPermituationScore[i]=ESsmallpermituation(rankedCorrolation,rankedEntrezIdsNULLreplacedWithZero,as.numeric(singleGeneSets[i,]),10)
#
# now calculate p-value


#for(i in 1:nrow(singleGeneSets)){
#for(i in 1:nrow(singleGeneSets)){
        #mainES[i]=ES(rankedCorrolation,rankedEntrezIdsNULLreplacedWithZero,as.numeric(singleGeneSets[i,]))[[2]]
#        mainPermituationScore[i]=ESsmallpermituation(rankedCorrolation,rankedEntrezIdsNULLreplacedWithZero,as.numeric(singleGeneSets[i,]),2500)
#}

fisherEnrichmentAnalysis=function(rankedGeneListsWithEntrezIds,lengthOfTopDEgenes, TotalGenesets,pvalueThreshold){
enrichment=list()
enrichmentPvalue=list()
for(j in 1:length(TotalGenesets)){
	singleGeneSets=TotalGenesets[[j]]
	enriched=c()
	for(i in 1:nrow(singleGeneSets)){

		enrichmentPvalue[i]=EnrichmentFisherTest(rankedGeneListsWithEntrezIds,lengthOfTopDEgenes,as.numeric(singleGeneSets[i,]))

		if(enrichmentPvalue[i]<.05){

			enriched=c(enriched, rownames(singleGeneSets[i,]))
			#enriched=c(enriched, singleGeneSets[i,])
		}

	
	}
	enrichment[j]=list(enriched)
	}
	
	return(enrichment)
}

#EnrichmentResult=data.frame(KEGG_enrichment=c(enrichment[[1]],rep(NA, max.row-length(enrichment[[1]]))),miRNA_enrichment=c(enrichment[[2]],rep(NA, max.row-length(enrichment[[2]]))),tf_enrichment=c(enrichment[[3]],rep(NA, max.row-length(enrichment[[3]]))),GO_bpEnrichment=c(enrichment[[4]],rep(NA, max.row-length(enrichment[[4]]))),GOCC_enrichment=c(enrichment[[5]],rep(NA, max.row-length(enrichment[[5]]))),GOMF_enrichment=c(enrichment[[6]],rep(NA, max.row-length(enrichment[[6]]))))

#write.table(EnrichmentResult,file="enrichmentResultTable.xsl", col.names="TRUE", sep="\t")




