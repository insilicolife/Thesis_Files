phynotipicCorolation=function(matrix, logratio,lowExpressionSampleGroup, highExpressionSampleGroup){

source("corCalculater.r")
corvalues=c()
i=1
	for(gene in names(logratio)){
	
	corvalues[i]=corolations(matrix[gene,lowExpressionSampleGroup[[gene]]],matrix[gene,highExpressionSampleGroup[[gene]]])

 		i=i+1

	}
	names(corvalues)=names(logratio)
	return(corvalues)
}


ES=function(rankedGeneListWithCorrolation,rankedEntrezIdList,geneSets){

	phit=c(rep(0, length(rankedEntrezIdList)))
	pmiss=c(rep(0, length(rankedEntrezIdList)))

	for(i in 1:length(rankedEntrezIdList)){
	
		if(!isEmpty(intersect(as.numeric(unlist(rankedEntrezIdList[[names(rankedEntrezIdList)[i]]])),geneSets))){
			
			#print(as.numeric(unlist(rankedEntrezIdList[[names(rankedEntrezIdList)[i]]])))
			if(i==1){
				phit[i]=phit[i]+rankedGeneListWithCorrolation[i]
                                pmiss[i]=pmiss[i]+0
			}else{
				phit[i]=phit[i-1]+rankedGeneListWithCorrolation[i]	
				pmiss[i]=pmiss[i-1]
				}
		}else{

			#print(geneSets)

				if(i==1){	

					pmiss[i]=1
                              
				}else{
					pmiss[i]=pmiss[i-1]+1
					phit[i]=phit[i-1]
				}
			}
	
	}
	#print(pmiss[1:10])
	#print(phit[1:10])
	if(sum(phit)>0){
		phit=phit/sum(phit)}
	pmiss=pmiss/(length(rankedGeneListWithCorrolation)-length(geneSets))

	Escore=phit-pmiss
	Es=max(abs(Escore))	
	return(list(Escore,Es))
	
}
ESpermituation=function(rankedGeneListWithCorrolation,rankedEntrezIdList,geneSets, perm=100){


	permMatrix=matrix(data=NA,ncol=length(rankedEntrezIdList),nrow=perm)
	permituatedES=c()
	for(i in 1:perm){
	
		permituate<-sample(1:length(rankedEntrezIdList),length(rankedEntrezIdList))
		newRankedEntrezIdList=rankedEntrezIdList[permituate]
		newRankedGeneListWithCorrolation=rankedGeneListWithCorrolation[permituate]
		permituatedES[i]=ES(newRankedGeneListWithCorrolation,newRankedEntrezIdList,geneSets)[[2]]
		#permMatrix[i,]=ES(newRankedGeneListWithCorrolation,newRankedEntrezIdList,geneSets)[[1]]
	}
	

}

indicator=function(condition) ifelse(condition,1,0)

ESsmallpermituation=function(rankedGeneListWithCorrolation,rankedEntrezIdList,geneSets, perm=10){

	xo=ES(rankedGeneListWithCorrolation,rankedEntrezIdList,geneSets)[[2]]
	N=perm
	yn=c()
	for(i in 1:N){
		
		permt=sample(1:length(rankedEntrezIdList),N)
		newRankedEntrezIdList=rankedEntrezIdList[permt]
                newRankedGeneListWithCorrolation=rankedGeneListWithCorrolation[permt]
		#print(permt)
		#print(newRankedEntrezIdList)
		#print(newRankedGeneListWithCorrolation)
		yn[i]=ES(newRankedGeneListWithCorrolation,newRankedEntrezIdList,geneSets)[[2]]
	}
	#print(yn)
	M=sum(indicator(abs(yn)>xo))
	print(M)
	if(M>10){
		#two-side testing p-value
		Pecdf=(1+mean(indicator(abs(yn)>=xo)))/length(yn)
		#print(Pecdf)
		return(list(Pecdf,xo))
	}else{
		ynSorted=sort(yn)
		Nexc=length(yn)-round(.05*length(yn))
		zi=ynSorted[which(abs(ynSorted)>xo)]-xo
		t=(ynSorted[Nexc]+ynSorted[Nexc+1])/2
		#ML estimate of the gpd parameters
		print(zi)	
		print(t)	

	if("fExtremes" %in% rownames(installed.packages()) == FALSE){
        
		install.packages("fExtremes")
	}
		
	if("ismev" %in% rownames(installed.packages()) == FALSE){

                install.packages("ismev")
	}
		library("fExtremes")
		library("ismev")
		if(length(zi)>2){
			modelFit=gev.fit(zi)
			modelFitParameters=modelFit$mle
			print(modelFitParameters)
			Pgpd=Nexc*(1-pgpd(xo-t,xi=modelFitParameters[2],mu=modelFitParameters[1],beta=modelFitParameters[3]))/length(ynSorted)
		}else{
			Pgpd=Nexc*(1-pgpd(xo-t))/length(ynSorted)
		}
		return(list(Pgpd[[1]][1],xo))
	}
	
}







EnrichmentFisherTest=function(orderedGeneList, cutoff, geneSet){

	geneSet=geneSet[!is.na(geneSet)]
	geneList_GenesetIntersection=length(intersect(as.numeric(unlist(orderedGeneList[1:cutoff])),geneSet))
	totalGeneset=length(geneSet)
	geneList=cutoff
	allGenes=length(orderedGeneList)
	fisher=matrix(c(geneList_GenesetIntersection,totalGeneset,geneList,allGenes),ncol=2)
	fisher=t(fisher)
	fisherPvalue=fisher.test(fisher)$p.value
	#print(fisher)
	#rowSum=apply(fisher,1,sum)
	#colSum=apply(fisher,2,sum)
	#rowFactorial=lapply(rowSum, lfactorial)
	#colFactorail=lapply(colSum, lfactorial)
	#print(rowFactorial)
	#print(colFactorail)
	#productOfRowFactorial=prod(unlist(rowFactorial))
	#productOfColFactorial=prod(unlist(colFactorail))
	#N=lapply(rowSum,sum)
	#Nfactorial=unlist(lfactorial(N))
	#factorialOfEachElement=lfactorial(fisher)
	#productOfFactorialOfEachElement=prod(exp(factorialOfEachElement))

	#fisherPvalue=productOfRowFactorial*productOfColFactorial/(productOfFactorialOfEachElement*Nfactorial)	
	#rbind(fisher,geneList,allGenes )
	return(fisherPvalue)
}









