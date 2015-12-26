
getPeomoterAndGeneSeq=function(topRnankedEnsembolGeneIds,ensembl){

        newGbmGeneInfo=getBM(attributes=c("ensembl_gene_id","chromosome_name","start_position","end_position","hgnc_symbol"), filters = "ensembl_gene_id", values=topRnankedEnsembolGeneIds, mart=ensembl)

        newGbmGeneInfofilteredByGeneSymbols=newGbmGeneInfo[which(newGbmGeneInfo$hgnc_symbol!=""|newGbmGeneInfo$hgnc_symbol!=NA),]
        newGbmGeneInfoGeneSymbols=newGbmGeneInfofilteredByGeneSymbols$hgnc_symbol
        newGbmGeneInfoWithOutGeneSymbole=newGbmGeneInfo[which(newGbmGeneInfo$hgnc_symbol==""|newGbmGeneInfo$hgnc_symbol==NA),]
        newGbmGeneInfoWithOutGeneSymboleEnsebolId=newGbmGeneInfoWithOutGeneSymbole$ensembl_gene_id

        chr=newGbmGeneInfoWithOutGeneSymbole$chromosome_name
        start=newGbmGeneInfoWithOutGeneSymbole$start_position
        end=newGbmGeneInfoWithOutGeneSymbole$end_position

        #retriving the upstream sequences using Gene symbol

        upstreamSeqBySymbol=getSequence(id=newGbmGeneInfoGeneSymbols,type="hgnc_symbol", seqType = "coding_gene_flank", upstream=2000, mart = ensembl)
        GeneSeq=getSequence(id=newGbmGeneInfoGeneSymbols,type="hgnc_symbol", seqType = "gene_exon_intron",mart = ensembl)

        return(list(upstreamSeqBySymbol,GeneSeq,newGbmGeneInfoWithOutGeneSymbole))

}
getBackgroundModel=function(geneSeq){


        notDuplicatedDNAseq=geneSeq[!duplicated(geneSeq$hgnc_symbol),1]
        names(notDuplicatedDNAseq)=geneSeq[!duplicated(geneSeq$hgnc_symbol),2]
        DNAStrings=DNAStringSet(notDuplicatedDNAseq)
        nucliotideCount=alphabetFrequency(DNAStrings,baseOnly=TRUE)[,1:4]
        rowSum=apply(nucliotideCount,1,sum)
        nucliotideFrequency=nucliotideCount/rowSum
        rownames(nucliotideFrequency)=names(notDuplicatedDNAseq)

        return(nucliotideFrequency)
}
getSingleBackground=function(seq){

        #library(Biostrings)
        DNAStrings=DNAStringSet(seq)
        nucliotideCount=alphabetFrequency(DNAStrings,baseOnly=TRUE)[,1:4]
        nucliotideSum=sum(nucliotideCount)
        nucliotideFrequency=nucliotideCount/nucliotideSum
        #rownames(nucliotideFrequency)=names(seq)

        return(nucliotideFrequency)
}

getSingleRNABackground=function(seq){

        #library(Biostrings)
        RNAStrings=RNAStringSet(seq)
        nucliotideCount=alphabetFrequency(RNAStrings,baseOnly=TRUE)[,1:4]
        nucliotideSum=sum(nucliotideCount)
        nucliotideFrequency=nucliotideCount/nucliotideSum
        #rownames(nucliotideFrequency)=names(seq)

        return(nucliotideFrequency)
}



getPFMhPDI=function(directoryName){
        library(stringr)
        filePath=directoryName
        file_list=list.files(filePath)
        #print(length(file_list))
        for(file in file_list){
                if (!exists("PFM")){

                rowPFM=read.table(paste(filePath,file, sep=""), header=FALSE, sep="\t")
                PFM=rowPFM[,-1]
                rownames(PFM)=rowPFM[,1]
                rm(rowPFM)
                PFM=list(PFM)
                }

                else if(exists("PFM")){

                row_temp_PFM=read.table(paste(filePath,file, sep=""), header=FALSE, sep="\t")
                temp_PFM=row_temp_PFM[,-1]
                rownames(temp_PFM)=row_temp_PFM[,1]
                temp_PFM=list(temp_PFM)
                PFM=c(PFM,temp_PFM)
                rm(temp_PFM)
                rm(row_temp_PFM)
                }else{
                        print("There might be some error")
                }
        }

        proteins=word(file_list,1, sep=fixed('.'))
        names(PFM)=proteins
        return(PFM)

}

getRNACpmpleteRBPpfm=function(directoryName){
	
	#load the RNAComplete Database PFM
	#XLConnect is the package that loads xlsx file
	if("XLConnect" %in% rownames(installed.packages()) == FALSE){
	install.packages("XLConnect", repos="http://ftp.acc.umu.se/mirror/CRAN/")
	}
	library("XLConnect")
	library(stringr)
        filePath=directoryName
        file_list=list.files(filePath)
	for(file in file_list){
                if (!exists("RNACompleteRBPpfm")){
		
                RNACompleteRBPpfm=list(t(read.table(paste(filePath,file, sep=""), header=TRUE, sep="")[,-1]))
                #RBPPFM=rowPFM
                #rownames(RBPPFM)=c("A","C","G","U")
                #RBPPFM=list(RBPPFM)
                #rm(rowPFM)
                }
                else if(exists("RNACompleteRBPpfm")){

                row_temp_PFM=list(t(read.table(paste(filePath,file, sep=""), header=TRUE, sep="")[,-1]))
                #temp_PFM=row_temp_PFM
                #rownames(temp_PFM)=c("A","C","G","U")
                #temp_PFM=list(temp_PFM)
                RNACompleteRBPpfm=c(RNACompleteRBPpfm,row_temp_PFM)
                #rm(temp_PFM)
                rm(row_temp_PFM)
                }else{
                        print("There might be some error")
                }
        }

	RNACompleteInfo=loadWorkbook("../OfficeFiles/Tools/RNACompletePFM/SupplementaryData1_RNAcompete_master_file.xlsx")
	setMissingValue(RNACompleteInfo, value = "NA")
	motifIdProteinMap=readWorksheet(RNACompleteInfo, sheet = "Master List--Plasmid info")
	names(RNACompleteRBPpfm)=motifIdProteinMap$Motif_ID
	motifIdProteinMapHS=motifIdProteinMap[motifIdProteinMap$Species=="Homo_sapiens",]
	RNACompleteRBPpfmHS=RNACompleteRBPpfm[motifIdProteinMapHS$Motif_ID]
	names(RNACompleteRBPpfmHS)=motifIdProteinMapHS$Protein_name
	return(RNACompleteRBPpfmHS)
}
getRBPpfm=function(directoryName){
	
	library(stringr)
        filePath=directoryName
        file_list=list.files(filePath)
         #print(length(file_list))
        for(file in file_list){
                if (!exists("RBPPFM")){
		
                rowPFM=read.table(paste(filePath,file, sep=""), header=FALSE, sep="", stringsAsFactors=FALSE)
                RBPPFM=rowPFM
                rownames(RBPPFM)=c("A","C","G","U")
                RBPPFM=list(RBPPFM)
		rm(rowPFM)
                }
                else if(exists("RBPPFM")){

                row_temp_PFM=read.table(paste(filePath,file, sep=""), header=FALSE, sep="", stringsAsFactors=FALSE)
                temp_PFM=row_temp_PFM
                rownames(temp_PFM)=c("A","C","G","U")
                temp_PFM=list(temp_PFM)
                RBPPFM=c(RBPPFM,temp_PFM)
                rm(temp_PFM)
                rm(row_temp_PFM)
                }else{
                        print("There might be some error")
                }
        }
	
	expId=word(file_list,1, sep=fixed('_'))
	exprIdToProteinIdTable=read.table("../OfficeFiles/Tools/matrices_human/RBPDB_v1.3.1_human_2012-11-21_TDT/RBPDB_v1.3.1_protExp_human_2012-11-21.tdt",header=FALSE, sep="\t")
	proteinAnnotationTable=read.delim("../OfficeFiles/Tools/matrices_human/RBPDB_v1.3.1_human_2012-11-21_TDT/RBPDB_v1.3.1_proteins_human_2012-11-21.tdt",header=FALSE, sep="\t")
	
	proteinIds=exprIdToProteinIdTable[which(exprIdToProteinIdTable$V2 %in% expId),1]
	proteinInfo=proteinAnnotationTable[match(proteinIds,proteinAnnotationTable$V1),]
	names(RBPPFM)=proteinInfo$V5
	return(RBPPFM)
}
getPWM=function(PFMs, backgroundModel){

        PWM=list()
        for(PFM in 1:length(PFMs)){
                #avaoiding a zero freqency by adding the psedocount
                correctedFrequency=t(apply(PFMs[[PFM]], 1, function(x) x+backgroundModel))
                correctedProbability=correctedFrequency/(ncol(PFMs[[PFM]])+ sum(backgroundModel))
                w=log2(correctedProbability/backgroundModel)
                PWM=c(PWM,list(w))
        }

        names(PWM)=names(PFMs)
        return(PWM)
}

getSinglePWM=function(singlePFM,singleBackground){

        correctedFrequency=singlePFM+singleBackground
        correctedProbability=correctedFrequency/(ncol(singlePFM)+ sum(singleBackground))
        w=log2(correctedProbability/singleBackground)
        return(w)
        }
getSingleRBPPWM=function(singlePFM, singleBackground){
	
 	if(max(singlePFM)<=1 & min(singlePFM)==0){
		
		correctedProbability=singlePFM+.0000001
		w=log2(correctedProbability/singleBackground)
        	return(w)

	}else if(max(singlePFM)>1 & min(singlePFM)==0){
		correctedFrequency=singlePFM+singleBackground
		correctedProbability=correctedFrequency/(ncol(singlePFM)+ sum(singleBackground))
		w=log2(correctedProbability/singleBackground)
		return(w)
	}else if(max(singlePFM)<=1 & min(singlePFM)>0){

		w=log2(singlePFM/singleBackground)
        	return(w)	
	}else{
		return(c("ERROR"))
	}
}

getWindow=function(seqToVector,windowLength){

        #seqToVector=strsplit(seq,"")[[1]]
        window=split(seqToVector, ceiling(seq_along(seqToVector)/windowLength))
        return(window)

}


getSlidingWindow=function(seq, windowLength){

        seqToVector=strsplit(seq,"")[[1]]
        slidingWindow=list()
        slidingWindow[[1]]=getWindow(seqToVector, windowLength)
        for(i in 2:windowLength){

                slidingWindow[[i]]=getWindow(seqToVector[-(seq(1,i-1,1))],windowLength)
        }
        return(slidingWindow)
}


getSinglePWMScore=function(seq,PWM){


        if(length(seq)==ncol(PWM)){

                posA=which(seq=="A")
                posC=which(seq=="C")
                posG=which(seq=="G")
                posT=which(seq=="T")
                scoreA=if(length(posA)>0) sum(PWM["A",posA]) else 0
                scoreC=if(length(posC)>0) sum(PWM["C",posC]) else 0
                scoreG=if(length(posG)>0) sum(PWM["G",posG]) else 0
                scoreT=if(length(posT)>0) sum(PWM["T",posT]) else 0
                totalScore=sum(scoreA,scoreC,scoreG,scoreT)
                #names(totalScore)=names(PWM)
                return(sum(scoreA,scoreC,scoreG,scoreT))
        }else return("NA")
}

getSingleRNAPWMScore=function(seq,PWM){


        if(length(seq)==ncol(PWM)){

                posA=which(seq=="A")
                posC=which(seq=="C")
                posG=which(seq=="G")
                posU=which(seq=="U")
                scoreA=if(length(posA)>0) sum(PWM["A",posA]) else 0
                scoreC=if(length(posC)>0) sum(PWM["C",posC]) else 0
                scoreG=if(length(posG)>0) sum(PWM["G",posG]) else 0
                scoreU=if(length(posU)>0) sum(PWM["U",posU]) else 0
                totalScore=sum(scoreA,scoreC,scoreG,scoreU)
                #names(totalScore)=names(PWM)
                return(sum(scoreA,scoreC,scoreG,scoreU))
        }else return("NA")
}


getAllGeneAllTFPWM<- function (seqVector,totalPWM) {
        cl<-makeCluster(23)
        registerDoParallel(cl)

        test=foreach(seq =seqVector, .combine =c) %:%

                foreach(TFwindowPWM=totalPWM, .combine=c,.export = c("getWindow","totalPWM","singlePWMScore")) %dopar% {

                                list(sapply(getWindow(unlist(seq),ncol(unlist(TFwindowPWM))),function(x) singlePWMScore(x,TFwindowPWM)))
                        }

        stopCluster(cl)
        return(test)
        }

getSingleGeneAllTFPWM=function(promoterSeq,totalPFM,geneSeq){

        cl<-makeCluster(10)
        registerDoParallel(cl)

        scor=foreach(TFPFM=totalPFM, .combine=c,.export = c("getWindow","getPWM","getSinglePWM","singlePWMScore","getSingleBackground")) %dopar% {
                                #print(TFPFM)
                                list(sapply(getWindow(unlist(promoterSeq),ncol(TFPFM)),function(x) singlePWMScore(unlist(x),getSinglePWM(TFPFM[[1]],getSingleBackground(geneSeq)))))
                                #getPWM(TFPFM,getSingleBackground(geneSeq))
                                #singlePWMScore(unlist(x),getSinglePWM(TFPFM,getSingleBackground(geneSeq)))
                                #list(sapply(getWindow(unlist(promoterSeq),ncol(TFPFM)),function(x) TFPFM))
                        #rownames(TFPFM[[1]])
                        #getSingleBackground(geneSeq)
                        #noquote(names(promoterSeq))
                }
        stopCluster(cl)
        rm(cl)
        return(scor)
}


PgetSingleGeneAllTFPWM=function(promoterSeq,totalPFM,geneSeq){

       # cl<-makeCluster(23)
        #registerDoParallel(cl)
        #paste("scanning for",names(promoterSeq))

        Map=vector("list", length=length(totalPFM))
        for(i in 1:length(totalPFM)){

                        Map[[i]]=getSlidingWindow(unlist(promoterSeq),ncol(totalPFM[[i]]))

                }

        #Reduce=mccollect(Map, wait=TRUE)
        myScore=list()
        for(i in 1:length(Map)){
		cat("going for", i, "th TF scanning\n")
                slidStartingFrom=list()
                for(j in 1:length(Map[[i]])){
			cat("going for",j,"th Orientation\n")
                        slidStartingFrom[[j]]=sapply(Map[[i]][[j]], function(x) getSinglePWMScore(unlist(x),getSinglePWM(totalPFM[[i]],getSingleBackground(geneSeq))))
                }
                myScore[[i]]=slidStartingFrom
        }
         #parallel:::mckill(parallel:::children(Map), tools::SIGTERM)
        #cat("TF scan killed\n")
        return(myScore)
}


PgetSingleRNAAllProteinsPWM=function(RNASeq,totalPFM){

       # cl<-makeCluster(23)
        #registerDoParallel(cl)
        #paste("scanning for",names(promoterSeq))

        Map=vector("list", length=length(totalPFM))
        for(i in 1:length(totalPFM)){

                        Map[[i]]=getSlidingWindow(unlist(RNASeq),ncol(totalPFM[[i]]))

                }

        #Reduce=mccollect(Map, wait=TRUE)
        myScore=list()
        for(i in 1:length(totalPFM)){
                cat("going for", i, "th RBP scanning\n")
                slidStartingFrom=list()
                for(j in 1:length(Map[[i]])){
                        cat("going for",j,"th Orientation\n")
                        slidStartingFrom[[j]]=sapply(Map[[i]][[j]], function(x) getSingleRNAPWMScore(unlist(x),getSinglePWM(totalPFM[[i]],getSingleRNABackground(unlist(RNASeq)))))
                }
                myScore[[i]]=slidStartingFrom
        }
         #parallel:::mckill(parallel:::children(Map), tools::SIGTERM)
        #cat("TF scan killed\n")
        return(myScore)
}




getACGTSeq=function(numberedVector){

	sapply(numberedVector,function(x){
		if(x==1)return("A")
		else if(x==2) return("C")
		else if(x==3) return("G")
		else if(x==4) return("T")
		else return("error") 
		})

}
getACGUSeq=function(nuberedVector){

	sapply(nuberedVector,function(x){
                if(x==1)return("A")
                else if(x==2) return("C")
                else if(x==3) return("G")
                else if(x==4) return("U")
                else return("error")
                })

}
singleBestMatchScore=function(Seq, singlePFM){
	
	PFM1=singlePFM
	PFM1[PFM1==0]<-NA
	bestMatch=apply(PFM1,2,which.min)
	bestMatchSeq=getACGTSeq(bestMatch)
	backgroundModel=getSingleBackground(Seq)
	bestMatchPWM=getSinglePWM(singlePFM,backgroundModel)
	bestMatchScore=getSinglePWMScore(bestMatchSeq,bestMatchPWM)
	return(bestMatchScore)
}
getSingleRNABestMuchScore=function(RNASeq, singlePFM){

	PFM1=singlePFM
	if(min(PFM1)==0){
        	PFM1[PFM1==0]=NA
	}else if(min(PFM1)>0){
		PFM1[PFM1==min(PFM1)]=NA
	}
        bestMatch=apply(PFM1,2,which.min)
        bestMatchSeq=getACGUSeq(bestMatch)
        backgroundModel=getSingleRNABackground(RNASeq)
        bestMatchPWM=getSinglePWM(singlePFM,backgroundModel)
        bestMatchScore=getSingleRNAPWMScore(bestMatchSeq,bestMatchPWM)
        return(bestMatchScore)
}
getLincRNABounderyExones=function(lincRNAInfo){

	baounderiesExon=lincRNAInfo[,paste0("Boudaries.Exon",1:10)]
	#rownames(baounderiesExon)=lincRNAInfo[,2]
	baounderiesExonList=list()
	for(i in 1:nrow(baounderiesExon)){

		baounderiesExonList[[i]]=baounderiesExon[i,!is.na(baounderiesExon[i,])]
		}
	names(baounderiesExonList)=lincRNAInfo[,2]
	
	#print(baounderiesExonList)
		
	exonCoordinate=list()
	nullExonecoordinate=c()
	#noExon=list()
	 #exon_start=c()
         #exon_end=c()
	for(i in 1:length(baounderiesExonList)){
		
		exon_start=c()
	        exon_end=c()


		if(length(baounderiesExonList[[i]])>0){
			for(j in 1:length(baounderiesExonList[[i]])){

				exon_start=c(exon_start,na.omit(unlist(strsplit(baounderiesExonList[[i]][[j]],"[[:punct:][:space:]]+"))[1]))
				exon_end=c(exon_end,na.omit(unlist(strsplit(baounderiesExonList[[i]][[j]],"[[:punct:][:space:]]+"))[2]))
			}
		
			exonCoordinate[[i]]=list(start=exon_start, end=exon_end)
		}else{
			nullExonecoordinate=c(nullExonecoordinate,names(baounderiesExonList[i]))
			cat(names(baounderiesExonList[i]),"does not have exons\n")
			exonCoordinate[[i]]=list(start=NULL, end=NULL)
		}

	}
	names(exonCoordinate)=lincRNAInfo[,2]
	#a=sapply(baounderiesExonList, function(x) sapply(strsplit(unlist(x), "-"), "["))
	exonCoordinate=lapply(exonCoordinate, function(x) x[!is.na(x)])
	return(list(exonCoordinate,nullExonecoordinate))
}
getLincRNAChromosome=function(lincRNAName){

	return(as.character(unlist(strsplit(lincRNAName, "-"))[2]))
}
getLincRNASeq=function(lincRNAExon,strnd){
	
#  	if("DASiRocLite" %in% rownames(installed.packages()) == FALSE){
#        	source("http://bioconductor.org/biocLite.R")
#        	biocLite("DASiRocLite")
#        }
#	if("GenomicRanges" %in% rownames(installed.packages()) == FALSE){
#                source("http://bioconductor.org/biocLite.R")
#                biocLite("GenomicRanges")
#       }
#	library(DASiR)
#	library(GenomicRanges)
	
	setDasServer("http://vega.sanger.ac.uk/das")
	source = "Homo_sapiens.VEGA51.reference"
	ranges=list()
	i=1
	for(linc in lincRNAExon){

        	index=which(as.numeric(linc$start)>as.numeric(linc$end))
        	temp=linc$start[index]
        	linc$start[index]=linc$end[index]
        	linc$end[index]=temp
        	#strand=rep("+",length(linc$start))
        	#strand[index]="-"
        	cat(names(lincRNAExon[i]),"\n")
        	ranges[i]=GRanges(getLincRNAChromosome(names(lincRNAExon[i])), IRanges(start=as.numeric(linc$start), end=as.numeric(linc$end)), strand=strnd[i])
        	i=i+1
        	#cat(as.numeric(linc$start), as.numeric(linc$end), "\n")
	}
	sequences =sapply(ranges, function(x) toString(RNAString(DNAString(paste(getDasSequence(source,x),collapse = '')))))
	RNASeq=lapply(sequences, function(x) strsplit(x,""))
	names(RNASeq)=names(lincRNAExon)
	lincRNAseq=sapply(RNASeq, function(x) if(length(unlist(x))>200) return(x))
	shortNcRNA=c()
	lincRNAseqFiltered=c()
	for(i in 1:length(lincRNAseq)){
        	if(length(unlist(lincRNAseq[i][[1]])>0)){
                	lincRNAseqFiltered=c(lincRNAseqFiltered,lincRNAseq[i])
                	#shortNcRNA=c(shortNcRNA,names(lincRNAseq[i]))
        	}else{
                	shortNcRNA=c(shortNcRNA,names(lincRNAseq[i]))
                	#lincRNAseqFiltered=c(lincRNAseqFiltered,lincRNAseq[i])
        	}
	}	
	shortNcRNAVect=RNASeq[shortNcRNA]
	shortNcRNASeq=lapply(shortNcRNAVect,function(x) paste(unlist(x),collapse=""))
	lincRNAseqFilteredToFasta=lapply(lincRNAseqFiltered,function(x) paste(unlist(x),collapse=""))
	library(seqinr)
	#write.fasta(sequences=shortNcRNASeq, names=names(shortNcRNASeq), file.out="sncRNA.fasta", open = "w", nbchar = 100, as.string = TRUE)
	write.fasta(sequences=lincRNAseqFilteredToFasta, names=names(lincRNAseqFilteredSeq), file.out="lincRNA.fasta", open = "w", nbchar = 100, as.string = TRUE)
	
	return(lincRNAseqFilteredToFasta)
}
getAllLincRNARegionsFromHG19=function(lincRNARegionCoordinates){

	totalLincRNARegionsSeq=list()
	for(i in 1:nrow(lincRNARegionCoordinates)){

		totalLincRNARegionsSeq[[i]]=getSingleLincRNARigionFromHg19(lincRNARegionCoordinates[i,2],lincRNARegionCoordinates[i,3],lincRNARegionCoordinates[i,4])	
		
	}
	return(totalLincRNARegionsSeq)
}
# serial computing
RNAProteinInteractionScanner=function(RNASeq,PFM1){
	
	startTime=Sys.time()
	#totalRNAScore=list()
        totalRNAScore=mclapply(RNASeq, function(x) PgetSingleRNAAllProteinsPWM(x,PFM1),mc.cores = getOption("mc.cores", 17L))
	#totalRNAScore2=mclapply(RNASeq[19:24], function(x) PgetSingleRNAAllProteinsPWM(x,PFM1),mc.cores = getOption("mc.cores", 6L))
        #totalRNAScore=c(totalRNAScore1,totalRNAScore2)
	endTime=Sys.time()
	processTime=difftime(endTime,startTime)
	
	#names(totalGeneScore)=c("forward", "reverse")
	#names(totalGeneScore[[1]])=names(promoterSeqNotDupfiltered)
	
	names(totalRNAScore)=names(RNASeq)

	for(i in 1:length(totalRNAScore)){
        	#names(totalRNAScore[[i]])=names(PFM1)
		
        	names(totalRNAScore[[i]])=names(PFM1)
	}
	cat(names(PFM1),"\n")
	cat("processing time=",processTime,"\n")
	return(totalRNAScore)
}

getRNABestMatch=function(totalRNAScore,RNASeq,PFM1){
	
	startTime=Sys.time()
		#bestRBPMatch=list()
        	RNAAndRBP=list()
        	for(j in 1:length(totalRNAScore)){
                	bindingSite=list()
                	for(i in 1:length(totalRNAScore[[j]])){
                        	bindingPosition=list()
                        	for(z in 1:length(totalRNAScore[[j]][[i]])){                               
				 	bindingPosition[[z]]=totalRNAScore[[j]][[i]][[z]][which(na.omit(totalRNAScore[[j]][[i]][[z]])>getSingleRNABestMuchScore(RNASeq[[names(totalRNAScore[j])]],PFM1[[names(totalRNAScore[[j]][i])]]))]
                        }
                        	bindingSite[[i]]=bindingPosition
                	}
                	names(bindingSite)=names(totalRNAScore[[j]])
                	RNAAndRBP[[j]]=bindingSite
        	}
        	names(RNAAndRBP)=names(totalRNAScore)
        	#bestRBPMatch=RNAAndRBP
	
	#names(bestMatch)=names(totalScore)
	endTime=Sys.time()
	processTime=difftime(endTime,startTime)
	cat("Processing Times:", processTime)
	return(RNAAndRBP)
}

DNAProteinInteractionScanner=function(bothStrandPromoterSeq,PFMs){
	startTime=Sys.time()
	totalGeneScore=list()
	for(i in 1:length(bothStrandPromoterSeq)){
        	if(i==1){
                	cat("going for forward strand scanning")
                	}
        	if(i==2){
                	cat("going for reverse strand scanning")
                }
        totalGeneScore[[i]]=mclapply(names(lapply(bothStrandPromoterSeq[[i]],names)), function(x) PgetSingleGeneAllTFPWM(bothStrandPromoterSeq[[i]][x],PFMs,bothStrandPromoterSeq[[i]][x]),mc.cores = getOption("mc.cores", 10L))
        }

	names(totalGeneScore)=c("forward", "reverse")
	#names(totalGeneScore[[1]])=names(promoterSeqNotDupfiltered)
	#names(totalGeneScore[[2]])=names(promoterSeqNotDupfiltered)
	for(i in 1:length(totalGeneScore[[1]])){
        	names(totalGeneScore[[1]][[i]])=names(PFMs)
        	names(totalGeneScore[[2]][[i]])=names(PFMs)
	}
	
	endTime=Sys.time()
        processTime=difftime(endTime,startTime)
	cat("Total time for scanning is", processTime,"\n")
	
	return(totalGeneScore)
}

DNAProteinBestMatch=function(totalGeneScore){
	bestMatch=list()
	for(k in 1:length(totalGeneScore)){
        	geneAndTF=list()
        	for(j in 1:length(totalGeneScore[[k]])){
                	bindingSite=list()
                	for(i in 1:length(totalGeneScore[[k]][[j]])){
                        	bindingPosition=list()
                        	for(z in 1:length(totalGeneScore[[k]][[j]][[i]])){
                                 	bindingPosition[[z]]=totalGeneScore[[k]][[j]][[i]][[z]][which(totalGeneScore[[k]][[j]][[i]][[z]]>singleBestMatchScore(promoterSeqNotDupfiltered[names(totalGeneScore[[k]][j])],PFM[[names(totalGeneScore[[k]][[j]][i])]]))]
                        }
                        bindingSite[[i]]=bindingPosition
                }
                names(bindingSite)=names(totalGeneScore[[k]][[j]])
                geneAndTF[[j]]=bindingSite
        }
        names(geneAndTF)=names(totalGeneScore[[k]])
        bestMatch[[k]]=geneAndTF
	}
        names(bestMatch)=names(totalGeneScore)
        #endTime=Sys.time()
        #processTime=difftime(endTime,startTime)
	#cat(processTime)
        return(bestMatch)
}
getSinglePvalue=function(genePWMScoresAllWindows, controlPWMScoreAllWindows){

	return(wilcox.test(as.numeric(genePWMScoresAllWindows),as.numeric(controlPWMScoreAllWindows))$p.value)
}

getTotalProteinPvalueForSingleGene=function(totalProteinPWMscore,totalControlPWMScore){

	totalPvalues=c()
	for(i in 1:length(totalProteinPWMscore)){
			
		totalPvalues[i]=getSinglePvalue(unlist(totalProteinPWMscore[[i]]),unlist(totalControlPWMScore[[i]]))
		
	}
	return(totalPvalues)
}
getTotalProteinPvalueForSingleRNA=function(totalProteinPWMscore,totalControlPWMScore){
	
	totalPvalues=c()
	for(i in 1:length(totalProteinPWMscore)){
		
		totalPvalues[i]=getSinglePvalue(unlist(totalProteinPWMscore[[i]]),unlist(totalControlPWMScore[[i]]))
			
	}

	return(totalPvalues)
}
getTotalGenePvalues=function(allGenesPWMscore, totalControl){

	forwardStrandPvalues=mclapply(allGenesPWMscore[[1]], function(x) getTotalProteinPvalueForSingleGene(x,totalControl[[1]][[1]]), mc.cores = getOption("mc.cores", 17L))	
	reverseStrandPvalues=mclapply(allGenesPWMscore[[2]], function(x) getTotalProteinPvalueForSingleGene(x, totalControl[[2]][[1]]),mc.cores = getOption("mc.cores", 17L))
	return(list(forwardStrandPvalues,reverseStrandPvalues))
}
getTotalRNAPvalues=function(allRNAPWMscore,controlPWMscore){
	
	#Pvalues=mclapply(allRNAPWMscore, function(x) getTotalProteinPvalueForSingleRNA(x,controlPWMscore[paste("control",names(x))]), mc.cores = getOption("mc.cores", 17L))
	pvalues=c()
	for(i in 1:length(allRNAPWMscore)){

		pvalues[i]=list(getTotalProteinPvalueForSingleRNA(allRNAPWMscore[[i]],controlPWMscore[[i]]))
	}	
	return(pvalues)
}
getSignificantDNABindingProtein=function(TotalPvlues,Gene,scanningDirection, thresholdValue){
	if(scanningDirection=="both"){
		
		df=data.frame(p_value=pValues[["Forward strand"]][[Gene]][intersect(names(which(pValues[["Forward strand"]][[Gene]]<thresholdValue)),names(which(pValues[["Reverse strand"]][[Gene]]<thresholdValue)))])
        }else{
			df=data.frame(p_value=TotalPvlues[[scanningDirection]][[Gene]][which(TotalPvlues[[scanningDirection]][[Gene]]<thresholdValue)])
	}
	return(df)
}
getSignificantRNABindingProtein=function(totalPvalues, threshold=.01){

	df=data.frame(RBP=names(totalPvalues[which(totalPvalues<threshold)]) , P_value=totalPvalues[which(totalPvalues<threshold)])

	return(df)
}

RBPBindingBothLincRNAAndDNA=function(RBPsInDNA,RNPsINlIncRNA){


DNARNABindingTotal=list()
for(i in 1:length(RNPsINlIncRNA)){


        DNARNABinding=list()
        lincRNANames=c()
        GeneName=c()
        RBPs=c()
        
        for(j in 1:length(RBPsInDNA)){
                if(length(intersect(RNPsINlIncRNA[[i]]$RBP, RBPsInDNA[[j]]))>0)
                {
                        cat("lincRNA Name=",names(RNPsINlIncRNA[i]),"\n DNA=",names(RBPsInDNA[j]),"\n RBP=",intersect(RNPsINlIncRNA[[i]]$RBP, RBPsInDNA[[j]]), "\n")
                        lincRNANames=c(lincRNANames,names(RNPsINlIncRNA[i]))
                        GeneName=c(GeneName,names(RBPsInDNA[j]))
                        RBPs=c(RBPs, list(intersect(RNPsINlIncRNA[[i]]$RBP, RBPsInDNA[[j]])))
                }
        }
        DNARNABindingTotalList[[i]]=c(list(GeneName), list(RBPs))
        names(DNARNABindingTotalList[[i]])=c("GeneNames","RBProteins")
}
names(DNARNABindingTotalList)=names(RNPsINlIncRNA)

return(DNARNABindingTotalList)

}

writeInteractionsToxslxFiles=function(InteractionDataFrame, fileName){
for(i in 1:length(InteractionDataFrame)){

        if(i==1){                write.xlsx(InteractionDataFrame[i], file=fileName, sheetName= names(InteractionDataFrame[i]),col.names=TRUE, row.names=TRUE, append=FALSE)
                }
        else{   
                        write.xlsx(InteractionDataFrame[i], file=fileName, sheetName= names(InteractionDataFrame[i]),col.names=TRUE, row.names=TRUE, append=TRUE)              
                }
        }

}

#packages
library(Biostrings) 
library(parallel)

if("biomaRt" %in% rownames(installed.packages()) == FALSE){
source("http://bioconductor.org/biocLite.R")
biocLite("biomaRt")
}
if("DASiRocLite" %in% rownames(installed.packages()) == FALSE){
      source("http://bioconductor.org/biocLite.R")
      biocLite("DASiRocLite")
        }
if("GenomicRanges" %in% rownames(installed.packages()) == FALSE){
      source("http://bioconductor.org/biocLite.R")
      biocLite("GenomicRanges")
        }
if("mirbase.db" %in% rownames(installed.packages()) == FALSE){
        source("http://bioconductor.org/biocLite.R")
        biocLite("mirbase.db")
}

#library(Biostrings)
#library(parallel)
library("biomaRt")
library(DASiR)
library(GenomicRanges)
library(mirbase.db)
#source("DEAnalyzer.r")

#get promoter and gene sequences

ensembl=useMart("ensembl")
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
filters = listFilters(ensembl)

top61PromoterandGeneSeq=getPeomoterAndGeneSeq(top61Genes,ensembl)
top61PromoterSeq=top61PromoterandGeneSeq[[1]]
promoterSeq=top61PromoterSeq[,1]
names(promoterSeq)=top61PromoterSeq[,2]
promoterSeqNotDup=promoterSeq[!duplicated(names(promoterSeq))]
promoterSeqNotDupfiltered=promoterSeqNotDup[promoterSeqNotDup!="Sequence unavailable"]
top61GeneSeq=top61PromoterandGeneSeq[[2]]
geneSeq=top61GeneSeq[,1]
names(geneSeq)=top61GeneSeq[,2]
geneSeqNotDup=geneSeq[!duplicated(names(geneSeq))]
unidentifiedGeneSeq=top61PromoterandGeneSeq[[3]]
promoterSeqNotDupOppositeStrand=sapply(promoterSeqNotDupfiltered, function(x) toString(complement(DNAString(x))))
bothStrandPromoterSeq=list(promoterSeqNotDupfiltered,promoterSeqNotDupOppositeStrand)


#get the LncRNA sequence
lincCoordinate=gbm53lincRNAInfo[,2]
filterBasedOnTranscriptionDirection=gbm53lincRNAInfo[which(gbm53lincRNAInfo$"Transcription.direction"!="NA"),]
filterBasedOnTranscriptionDirection=filterBasedOnTranscriptionDirection[filterBasedOnTranscriptionDirection$"Transcription.direction"!="?",]
filterBasedOnTranscriptionDirection=filterBasedOnTranscriptionDirection[filterBasedOnTranscriptionDirection$"Transcription.direction"!="+ (?)",]
filterBasedOnTranscriptionDirection=filterBasedOnTranscriptionDirection[filterBasedOnTranscriptionDirection$"Transcription.direction"!="-?",]
filterBasedOnTranscriptionDirection=filterBasedOnTranscriptionDirection[filterBasedOnTranscriptionDirection$"Gene"!="Gene",]
lincRNAExons=getLincRNABounderyExones(filterBasedOnTranscriptionDirection)
filteringIndex=match(lincRNAExons[[2]],names(lincRNAExons[[1]]))
lincRNAExon=lincRNAExons[[1]][-filteringIndex]
lincRNAExon$"TCGA_gbm-6-99298751"$start=lincRNAExon$"TCGA_gbm-6-99298751"$start[-3]
lincRNAExon$"TCGA_gbm-6-23626001"$start=lincRNAExon$"TCGA_gbm-6-23626001"$start[-8]
lincRNAExon$"TCGA_gbm-6-23626001"$end=lincRNAExon$"TCGA_gbm-6-23626001"$end[-8]
index=match("TCGA_gbm-7-54802501",names(lincRNAExon))
lincRNAExon=lincRNAExon[-index]
strand=gbm53lincRNAInfo[match(names(lincRNAExon), gbm53lincRNAInfo$Gene),"Transcription.direction"]
strand[17]="+"
lincRNASequences=getLincRNASeq(lincRNAExon,strand)
lincRNASequencesDE=lincRNASequences[intersect(intersect(rownames(top61TranscriptExpressionMatrix),gbm53lincRNAInfo$Gene),names(lincRNASequences))]

#quering the lincRNA sequence from HG19 genome
#read genomic coordinate and score value from gene identfying algorithm
#bedtools getfasta -fi ../../../data/csb/organisms/homo_sapiens/hg19.fa -bed exon.bed -s -name -fo LincRNASeq.fa.out
#bedtools getfasta -fi ../../../data/csb/organisms/homo_sapiens/hg19.fa -bed LincRNARegion.bed -s -name -fo LincRNARegionSeq.fa.out
library("seqinr")
lincRNARegionSeq=read.fasta("LincRNARegionSeq.fa.out",as.string=1)
lincRNARegionSeq=unlist(lincRNARegionSeq)
LincRNASeq=read.fasta("LincRNASeq.fa.out", as.string=1)
LincRNASeq=unlist(LincRNASeq)

LincRNASeqDNAStringObj=sapply(LincRNASeq, function(x) DNAString(x))
LincRNASeqRNAStringObj=sapply(LincRNASeqDNAStringObj, function(x) RNAString(complement(x)))
# add Poly A tail based on the janne's analysis result
lincRNAAdditionalInfo=read.table("significantLincRNAseq.txt",header=TRUE, sep="\t",stringsAsFactors = FALSE)
names(LincRNASeqDNAStringObj)=lincRNAAdditionalInfo[,1]
names(LincRNASeqRNAStringObj)=lincRNAAdditionalInfo[,1]

namesOfLincRNAWithPoliAtail=lincRNAAdditionalInfo[lincRNAAdditionalInfo$PolyA_site>0,1]
LincRNASeqRNAStringObj[namesOfLincRNAWithPoliAtail]=sapply(LincRNASeqRNAStringObj[namesOfLincRNAWithPoliAtail], function(x) RNAString(paste0(toString(x),"AAAAAAAAAAAAAAA")))
LincRNASeqHg19=sapply(LincRNASeqRNAStringObj,function(x) toString(x))
write.fasta(as.list(LincRNASeqHg19), names(LincRNASeqHg19), nbchar = 300, "LincRNASeqHg19.fasta", open = "w")
#add the ploy A tail accourding to the result of gene identiying algorithum
#merging exones of "TCGA_gbm-2-104096501" 
#LincRNASeq[["TCGA_gbm-2-104096501Ex1"]]=c(LincRNASeq[["TCGA_gbm-2-104096501Ex1"]],LincRNASeq[["TCGA_gbm-2-104096501Ex2"]],LincRNASeq[["TCGA_gbm-2-104096501Ex3"]])
#LincRNASeq=LincRNASeq[-c(4,5)]
#eddit the fasta file for future use
#write.fasta(sequences=LincRNASeq,names=names(LincRNASeq),file.out="LincRNASeq.fa.out")

#Now we got both the region in which the gene identifing algorithm calculated  the coding potential and the LincRNA sequences from hg19 genome assembly.
#Taking the whole region sequences in search of UTR regions

#chr=sapply(names(lincRNASequencesDE), function(x) getLincRNAChromosome(x))
#coordinate=gbm53lincRNAInfo[match(names(lincRNASequencesDE),gbm53lincRNAInfo$Gene),c(18,19,23)]
#names(coordinate)=names(lincRNASequencesDE)

#setDasServer("http://vega.sanger.ac.uk/das")
#source = "Homo_sapiens.VEGA51.reference"

#totalRanges=GRanges(chr, IRanges(start=as.numeric(coordinate[,1]), end=as.numeric(coordinate[,2])), strand=Rle(strand(coordinate[,3]), rep(1,17)))
#bulkLincRNARefDNA=sequences =sapply(totalRanges, function(x) toString(RNAString(DNAString(paste(getDasSequence(source,x),collapse = '')))))


#get the miRNA sequences
totalmiRNASeqs=mirbaseSEQUENCE
miRNAsequences=as.list(totalmiRNASeqs[tolower(names(as.list(totalmiRNASeqs))) %in% tolower(names(top61SignificantmiRNA))])
foundIndex=match(names(miRNAsequences),tolower(names(top61SignificantmiRNA)))
write.fasta(miRNAsequences,names(miRNAsequences),file.out="miRNASeqmiRbaseDB.fasta", open="w",nbchar =60)
#sequences not found in the miRbase database
notfoundmiRNA=tolower(names(top61SignificantmiRNA))[-foundIndex]

relatedNotFoundmiRNA=lapply(notfoundmiRNA[-3], function(x) as.list(totalmiRNASeqs[grep(paste0(paste(unlist(strsplit(tolower(x),"-"))[1:3],collapse="-"),"+[:punct:]*[:alpha:]*"),tolower(names(as.list(totalmiRNASeqs))))]))




#Get PFM fro the DNA binding proteins
	
#PFM=getPFMhPDI("hPDI/all_pwm_hPDI/")

# run DNA protein intercation scanner in the promoter region of the DE genes
#be carefull this will take about 2 and half days using 10 core parallaly
DNAProteinInteractionMotifs=DNAProteinInteractionScanner(bothStrandPromoterSeq,PFM)
bestMatches=DNAProteinBestMatch(DNAProteinInteractionMotifs)

#Check the significance of the PWM score by taking random intergenic sequence of length 2000 bp

#lets take an intergenic genomic region: chromosome 4, start coordinate= 54985200, end coordinate= 54987200 as the control and calculate the significance of each window's PWM score
setDasServer("http://vega.sanger.ac.uk/das")
source = "Homo_sapiens.VEGA51.reference"
range=GRanges(4, IRanges(start=54985200, end=54987200), strand="+")
controlSeq=toString(DNAString(paste(getDasSequence(source,range),collapse = '')))
names(controlSeq)="intergenicControlPosetiveStrand"
controlSeqOppositStrand=toString(complement(DNAString(controlSeq)))
names(controlSeqOppositStrand)="intergenicControlNegetiveStrand"
control=list(controlSeq,controlSeqOppositStrand)

#now lets scan the control sequence looking for protein binding motif
controlDNAProteinInteractionMotifs=DNAProteinInteractionScanner(control,PFM)

# Now lets calculate the p-value of the PWM scores with the controle one.

#since we are not sure of the distibutions of the PWM scores I just would like to use many-witny(walcoxon test)
pValues=getTotalGenePvalues(DNAProteinInteractionMotifs,controlDNAProteinInteractionMotifs)
names(pValues)=c("Forward strand","Reverse strand")
names(pValues[[1]])=names(promoterSeqNotDupfiltered)
names(pValues[[2]])=names(promoterSeqNotDupfiltered)
for(i in 1:length(pValues[[1]])){

	names(pValues[[1]][[i]])=names(PFM)
	names(pValues[[2]][[i]])=names(PFM)

}
#significant DNA-protein  interactions

forwardStrandBinding=lapply(names(pValues[[1]]), function(x) getSignificantDNABindingProtein(pValues,x,"Forward strand", .01))
reverseStrandBinding=lapply(names(pValues[[1]]), function(x) getSignificantDNABindingProtein(pValues,x,"Reverse strand", .01))
proteinsBindingInBothStrand=lapply(names(pValues[[1]]), function(x) getSignificantDNABindingProtein(pValues,x,"both", .01))

names(forwardStrandBinding)=names(promoterSeqNotDupfiltered)
names(reverseStrandBinding)=names(promoterSeqNotDupfiltered)
names(proteinsBindingInBothStrand)=names(promoterSeqNotDupfiltered)

#protein types used in the PWM scanings
library(stringi)
TFs=read.table("DNABindingProteinClass/TFs.txt", header=TRUE, sep="\t")
TFs$Protein=sapply(TFs$Protein,function(x)(stri_trim(as.character(x), "right")))
RBPs=read.table("DNABindingProteinClass/RBP.txt", header=TRUE, sep="\t")
RBPs$Protein=sapply(RBPs$Protein,function(x)(stri_trim(as.character(x), "right")))
kainese_protein=read.table("DNABindingProteinClass/kainese_proteins.txt", header=TRUE, sep="\t")
kainese_protein$Protein=sapply(kainese_protein$Protein,function(x)(stri_trim(as.character(x), "right")))
chromatin_associated_protein=read.table("DNABindingProteinClass/kromatin_associated_prorein.txt", header=TRUE, sep="\t")
chromatin_associated_protein$Protein=sapply(chromatin_associated_protein$Protein,function(x)(stri_trim(as.character(x), "right")))
mitokondorial_protein=read.table("DNABindingProteinClass/mitokondorial_protein.txt", header=TRUE, sep="\t")
mitokondorial_protein$Protein=sapply(mitokondorial_protein$Protein,function(x)(stri_trim(as.character(x), "right")))
niclicAcid_binding_protein=read.table("DNABindingProteinClass/nucleic_acide_binding_protein.txt", header=TRUE, sep="\t")
niclicAcid_binding_protein$Protein=sapply(niclicAcid_binding_protein$Protein,function(x)(stri_trim(as.character(x), "right")))

#export the results to excel sheets

forwardStrandBindingRBPs=sapply(forwardStrandBinding, function(x) intersect(RBPs$Protein,rownames(x)))
reverseStrandBindingRBPs=sapply(reverseStrandBinding, function(x) intersect(RBPs$Protein,rownames(x)))
proteinsBindingInBothStrandRBPs=sapply(proteinsBindingInBothStrand, function(x) intersect(RBPs$Protein,rownames(x)))

############################### Protein-RNA interaction ###############################################

#get PFM for  RBP 
	#PFM from RBPDB database
	RBPDBpfm=getRBPpfm("../OfficeFiles/Tools/matrices_human/PFMDir/")
	#PFM from RNAComplete database
	RNACompletepfm=getRNACpmpleteRBPpfm("../OfficeFiles/Tools/RNACompletePFM/top10align_motifs/")

#totalRBPpfm=c(RBPDBpfm,RNACompletepfm)
#run the RNA protein intercation scanner to the entaier region of the lincRNA
#RNAProteinInteractionMotifs=RNAProteinInteractionScanner(lincRNASequences,RBPDBpfm)
#RNAProteinInteractionMotifs1=RNAProteinInteractionScanner(lincRNASequences,RBPDBpfm)
#RNAProteinInteractionMotifs1DE=RNAProteinInteractionMotifs1[lincRNASequencesDE]
#RNAProteinInteractionMotifsRNAComplete=RNAProteinInteractionScanner(lincRNASequences,RNACompletepfm)
#RNAProteinInteractionMotifsRNACompleteDE=RNAProteinInteractionScanner(lincRNASequencesDE,RNACompletepfm)
#RNAProteinInteractionMotifsTotal=c(RNAProteinInteractionMotifs,RNAProteinInteractionMotifsRNAComplete)

#Running RNA-Protein PWM score calculator for lincRNA extracted from Hg19

RNAProteinInteractionMotifsRBPDBpfmHg19=RNAProteinInteractionScanner(LincRNASeqHg19,RBPDBpfm)
RNAProteinInteractionMotifsRNACompletepfmHg19=RNAProteinInteractionScanner(LincRNASeqHg19,RNACompletepfm)

#calculate the PWM score of the the control sequence from random exonic sequence  of hg19 genome.

controlExonicSeq=read.fasta("lincRNAControlSeq.fa",as.string=1)
controlExonicSeq=unlist(controlExonicSeq)
#convert the exonic DNA sequence to RNA sequence
controlExonicSeqDNAObj=DNAString(controlExonicSeq)

controlExonicSeqRNAObj=RNAString(complement(controlExonicSeqDNAObj))
controlExonicSeqString=toString(controlExonicSeqRNAObj)
#get the length of each lincRNA sequences
LincRNASeqHg19Length=sapply(LincRNASeqHg19, function(x) nchar(x))
controlExonicSeqStrings=sapply(LincRNASeqHg19Length, function(x) substr(controlExonicSeqString,1,x))
names(controlExonicSeqStrings)=paste("control",names(controlExonicSeqStrings))
#calculate the PWM score for the control RNA sequence
controlRNAProteinInteractionMotifsRBPDBpfmHg19=RNAProteinInteractionScanner(controlExonicSeqStrings,RBPDBpfm)
controlRNAProteinInteractionMotifsRNACompletepfmHg19=RNAProteinInteractionScanner(controlExonicSeqStrings,RNACompletepfm)

# calculate the significance p-values

#RNAProteinInteractionPvalues=getTotalRNAPvalues(RNAProteinInteractionMotifsRNACompletepfmHg19,controlRNAProteinInteractionMotifsRNACompletepfmHg19)

totalPvaluesForRBPDBlincRNAs=getTotalRNAPvalues(RNAProteinInteractionMotifsRBPDBpfmHg19,controlRNAProteinInteractionMotifsRBPDBpfmHg19)
totalPvaluesForRNACompletelincRNAs=getTotalRNAPvalues(RNAProteinInteractionMotifsRNACompletepfmHg19,controlRNAProteinInteractionMotifsRNACompletepfmHg19)

for(i in 1:length(totalPvaluesForRBPDBlincRNAs)){

	names(totalPvaluesForRBPDBlincRNAs[[i]])=names(RBPDBpfm) 
}
for(i in 1:length(totalPvaluesForRNACompletelincRNAs)){

	names(totalPvaluesForRNACompletelincRNAs[[i]])=names(RNACompletepfm)
}

names(totalPvaluesForRBPDBlincRNAs)=names(RNAProteinInteractionMotifsRNACompletepfmHg19)
names(totalPvaluesForRNACompletelincRNAs)=names(RNAProteinInteractionMotifsRNACompletepfmHg19)

#totalPvaluesForRBPDBlincRNAsDataFrame=lapply(totalPvaluesForRBPDBlincRNAs,function(x) data.frame(t(sapply(x,c))))
#totalPvaluesForRNACompletelincRNAsDataFrame=lapply(totalPvaluesForRNACompletelincRNAs, function(x) data.frame(t(sapply(x,c))))

# significant RNA protein interactions 
significantRNABindingProteinFromRBPDB=lapply(totalPvaluesForRBPDBlincRNAs, function(x) getSignificantRNABindingProtein(x))
significantRNABindingProteinFromRNAComplet=lapply(totalPvaluesForRNACompletelincRNAs, function(x) getSignificantRNABindingProtein(x))
#write to the xlsx file
writeInteractionsToxslxFiles(significantRNABindingProteinFromRBPDB, "significantRNABindingProteinFromRBPDB.xlsx")
writeInteractionsToxslxFiles(significantRNABindingProteinFromRNAComplet, "significantRNABindingProteinFromRNAComplet.xlsx")

#bestRNAProteinMatch=getRNABestMatch(RNAProteinInteractionMotifsTotal,lincRNASequences,totalRBPpfm)
#find minimum best match score
#bestRNAProteinMatchRBPDB=getRNABestMatch(RNAProteinInteractionMotifs1DE,lincRNASequencesDE,RBPDBpfm)
#bestRNAProteinMatchRNAComplet=getRNABestMatch(RNAProteinInteractionMotifsRNACompleteDE,lincRNASequencesDE,RNACompletepfm)


#RBPs that bind both in LincRNA and DNA
forwardStrandBindingRBPs=sapply(forwardStrandBinding, function(x) intersect(RBPs$Protein,rownames(x)))
reverseStrandBindingRBPs=sapply(reverseStrandBinding, function(x) intersect(RBPs$Protein,rownames(x)))
proteinsBindingInBothStrandRBPs=sapply(proteinsBindingInBothStrand, function(x) intersect(RBPs$Protein,rownames(x)))
#RBP from RBPDB database
RBPInBothLincRNAandForwadStrandDNA=RBPBindingBothLincRNAAndDNA(forwardStrandBindingRBPs,significantRNABindingProteinFromRBPDB)
RBPInBothLincRNAandReverseStrandDNA=RBPBindingBothLincRNAAndDNA(reverseStrandBindingRBPs,significantRNABindingProteinFromRBPDB)
RBPInBothLincRNAandBothStrandDNA=RBPBindingBothLincRNAAndDNA(proteinsBindingInBothStrandRBPs,significantRNABindingProteinFromRBPDB)
#RBP from RNAComplete Databases
RBPInBothLincRNAandForwadStrandDNARNACompleteDB=RBPBindingBothLincRNAAndDNA(forwardStrandBindingRBPs,significantRNABindingProteinFromRNAComplet)
RBPInBothLincRNAandReverseStrandDNARNACompleteDB=RBPBindingBothLincRNAAndDNA(reverseStrandBindingRBPs,significantRNABindingProteinFromRNAComplet)
RBPInBothLincRNAandBothStrandDNARNACompleteDB=RBPBindingBothLincRNAAndDNA(proteinsBindingInBothStrandRBPs,significantRNABindingProteinFromRNAComplet)


#storing interactions in a data frame RBPDB database
forwardInteractionInDataFrameRBPDB=lapply(RBPInBothLincRNAandForwadStrandDNA, function(x) data.frame(t(sapply(x,c))))
reverseInteractionInDataFrameRBPDB=lapply(RBPInBothLincRNAandReverseStrandDNA, function(x) data.frame(t(sapply(x,c))))
bothStrandInteractionInDataFrameRBPDB=lapply(RBPInBothLincRNAandBothStrandDNA, function(x) data.frame(t(sapply(x,c))))
#storing interactions in a data frame RNAComplete database
forwardInteractionInDataFrame=lapply(RBPInBothLincRNAandForwadStrandDNARNACompleteDB, function(x) data.frame(t(sapply(x,c))))
reverseInteractionInDataFrame=lapply(RBPInBothLincRNAandReverseStrandDNARNACompleteDB, function(x) data.frame(t(sapply(x,c))))
bothStrandInteractionInDataFrame=lapply(RBPInBothLincRNAandBothStrandDNARNACompleteDB, function(x) data.frame(t(sapply(x,c))))
#adjust RBPDB interactions the data frame to write to xslx file
forwardInteractionInDataFrameRBPDBToxslx=lapply(forwardInteractionInDataFrameRBPDB, function(x) data.frame(lapply(x, as.character), stringsAsFactors=FALSE))
reverseInteractionInDataFrameRBPDBToxslx=lapply(reverseInteractionInDataFrameRBPDB, function(x) data.frame(lapply(x, as.character), stringsAsFactors=FALSE))
bothStrandInteractionInDataFrameRBPDBToxslx=lapply(bothStrandInteractionInDataFrameRBPDB, function(x) data.frame(lapply(x, as.character), stringsAsFactors=FALSE))
#adjusting the RNAComplete DB interaction data frames to write to xslx file
forwardInteractionInDataFrameToxslx=lapply(forwardInteractionInDataFrame, function(x) data.frame(lapply(x, as.character), stringsAsFactors=FALSE))
reverseInteractionInDataFrameToxslx=lapply(reverseInteractionInDataFrame, function(x) data.frame(lapply(x, as.character), stringsAsFactors=FALSE))
bothStrandInteractionInDataFrameToxslx=lapply(bothStrandInteractionInDataFrame, function(x) data.frame(lapply(x, as.character), stringsAsFactors=FALSE))


install.packages("xlsx")
library("xlsx")
#write the interaction results from RBPDB to the Xslx file one sheet for each lincRNA that contains information of its RBPs and genes that same RBP interact with.
writeInteractionsToxslxFiles(forwardInteractionInDataFrameRBPDBToxslx, "forwardInteractionInDataFrameRBPDBToxslx.xlsx")

writeInteractionsToxslxFiles(reverseInteractionInDataFrameRBPDBToxslx, "reverseInteractionInDataFrameRBPDBToxslx.xlsx")

writeInteractionsToxslxFiles(bothStrandInteractionInDataFrameRBPDBToxslx, "bothStrandInteractionInDataFrameRBPDBToxslx.xlsx")
#write the interaction results from RBPDB to the Xslx file one sheet for each lincRNA that contains information of its RBPs and genes that same RBP interact with.
writeInteractionsToxslxFiles(forwardInteractionInDataFrameToxslx, "forwardInteractionInDataFrameToxslx.xlsx")

writeInteractionsToxslxFiles(reverseInteractionInDataFrameToxslx, "reverseInteractionInDataFrameToxslx.xlsx")

writeInteractionsToxslxFiles(bothStrandInteractionInDataFrameToxslx, "bothStrandInteractionInDataFrameToxslx.xlsx")

#Lets be selective on those lincRNAs that are with top 5 differently expressed and study them
Top5DELincRNAs=data.frame(FoldChange=subset(top61SignificantTranscript[names(LincRNASeqHg19)],top61SignificantTranscript[names(LincRNASeqHg19)]>6))




