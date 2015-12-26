corolations=function(vect1,vect2){


	 if(length(vect1)==length(vect2)){
	
		corl=cor(as.numeric(vect1),as.numeric(vect2))
	}

	return(corl)
}
corAssociation=function(mat1,mat2){

	if(nrow(mat1)==nrow(mat2)){
		
		corMatrix=matrix(0, nrow=nrow(mat1), ncol=nrow(mat2))
		 for(i in 1:nrow(mat1)){

                        for(j in 1:nrow(mat2))

                                corMatrix[i,j]=corolations(mat1[i,],mat2[j,])

                }
	
	}
	rownames(corMatrix)=rownames(mat1)
	colnames(corMatrix)=rownames(mat2)
	return(corMatrix)

}
#visualize the cirrolation with igraph package
igrapher=function(geneTranscriptNetwork,origColor,targColor){
        net=as.data.frame(which(geneTranscriptNetwork==TRUE, arr.ind=TRUE))
        origen=rownames(net)
        target=colnames(geneTranscriptNetwork[,net$col])
#origen=rep(rownames(geneTranscriptNetwork), ncol(geneTranscriptNetwork))
#target=c(sapply(colnames(geneTranscriptNetwork),function(x) rep(x, nrow(geneTranscriptNetwork))))
        directed=rep("TRUE",length(origen))
        geneTranscriptNetworkGraph=data.frame(origen,target,directed)
        geneTranscriptNetworkGraph.net=graph.data.frame(geneTranscriptNetworkGraph,directed=F)
        V(geneTranscriptNetworkGraph.net)$color=origColor
        V(geneTranscriptNetworkGraph.net)[match(target,V(geneTranscriptNetworkGraph.net)$name)]$color=targColor
#V(geneTranscriptNetworkGraph.net)$size=degree(geneTranscriptNetworkGraph.net)/10

        return(geneTranscriptNetworkGraph.net)
}



