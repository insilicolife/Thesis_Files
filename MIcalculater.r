if("entropy" %in% rownames(installed.packages()) == FALSE){
        biocLite("entropy")
        }
library("entropy")
MIs=function(vect1,vect2){

	if(length(vect1)==length(vect2)){
		

		vect1=cut(rank(vect1),5,labels=F)
		vect2=cut(rank(vect2),5,labels=F)
                ttable=table(vect1,vect2)
               	pofVect1andvect2=ttable/sum(ttable)
                MIs=mi.plugin(pofVect1andvect2,unit="log2")
        }
		
	return(MIs)
}

MIAssociation=function(mat1,mat2){
	if(nrow(mat1)==nrow(mat2)){
		MImatrix=matrix(0, nrow=nrow(mat1), ncol=nrow(mat2))
		for(i in 1:nrow(mat1)){

			for(j in 1:nrow(mat2))
				
				MImatrix[i,j]=MIs(mat1[i,],mat2[j,])

		}
	}
	rownames(MImatrix)=rownames(mat1)
	colnames(MImatrix)=rownames(mat2)
	return(MImatrix)
}


