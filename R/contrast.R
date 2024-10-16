#' Compute the treatment effects
#'
#' @param sele a list specifying the positions of
#' interactions terms that will be used to compute the treatment effects
#' @param coeVec the estimates of fixed effects of the linear mixed-effects model
#' @param varMatrix the variance matrix of the estimates of the fixed effects
#' @return  a list for which the first one is from model while the second one is from origin approach
contrast=function(sele,coeVec,varMatrix){
## sele is a list of length 2 each of whose components is a vector specifying which entries have to be selected
## coeVec gives the estimated coefficient with the first element representing
##  reference group and the rest elements representing differences between other groups and representing group
##varMatrix gives the covariance matrix for coeVec
  if (any(names(coeVec)!=rownames(varMatrix)) | any(names(coeVec)!= colnames(varMatrix)) ){
    stop("coeVec and varMatrix should have the same names!")
  }
  modelname=names(coeVec)
finalV=matrix(nrow=length(sele[[1]]),ncol = length(sele[[2]]))
finalP=matrix(nrow=length(sele[[1]]),ncol = length(sele[[2]]))
finalT=matrix(nrow=length(sele[[1]]),ncol = length(sele[[2]]))

rownames(finalV)=names(sele[[1]])
colnames(finalV)=names(sele[[2]])
  coev1=coeVec[sele[[1]][-1]]
  coev2=coeVec[sele[[2]][-1]]

  finalE1= coev1%*%t(rep(1,length(coev2)))+rep(1,length(coev1))%*%t(coev2)+coeVec[2]
  if (length(sele)==3){
  coev3=coeVec[sele[[3]][-1]]
  finalE1=finalE1+matrix(coev3,nrow = length(sele[[1]])-1)
  }

finalE2=cbind(coev1+coeVec[2],finalE1)
finalE3=c(coeVec[2],coeVec[2]+coev2)
finalE=rbind(finalE3,finalE2)


rownames(finalE)=names(sele[[1]])
colnames(finalE)=names(sele[[2]])
rownames(finalT)=names(sele[[1]])
colnames(finalT)=names(sele[[2]])
rownames(finalP)=names(sele[[1]])
colnames(finalP)=names(sele[[2]])

if (length(sele)==2){
  seleLong=c(2,c(sele[[1]][-1],sele[[2]][-1]))
  varm=varMatrix[seleLong,seleLong]
finalV[1,1]=varm[1,1]
for (i in 2:length(sele[[1]])) {
finalV[i,1]=t(c(1,1))%*%varm[c(1,i),c(1,i)]%*%c(1,1)
}

for (i in 2:length(sele[[2]])) {
  finalV[1,i]=t(c(1,1))%*%varm[c(1,length(sele[[1]])-1+i-1),c(1,length(sele[[1]])-1+i-1)]%*%c(1,1)
}


for (i in 2:length(sele[[1]])) {
  for (j in 2:length(sele[[2]])){

   finalV[i,j]= t(c(1,1,1))%*%varm[c(1,i,length(sele[[1]])+j-1),c(1,i, length(sele[[1]])+j-1)]%*%c(1,1,1)
  }
}
}


if (length(sele)==3){
  seleLong=c(2,c(sele[[1]][-1],sele[[2]][-1]),sele[[3]][-1])
  varm=varMatrix[seleLong,seleLong]
  finalV[1,1]=varm[1,1]
  for (i in 2:length(sele[[1]])) {
    finalV[i,1]=t(c(1,1))%*%varm[c(1,i),c(1,i)]%*%c(1,1)
  }

  for (i in 2:length(sele[[2]])) {
    finalV[1,i]=t(c(1,1))%*%varm[c(1,length(sele[[1]])+i-1),c(1,length(sele[[1]])+i-1)]%*%c(1,1)
  }


  for (j in 2:length(sele[[2]])) {
    for (i in 2:length(sele[[1]])){
      finalV[i,j]= t(c(1,1,1,1))%*%varm[c(1,i,length(sele[[1]])+j-1,length(sele[[1]])+length(sele[[2]])-1+(length(sele[[1]])-1)*(j-2)+(i-1)),
                                        c(1,i, length(sele[[1]])+j-1,length(sele[[1]])+length(sele[[2]])-1+(length(sele[[1]])-1)*(j-2)+(i-1))]%*%c(1,1,1,1)
  }
}



}


for (i in 1:nrow(finalV)) {
  for (j in 1:ncol(finalV)){
finalT[i,j]=abs(finalE[i,j])/(finalV[i,j])^0.5
    finalP[i,j]=pnorm(finalT[i,j],lower.tail = F)*2
  }
}


return(list(effect=finalE,variance=finalV,tvalue=finalT,pvalue=finalP))
  }





