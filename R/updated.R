

danovir1=function(data,lambda,covariate,...){
  fullnames=colnames(data)
  if (!all(c("response","subjectid","treatment","assay","antigen","antibody",covariate)
           %in% fullnames)){
    stop("Please re-specify the column names of data which should include variables **response**,
         subjectid,treatment,assay,antigen,antibody and all the variables include in
         covariate")
  }


  data$treatment=as.factor(data$treatment)
  data$antigen=as.factor(data$antigen)
  data$antibody=as.factor(data$antibody)
  for (i in 1:length(covariate)) {
    if (is.character(data[,covariate[i]])){
      data[,covariate[i]]=as.factor(data[,covariate[i]])
    }
  }
  attach(data)


  if (length(levels(treatment))!=2){
    stop("Treatment should have two levels!")
  }
  n1=length(levels(antigen))
  n2=length(levels(antibody))
  cn=length(covariate)
  n3=0
  for (i in 1:cn) {
    n3=ifelse(is.factor(data[,covariate[i]]),n3+length(levels(data[,covariate[i]]))-1,n3+1)
  }
  ## nn represent the number of terms before second order terms
  nn=(n1-1)+ #antigen
    (n2-1)+  # antibody
    n3+ #covariate
    1+ # treatment
    1  #intercept
  ## nnn represent the number of terms before third order terms
  nnn=1+ # intercept
    1+ # treatment
    n1-1+ n2-1+ # antigen and antibody
    1*(n1-1)+1*(n2-1)+(n1-1)*(n2-1) + # second order terms
    n3 # covariate
  data=groupedData(response~treatment|subjectid,data=data)
  if (!all(type%in% c(0,1))){stop("Please specify the model in correct way!")}


  data$feature=paste0(data$assay,data$antigen,data$detect_reagent,sep="")
  ff=as.formula(paste("response~treatment+ antigen+ antibody+ ",paste(covariate,sep="+"),"+
                    treatment:antigen+treatment:antibody+antigen:antibody+treatment:antigen:antibody"))
  mm=model.matrix(lm(ff,data=data))
   net=glmnet(x=mm,y=response,lambda = lambda)
  results=lme(fixed=ff,
              random=pdDiag(~1),weights= varIdent(form=~1|feature),
              data = data,method = "ML",
              control = control)
  sele=list(antigen=c(2,(nn+1):(nn+n1-1)),antibody=c(2,(nn+n1):(nn+n1-1+n2-1)),
            three=c(2,(nnn+1):(nnn+(n1-1)*(n2-1))))



  aa=summary(results)
  names(sele$antigen)=levels(data$antigen)
  names(sele$antibody)=levels(data$antibody)
  coeVec=aa$coefficients$fixed
  varMatrix=aa$varFix
  bb=contrast(sele=sele,coeVec = coeVec,varMatrix = varMatrix)


  antigen=rownames(bb$effect)
  antibody=colnames(bb$effect)
  basic=matrix(nrow=length(antigen),ncol=length(antibody))
  pp= matrix(nrow=length(antigen),ncol=length(antibody))
  vv= matrix(nrow=length(antigen),ncol=length(antibody))
  for (i in 1:length(antigen)) {
    for (j in 1:length(antibody)) {
      index=which(data$antigen==antigen[i] & data$antibody==antibody[j])
      dd=data[index,]
      x=dd$response[which(dd$treatment==levels(data$treatment)[1])]
      y=dd$response[which(dd$treatment==levels(data$treatment)[2])]
      basic[i,j]=t.test(y,x)$estimate[1]-t.test(y,x)$estimate[2]
      pp[i,j]=t.test(y,x)$p.value
      vv[i,j]=var(x)/length(x)+var(y)/length(y)
    }
  }

  seleFeature1=data.frame(antigen=rep(rownames(bb$effect),ncol(bb$effect)),
                          antibody=rep(colnames(bb$effect),rep(nrow(bb$effect),
                                                               ncol(bb$effect))),
                          effects=c(bb$effect),variance=c(bb$variance),pvalue=c(bb$pvalue))


  seleFeature2=data.frame(antigen=rep(rownames(bb$effect),ncol(bb$effect)),
                          antibody=rep(colnames(bb$effect),rep(nrow(bb$effect),
                                                               ncol(bb$effect))),
                          effects=c(basic),variance=c(vv),pvalue=c(pp))


  seleFeature=list(lmeresults=results,selected=seleFeature1,benchmark=seleFeature2)
  return(seleFeature)
}











