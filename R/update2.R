

danovir=function(treatment,antigen,antibody,covariate=NULL,type=1,
                  nfold=10,nlambda=100,data){

  fullnames=colnames(data)
  if (!all(c("response","subjectid","treatment","antigen","antibody","feature",covariate)
           %in% fullnames)){
    stop("Please re-specify the column names of data which should include variables **response**,
         subjectid,treatment,assay,antigen,antibody and all the variables include in
         covariate")
  }

  if (!is.null(covariate)){
  for (i in 1:length(covariate)) {
    if (is.character(data[,covariate[i]])){
      data[,covariate[i]]=as.factor(data[,covariate[i]])
    }
  }
}


    data$feature=paste0(data$antigen,data$antibody,sep="")
  ff=as.formula(paste("response~treatment+ antigen+ antibody+ ",paste(covariate,sep="+"),"+
                    treatment:antigen+treatment:antibody+antigen:antibody+treatment:antigen:antibody"))

fff=as.formula(paste("response~treatment+ feature+ ",paste(covariate,sep="+"),"+
                    treatment:feature"))



  if (length(levels(data$treatment))!=2){
    stop("Treatment should have two levels!")
  }
  n1=length(levels(data$antigen))
  n2=length(levels(data$antibody))
  cn=length(covariate)
  n3=0
  if (cn>0){
  for (i in 1:cn) {
    n3=ifelse(is.factor(data[,covariate[i]]),n3+length(levels(data[,covariate[i]]))-1,n3+1)
  }
  }
  ## nn represent the number of terms before second order terms
  nn=(n1-1)+ #antigen
    (n2-1)+  # antibody
    n3+ #covariate
    1+ # treatment
    1  #intercept
  ## nnn represent the number of terms before third order terms
  # nnn=1+ # intercept
  #   1+ # treatment
  #   n1-1+ n2-1+ # antigen and antibody
  #   1*(n1-1)+1*(n2-1)+(n1-1)*(n2-1) + # second order terms
  #   n3 # covariate
  nnn=nn+
    1*(n1-1)+1*(n2-1)+(n1-1)*(n2-1)  # second order terms
  data=groupedData(response~treatment|subjectid,data=data)
  if (!all(type%in% c(1,2,3,4))){stop("Please specify the model in correct way!")}

  control$maxIter=10^10

  if (type==1){
  results=lme(fixed=ff,
              random=pdDiag(~1),
              weights= varIdent(form=~1|feature),
              data = data,method = "ML",
              control = control)
  aa=summary(results)
  coeVec=aa$coefficients$fixed
  varMatrix=aa$varFix
  }


  if (type==2){
    results=lme(fixed=ff,
                random=pdDiag(~antigen),
                data = data,method = "ML",
                control = control)
    aa=summary(results)
    coeVec=aa$coefficients$fixed
    varMatrix=aa$varFix
  }


  if (type==3){
    #gw=gls(model=ff,data = data, weights = varIdent(form=~1|feature),method = "ML")
    gw=gls(model=ff,data = data,method = "ML")
    data$res=gw$residuals
    feature_uniq=unique(data$feature)
    weigh=c()
    for (i in 1:length(feature_uniq)) {
      index=which(data$feature==feature_uniq[i])
      weigh=c(weigh,1/var(data$res[index]))
    }
    wei=data.frame(feature=feature_uniq,weigh=weigh)
    data=merge(data,wei)

    mm=model.matrix(lm(ff,data=data))
    response=data$response
    weigh=data$weigh
    results=cv.glmnet(x=mm,y=response,weights=weigh,nfold=nfold,nlambda = nlambda)
    #results=cv.glmnet(x=mm,y=response,nlambda = nlambda)
    #sslambda=ncol(results$glmnet.fit$beta)
    #print(sslambda)
    #print(results$lambda)
    #coeVec=coef(results$glmnet.fit)[-2,sslambda]
    coeVec=coef(results$glmnet.fit)[-2,results$index[1]]
    varMatrix=diag(1,length(coeVec))
  }

  # if (type==4){
  #   #gw=gls(model=ff,data = data, weights = varIdent(form=~1|feature),method = "ML")
  #   mm=model.matrix(lm(ff,data=data))
  #   response=data$response
  #   results=glmmLasso(fixed=ff, rnd= ,nfold=nfold,nlambda = nlambda)
  #   #results=cv.glmnet(x=mm,y=response,nlambda = nlambda)
  #   #sslambda=ncol(results$glmnet.fit$beta)
  #   #print(sslambda)
  #   #print(results$lambda)
  #   #coeVec=coef(results$glmnet.fit)[-2,sslambda]
  #   coeVec=coef(results$glmnet.fit)[-2,results$index[1]]
  #   varMatrix=diag(1,length(coeVec))
  # }



  sele=list(antigen=c(2,(nn+1):(nn+n1-1)),antibody=c(2,(nn+n1):(nn+n1-1+n2-1)),
            three=c(2,(nnn+1):(nnn+(n1-1)*(n2-1))))
  names(sele$antigen)=levels(data$antigen)
  names(sele$antibody)=levels(data$antibody)
  bb=contrast(sele=sele,coeVec = coeVec,varMatrix = varMatrix)


  mainEffect2=gls(model=ff,data = data,
                  weights = varIdent(form=~1|feature),method = "ML")
  aa=summary(mainEffect2)
  coeVec1=aa$coefficients
  varMatrix1=aa$varBeta
  bbb=contrast(sele=sele,coeVec = coeVec1,varMatrix = varMatrix1)

  # antigen=rownames(bb$effect)
  # antibody=colnames(bb$effect)
  # basic=matrix(nrow=length(antigen),ncol=length(antibody))
  # pp= matrix(nrow=length(antigen),ncol=length(antibody))
  # vv= matrix(nrow=length(antigen),ncol=length(antibody))
  # for (i in 1:length(antigen)) {
  #   for (j in 1:length(antibody)) {
  #     index=which(data$antigen==antigen[i] & data$antibody==antibody[j])
  #     dd=data[index,]
  #     x=dd$response[which(dd$treatment==levels(data$treatment)[1])]
  #     y=dd$response[which(dd$treatment==levels(data$treatment)[2])]
  #     basic[i,j]=t.test(y,x)$estimate[1]-t.test(y,x)$estimate[2]
  #     pp[i,j]=t.test(y,x)$p.value
  #     vv[i,j]=var(x)/length(x)+var(y)/length(y)
  #   }
  # }

  seleFeature1=data.frame(antigen=rep(rownames(bb$effect),ncol(bb$effect)),
                          antibody=rep(colnames(bb$effect),rep(nrow(bb$effect),
                                                               ncol(bb$effect))),
                          effects=c(bb$effect),variance=c(bb$variance),pvalue=c(bb$pvalue))


  seleFeature2=data.frame(antigen=rep(rownames(bb$effect),ncol(bb$effect)),
                          antibody=rep(colnames(bb$effect),rep(nrow(bb$effect),
                                                               ncol(bb$effect))),
                          effects=c(bbb$effect),variance=c(bbb$variance),pvalue=c(bbb$pvalue))

fixedStandard=list(coeVec=coeVec1,varMat=varMatrix1)
fixedDecomposed=list(coeVec=coeVec,varMat=varMatrix)
  results=list(fixedStandard=fixedStandard,fixedDecomposed=fixedDecomposed,formula=ff,
                   selected=seleFeature1,benchmark=seleFeature2)
  return(results)
}











