

danovir=function(response,treatment,antigen,antibody,subjectid,covariate=NULL,
                 corr=0,denoise=FALSE,nfold=10,nlambda=100){

  if (!is.null(covariate)){
    for (i in 1:length(covariate)) {
      if (is.character(covariate[i])){
        covariate[i]=as.factor(covariate[i])
      }
    }
  }


  feature=paste0(antigen,antibody,sep="")
  ff=as.formula(paste("response~treatment+ antigen+ antibody+ ",paste(covariate,sep="+"),"+
                    treatment:antigen+treatment:antibody+antigen:antibody+treatment:antigen:antibody"))
#
#   fff=as.formula(paste("response~treatment+ feature+ ",paste(covariate,sep="+"),"+
#                     treatment:feature"))



  if (length(unique(treatment))!=2){
    stop("Treatment should have two levels!")
  }
  n1=length(levels(antigen))
  n2=length(levels(antibody))
  cn=length(covariate)
  n3=0
  if (cn>0){
    for (i in 1:cn) {
      if (is.numeric(covariate[i])) {
        n3=n3+1
      }else{
        covariate[i]=as.factor(covariate[i])
        n3=n3+length(levels(covariate[i]))-1
      }
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

  data=cbind.data.frame(response=response,treatment=treatment,antigen=antigen,antibody=antibody,
                  subjectid=subjectid,covariate)
  data=groupedData(response~treatment|subjectid,data=data)


  control$maxIter=10^10

  if (denoise==F & corr==0){
  results=gls(value_reported~arm+ antigen+detect_reagent+assay+ arm:antigen+arm:detect_reagent+antigen:detect_reagent+arm:antigen:detect_reagent,data = data15 ,weights = varIdent(form=~1|feature),method = "ML")
  aa=summary(results)
  coeVec=aa$coefficients
  varMatrix=aa$varBeta
}


  if (denoise==F & corr=="uniform"){
    results=lme(fixed=ff,
                random=pdDiag(~1),
                weights= varIdent(form=~1|feature),
                data = data,method = "ML",
                control = control)
    aa=summary(results)
    coeVec=aa$coefficients$fixed
    varMatrix=aa$varFix
  }


  if (denoise==F & corr=="antigen"){
    results=lme(fixed=ff,
                random=pdDiag(~antigen),
                weights= varIdent(form=~1|feature),
                data = data,method = "ML",
                control = control)
    aa=summary(results)
    coeVec=aa$coefficients$fixed
    varMatrix=aa$varFix
  }

  if (denoise==F & corr=="antibody"){
    results=lme(fixed=ff,
                random=pdDiag(~antibody),
                weights= varIdent(form=~1|feature),
                data = data,method = "ML",
                control = control)
    aa=summary(results)
    coeVec=aa$coefficients$fixed
    varMatrix=aa$varFix
  }



  # if (type==3){
  #   #gw=gls(model=ff,data = data, weights = varIdent(form=~1|feature),method = "ML")
  #   gw=gls(model=ff,data = data,method = "ML")
  #   data$res=gw$residuals
  #   feature_uniq=unique(data$feature)
  #   weigh=c()
  #   for (i in 1:length(feature_uniq)) {
  #     index=which(data$feature==feature_uniq[i])
  #     weigh=c(weigh,1/var(data$res[index]))
  #   }
  #   wei=data.frame(feature=feature_uniq,weigh=weigh)
  #   data=merge(data,wei)
  #
  #   mm=model.matrix(lm(ff,data=data))
  #   response=data$response
  #   weigh=data$weigh
  #   results=cv.glmnet(x=mm,y=response,weights=weigh,nfold=nfold,nlambda = nlambda)
  #   #results=cv.glmnet(x=mm,y=response,nlambda = nlambda)
  #   #sslambda=ncol(results$glmnet.fit$beta)
  #   #print(sslambda)
  #   #print(results$lambda)
  #   #coeVec=coef(results$glmnet.fit)[-2,sslambda]
  #   coeVec=coef(results$glmnet.fit)[-2,results$index[1]]
  #   varMatrix=diag(1,length(coeVec))
  # }
  #
  # if (type==4){
  #   #gw=gls(model=ff,data = data, weights = varIdent(form=~1|feature),method = "ML")
  #   mm=model.matrix(lm(ff,data=data))
  #   response=data$response
  #   results=glmmLasso(fix=ff, rnd=list(subjectid~1+antigen),data=data,lambda = lambda)
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











