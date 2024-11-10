
#' Denoised analysis of vaccine-induced immune response
#' @description This is the main function of the package. This function outputs the
#'   denoised estimate of the treatment effects and their variances. The traditional estimates and their variances are also provided which
#'   can be used to determine how the model based method affects the estimates.
#'
#' @param data data set that include all the relevant variables
#' @param type a vector of length 4 for which the first component specify
#' if \code{antigen},\code{antibody} interaction should be include in the model,
#' the second component specifies whether three factor interaction
#'  \code{treatment-antigen-antibody} should be included in the model,
#'  the third component specifies whether \code{antigen} should be included
#'   in the random effects,  and the fourth component specify whether
#'    \code{antibody} should be included in the random effects.
#' @param covariate variables that need to be included in the model
#' @param reproduce If true, then model include a random intercept and
#' variance of error term is feature wise constant. If false, then model include
#' random intercept, antigen effect and antibody effect. Variance of errow term is constant.
#' @return A \code{danovir} object which is a list of length 3. The first
#' component  gives the fitting results of the mixed-effects model;
#' the second component gives the model-based estimates of treatment effects,
#' while the third component gives the results from traditional approach.
#' @export
#' @details The underlying basic model in the package is
#'
#'
#' \code{response_{ijk}=mu+alpha_j+beta_k+lambda_l+gamma_{jl}+gamma_{kl}+u_i+e_{ijk}}
#'
#'
#' Here the meanings of the terms are, respectively,
#'
#'
#'  - \code{mu}, the intercept
#'
#'
#'   - \code{alpha_j} , the main antigen effects
#'
#'
#'   - \code{beta_k} , the main antibody effects
#'
#'
#'   - \code{lambda_l} , the main treatment effects
#'
#'
#'   - \code{gamma_{jl}} , the interaction between treatment and antigen
#'
#'
#'   - \code{gamma_{kl}} , the interaction between treatment and antibody
#'
#'
#'  - \code{u_i} , the random intercept for subject
#'
#'
#'   - \code{e_{ijk}} , the error term
#'
#'
#' There are four additional terms that are optional which
#' should be chosen by users. These choices
#'  should be made based on the characteristics
#'   of the data at hand. These four terms are, respectively
#'
#'
#'  - \code{gamma_{jk}}, interaction between antigen and antibody
#'
#'
#'  - \code{gamma_{jkl}}, interaction between treatment, antigen and antibody
#'
#'
#'  -  \code{u_{ij}}, random antigen effects
#'
#'
#'  -  \code{u_{ik}}, random antibody effects


danovir=function(data,type,covariate,reproduce=FALSE,...){
  fullnames=colnames(data)
  if (!all(c("response","subjectid","treatment","assay","antigen","antibody",covariate)
           %in% fullnames)){
    stop("Please re-specify the column names of data which should include variables response,
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
    if (all(type==c(0,0,0,0))){
      ff=as.formula(paste("response~treatment+ antigen+ antibody+ ",paste(covariate,sep="+"),"+
                    treatment:antigen+treatment:antibody"))
      results=lme(fixed=ff,
                  random=pdDiag(~1),data = data,method = "ML",
                  control = control)
      sele=list(antigen=c(2,(nn+1):(nn+n1-1)),antibody=c(2,(nn+n1):(nn+n1-1+n2-1)))
    }
    if (all(type==c(0,0,0,1))){
      ff=as.formula(paste("response~treatment+ antigen+ antibody+ ",paste(covariate,sep="+"),"+
                    treatment:antigen+treatment:antibody"))
      results=lme(fixed=ff,
                  random=pdDiag(~antibody),data = data,method = "ML",
                  control = control)
      sele=list(antigen=c(2,(nn+1):(nn+n1-1)),antibody=c(2,(nn+n1):(nn+n1-1+n2-1)))
    }
    if (all(type==c(0,0,1,0))){
      ff=as.formula(paste("response~treatment+ antigen+ antibody+ ",paste(covariate,sep="+"),"+
                    treatment:antigen+treatment:antibody"))
      results=lme(fixed=ff,
                  random=pdDiag(~antigen),data = data,method = "ML",
                  control = control)
      sele=list(antigen=c(2,(nn+1):(nn+n1-1)),antibody=c(2,(nn+n1):(nn+n1-1+n2-1)))
    }
    if (all(type==c(0,0,1,1))){
      ff=as.formula(paste("response~treatment+ antigen+ antibody+ ",paste(covariate,sep="+"),"+
                    treatment:antigen+treatment:antibody"))
      results=lme(fixed=ff,
                  random=pdDiag(~antigen+antibody),data = data,method = "ML",
                  control = control)
      sele=list(antigen=c(2,(nn+1):(nn+n1-1)),antibody=c(2,(nn+n1):(nn+n1-1+n2-1)))
    }




    if (all(type==c(0,1,0,0))){
      ff=as.formula(paste("response~treatment+ antigen+ antibody+ ",paste(covariate,sep="+"),"+
                    treatment:antigen+treatment:antibody+treatment:antigen:antibody"))
      results=lme(fixed=ff,
                  random=pdDiag(~1),data = data,method = "ML",
                  control = control)
      sele=list(antigen=c(2,(nn+1):(nn+n1-1)),antibody=c(2,(nn+n1):(nn+n1-1+n2-1)),
                three=c(2,(nnn+1):(nnn+(n1-1)*(n2-1))))
    }
    if (all(type==c(0,1,0,1))){
      ff=as.formula(paste("response~treatment+ antigen+ antibody+ ",paste(covariate,sep="+"),"+
                    treatment:antigen+treatment:antibody+treatment:antigen:antibody"))
      results=lme(fixed=ff,
                  random=pdDiag(~antibody),data = data,method = "ML",
                  control = control)
      sele=list(antigen=c(2,(nn+1):(nn+n1-1)),antibody=c(2,(nn+n1):(nn+n1-1+n2-1)),
                three=c(2,(nnn+1):(nnn+(n1-1)*(n2-1))))
    }
    if (all(type==c(0,1,1,0))){
      ff=as.formula(paste("response~treatment+ antigen+ antibody+ ",paste(covariate,sep="+"),"+
                    treatment:antigen+treatment:antibody+treatment:antigen:antibody"))
      results=lme(fixed=ff,
                  random=pdDiag(~antigen),data = data,method = "ML",
                  control = control)
      sele=list(antigen=c(2,(nn+1):(nn+n1-1)),antibody=c(2,(nn+n1):(nn+n1-1+n2-1)),
                three=c(2,(nnn+1):(nnn+(n1-1)*(n2-1))))
    }
    if (all(type==c(0,1,1,1))){
      ff=as.formula(paste("response~treatment+ antigen+ antibody+ ",paste(covariate,sep="+"),"+
                    treatment:antigen+treatment:antibody+treatment:antigen:antibody"))
      results=lme(fixed=ff,
                  random=pdDiag(~antigen+ antibody),data = data,method = "ML",
                  control = control)
      sele=list(antigen=c(2,(nn+1):(nn+n1-1)),antibody=c(2,(nn+n1):(nn+n1-1+n2-1)),
                three=c(2,(nnn+1):(nnn+(n1-1)*(n2-1))))
    }






    if (all(type==c(1,0,0,0))){
      ff=as.formula(paste("response~treatment+ antigen+ antibody+ ",paste(covariate,sep="+"),"+
                    treatment:antigen+treatment:antibody+antigen:antibody"))
      results=lme(fixed=ff,
                  random=pdDiag(~1),data = data,method = "ML",
                  control = control)
      sele=list(antigen=c(2,(nn+1):(nn+n1-1)),antibody=c(2,(nn+n1):(nn+n1-1+n2-1)))
    }
    if (all(type==c(1,0,0,1))){
      ff=as.formula(paste("response~treatment+ antigen+ antibody+ ",paste(covariate,sep="+"),"+
                    treatment:antigen+treatment:antibody+antigen:antibody"))
      results=lme(fixed=ff,
                  random=pdDiag(~antibody),data = data,method = "ML",
                  control = control)
      sele=list(antigen=c(2,(nn+1):(nn+n1-1)),antibody=c(2,(nn+n1):(nn+n1-1+n2-1)))
    }
    if (all(type==c(1,0,1,0))){
      ff=as.formula(paste("response~treatment+ antigen+ antibody+ ",paste(covariate,sep="+"),"+
                    treatment:antigen+treatment:antibody+antigen:antibody"))
      results=lme(fixed=ff,
                  random=pdDiag(~antigen),data = data,method = "ML",
                  control = control)
      sele=list(antigen=c(2,(nn+1):(nn+n1-1)),antibody=c(2,(nn+n1):(nn+n1-1+n2-1)))
    }
    if (all(type==c(1,0,1,1))){
      ff=as.formula(paste("response~treatment+ antigen+ antibody+ ",paste(covariate,sep="+"),"+
                    treatment:antigen+treatment:antibody+antigen:antibody"))
      results=lme(fixed=ff,
                  random=pdDiag(~antigen+antibody),data = data,method = "ML",
                  control = control)
      sele=list(antigen=c(2,(nn+1):(nn+n1-1)),antibody=c(2,(nn+n1):(nn+n1-1+n2-1)))
    }




    if (all(type==c(1,1,0,0))){
      ff=as.formula(paste("response~treatment+ antigen+ antibody+ ",paste(covariate,sep="+"),"+
                    treatment:antigen+treatment:antibody+antigen:antibody+treatment:antigen:antibody"))
      results=lme(fixed=ff,
                  random=pdDiag(~1),data = data,method = "ML",
                  control = control)
      sele=list(antigen=c(2,(nn+1):(nn+n1-1)),antibody=c(2,(nn+n1):(nn+n1-1+n2-1)),
                three=c(2,(nnn+1):(nnn+(n1-1)*(n2-1))))
    }
    if (all(type==c(1,1,0,1))){
      ff=as.formula(paste("response~treatment+ antigen+ antibody+ ",paste(covariate,sep="+"),"+
                    treatment:antigen+treatment:antibody+antigen:antibody+treatment:antigen:antibody"))
      results=lme(fixed=ff,
                  random=pdDiag(~antibody),data = data,method = "ML",
                  control = control)
      sele=list(antigen=c(2,(nn+1):(nn+n1-1)),antibody=c(2,(nn+n1):(nn+n1-1+n2-1)),
                three=c(2,(nnn+1):(nnn+(n1-1)*(n2-1))))
    }
    if (all(type==c(1,1,1,0))){
      ff=as.formula(paste("response~treatment+ antigen+ antibody+ ",paste(covariate,sep="+"),"+
                    treatment:antigen+treatment:antibody+antigen:antibody+treatment:antigen:antibody"))
      results=lme(fixed=ff,
                  random=pdDiag(~antigen),data = data,method = "ML",
                  control = control)
      sele=list(antigen=c(2,(nn+1):(nn+n1-1)),antibody=c(2,(nn+n1):(nn+n1-1+n2-1)),
                three=c(2,(nnn+1):(nnn+(n1-1)*(n2-1))))
    }
    if (all(type==c(1,1,1,1)) & reproduce==FALSE){
      ff=as.formula(paste("response~treatment+ antigen+ antibody+ ",paste(covariate,sep="+"),"+
                    treatment:antigen+treatment:antibody+antigen:antibody+treatment:antigen:antibody"))
      results=lme(fixed=ff,
                  random=pdDiag(~antigen+antibody),data = data,method = "ML",
                  control = control)
      sele=list(antigen=c(2,(nn+1):(nn+n1-1)),antibody=c(2,(nn+n1):(nn+n1-1+n2-1)),
                three=c(2,(nnn+1):(nnn+(n1-1)*(n2-1))))
    }

      if (all(type==c(1,1,1,1)) & reproduce==TRUE){
        data$feature=paste0(data$assay,data$antigen,data$detect_reagent,sep="")
        ff=as.formula(paste("response~treatment+ antigen+ antibody+ ",paste(covariate,sep="+"),"+
                    treatment:antigen+treatment:antibody+antigen:antibody+treatment:antigen:antibody"))
        results=lme(fixed=ff,
                    random=pdDiag(~1),weights= varIdent(form=~1|feature),
                    data = data,method = "ML",
                    control = control)
        sele=list(antigen=c(2,(nn+1):(nn+n1-1)),antibody=c(2,(nn+n1):(nn+n1-1+n2-1)),
                  three=c(2,(nnn+1):(nnn+(n1-1)*(n2-1))))
      }
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
