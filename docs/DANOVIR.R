## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----echo=FALSE---------------------------------------------------------------
library(nlme)
library(ggplot2)
options(useFancyQuotes = FALSE)
par(mar=c(1,1,1,1))

## ----echo=FALSE---------------------------------------------------------------
library(RColorBrewer)
# colors=c("#E04C5C","#7DAF4C","#23AECE", "#FB894B", "#E7DA36",  "#187A51",
#          "#5EA4A2",  "#3D3C4E", "#4D1836", "#C51B7D",
#          "#E9A3C9",  "#B35806", "#F1A340", "#FEE08B", "#D9EF8B",
#          "#91CF60", "#C7EAE5", "#5AB4AC", "#01665E", "#E7D4E8",
#          "#AF8DC3", "#762A83","#FC0FC0","#F9C7DE","#f3a0c4")
# names(colors)=c("IgA2",      "IgA1",      "FcaR",      "FcgR2A131", "FcgR3A158", "IgG" ,      "ELISA IgG", "IgG4",      "IgG2",      "IgA"  ,     "IgG3"  ,    "IgM"    ,   "FcgR1A"  ,  "MN"    ,    "FcgR2b"  ,  "IgG1" ,
#  "HAI",       "NAI"  ,     "FcgR3b","ADCC","ADCP","ADCD")
colors = c(IgG = "#1e90ff",

                   IgA = "#3f007d",

                   IgA1 = "#C7EAE5",

                   IgA2 = "#9e9ac8",

                   IgD = "#662506",

                   IgM = "#8b7765",

                   IgG1 = "#C6E2FF",

                   IgG2 = "#CC3232",

                   IgG3 = "#4F94CD",

                   IgG4 = "#1D3700",

                   FcgR2A131 = "#157B74",

                   FcgR2b = "#CDCD00",

                   FcgR3A158 = "#722841",

                   FcgR3b = "#23AECE",

                   FcaR = "#bcbddc",
            FcRn = "#7B5D56",

                   Func="#000000",
           
           ADCC= "#80CBC4",
           ADCP="#A5D6A7",
           ADCD="#2166AC",
           MN= "#23AECE",
           HAI=  "#C7EAE5",
           NAI=  "#5AB4AC",
          "ELISA IgG"= "#7DAF4C",
           CA07="#1B4F72",
            FcgR2B = "#CDCD00",
           FcgR1A="#AA4371",

                   FcgR3A = "#722841",

                   FcgR3B = "#23AECE",
           FcgR2A = "#458B74",
          
 "ELISA"= "#7DAF4C"
           )
shapes=c(15,17,18)
names(shapes)=c("HASK",        "NA",         "HAFL")

## ----setup--------------------------------------------------------------------
library(danovir)

## ----echo=FALSE---------------------------------------------------------------
colnames(dataset)[5]="treatment"
colnames(dataset)[6]="subjectid"
colnames(dataset)[13]="response"
colnames(dataset)[15]="antibody"
dataset$treatment=as.factor(dataset$treatment)

## ----message=FALSE------------------------------------------------------------
results=danovir(data=dataset,covariate=c("assay"),type = c(1,0,1,0),control = control)

## ----echo=FALSE,message=FALSE-------------------------------------------------
volcanoData=results$selected
volcanoData$fdr=p.adjust(volcanoData$pvalue,method = "BH")
volcanoData=volcanoData[order(volcanoData$pvalue),]


cutoff <- ifelse(volcanoData$fdr<0.2,0.2,volcanoData$fdr)
index=max(which(volcanoData$fdr<0.2))
h=-log10(volcanoData$pvalue[index])
xx=max(abs(volcanoData$effects))
rownames(volcanoData)=c()
volcanoData
pp=ggplot(data=volcanoData,aes(x=effects,y=-log10(pvalue)))+
    geom_point(aes(shape=antigen,color=antibody),size=4)+
     scale_shape_manual(values = shapes)  +
    scale_color_manual(values = colors)+
   theme_classic()+
  theme(legend.key=element_blank())+
  geom_hline(yintercept=h,linetype="dashed")+
  geom_vline(xintercept=0,linetype="dashed")+
  labs(x=bquote("Fold Change"), y=bquote(-Log[10](p-value)))+
  #ggtitle( paste0("Day0: NP vs P15(R^2=",round(rsqr,2),")" ) )+
   theme(legend.position = "right",
          legend.box = "vertical")+
   guides(color = guide_legend(nrow = 8,title="Detection Reagent"))+
    geom_text(aes(xx/2,h,label = paste("FDR=", 0.2, sep=""), vjust = -1),show.legend = FALSE)
print(pp)

## ----echo=FALSE,message=FALSE-------------------------------------------------
volcanoData=results$benchmark
volcanoData$fdr=p.adjust(volcanoData$pvalue,method = "BH")
volcanoData=volcanoData[order(volcanoData$pvalue),]


cutoff <- ifelse(volcanoData$fdr<0.2,0.2,volcanoData$fdr)
index=max(which(volcanoData$fdr<0.2))
h=-log10(volcanoData$pvalue[index])
xx=max(abs(volcanoData$effects))
rownames(volcanoData)=c()
volcanoData
ggplot(data=volcanoData,aes(x=effects,y=-log10(pvalue)))+
    geom_point(aes(shape=antigen,color=antibody),size=4)+
     scale_shape_manual(values = shapes)  +
    scale_color_manual(values = colors)+
   theme_classic()+
  #xlim(-1,1)+
  theme(legend.key=element_blank())+
  geom_hline(yintercept=h,linetype="dashed")+
  geom_vline(xintercept=0,linetype="dashed")+
  labs(x=bquote("Fold Change"), y=bquote(-Log[10](p-value)))+
    #ggtitle( paste0("Day0: P15 vs P30(R^2=",round(rsqr,2),")" ) )+
   theme(legend.position = "right",
          legend.box = "vertical")+
   guides(color = guide_legend(nrow = 8,title="Detection Reagent"))+
    geom_text(aes(xx/2,h,label = paste("FDR=", 0.2, sep=""), vjust = -1),show.legend = FALSE)

## -----------------------------------------------------------------------------
plot(c(results$selected$effects),c(results$benchmark$effects),xlab = "Model-based effects",ylab = "Original effects")
abline(coef = c(0,1))
cor(c(results$selected$effects),c(results$benchmark$effects))
plot(c(results$selected$variance),c(results$benchmark$variance),
     xlab = "Model-based variances",
     ylab = "Original variances")
abline(coef = c(0,1))

