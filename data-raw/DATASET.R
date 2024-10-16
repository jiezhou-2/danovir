
library(nlme)
source("./data-raw/data_cleaning_no_buffer.R")
## exclude non Lumbaa assay
device=do.call(rbind,strsplit(combined_experiments$detect_method,split="_"))
combined_experiments=cbind(combined_experiments,assay=device[,1])
index=which(device[,1]=="LumBAA" | device[,1]=="NAI" | device[,1]=="MN" | device[,1]=="HAI" | device[,1]=="ELISA"
            | device[,1]=="ADCC" | device[,1]=="ADCP" | (device[,1]=="ADCD" & device[,3]=="dil3200" & device[,2]=="set1" ))

combined_experiments=combined_experiments[index,]
index=which(combined_experiments$detect_method=="ADCD_set1_dil3200")
combined_experiments$detect_method[index]="ADCD"



## exclude non ca data
index=grep("CA",combined_experiments$analyte_reported)
combined_experiments=combined_experiments[index,]
## exclude HA_dtm
index=which(combined_experiments$analyte_type=="HA_dTM")
combined_experiments=combined_experiments[-index,]
## create features
temp = paste(combined_experiments$assay, combined_experiments$analyte_reported, sep="_")
temp=ifelse(temp=="NAI_NAI_CA07","NAI_CA07", ifelse(temp=="MN_MN_CA07","MN_CA07",ifelse(temp=="HAI_HAI_CA07","HAI_CA07",
                                                                                        ifelse(temp=="ADCC_ADCC_CA07","ADCC_CA07",ifelse(temp=="ADCP_ADCP_CA07","ADCP_CA07",
                                                                                                                                         ifelse(temp=="ADCD_ADCD_CA07_HA","ADCD_CA07_HA",ifelse(temp=="ADCD_ADCD_CA07_HA_FL","ADCD_CA07_HA_FL",ifelse(temp=="ADCD_ADCD_CA07_HA_stalk","ADCD_CA07_HA_stalk",ifelse(temp=="ADCD_ADCD_CA07_NA","ADCD_CA07_NA",temp)))))))))
combined_experiments$temp=temp
## for each assay, assign a shape number for plotting in the following.
#assay=unique(combined_experiments$assay)
shapes=c(15,17,18)
names(shapes)=c("HASK",        "NA",         "HAFL")





index=which(combined_experiments$arm_name=="Pregnant 30 mcg" | combined_experiments$visit_name=="Delivery Visit" | combined_experiments$visit_name=="Visit 2")
data15=combined_experiments[-index,]
data15$arm=as.factor(ifelse(data15$arm_name=="Pregnant 15 mg","P15","NP15"))
data15$antigen= ifelse(data15$analyte_type=="HA","HA",ifelse(data15$analyte_type=="HA_FL","HAFL", ifelse(data15$analyte_type=="HA_stalk","HASK",ifelse(data15$analyte_type=="live_virus", "LV","NA"))))
index=which(data15$antigen=="LV")
data15=data15[-index,]
data15$antigen=factor(data15$antigen,levels = c("NA","HA","HAFL","HASK"))
data15$value_reported=log2(data15$value_reported+1)
data15$detect_reagent=factor(data15$detect_reagent,levels=c("IgG", "FcaR",       "FcgR1A", "FcgR2A131", "FcgR2b", "FcgR3A158", "FcgR3b", "IgA", "IgA1", "IgA2", "IgG1", "IgG2", "IgG3", "IgG4", "IgM"))
data15$assay=as.factor(data15$assay)
index=which(data15$antigen=="HA")
data15=data15[-index,]
data15$antigen=droplevels(data15$antigen)
data15$sub001=ifelse(data15$subject_accession=="SUB8888001",1,0)
control=lmeControl(maxIter = 100,msMaxIter = 100,niterEM = 100)
cohort1=unique(data15$subject_accession[which(data15$arm=="NP15")])[1:10]
cohort2=unique(data15$subject_accession[which(data15$arm=="P15")])[1:10]
cohort=c(cohort1,cohort2)
index=which(data15$subject_accession %in% cohort)
dataset=data15[index,]
usethis::use_data(dataset, overwrite = TRUE)
usethis::use_data(control, overwrite = TRUE)
