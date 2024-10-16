
#source("./code/connection.R", local=TRUE)

library(DBI)
library(DT)
library(RPostgres)

## Read sample data
conn <- dbConnect(RPostgres::Postgres(), dbname = Sys.getenv("db"),
                  host=Sys.getenv("db_host"), port=Sys.getenv("db_port"), user=Sys.getenv("db_userid_x"),
                  password=Sys.getenv("db_pwd_x"), options="-c search_path=madi_dat")

#selected_study <- paste0("SDY8001")
#safe_study <- dbQuoteLiteral(conn, selected_study)

### MADI combined set of all experimental data for the SDY8001 project ###
### truncated set of columns to make more manageable ###
select_query <- "SELECT study_accession, experiment_accession, experiment_name, measurement_technique, arm_name, subject_accession, gender, biosample_accession, biosample_type, visit_name, unit_reported, analyte_reported, value_reported, detect_method, detect_reagent, analyte_type, virus_strain_abrev, virus_strain
,CASE
WHEN analyte_type IN ('HA_FL','HA_stalk','HA','HA1') AND virus_strain_abrev = 'CA07' THEN 'primary vaccine antigen'
WHEN analyte_type IN ('NA')  AND virus_strain_abrev = 'CA07' THEN 'secondary vaccine antigen'
ELSE 'not vaccine related/not focus' END AS antigen_class
FROM madi_results.mv_sdy8001_ext888011;"

combined_experiments <- DBI::dbGetQuery(conn,select_query)
temp <- paste(combined_experiments$detect_method, combined_experiments$analyte_reported, sep="_")
combined_experiments=cbind(combined_experiments,temp)
#analyte_sample=unique(combined_experiments$analyte_reported)
combined_experiments=unique(combined_experiments)


select_query <- "SELECT control_sample_accession, CASE WHEN position('_CB_' IN control_sample.assay_id) > 0 THEN 'cord blood' ELSE 'serum' END AS sample_type, dilution_factor
, result_schema, source,  source_type, mbaa_result.analyte_reported, mbaa_result.assay_id, mbaa_result.mfi, j.*
  FROM madi_dat.control_sample
INNER JOIN madi_dat.mbaa_result ON mbaa_result.source_accession = control_sample.control_sample_accession
LEFT OUTER JOIN madi_results.vw_analyte_annotations AS j ON j.revised_analyte = mbaa_result.analyte_reported
WHERE control_sample.experiment_accession = 'EXP888010'"


control_sample_exp888010  <- DBI::dbGetQuery(conn,select_query)
index=which(control_sample_exp888010$source=="Buffer")
control_sample_exp888010=control_sample_exp888010[index,]
value_reported=control_sample_exp888010$mfi
control_sample_exp888010=cbind(control_sample_exp888010,value_reported)
control_sample_exp888010$value_reported=as.numeric(control_sample_exp888010$value_reported)
index=which(control_sample_exp888010$value_reported<=1500)
control_sample_exp888010=control_sample_exp888010[index,]

select_query_ADCD <- "SELECT control_sample.control_sample_accession, 'ADCD' AS detect_method, assay_group_id, dilution_factor, control_sample.experiment_accession, source, virus_strain_reported AS analyte_reported, value_reported
	FROM madi_dat.control_sample
	INNER JOIN madi_dat.hai_control ON hai_control.control_sample_accession = control_sample.control_sample_accession
	WHERE control_sample.experiment_accession = 'EXP888042';"
control_ADCD= DBI::dbGetQuery(conn,select_query_ADCD)
index=which(control_ADCD$source=="Buffer")
control_ADCD=control_ADCD[index,]


buffername=unique(control_sample_exp888010$analyte_reported)
bufferaverage=c()
for (i in 1:length(buffername)) {
  index=which(control_sample_exp888010$analyte_reported==buffername[i])
  meani=mean(as.numeric(control_sample_exp888010$value_reported[index]),na.rm = T)
  #print(meani)
  #sdi=ifelse(length(index)==1,0,sd(as.numeric(control_sample_exp888010$value_reported[index]),na.rm = T))
  bufferaverage=c(bufferaverage, meani)

}
names(bufferaverage)=buffername



buffername_ADCD=unique(control_ADCD$analyte_reported)
bufferaverage_ADCD=c()
for (i in 1:length(buffername_ADCD)) {
  index=which(control_ADCD$analyte_reported==buffername_ADCD[i])

  meani_ADCD=mean(as.numeric(control_ADCD$value_reported[index]),na.rm = T)
  #print(meani_ADCD)
  bufferaverage_ADCD=c(bufferaverage_ADCD,meani_ADCD)
}
names(bufferaverage_ADCD)=buffername_ADCD


# select_delivery_set <- "SELECT * ,
# CASE
# WHEN antigen_type IN ('HA_FL','HA_stalk','HA','HA1') AND virus_strain_abrev = 'CA07' THEN 'primary vaccine antigen'
# WHEN antigen_type IN ('NA')  AND virus_strain_abrev = 'CA07' THEN 'secondary vaccine antigen'
# ELSE 'not vaccine related/not focus' END AS antigen_class
# FROM madi_results.mv_sdy8001_ext888053;"

select_delivery_set <- "SELECT * ,
CASE
WHEN antigen_type IN ('HA_FL','HA_stalk','HA','HA1','HAFull', 'HAstalk') AND virus_strain_abrev = 'CA07' THEN 'primary vaccine antigen'
WHEN antigen_type IN ('NA')  AND virus_strain_abrev = 'CA07' THEN 'secondary vaccine antigen'
ELSE 'not vaccine related/not focus' END AS antigen_class,
antigen_type AS analyte_type
FROM madi_results.mv_sdy8001_ext888053;"



combined_experiments_single <- DBI::dbGetQuery(conn,select_delivery_set)

combined_experiments_single <- DBI::dbGetQuery(conn,select_delivery_set)
combined_experiments_single$antigen_type <- ifelse(combined_experiments_single$antigen_type == "HAstalk","HA_stalk",combined_experiments_single$antigen_type)
combined_experiments_single$antigen_type <- ifelse(combined_experiments_single$antigen_type == "HAFull","HA_FL",combined_experiments_single$antigen_type)




select_delivery_set_control <- "SELECT control_sample_accession, dilution_factor
, result_schema, source,  source_type, mbaa_result.analyte_reported, mbaa_result.assay_id, mbaa_result.mfi
	FROM madi_dat.control_sample
	INNER JOIN madi_dat.mbaa_result ON mbaa_result.source_accession = control_sample.control_sample_accession
	WHERE control_sample.experiment_accession = 'EXP888053';"

combined_experiments_single_control <- DBI::dbGetQuery(conn,select_delivery_set_control)


rm(list=setdiff(ls(),c("combined_experiments","bufferaverage","bufferaverage_ADCD",
                       "combined_experiments_single","combined_experiments_single_control")))






