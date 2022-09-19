R
library(reshape)
library(plyr)
library(data.table)
library(tidyverse)
library(stringr)
library(dplyr)
library(PheWAS)

bd<-fread("/ukb39726.tab.gz", header=TRUE, sep="\t")
extract<-bd[,c(1,11032,12,10968,10976,11037,11040:11049)]
colnames(extract)<-c("eid","sex","birth_year","age_assessment","age_recruitment","genetic_ethnicity","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")
qc<-read.table("/ukb_sqc_v2.txt.gz", as.is=TRUE, header=FALSE, stringsAsFactors=FALSE)
fam<-read.table("/ukb_cal_v2.fam", as.is=TRUE, header=FALSE, stringsAsFactors=FALSE)

qc1<-cbind(fam,qc)
qc1<-qc1[,c(1,9)]
qc1$bileve<-ifelse(qc1$V3=="UKBL",1,0)

extract2<-merge(extract,qc1[,c(1,3)],by.x="eid",by.y="V1")
extract2[is.na(extract2)] <- "NA"

write.table(extract2,file="/UKB.covar.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote = FALSE)


df1<- fread("./hesin_diag_24Jun2021.txt.gz",header=TRUE, stringsAsFactors=FALSE)
colnames(df1)
df1 <- df1[, -c(2:6, 8)]
df1 <- df1 %>% count(eid, diag_icd10)
colnames(df1)[1] <- "id"
colnames(df1)[2] <- "code"
colnames(df1)[3] <- "index"
df1 <- df1[as.character(df1$code)!= "" ,] #remove entries with no icd10
df1$vocabulary_id <- "ICD10CM"
df1 <- df1[,c(1,4,2,3)]
df1$code <- sub("((?:^a.|^[^a]).{2})(.+)", "\\1.\\2", df1$code)

samples<-read.table("/ukb_imp_v3.sample", as.is=TRUE,header=FALSE, stringsAsFactors=FALSE,skip=2)
pheno<-createPhenotypes(df1, min.code.count=1,add.phecode.exclusions=F, translate=T,full.population.ids=unique(samples[[1]]),aggregate.fun=PheWAS:::default_code_agg,vocabulary.map=PheWAS::phecode_map,rollup.map=PheWAS::phecode_rollup_map,exclusion.map=PheWAS::phecode_exclude)

write.table(pheno,file="/HES_icd10_pheno.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote = FALSE)


bd<-fread("/ukb39726.tab.gz", header=TRUE, sep="\t")

extract<-bd[,c(1,412,1595,1603,6942,7142,11186,11187,11188,11193,11194,11195,11196,11204,11205,11210,11211,12475,12478,12481,12484,12705,12707,12709,12711,12713,12715,12717,12769,12771,12773,12775,12865,12867,12870,12871,13480,13493,13506,13519,13532,13545,13558,13571,13584,13597,13610,13623,13636,13649,13662,13675,13688,13701,13714,13727,13740,13753,13766,13779,13792,13805,13818,13831,13844,13857,13870,13901,13915,13929,13943,13957,13971,13985,13999,14013,14027,14041,14055,14069,14083,14097,14111,14123,14137,14151,14165,14179,14193,14207,14221,14235,14249,14263,14277,14291,14305)]

bd<-fread("/ukb41208.tab.gz", header=TRUE, sep="\t")
extract<-merge(extract,bd,by="f.eid",all.x=TRUE)
bd<-fread("/ukb44425.tab.gz", header=TRUE, sep="\t")
extract<-merge(extract,bd[,c(1,2,6,9,13)],by="f.eid",all.x=TRUE)

colnames(extract)<-c("eid","Pulse_rate","DBP","SBP","Rx.1","Rx.2","Liver_Fe","LIF","PDFF","VAT","ASAT","Total_thigh_muscle_volume","Total_trunk_fat_volume","Total_adipose_tissue_volume","Total_lean_tissue_volume","Cardiac_output","Cardiac_index","Mean_carotid_IMT_120d","Mean_carotid_IMT_150d","Mean_carotid_IMT_210d","Mean_carotid_IMT_240d","Body_fat_percentage","Whole_body_fat_mass","Whole_body_fat_free_mass","Whole_body_water_mass","BMI","BMR","Whole_body_Impedance","Trunk_fat_percentage","Trunk_fat_mass","Trunk_fat_free_mass","Trunk_predicted_mass","Total_tissue_fat_percentage","Total_mass","Trunk_tissue_fat_percentage","Trunk_total_mass","White_blood_cell_leukocyte_count","Red_blood_cell_erythrocyte_count","Haemoglobin_concentration","Haematocrit_percentage","Mean_corpuscular_volume","Mean_corpuscular_haemoglobin","Mean_corpuscular_haemoglobin_concentration","Red_blood_cell_erythrocyte_distribution_width","Platelet_count","Platelet_crit","Mean_platelet_thrombocyte_volume","Platelet_distribution_width","Lymphocyte_count","Monocyte_count","Neutrophill_count","Eosinophill_count","Basophill_count","Nucleated_red_blood_cell_count","Lymphocyte_percentage","Monocyte_percentage","Neutrophill_percentage","Eosinophill_percentage","Basophill_percentage","Nucleated_red_blood_cell_percentage","Reticulocyte_percentage","Reticulocyte_count","Mean_reticulocyte_volume","Mean_sphered_cell_volume","Immature_reticulocyte_fraction","High_light_scatter_reticulocyte_percentage","High_light_scatter_reticulocyte_count","Albumin","Alkaline_phosphatase","Alanine_aminotransferase","ApoA","ApoB","Aspartate_aminotransferase","Direct_bilirubin","Urea","Calcium","Cholesterol","Creatinine","CRP","Cystatin_C","Gamma_glutamyltransferase","Glucose","HbA1c","HDL_cholesterol","IGF1","LDL_direct","Lipoprotein_A","Oestradiol","Phosphate","Rheumatoid_factor","SHBG","Total_bilirubin","Testosterone","Total_protein","Triglycerides","Urate","VitaminD","cT1","Urine_Microalbumin","Urine_Creatinine","Urine_Potassium","Urine_Sodium")

extract$Rx.1<-as.numeric(extract$Rx.1)
extract$Rx.2<-as.numeric(extract$Rx.2)

extract$Statins<-ifelse(extract$Rx.1==1 | extract$Rx.2==1, "1", "0")
extract[is.na(extract)] <- "NA"
extract$Statins<-gsub("NA", 0, extract$Statins)

write.table(extract,file="/UKB.lab.pheno.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote = FALSE)


###########
#Run the HES PheWAS
library(data.table)
library(GenABEL)
library(plyr)
library(PheWAS)

pheno<-fread("/HES_icd10_pheno.txt",header=TRUE, stringsAsFactors=FALSE)
geno<-fread("/UKB.GLGC.PRS.txt",header=TRUE, stringsAsFactors=FALSE)
covar<-fread("/UKB.covar.txt",header=TRUE, stringsAsFactors=FALSE)

data<-merge(pheno,geno,by.x="id",by.y="IID")
data<-merge(data,covar,by.x="id",by.y="eid")

data<-subset(data,genetic_ethnicity=="1")

data$n.PRS.LDL <- rntransform(data$PRS.LDL)
data$n.PRS.HDL <- rntransform(data$PRS.HDL)
data$n.PRS.logTG <- rntransform(data$PRS.logTG)
data$n.PRS.TC <- rntransform(data$PRS.TC)
data$n.PRS.nonHDL <- rntransform(data$PRS.nonHDL)


set.seed(1)

results.LDL<-phewas(data=data,phenotypes=names(data[,c(2:1581)]),predictors=names(data[,c(1603)]),covariates=names(data[,c(1587,1590,1592:1602)]), cores=10,additive.genotypes = FALSE,significance.threshold="bonferroni")
results_d=addPhecodeInfo(results.LDL)
write.table(results_d,file="/PheWAS_results_HES_LDL.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote = FALSE)
png(file="/PheWAS.HES.LDL.png",width = 14, height = 8, units = 'in', res = 300)
phewasManhattan(results.LDL,title="UKB EUR HES ICD10 - LDL PRS PheWAS",OR.direction=TRUE,annotate.size=3,point.size=1,y.axis.interval=50)
dev.off()


#Run the Lab PheWAS
library(data.table)
library(GenABEL)
library(plyr)
library(PheWAS)

pheno<-fread("/UKB.lab.pheno.txt",header=TRUE, stringsAsFactors=FALSE,colClasses=list(numeric=2:98))
geno<-fread("/UKB.GLGC.PRS.txt",header=TRUE, stringsAsFactors=FALSE)
covar<-fread("/UKB.covar.txt",header=TRUE, stringsAsFactors=FALSE)
description<-fread("/UKB_lab_description.txt",header=TRUE, stringsAsFactors=FALSE)

pheno$LDL_cor<-ifelse(pheno$Statins=="1",pheno$LDL_direct/0.7,pheno$LDL_direct)
pheno$TC_cor<-ifelse(pheno$Statins!="1",pheno$Cholesterol/0.8,pheno$Cholesterol)

data<-merge(pheno,geno,by.x="eid",by.y="IID")
data<-merge(data,covar,by="eid")

data<-subset(data,genetic_ethnicity=="1")

data$n.PRS.LDL <- rntransform(data$PRS.LDL)
data$n.PRS.HDL <- rntransform(data$PRS.HDL)
data$n.PRS.logTG <- rntransform(data$PRS.logTG)
data$n.PRS.TC <- rntransform(data$PRS.TC)
data$n.PRS.nonHDL <- rntransform(data$PRS.nonHDL)


set.seed(1)

results.LDL<-phewas(data=data,phenotypes=names(data[,c(2:4,7:76,78:85,87:102,104:105)]),predictors=names(data[,c(127)]),covariates=names(data[,c(111,114,116:126)]), cores=1,additive.genotypes = FALSE,significance.threshold="bonferroni")
results.LDL<-merge(description,results.LDL,by="phenotype",all.y=TRUE)
write.table(results.LDL,file="/PheWAS_results_Lab_LDL.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote = FALSE)
png(file="/PheWAS.Lab.LDL.png",width = 14, height = 8, units = 'in', res = 300)
phenotypeManhattan(results.LDL,title="UKB EUR Lab - LDL PRS PheWAS",annotate.size=3,point.size=1,y.axis.interval=50,x.group.labels=T,annotate.phenotype=T)
dev.off()


####MVP
library(plyr)
library(PheWAS)
library(tidyverse)

lipid<-c("LDL","HDL","TC","logTG","nonHDL")

for (i in 1:length(lipid)) {

mvp<-read.table(paste0("/MVP_EUR_",lipid[i],"_PheWAS.txt"),header=TRUE, stringsAsFactors=FALSE, as.is=TRUE,sep="\t",fill=TRUE,quote="\"")

mvp$beta<-as.numeric(mvp$beta)
mvp$OR<-exp(mvp$beta)
names(mvp)[names(mvp)=="con"] <- "phenotype"
names(mvp)[names(mvp)=="pval"] <- "p"
mvp$phenotype<-as.character(mvp$phenotype)
mvp$p<-as.numeric(mvp$p)
mvp<-subset(mvp,grp!="Lab")


mvp$phenotype<-as.numeric(mvp$phenotype)

ukb<-read.table(paste0("/PheWAS_results_HES_",lipid[i],".txt"),header=TRUE, stringsAsFactors=FALSE, as.is=TRUE,sep="\t",quote="\"")
ukb$phenotype<-as.character(ukb$phenotype)

mvp$bonferroni<-ifelse(mvp$p <= 0.05/(length(which(mvp$p!="NA"))),"TRUE","FALSE")

combine<-merge(ukb,mvp,by="phenotype",all.x=TRUE,all.y=TRUE)
combine[is.na(combine)] <- "NA"
combine$description<-ifelse(combine$description=="NA",combine$info,combine$description)
combine$group<-ifelse(combine$group=="NA",combine$grp,combine$group)
combine<-combine[,c("phenotype","type","description","group","beta.x","SE","OR.x","p.x","n_total","n_cases","n_controls","bonferroni.x","beta.y","se","OR.y","p.y","n","nca","nco","bonferroni.y")]
colnames(combine)<-c("phenotype","type","description","group","beta.ukb","SE.ukb","OR.ukb","p.ukb","n_total.ukb","n_cases.ukb","n_controls.ukb","bonferroni.ukb","beta.mvp","SE.mvp","OR.mvp","p.mvp","n_total.mvp","n_cases.mvp","n_controls.mvp","bonferroni.mvp")

combine$beta.ukb<-as.numeric(combine$beta.ukb)
combine$beta.mvp<-as.numeric(combine$beta.mvp)
combine$color<-ifelse(combine$beta.ukb>0 & combine$beta.mvp>0 | combine$beta.ukb<0 & combine$beta.mvp<0,"black","red")

png(file=paste0("/PheWAS.",lipid[i],".PRS.EUR.UKBvsMVP.BETAS.ALL.png"),width = 8, height = 8, units = 'in', res = 300)
plot(combine$beta.ukb,combine$beta.mvp,pch=20,col=combine$color,cex =0.4,main=c(paste0("PheWAS BETAS ALL ",lipid[i],"-PRS EUR")),xlab=c("UKB betas"),ylab=c("MVP betas"),abline(0,1,col="darkblue"))
dev.off()

select<-subset(combine,bonferroni.ukb=="TRUE" | bonferroni.mvp=="TRUE")
png(file=paste0("/GLGC_PheWAS/PheWAS.",lipid[i],".PRS.EUR.UKBvsMVP.BETAS.BONF.png"),width = 8, height = 8, units = 'in', res = 300)
plot(select$beta.ukb,select$beta.mvp,pch=20,col=select$color,cex =0.6,main=c(paste0("PheWAS BETAS BONFERRONI ",lipid[i],"-PRS EUR")),xlab=c("UKB betas"),ylab=c("MVP betas"),abline(0,1,col="darkblue"))
dev.off()
}

library(data.table)
library(GenABEL)
library(plyr)
library(PheWAS)
library(dplyr) 

lipid<-c("LDL", "HDL", "TC", "nonHDL", "logTG")

for (i in 1:length(lipid)) {

ukb1<-read.table(paste0("/PheWAS_results_HES_",lipid[i],".txt"),header=TRUE, stringsAsFactors=FALSE, as.is=TRUE,sep="\t",quote="\"")
ukb1<-subset(ukb1,n_cases>=100)
ukb1<-ukb1[,c(1,2,4:11,17,18)]
colnames(ukb1)[1] <- "HES_code"
colnames(ukb1)[11] <- "phenotype"

ukb2<-read.table(paste0("/PheWAS_results_Lab_",lipid[i],".txt"),header=TRUE, stringsAsFactors=FALSE, as.is=TRUE,sep="\t",fill=TRUE,quote="\"")
ukb2<-subset(ukb2,n_total>=100)
ukb2<-ukb2[,c(1,2,5,7:14)]
ukb2$OR<- exp(ukb2$beta)
colnames(ukb2)[2] <- "HES_code"
ukb2$group<-c("Biomarkers & Measurements")
ukb1<-ukb1[,c(11,1:10,12)]

ukb<-rbind(ukb1,ukb2)
bonf<-nrow(ukb)

ukb$bonferroni<-ifelse(ukb$p > 0.05/bonf, "FALSE", "TRUE")
ukb$groupnum<-as.character(ukb$phenotype)

ukb<-ukb %>% mutate(groupnum = as.integer(factor(group, levels = unique(group))) )
ukb<-as.data.frame(ukb)

ukb$groupnum<-as.factor(ukb$groupnum)
png(file=paste0("/PheWAS.ALL.UKB.EUR.",lipid[i],".PRS.png"),width = 14, height = 8, units = 'in', res = 300)
phenotypeManhattan(ukb,title=paste0("PheWAS ",lipid[i]," PRS - UKBB EUR"),annotate.size=3,point.size=1,y.axis.interval=50,x.group.labels=T,annotate.phenotype=T,color.palette=topo.colors(18),OR.direction=TRUE,suggestive.line=0.05/bonf,significant.line=5e-08,annotate.level=5e-08)
dev.off()

write.table(ukb,file=paste0("/PheWAS.ALL.UKB.EUR.",lipid[i],".PRS.qced.txt"), sep="\t", row.names=FALSE, col.names=TRUE, quote = FALSE)

}

for (i in 1:length(lipid)) {

mvp<-read.table(paste0("/MVP_EUR_",lipid[i],"_PheWAS.txt"),header=TRUE, stringsAsFactors=FALSE, as.is=TRUE,sep="\t",fill=TRUE,quote="\"")
mvp[is.na(mvp)] <- "NA"
colnames(mvp)<-c("HES_code", "phenotype","group","p","beta","SE","n_total","n_controls","n_cases")
mvp$beta<-as.numeric(mvp$beta)
mvp$OR<-exp(mvp$beta)
mvp$p<-as.numeric(mvp$p)
mvp<-subset(mvp,p!="NA")
mvp$n_cases<-as.numeric(mvp$n_cases)
mvp$type<-ifelse(mvp$group=="Lab" | mvp$group=="Lab at enrollment" | mvp$group=="Lifestyle" | mvp$group=="Vital at enrollment","linear","logistic")
mvp$predictor<-paste0("n.PRS.",lipid[i])
mvp1<-subset(mvp,type=="logistic")
mvp1<-subset(mvp1,n_cases>=100)
mvp2<-subset(mvp,type=="linear")
mvp2<-subset(mvp2,n_total>=100)
mvp<-rbind(mvp1,mvp2)
bonf<-nrow(mvp)

mvp$bonferroni<-ifelse(mvp$p > 0.05/bonf, "FALSE", "TRUE")
mvp$groupnum<-as.character(mvp$phenotype)

mvp<-mvp %>% mutate(groupnum = as.integer(factor(group, levels = unique(group))) )
mvp<-as.data.frame(mvp)

mvp$groupnum<-as.factor(mvp$groupnum)

mvp<-mvp[,c("phenotype","HES_code","predictor","beta","SE","OR","p","type","n_total","n_cases","n_controls","group","bonferroni","groupnum")] 

png(file=paste0("/PheWAS.ALL.MVP.EUR.",lipid[i],".PRS.png"),width = 14, height = 8, units = 'in', res = 300)
phenotypeManhattan(mvp,title=paste0("PheWAS ",lipid[i]," PRS - MVP EUR"),annotate.size=3,point.size=1,y.axis.interval=50,x.group.labels=T,annotate.phenotype=T,color.palette=topo.colors(21),OR.direction=TRUE,suggestive.line=0.05/bonf,significant.line=5e-08,annotate.level=5e-08)
dev.off()

write.table(mvp,file=paste0("/PheWAS.ALL.MVP.EUR.",lipid[i],".PRS.qced.txt"), sep="\t", row.names=FALSE, col.names=TRUE, quote = FALSE)

}


bridge<-read.table("/lab_ukb_mvp_bridge.txt",header=TRUE, stringsAsFactors=FALSE, as.is=TRUE)

for (i in 1:length(lipid)) {
ukb<-read.table(paste0("/PheWAS.ALL.UKB.EUR.",lipid[i],".PRS.qced.txt"),header=TRUE, stringsAsFactors=FALSE, as.is=TRUE,sep="\t",fill=TRUE,quote="\"")

mvp<-read.table(paste0("/PheWAS.ALL.MVP.EUR.",lipid[i],".PRS.qced.txt"),header=TRUE, stringsAsFactors=FALSE, as.is=TRUE,sep="\t",fill=TRUE,quote="\"")
mvp$group<-ifelse(mvp$group=="Lab" | mvp$group=="Lab at enrollment" | mvp$group=="Lifestyle" | mvp$group=="Vital at enrollment","Biomarkers & Measurements",mvp$group)
mvp<-merge(mvp,bridge,by.x="phenotype",by.y="con",all.x=TRUE)
mvp[is.na(mvp)] <- "NA"
mvp$phenotype<-ifelse(mvp$group!="Biomarkers & Measurements",mvp$phenotype,ifelse(mvp$group=="Biomarkers & Measurements" & mvp$phenotype.y=="NA",mvp$phenotype,mvp$phenotype.y))
mvp$HES_code<-ifelse(mvp$group!="Biomarkers & Measurements",mvp$HES_code,ifelse(mvp$group=="Biomarkers & Measurements" & mvp$phenotype.y=="NA",mvp$phenotype,mvp$phenotype.y))

ukb$study<-"ukb"
mvp$study<-"mvp"
results<-rbind(ukb,mvp[,c(1:14,16)])
results$snp<-"PRS"
results$adjustment<-"NA"
results$n_cases<-as.numeric(results$n_cases)
results$n_controls<-as.numeric(results$n_controls)
results$n_total<-as.numeric(results$n_total)

results.meta=phewasMeta(results, fixed=TRUE, keep.both=TRUE)

meta<-as.data.frame(results.meta)
meta<-subset(meta,k_studies==2)
bonf<-nrow(meta)

meta$OR.fixed<-ifelse(meta$type=="linear",exp(meta$beta.fixed),meta$OR.fixed)
meta$OR.random<-ifelse(meta$type=="linear",exp(meta$beta.random),meta$OR.random)
meta$p<-ifelse(meta$Q.p > 0.05/bonf,meta$p.fixed, meta$p.random)
meta$beta<-ifelse(meta$Q.p > 0.05/bonf,meta$beta.fixed, meta$beta.random)
meta$SE<-ifelse(meta$Q.p > 0.05/bonf,meta$SE.fixed, meta$SE.random)
meta$OR<-ifelse(meta$Q.p > 0.05/bonf,meta$OR.fixed, meta$OR.random)
meta$bonferroni<-ifelse(meta$p > 0.05/bonf, "FALSE", "TRUE")

meta<-merge(meta,ukb[,c("phenotype","group","HES_code","beta","SE","p","n_total","n_cases","n_controls")],by="phenotype",all.x=TRUE)
meta<-merge(meta,mvp[,c("phenotype","beta","SE","p","n_total","n_cases","n_controls")],by="phenotype",all.x=TRUE)
meta$phenotype<-gsub("_cor", "_cholesterol", meta$phenotype)
colnames(meta)<-c("phenotype","snp","adjustment","beta","OR","SE","p","type","n_total","n_cases","n_controls","HWE_p.min","allele_freq","n_no_snp","k_studies","tau2","I2.percent","Q","Q.df","Q.p","beta.fixed","OR.fixed","SE.fixed","p.fixed","beta.random","OR.random","SE.random","p.random","bonferroni","group","HES_code","beta.ukb","SE.ukb","p.ukb","n_total.ukb","n_cases.ukb","n_controls.ukb","beta.mvp","SE.mvp","p.mvp","n_total.mvp","n_cases.mvp","n_controls.mvp")
meta[is.na(meta)] <- "NA"
meta$n_total<-ifelse(meta$n_total=="NA",meta$n_total.ukb + meta$n_total.mvp, meta$n_total)
meta$n_total<-as.numeric(meta$n_total)
meta$n_cases<-as.numeric(meta$n_cases)
meta1<-subset(meta,type=="logistic")
meta1<-subset(meta1,n_cases>=500)
meta1 <- meta1[order(meta1$HES_code),]
meta2<-subset(meta,type=="linear")
meta2<-subset(meta2,n_total>=500)
meta<-rbind(meta1,meta2)
meta$groupnum<-as.character(meta$phenotype)
meta<-meta %>% mutate(groupnum = as.integer(factor(group, levels = unique(group))) )
meta<-as.data.frame(meta)

meta$groupnum<-as.factor(meta$groupnum)
png(file=paste0("/PheWAS.UKB.MVP.EUR.META.ALL.",lipid[i],".PRS.v2.png"),width = 14, height = 8, units = 'in', res = 300)
phenotypeManhattan(meta,title=paste0("PheWAS ",lipid[i]," PRS - Meta-analysis UKB-MVP (EUR)"),annotate.size=3,point.size=1,y.axis.interval=50,x.group.labels=T,annotate.phenotype=T,color.palette=topo.colors(17),OR.direction=TRUE,suggestive.line=0.05/bonf,significant.line=5e-08,annotate.level=0.05/bonf)
dev.off()


}


#############################################
#Run PheWAS for index variants
#############################################

#Run the HES PheWAS
library(data.table)
library(plyr)
library(PheWAS)

pheno<-fread("/HES_icd10_pheno.txt",header=TRUE, stringsAsFactors=FALSE)
geno<-fread("/UKB.GLGC.TRANSETHNIC.INDEX.txt",header=TRUE, stringsAsFactors=FALSE)
covar<-fread("/UKB.covar.txt",header=TRUE, stringsAsFactors=FALSE)

data<-merge(pheno,covar,by.x="id",by.y="eid")
data<-merge(data,geno,by.x="id",by.y="IID")

data<-subset(data,genetic_ethnicity=="1")


set.seed(1)

results<-phewas(data=data,phenotypes=names(data[,c(2:1581)]),genotypes=names(data[,c(1599:3417)]),covariates=names(data[,c(1582,1585,1587:1597)]), cores=1,additive.genotypes = TRUE,significance.threshold="bonferroni",min.records=500)
results_d=addPhecodeInfo(results)
write.table(results_d,file="/PheWAS_results_HES_Index_Transethnic.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote = FALSE)


pheno<-fread("/UKB.lab.pheno.txt",header=TRUE, stringsAsFactors=FALSE,colClasses=list(numeric=2:98))
description<-fread("/UKB_lab_description.txt",header=TRUE, stringsAsFactors=FALSE)

pheno$LDL_cor<-ifelse(pheno$Statins=="1",pheno$LDL_direct/0.7,pheno$LDL_direct)
pheno$TC_cor<-ifelse(pheno$Statins!="1",pheno$Cholesterol/0.8,pheno$Cholesterol)

data<-merge(pheno,geno,by.x="eid",by.y="IID")
data<-merge(data,covar,by="eid")

data<-subset(data,genetic_ethnicity=="1")


set.seed(1)

results<-phewas(data=data,phenotypes=names(data[,c(2:4,7:76,78:85,87:102,104:105)]),genotypes=names(data[,c(107:1925)]),covariates=names(data[,c(1906,1929,1931:1941)]), cores=1,additive.genotypes = TRUE,significance.threshold="bonferroni",min.records=500)
results<-merge(description,results,by="phenotype",all.y=TRUE)
write.table(results,file="/PheWAS_results_LAB_Index_Transethnic.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote = FALSE)

rm(list=ls(all=TRUE))


results<-read.table("/PheWAS_results_HES_Index_Transethnic.txt", header=TRUE, as.is=TRUE,stringsAsFactors=FALSE,sep="\t",fill=TRUE)
phewas<-read.table("/PheWAS_results_LAB_Index_Transethnic.txt", header=TRUE, as.is=TRUE,stringsAsFactors=FALSE,sep="\t",fill=TRUE)

colnames(results)[1]<-"code"
colnames(results)[17]<-"phenotype"

phewas$group<-ifelse(phewas$group=="Blood Count","Blood Biomarkers",phewas$group)
phewas$code<-"lab"
phewas<-phewas[,c(20,5:19,1,2)]

combine<-rbind(results,phewas)
combine$snp<-gsub("`","",combine$snp)
combine$snp<-gsub(":","_",combine$snp)
temp<-colsplit(combine$snp,"_",c("chr","pos","A1","A2","refAllele"))
temp$snp<-paste(temp$chr,":",temp$pos,"_",temp$A1,"_",temp$A2, sep="")
combine<-cbind(combine,temp[,c(5,6)])
combine<-combine[,c(1,20,19,3:18)]

combine<-subset(combine,p!="NA")

write.table(combine,file="/PheWAS_results_UKB_EUR_Index_Transethnic_GLGC.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote = FALSE)

combine<-subset(combine,p<=5e-08)

write.table(combine,file="/PheWAS_results_UKB_EUR_Index_Transethnic_GLGC_GWS.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote = FALSE)


mvp$var<-gsub("-","_",mvp$var)
temp<-colsplit(mvp$var,"_",c("chr","pos","A1","A2","refAllele"))
temp$snpid<-paste0(temp$chr,":",temp$pos,"_",temp$A1,"_",temp$A2)
temp$snp<-paste0(temp$chr,":",temp$pos)
mvp<-cbind(mvp,temp[,c(7,5,6)])

colnames(ukb)[2]<-"snpid"
temp<-colsplit(ukb$snpid,"_",c("snp","A1","A2"))
ukb<-cbind(ukb,temp[,c(1)])
colnames(ukb)[20]<-"snp"

temp<-ukb[,c(20,3)]
temp<-temp %>% distinct(snp, .keep_all = TRUE)
mvp2<-merge(mvp,temp,by="snp",all.x=TRUE)
mvp2$refAllele.x<-as.character(mvp2$refAllele.x)
mvp2$refAllele.y<-as.character(mvp2$refAllele.y)
mvp2[is.na(mvp2)] <- "NA"
mvp2$beta<-as.numeric(mvp2$beta)
mvp2$beta_new<-ifelse(mvp2$refAllele.y=="NA",mvp2$beta,ifelse(mvp2$refAllele.x==mvp2$refAllele.y,mvp2$beta, - mvp2$beta))
mvp2$refAllele<-ifelse(mvp2$refAllele.y=="NA",mvp2$refAllele.x,ifelse(mvp2$refAllele.x==mvp2$refAllele.y,mvp2$refAllele.x, mvp2$refAllele.y))

bridge1<-read.table("/lab_ukb_mvp_bridge.txt",header=TRUE, stringsAsFactors=FALSE, as.is=TRUE)
bridge2<-read.table("/UKB_lab_description.txt",header=TRUE, stringsAsFactors=FALSE, as.is=TRUE,sep="\t",fill=TRUE,quote="\"")
bridge<-merge(bridge1,bridge2,by="phenotype")
bridge$group<-ifelse(bridge$group=="Blood Count","Blood Biomarkers",bridge$group)

mvp2<-merge(mvp2,bridge,by="con",all.x=TRUE)
mvp2<-merge(mvp2,mvp_old[,c(1,3)],by="con",all.x=TRUE)

mvp2$phenotype<-gsub("LDL_cor","LDL",mvp2$phenotype)
mvp2$phenotype<-gsub("TC_cor","TC",mvp2$phenotype)

mvp2$phenotype<-ifelse(mvp2$grp.x=="Lab at enrollment (statin adjustment)" & mvp2$info=="LDL-C","LDL_cor",ifelse(mvp2$grp.x=="Lab at enrollment (statin adjustment)" & mvp2$info=="Total cholesterol","TC_cor",mvp2$phenotype))

mvp2$group<-ifelse(mvp2$grp.x=="Phecode", mvp2$grp.y, ifelse(mvp2$grp.x=="Lab", mvp2$group, ifelse(mvp2$phenotype=="LDL_cor" | mvp2$phenotype=="TC_cor", "Blood Biomarkers", ifelse(mvp2$phenotype=="BMI","Physical Measures","NA"))))

mvp2[is.na(mvp2)] <- "NA"

mvp2$phenotype<-ifelse(mvp2$grp.x=="Phecode",mvp2$info,ifelse(mvp2$grp.x!="Phecode" & mvp2$phenotype=="NA",mvp2$info,mvp2$phenotype))

mvp2$type<-ifelse(mvp2$grp.x=="Phecode","logistic","linear")
ukb$study<-"ukb"
mvp2$study<-"mvp"

mvp3<-mvp2[,c(17,18,2,7,8,6,9,11,10,22,23)]
ukb2<-ukb[,c(18,19,20,5,6,8,10,11,12,9,21)]
colnames(mvp3)<-c("phenotype", "group", "snp", "beta", "SE", "p", "n_total", "n_cases", "n_controls", "type","study")
colnames(ukb2)<-c("phenotype", "group", "snp", "beta", "SE", "p", "n_total", "n_cases", "n_controls", "type","study")


merged<-merge(ukb2,mvp3,by=c("phenotype","group","snp"))
mvp3<-merged[,c(1:3,12:19)]
ukb2<-merged[,c(1:11)]
colnames(mvp3)<-c("phenotype", "group", "snp", "beta", "SE", "p", "n_total", "n_cases", "n_controls", "type","study")
colnames(ukb2)<-c("phenotype", "group", "snp", "beta", "SE", "p", "n_total", "n_cases", "n_controls", "type","study")

results<-rbind(ukb2,mvp3)

results$adjustment<-"NA"
results$n_cases<-as.numeric(results$n_cases)
results$n_controls<-as.numeric(results$n_controls)
results$n_total<-as.numeric(results$n_total)
results$beta<-as.numeric(results$beta)
results$SE<-as.numeric(results$SE)
results$p<-as.numeric(results$p)

results.meta=phewasMeta(results, fixed=TRUE, keep.both=TRUE)

meta<-as.data.frame(results.meta)


meta<-read.table("/PheWAS.ALL.META.UKB.MVP.EUR.TRANS.INDEX.txt", header=TRUE, as.is=TRUE,stringsAsFactors=FALSE,sep="\t",fill=TRUE)


meta<-subset(meta,k_studies==2)

bonf<-nrow(meta)
meta$OR.fixed<-ifelse(meta$type=="linear",exp(meta$beta.fixed),meta$OR.fixed)
meta$OR.random<-ifelse(meta$type=="linear",exp(meta$beta.random),meta$OR.random)
meta$p<-ifelse(meta$Q.p > 0.05/bonf,meta$p.fixed, meta$p.random)
meta$beta<-ifelse(meta$Q.p > 0.05/bonf,meta$beta.fixed, meta$beta.random)
meta$SE<-ifelse(meta$Q.p > 0.05/bonf,meta$SE.fixed, meta$SE.random)
meta$OR<-ifelse(meta$Q.p > 0.05/bonf,meta$OR.fixed, meta$OR.random)
meta$bonferroni<-ifelse(meta$p > 0.05/bonf, "FALSE", "TRUE")

ukb2<-ukb[,c(18,19,20,5,6,8,10,11,12,9,21,3)]
colnames(ukb2)<-c("phenotype", "group", "snp", "beta", "SE", "p", "n_total", "n_cases", "n_controls", "type","study","reference_allele")

meta<-merge(meta,ukb2[,c("phenotype","group","snp","beta","SE","p","n_total","n_cases","n_controls","reference_allele")],by=c("phenotype","snp"),all.x=TRUE)
meta<-merge(meta,mvp3[,c("phenotype","snp","beta","SE","p","n_total","n_cases","n_controls")],by=c("phenotype","snp"),all.x=TRUE)
colnames(meta)<-c("phenotype","snp","adjustment","beta","OR","SE","p","type","n_total","n_cases","n_controls","HWE_p.min","allele_freq","n_no_snp","k_studies","tau2","I2.percent","Q","Q.df","Q.p","beta.fixed","OR.fixed","SE.fixed","p.fixed","beta.random","OR.random","SE.random","p.random","bonferroni","group","beta.ukb","SE.ukb","p.ukb","n_total.ukb","n_cases.ukb","n_controls.ukb","reference_allele","beta.mvp","SE.mvp","p.mvp","n_total.mvp","n_cases.mvp","n_controls.mvp")
meta[is.na(meta)] <- "NA"
meta$n_total<-ifelse(meta$n_total=="NA",meta$n_total.ukb + meta$n_total.mvp, meta$n_total)
meta$n_total<-as.numeric(meta$n_total)
meta$n_cases<-as.numeric(meta$n_cases)
meta1<-subset(meta,type=="logistic")
meta1<-subset(meta1,n_cases>=500)
meta2<-subset(meta,type=="linear")
meta2<-subset(meta2,n_total>=500)
meta<-rbind(meta1,meta2)

index<-read.table("/Trans_Ancestry_Gene_Prioritizations_Best_Guess.txt",header=TRUE, stringsAsFactors=FALSE, as.is=TRUE,sep="\t",fill=TRUE,quote="\"")
index$V1<-paste0(index$chr,":",index$pos)
index<-index[,c(1,2,7,8)]
index<-index[order(index$V1,index$trait),]
index1<- index %>% group_by(rsid_dbSNP150,V1) %>% mutate(trait = paste(trait,collapse =","))

index1<-as.data.frame(index1)

index1<-distinct(index1)


combine<-merge(meta,index1,by.x="snp",by.y="V1",all.x=TRUE,sort=FALSE)
combine[is.na(combine)] <- "NA"
  
write.table(combine,file="/PheWAS.ALL.META.UKB.MVP.EUR.TRANS.INDEX.WITH_GENES.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote = FALSE)

rm(list=ls(all=TRUE))
