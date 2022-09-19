R

library(tidyverse)
library(data.table)
library(reshape)

index_genes<-read.table("/trans_ethnic_index_variants_biomart_genes.txt", header=TRUE, stringsAsFactors=FALSE, as.is=TRUE)
mouse<-read.table("/LIPID_genes_lipid_annotated.tsv.gz", header=TRUE, stringsAsFactors=FALSE, sep="\t", as.is=TRUE,quote="\"",fill=TRUE,comment.char="",na.strings=c("","NA"))
pheno<-read.table("/mp-names.tsv.gz", header=TRUE, stringsAsFactors=FALSE, sep="\t", as.is=TRUE,quote="\"")
priority<-read.table("/phenotypes_gradesv2.txt", header=TRUE, stringsAsFactors=FALSE, as.is=TRUE,quote="\"",sep="\t")
index_variants<-read.table("/trans_ethnic_index_variants.txt",as.is=TRUE,header=TRUE, stringsAsFactors=FALSE)

mouse<-subset(mouse,num_phenotypes>0)


temp1<-subset(mouse,phenotypes_growth_body>0)
temp2<-subset(mouse,phenotypes_homeostasis>0)
temp3<-subset(mouse,phenotypes_lipid_homeostasis>0)
temp4<-subset(mouse,phenotypes_cholesterol_homeostasis>0)
temp5<-subset(mouse,phenotypes_metabolism>0)

temp1<-separate_rows(mouse, phenotypes_growth_body,sep=";", convert = TRUE)
temp2<-separate_rows(mouse, phenotypes_homeostasis,sep=";", convert = TRUE)
temp3<-separate_rows(mouse, phenotypes_lipid_homeostasis,sep=";", convert = TRUE)
temp4<-separate_rows(mouse, phenotypes_cholesterol_homeostasis,sep=";", convert = TRUE)
temp5<-separate_rows(mouse, phenotypes_metabolism,sep=";", convert = TRUE)

temp1<-as.data.frame(temp1[,c(1:8,11)])
temp1$phenotype_category<-"phenotypes_growth_body"
temp2<-as.data.frame(temp2[,c(1:8,13)])
temp2$phenotype_category<-"phenotypes_homeostasis"
temp3<-as.data.frame(temp3[,c(1:8,15)])
temp3$phenotype_category<-"phenotypes_lipid_homeostasis"
temp4<-as.data.frame(temp4[,c(1:8,17)])
temp4$phenotype_category<-"phenotypes_cholesterol_homeostasis"
temp5<-as.data.frame(temp5[,c(1:8,19)])
temp5$phenotype_category<-"phenotypes_metabolism"

names(temp1) <- gsub("_growth_body", "", names(temp1))
names(temp2) <- gsub("_homeostasis", "", names(temp2))
names(temp3) <- gsub("_lipid_homeostasis", "", names(temp3))
names(temp4) <- gsub("_cholesterol_homeostasis", "", names(temp4))
names(temp5) <- gsub("_metabolism", "", names(temp5))

mouse<-rbind(temp1,temp2,temp3,temp4,temp5)

query<-merge(mouse,index_genes,by="GeneID")
query<-merge(query,pheno,by.x="phenotypes",by.y="id")
query<-merge(query,priority,by.x="name",by.y="phenotype_name")

query[is.na(query)] <- "NA"

select<-query %>% group_by(rsid) %>% top_n(1, level_of_relevance)
select<-as.data.frame(select)

select2<-select[,c(1:9,12:14)]
select2<-distinct(select2)

select2<-merge(select2,index_variants,by="rsid")

select3<-select2 %>% 
     group_by(rsid,Symbol.y,trait) %>% 
     mutate(name2 = paste0(name, collapse = ";"), phenotypes2= paste0(phenotypes, collapse = ";"),Model_Description= paste0(model_description, collapse = ";")) 

select3<-as.data.frame(select3)

select4<-select3[,c(1,13:15,17,18,4:9,19,12)]	

select4<-distinct(select4)

select5<-select4 %>% 
     group_by(rsid,Symbol.x,trait) %>% 
     mutate(Model_ID= paste0(model_id, collapse = ";"))

select5<-as.data.frame(select5)

select5<-select5[,c(1:10,15,12:14)]	

select5<-distinct(select5)

final<-select5 %>% 
     group_by(rsid,Symbol.x,trait) %>% 
     mutate(Gene_ID= paste0(GeneID, collapse = ";"))

final<-as.data.frame(final)

final<-final[,c(1:6,15,8:14)]	

final<-distinct(final)
colnames(final)<-c("rsid","trait","chr","pos","name","phenotypes","Gene_ID","Symbol","HGNC_id","MGI_id","Model_ID","type_of_gene","model_description","priority_level")

final$Symbol<-sQuote(final$Symbol)

write.table(final,file="/trans_ethnic_mouse_prioritised_genes.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote = FALSE)


rm(list=ls(all=TRUE))
