rm(list=ls())
library(readxl)
library(Map2NCBI)
library(tidyverse)
trans_index <- read_excel("GLGC Gene prioritization.xlsx", 
                          sheet = "Trans-ancestry index variants")
specific_index <- read_excel("GLGC Gene prioritization.xlsx", 
                             sheet = "Ancestry-specific index variant")
###########
genes <- read_excel("Ensembl_biomart_genes_build37.xlsx")
genes <- data.frame(cbind(genes$`Gene name`,genes$`Chromosome/scaffold name`,genes$`Gene start (bp)`,genes$`Gene end (bp)`))
colnames(genes) <- c("FeatureName","chromosome","start","end")
genes <- genes[genes$chromosome %in% c(1:22),]
genes$chromosome <- as.numeric(genes$chromosome)
genes$start <- as.numeric(genes$start)
genes$end <- as.numeric(genes$end)
###########
trans <- data.frame(cbind(trans_index$rsid_dbSNP150,trans_index$chr,trans_index$`pos (build 37)`))
colnames(trans) <- c("Marker","chromosome","position")
trans$chromosome <- as.numeric(trans$chromosome)
trans$position <- as.numeric(trans$position)
trans <- unique(trans)

trans_closest <- MapMarkers(genes, trans, nAut=22, other = F, savefiles = T,destfile = getwd())
trans_closest$Distance <- as.numeric(trans_closest$Distance)
trans_closest <- trans_closest[trans_closest$Distance <= 500000,]
trans_closest <- trans_closest[,c(1:4,8)]
colnames(trans_closest) <- c("rsid_dbSNP150","chr","pos (build 37)","closest gene","distance")
###################
specific <- data.frame(cbind(specific_index$dbsnp_150,specific_index$chr,specific_index$`pos (build 37)`))
colnames(specific) <- c("Marker","chromosome","position")
specific$chromosome <- as.numeric(specific$chromosome)
specific$position <- as.numeric(specific$position)
specific <- unique(specific)

specific_closest <- MapMarkers(genes, specific, nAut=22, other = F, savefiles = T,destfile = getwd())
specific_closest$Distance <- as.numeric(specific_closest$Distance)
specific_closest <- specific_closest[specific_closest$Distance <= 500000,]
specific_closest <- specific_closest[,c(1:4,8)]
colnames(specific_closest) <- c("dbsnp_150","chr","pos (build 37)","closest gene","distance")
###################

write.table(trans_closest,file = "Trans-ancestry index variants closest gene.txt",row.names = F,col.names = T, quote = F)
write.table(specific_closest,file = 'Ancestry-specific index variant closest gene.txt',row.names = F,col.names = T,quote = F)

library("xlsx")
write.xlsx(fulldata,"Ancestry-specific index variant.xlsx",row.names = F)


###
trans_index_bed <- trans_index %>% 
  select(chr,`pos (build 37)`,rsid_dbSNP150) %>% 
  mutate(end = `pos (build 37)`) %>% 
  select(chr,`pos (build 37)`,end, rsid_dbSNP150) %>% 
  arrange(chr,`pos (build 37)`) %>% distinct()
trans_index_bed$chr <- paste("chr", trans_index_bed$chr, sep="")

genes <- genes[genes$`Chromosome/scaffold name` %in% c(1:22),]
genes_all_bed <- genes %>% 
  select(`Chromosome/scaffold name`,`Gene start (bp)`,`Gene end (bp)`,`Gene name`,`Gene type`,`Gene stable ID`) %>% 
  arrange(`Chromosome/scaffold name`,`Gene start (bp)`)
genes_all_bed$`Chromosome/scaffold name` <- paste("chr", genes_all_bed$`Chromosome/scaffold name`, sep="")

write.table(trans_index_bed, "trans_index.bed", row.names = F,col.names = F,quote = F,sep = '\t')
write.table(genes_all_bed, "genes_all.bed", row.names = F,col.names = F,quote = F,sep = '\t')


test <- read.table(file = "test_1000.bed")
test <- distinct(test)
test <- test[test$V10 <= 500000,]
trans_index_genes_all <- test %>% group_by(V1,V2) %>% slice(1:5) %>% select(V4,V1,V2,V8,V6,V7,V9,V10)
colnames(trans_index_genes_all) <- c("rsid_dbSNP150","chr","pos (build 37)","closest gene","gene start","gene end","gene type","distance")

test <- test[test$V9 =="protein_coding",]
trans_index_genes_coding <- test %>% group_by(V1,V2) %>% select(V4,V1,V2,V8,V6,V7,V9,V10)
colnames(trans_index_genes_coding) <- c("rsid_dbSNP150","chr","pos (build 37)","closest_gene","gene_start","gene_end","gene_type","distance")

write.table(trans_index_genes_all,"trans_index_all_genes_top5.txt",row.names = F,col.names = T,quote = F)
write.table(trans_index_genes_coding,"trans_index_protein_coding_genes_all.txt",row.names = F,col.names = T,quote = F,sep = "\t")

write.csv(trans_index_genes_coding,"trans_index_protein_coding_genes_all.csv",row.names = F,col.names = T,quote = F)


trans_index_genes_coding <- trans_index_genes_coding %>%
  group_by(rsid_dbSNP150,chr,`pos (build 37)`) %>%
  mutate(all_closest_gene = paste(closest_gene, collapse = ";"))

trans_index_genes_coding <- trans_index_genes_coding[,c(1,9)]
trans_index_genes_coding <- trans_index_genes_coding[!duplicated(trans_index_genes_coding),]
write.table(trans_index_genes_coding,"trans_index_protein_coding_genes_all.txt",row.names = F,col.names = T,quote = F)

