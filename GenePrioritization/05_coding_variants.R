credible_set <-read.delim("99_percent_credible_set_variants_with_PP_updated_summary_by_transancestry_locus.txt",header = T)
credible_set <- credible_set[,c(1,2,3,11)]
credible_set <- credible_set[credible_set$anno_variants!="",]

seperate <- function(x){
  paste2(unique(unlist(lapply(strsplit(unlist(strsplit(unlist(lapply(strsplit(unlist(strsplit(x[4],",")),"anno:"),`[[`, 2)),"\\|")),"\\("),`[[`, 1))),collapse=";")}
credible_set$credible_set <-apply(credible_set, 1, seperate)
credible_set <- credible_set[,c(1,2,3,5)]
colnames(credible_set) <-c("trait", "chr","pos (build 37)","credible_set")

write.csv(credible_set,file = "credible_set.csv",quote = F,row.names = F)

chrall_hdl <- read.delim("chrALL_HDL_MR-MEGA.out.result_rsID_pval_lt_5e-8_missense_stop_splice_frameshift_anno_new.txt",header = T)
chrall_hdl <- chrall_hdl[chrall_hdl$LD_r2_with_nearest_index != "unable_to_lookup_LD",]
chrall_hdl <- chrall_hdl[complete.cases(chrall_hdl[ ,43]),]
chrall_hdl$LD_r2_with_nearest_index <- as.numeric(chrall_hdl$LD_r2_with_nearest_index)
chrall_hdl <- chrall_hdl[chrall_hdl$LD_r2_with_nearest_index > 0.8,]
chrall_hdl <- chrall_hdl[,c(40:41)]
chrall_hdl <- chrall_hdl[!duplicated(chrall_hdl), ]
chrall_hdl <- aggregate(SnpEff_ensembl_summary ~ closest_index_variant,data=chrall_hdl,paste,collapse=';')

seperate <- function(x){
  paste2(unique(unlist(lapply(strsplit(unlist(strsplit(unlist(strsplit(x[2],";")),"\\|")),"\\("), `[[`, 1))),collapse=";")}

chrall_hdl$SnpEff_ensembl_summary <-apply(chrall_hdl, 1, seperate)
chrall_hdl$trait <- "HDL"
chrall_hdl$trait_pos <- paste(chrall_hdl$trait, chrall_hdl$closest_index_variant, sep=":")
chrall_hdl <- chrall_hdl[,c(2,4)]



chrall_ldl <- read.delim("chrALL_LDL_MR-MEGA.out.result_rsID_pval_lt_5e-8_missense_stop_splice_frameshift_anno_new.txt",header = T)
chrall_ldl <- chrall_ldl[chrall_ldl$LD_r2_with_nearest_index != "unable_to_lookup_LD",]
chrall_ldl <- chrall_ldl[complete.cases(chrall_ldl[ ,43]),]
chrall_ldl$LD_r2_with_nearest_index <- as.numeric(chrall_ldl$LD_r2_with_nearest_index)
chrall_ldl <- chrall_ldl[chrall_ldl$LD_r2_with_nearest_index > 0.8,]
chrall_ldl <- chrall_ldl[,c(40:41)]
chrall_ldl <- chrall_ldl[!duplicated(chrall_ldl), ]
chrall_ldl <- aggregate(SnpEff_ensembl_summary ~ closest_index_variant,data=chrall_ldl,paste,collapse=';')
seperate <- function(x){
  paste2(unique(unlist(lapply(strsplit(unlist(strsplit(unlist(strsplit(x[2],";")),"\\|")),"\\("), `[[`, 1))),collapse=";")}
chrall_ldl$SnpEff_ensembl_summary <-apply(chrall_ldl, 1, seperate)
chrall_ldl$trait <- "LDL"
chrall_ldl$trait_pos <- paste(chrall_ldl$trait, chrall_ldl$closest_index_variant, sep=":")
chrall_ldl <- chrall_ldl[,c(2,4)]





chrall_logTG <- read.delim("chrALL_logTG_MR-MEGA.out.result_rsID_pval_lt_5e-8_missense_stop_splice_frameshift_anno_new.txt",header = T)
chrall_logTG <- chrall_logTG[chrall_logTG$LD_r2_with_nearest_index != "unable_to_lookup_LD",]
chrall_logTG <- chrall_logTG[complete.cases(chrall_logTG[ ,43]),]
chrall_logTG$LD_r2_with_nearest_index <- as.numeric(chrall_logTG$LD_r2_with_nearest_index)
chrall_logTG <- chrall_logTG[chrall_logTG$LD_r2_with_nearest_index > 0.8,]
chrall_logTG <- chrall_logTG[,c(40:41)]
chrall_logTG <- chrall_logTG[!duplicated(chrall_logTG), ]
chrall_logTG <- aggregate(SnpEff_ensembl_summary ~ closest_index_variant,data=chrall_logTG,paste,collapse=';')
seperate <- function(x){
  paste2(unique(unlist(lapply(strsplit(unlist(strsplit(unlist(strsplit(x[2],";")),"\\|")),"\\("), `[[`, 1))),collapse=";")}
chrall_logTG$SnpEff_ensembl_summary <-apply(chrall_logTG, 1, seperate)
chrall_logTG$trait <- "logTG"
chrall_logTG$trait_pos <- paste(chrall_logTG$trait, chrall_logTG$closest_index_variant, sep=":")
chrall_logTG <- chrall_logTG[,c(2,4)]


chrall_nonHDL<- read.delim("chrALL_nonHDL_MR-MEGA.out.result_rsID_pval_lt_5e-8_missense_stop_splice_frameshift_anno_new.txt",header = T)
chrall_nonHDL<- chrall_nonHDL[chrall_nonHDL$LD_r2_with_nearest_index != "unable_to_lookup_LD",]
chrall_nonHDL<- chrall_nonHDL[complete.cases(chrall_nonHDL[ ,43]),]
chrall_nonHDL$LD_r2_with_nearest_index <- as.numeric(chrall_nonHDL$LD_r2_with_nearest_index)
chrall_nonHDL<- chrall_nonHDL[chrall_nonHDL$LD_r2_with_nearest_index > 0.8,]
chrall_nonHDL<- chrall_nonHDL[,c(40:41)]
chrall_nonHDL<- chrall_nonHDL[!duplicated(chrall_nonHDL), ]
chrall_nonHDL<- aggregate(SnpEff_ensembl_summary ~ closest_index_variant,data=chrall_nonHDL,paste,collapse=';')
seperate <- function(x){
  paste2(unique(unlist(lapply(strsplit(unlist(strsplit(unlist(strsplit(x[2],";")),"\\|")),"\\("), `[[`, 1))),collapse=";")}
chrall_nonHDL$SnpEff_ensembl_summary <-apply(chrall_nonHDL, 1, seperate)
chrall_nonHDL$trait <- "nonHDL"
chrall_nonHDL$trait_pos <- paste(chrall_nonHDL$trait, chrall_nonHDL$closest_index_variant, sep=":")
chrall_nonHDL<- chrall_nonHDL[,c(2,4)]


chrall_TC<- read.delim("chrALL_TC_MR-MEGA.out.result_rsID_pval_lt_5e-8_missense_stop_splice_frameshift_anno_new.txt",header = T)
chrall_TC<- chrall_TC[chrall_TC$LD_r2_with_nearest_index != "unable_to_lookup_LD",]
chrall_TC<- chrall_TC[complete.cases(chrall_TC[ ,43]),]
chrall_TC$LD_r2_with_nearest_index <- as.numeric(chrall_TC$LD_r2_with_nearest_index)
chrall_TC<- chrall_TC[chrall_TC$LD_r2_with_nearest_index > 0.8,]
chrall_TC<- chrall_TC[,c(40:41)]
chrall_TC<- chrall_TC[!duplicated(chrall_TC), ]
chrall_TC<- aggregate(SnpEff_ensembl_summary ~ closest_index_variant,data=chrall_TC,paste,collapse=';')
seperate <- function(x){
  paste2(unique(unlist(lapply(strsplit(unlist(strsplit(unlist(strsplit(x[2],";")),"\\|")),"\\("), `[[`, 1))),collapse=";")}
chrall_TC$SnpEff_ensembl_summary <-apply(chrall_TC, 1, seperate)
chrall_TC$trait <- "TC"
chrall_TC$trait_pos <- paste(chrall_TC$trait, chrall_TC$closest_index_variant, sep=":")
chrall_TC<- chrall_TC[,c(2,4)]

ld_coding <- rbind(chrall_hdl,chrall_ldl,chrall_logTG,chrall_nonHDL,chrall_TC)
ld_0.8 <- ld_coding %>% tidyr::separate(trait_pos, 
                                        c("trait", "chr","pos (build 37)", "EA","NEA"))
ld_0.8 <- ld_0.8[,c(2,3,4,1)]
colnames(ld_0.8)[4] <- "LD_coding_variants_0.8"

write.csv(ld_0.8,file = "ld08.csv",quote = F,row.names = F)


#####################VCF for ANNOVAR#################################
library(tidyverse)
paste2 <- function(...,sep=", ") {
  L <- list(...)
  L <- lapply(L,function(x) {x[is.na(x)] <- ""; x})
  ret <-gsub(paste0("(^",sep,"|",sep,"$)"),"",
             gsub(paste0(sep,sep),sep,
                  do.call(paste,c(L,list(sep=sep)))))
  is.na(ret) <- ret==""
  ret
}

credible_set <-read.delim("99_percent_credible_set_variants_with_PP_updated_summary_by_transancestry_locus.txt",header = T)
credible_set <- credible_set[,c(1,2,3,11)]
credible_set <- credible_set[credible_set$anno_variants!="",]


credible_set$variants <-apply(credible_set, 1, function(x){
  paste2(unlist((lapply(strsplit(unlist(strsplit(x[4],",")),"_"),`[[`, 1))),collapse=";")})


credible_set$genes <-apply(credible_set, 1, function(x){
  paste2(unique(unlist(lapply(strsplit(unlist(strsplit(unlist(lapply(strsplit(unlist(strsplit(x[4],",")),"anno:"),`[[`, 2)),"\\|")),"\\("),`[[`, 1))),collapse=";")})

credible_set <- separate_rows(credible_set,variants,sep=";")
credible_set <- separate(credible_set,variants,c("chr", "pos","ea","nea"),sep=":")
credible_set$end <- credible_set$pos
credible_set <- credible_set[,c(5,6,10,7,8,1,2,3,9)]

chrall_hdl <- read.delim("chrALL_HDL_MR-MEGA.out.result_rsID_pval_lt_5e-8_missense_stop_splice_frameshift_anno_new.txt",header = T)
chrall_hdl <- chrall_hdl[complete.cases(chrall_hdl[ ,43]),]
chrall_hdl <- chrall_hdl[chrall_hdl$LD_r2_with_nearest_index != "unable_to_lookup_LD",]
chrall_hdl$LD_r2_with_nearest_index <- as.numeric(chrall_hdl$LD_r2_with_nearest_index)
chrall_hdl <- chrall_hdl[chrall_hdl$LD_r2_with_nearest_index > 0.5,]
chrall_hdl$genes <-apply(chrall_hdl, 1, function(x){
  paste2(unique(unlist(lapply(strsplit(unlist(strsplit(unlist(strsplit(x[40],";")),"\\|")),"\\("), `[[`, 1))),collapse=";")})
chrall_hdl$trait <- "HDL"
chrall_hdl$end <- chrall_hdl$Position
chrall_hdl <- separate(chrall_hdl,closest_index_variant,c("index_chr", "index_pos","index_ea","index_nea"),sep=":")
chrall_hdl <- chrall_hdl[,c(2,3,50,5,4,49,41,42,48)]
colnames(chrall_hdl) <- colnames(credible_set)

chrall_ldl <- read.delim("chrALL_LDL_MR-MEGA.out.result_rsID_pval_lt_5e-8_missense_stop_splice_frameshift_anno_new.txt",header = T)
chrall_ldl <- chrall_ldl[complete.cases(chrall_ldl[ ,43]),]
chrall_ldl <- chrall_ldl[chrall_ldl$LD_r2_with_nearest_index != "unable_to_lookup_LD",]
chrall_ldl$LD_r2_with_nearest_index <- as.numeric(chrall_ldl$LD_r2_with_nearest_index)
chrall_ldl <- chrall_ldl[chrall_ldl$LD_r2_with_nearest_index > 0.5,]
chrall_ldl$genes <-apply(chrall_ldl, 1, function(x){
  paste2(unique(unlist(lapply(strsplit(unlist(strsplit(unlist(strsplit(x[40],";")),"\\|")),"\\("), `[[`, 1))),collapse=";")})
chrall_ldl$trait <- "LDL"
chrall_ldl$end <- chrall_ldl$Position
chrall_ldl <- separate(chrall_ldl,closest_index_variant,c("index_chr", "index_pos","index_ea","index_nea"),sep=":")
chrall_ldl <- chrall_ldl[,c(2,3,50,5,4,49,41,42,48)]
colnames(chrall_ldl) <- colnames(credible_set)

chrall_logTG <- read.delim("chrALL_logTG_MR-MEGA.out.result_rsID_pval_lt_5e-8_missense_stop_splice_frameshift_anno_new.txt",header = T)
chrall_logTG <- chrall_logTG[complete.cases(chrall_logTG[ ,43]),]
chrall_logTG <- chrall_logTG[chrall_logTG$LD_r2_with_nearest_index != "unable_to_lookup_LD",]
chrall_logTG$LD_r2_with_nearest_index <- as.numeric(chrall_logTG$LD_r2_with_nearest_index)
chrall_logTG <- chrall_logTG[chrall_logTG$LD_r2_with_nearest_index > 0.5,]
chrall_logTG$genes <-apply(chrall_logTG, 1, function(x){
  paste2(unique(unlist(lapply(strsplit(unlist(strsplit(unlist(strsplit(x[40],";")),"\\|")),"\\("), `[[`, 1))),collapse=";")})
chrall_logTG$trait <- "logTG"
chrall_logTG$end <- chrall_logTG$Position
chrall_logTG <- separate(chrall_logTG,closest_index_variant,c("index_chr", "index_pos","index_ea","index_nea"),sep=":")
chrall_logTG <- chrall_logTG[,c(2,3,50,5,4,49,41,42,48)]
colnames(chrall_logTG) <- colnames(credible_set)

chrall_nonHDL <- read.delim("chrALL_nonHDL_MR-MEGA.out.result_rsID_pval_lt_5e-8_missense_stop_splice_frameshift_anno_new.txt",header = T)
chrall_nonHDL <- chrall_nonHDL[complete.cases(chrall_nonHDL[ ,43]),]
chrall_nonHDL <- chrall_nonHDL[chrall_nonHDL$LD_r2_with_nearest_index != "unable_to_lookup_LD",]
chrall_nonHDL$LD_r2_with_nearest_index <- as.numeric(chrall_nonHDL$LD_r2_with_nearest_index)
chrall_nonHDL <- chrall_nonHDL[chrall_nonHDL$LD_r2_with_nearest_index > 0.5,]
chrall_nonHDL$genes <-apply(chrall_nonHDL, 1, function(x){
  paste2(unique(unlist(lapply(strsplit(unlist(strsplit(unlist(strsplit(x[40],";")),"\\|")),"\\("), `[[`, 1))),collapse=";")})
chrall_nonHDL$trait <- "nonHDL"
chrall_nonHDL$end <- chrall_nonHDL$Position
chrall_nonHDL <- separate(chrall_nonHDL,closest_index_variant,c("index_chr", "index_pos","index_ea","index_nea"),sep=":")
chrall_nonHDL <- chrall_nonHDL[,c(2,3,50,5,4,49,41,42,48)]
colnames(chrall_nonHDL) <- colnames(credible_set)

chrall_TC <- read.delim("chrALL_TC_MR-MEGA.out.result_rsID_pval_lt_5e-8_missense_stop_splice_frameshift_anno_new.txt",header = T)
chrall_TC <- chrall_TC[complete.cases(chrall_TC[ ,43]),]
chrall_TC <- chrall_TC[chrall_TC$LD_r2_with_nearest_index != "unable_to_lookup_LD",]
chrall_TC$LD_r2_with_nearest_index <- as.numeric(chrall_TC$LD_r2_with_nearest_index)
chrall_TC <- chrall_TC[chrall_TC$LD_r2_with_nearest_index > 0.5,]
chrall_TC$genes <-apply(chrall_TC, 1, function(x){
  paste2(unique(unlist(lapply(strsplit(unlist(strsplit(unlist(strsplit(x[40],";")),"\\|")),"\\("), `[[`, 1))),collapse=";")})
chrall_TC$trait <- "TC"
chrall_TC$end <- chrall_TC$Position
chrall_TC <- separate(chrall_TC,closest_index_variant,c("index_chr", "index_pos","index_ea","index_nea"),sep=":")
chrall_TC <- chrall_TC[,c(2,3,50,5,4,49,41,42,48)]
colnames(chrall_TC) <- colnames(credible_set)

demage_variant <- rbind(chrall_hdl,chrall_ldl,chrall_logTG,chrall_nonHDL,chrall_TC)
demage_variant <- arrange(demage_variant,chr,pos)

write.table(demage_variant, "coding.avinput", row.names = F,col.names = F,quote = F,sep = '\t')

### Gene based test
library(splitstackshape)
hdl_gene_based <- read.table("gene_based_results_2019Sep16/Results_HDL_ALL.txt",header = T)
hdl_gene_based <- hdl_gene_based[hdl_gene_based$p.value <0.001,]
hdl_gene_based <- hdl_gene_based[,c(1,17)]
hdl_gene_based <- hdl_gene_based[!duplicated(hdl_gene_based),]

ldl_gene_based <- read.table("lipid_exomes/gene_based_results_2019Sep16/Results_LDL_ADJ_ALL.txt",header = T)
ldl_gene_based <- ldl_gene_based[ldl_gene_based$p.value <0.001,]
ldl_gene_based <- ldl_gene_based[,c(1,17)]
ldl_gene_based <- ldl_gene_based[!duplicated(ldl_gene_based),]

nonHDL_gene_based <- read.table("lipid_exomes/gene_based_results_2019Sep16/Results_NON_HDL_ALL.txt",header = T)
nonHDL_gene_based <- nonHDL_gene_based[nonHDL_gene_based$p.value <0.001,]
nonHDL_gene_based <- nonHDL_gene_based[,c(1,17)]
nonHDL_gene_based <- nonHDL_gene_based[!duplicated(nonHDL_gene_based),]

logTG_gene_based <- read.table("lipid_exomes/gene_based_results_2019Sep16/Results_TG_ALL.txt",header = T)
logTG_gene_based <- logTG_gene_based[logTG_gene_based$p.value <0.001,]
logTG_gene_based <- logTG_gene_based[,c(1,17)]
logTG_gene_based <- logTG_gene_based[!duplicated(logTG_gene_based),]

tc_gene_based <- read.table("lipid_exomes/gene_based_results_2019Sep16/Results_TOTAL_ADJ_ALL.txt",header = T)
tc_gene_based <- tc_gene_based[tc_gene_based$p.value <0.001,]
tc_gene_based <- tc_gene_based[,c(1,17)]
tc_gene_based <- tc_gene_based[!duplicated(tc_gene_based),]

# #######################################################################################
# hdl <- read.table("lipid_exomes/gene_based_results_2019Sep16/Results_HDL_ALL.txt",header = T)
# hdl <- hdl[hdl$p.value <0.001,]
# hdl <- hdl[,c(1,17)]
# hdl <- hdl[!duplicated(hdl),]
# 
# datalist = list()
# 
# for (i in 1:nrow(hdl)) {
#   # ... make some data
#   dat <- data.frame(sapply(strsplit(unlist(strsplit(hdl[i,2],",")),"\\/"), `[[`, 1))
#   colnames(dat)[1] <- "pos"
#   dat$gene <- hdl[i,1]  # maybe you want to keep track of which iteration produced it?
#   datalist[[i]] <- dat # add it to your list
# }
# 
# hdl_data = do.call(rbind, datalist)
# hdl_data$trait <- "HDL"
# 
# 
# ldl <- read.table("lipid_exomes/gene_based_results_2019Sep16/Results_LDL_ADJ_ALL.txt",header = T)
# ldl <- ldl[ldl$p.value <0.001,]
# ldl <- ldl[,c(1,17)]
# ldl <- ldl[!duplicated(ldl),]
# 
# datalist = list()
# 
# for (i in 1:nrow(ldl)) {
#   # ... make some data
#   dat <- data.frame(sapply(strsplit(unlist(strsplit(ldl[i,2],",")),"\\/"), `[[`, 1))
#   colnames(dat)[1] <- "pos"
#   dat$gene <- ldl[i,1]  # maybe you want to keep track of which iteration produced it?
#   datalist[[i]] <- dat # add it to your list
# }
# 
# ldl_data = do.call(rbind, datalist)
# ldl_data$trait <- "LDL"
# 
# nonHDL <- read.table("lipid_exomes/gene_based_results_2019Sep16/Results_NON_HDL_ALL.txt",header = T)
# nonHDL <- nonHDL[nonHDL$p.value <0.001,]
# nonHDL <- nonHDL[,c(1,17)]
# nonHDL <- nonHDL[!duplicated(nonHDL),]
# 
# datalist = list()
# 
# for (i in 1:nrow(nonHDL)) {
#   # ... make some data
#   dat <- data.frame(sapply(strsplit(unlist(strsplit(nonHDL[i,2],",")),"\\/"), `[[`, 1))
#   colnames(dat)[1] <- "pos"
#   dat$gene <- nonHDL[i,1]  # maybe you want to keep track of which iteration produced it?
#   datalist[[i]] <- dat # add it to your list
# }
# 
# nonHDL_data = do.call(rbind, datalist)
# nonHDL_data$trait <- "nonHDL"
# 
# 
# logTG <- read.table("lipid_exomes/gene_based_results_2019Sep16/Results_TG_ALL.txt",header = T)
# logTG <- logTG[logTG$p.value <0.001,]
# logTG <- logTG[,c(1,17)]
# logTG <- logTG[!duplicated(logTG),]
# 
# datalist = list()
# 
# for (i in 1:nrow(logTG)) {
#   # ... make some data
#   dat <- data.frame(sapply(strsplit(unlist(strsplit(logTG[i,2],",")),"\\/"), `[[`, 1))
#   colnames(dat)[1] <- "pos"
#   dat$gene <- logTG[i,1]  # maybe you want to keep track of which iteration produced it?
#   datalist[[i]] <- dat # add it to your list
# }
# 
# logTG_data = do.call(rbind, datalist)
# logTG_data$trait <- "logTG"
# 
# 
# tc <- read.table("lipid_exomes/gene_based_results_2019Sep16/Results_TOTAL_ADJ_ALL.txt",header = T)
# tc <- tc[tc$p.value <0.001,]
# tc <- tc[,c(1,17)]
# tc <- tc[!duplicated(tc),]
# 
# datalist = list()
# 
# for (i in 1:nrow(tc)) {
#   # ... make some data
#   dat <- data.frame(sapply(strsplit(unlist(strsplit(tc[i,2],",")),"\\/"), `[[`, 1))
#   colnames(dat)[1] <- "pos"
#   dat$gene <- tc[i,1]  # maybe you want to keep track of which iteration produced it?
#   datalist[[i]] <- dat # add it to your list
# }
# 
# tc_data = do.call(rbind, datalist)
# tc_data$trait <- "TC"
# 
# # all_traits <- rbind(hdl_data,ldl_data,logTG_data,nonHDL_data,tc_data)
# # all_traits$chr <- "chr"
# # all_traits$tmp <- paste0(all_traits$chr,all_traits$pos)
# # all_traits$final <- paste0(all_traits$tmp,"-",all_traits$pos_2)
# # 
# # txt <- all_traits[,c(8)]
# # write.table(txt,"pos_build_38.txt",quote = F,row.names = F,col.names = F)
# # 
# # pos_build_37 <- read.table("pos_build_37.txt",header = F)
# # all_traits_37 <- cbind(all_traits[,c(2,3)],pos_build_37)
# # all_traits_37 <- cSplit(all_traits_37, "V1", sep="chr",drop = T)
# # all_traits_37 <- cSplit(all_traits_37, "V1_2", sep=":",drop = T)
# # all_traits_37 <- cSplit(all_traits_37, "V1_2_2", sep="-",drop = T)
# # all_traits_37 <- all_traits_37[,c(1,2,4,5)]
# # 
# # colnames(all_traits_37) <- c("Gene_Based", "trait","chr","pos (build 37)")
# # all_traits_37 <- all_traits_37[with(all_traits_37,order(chr,`pos (build 37)`)),]
# # write.csv(all_traits_37,"gene_based_build_37.csv",row.names = F,quote = F)
# # 
