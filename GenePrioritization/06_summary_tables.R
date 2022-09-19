paste2 <- function(...,sep=", ") {
  L <- list(...)
  L <- lapply(L,function(x) {x[is.na(x)] <- ""; x})
  ret <-gsub(paste0("(^",sep,"|",sep,"$)"),"",
             gsub(paste0(sep,sep),sep,
                  do.call(paste,c(L,list(sep=sep)))))
  is.na(ret) <- ret==""
  ret
}


top2_closest_protein_coding <- read.table("trans_index_protein_coding_genes_top2.txt",header = T)
tmp <- read_csv("trans-ancestry-02-16.csv", na = c("NA","#N/A"))
tmp <- arrange(tmp,chr,`pos (build 37)`,trait)

full_data <- read_csv("rans-ancestry.csv", na = c("NA","#N/A"))
full_data <- arrange(full_data,chr,`pos (build 37)`,trait)

full_data$TWAS <- tmp$TWAS
full_data$`Nearest_2_protein_coding_genes_(Ensembl)` <- tmp$`Nearest_2_protein_coding_genes_(Ensembl)`


trans_ancestry <- full_data%>% 
  select(trait,                                       # 1
         rsid_dbSNP150,                               # 2
         chr,                                         # 3
         `pos (build 37)`,                            # 4
         overlap,                                     # 5
         VEP_ensembl_summary,                         # 6
         Mendelian_gene,                              # 7
         eQTL,                                        # 8
         TWAS,                                        # 9
         PoPS_Top_2_Prioritized_Genes,                # 10
         `Nearest_2_protein_coding_genes_(Ensembl)`,  # 11
         `DEPICT (FDR < 0.05)`)                       # 12

# PoPS Top 20%  Top 1 #10

trans_ancestry[trans_ancestry$trait=="HDL",]$PoPS_Top_2_Prioritized_Genes <- as.character(apply(trans_ancestry[trans_ancestry$trait=="HDL",], 1, function(x){
  pops <- c(unlist(strsplit(unlist(strsplit(x[10],":")),";"))[1])
  result <-  intersect(pops,hdl_genes_pops$Genes)
  paste2(unlist(result),collapse=";")
}))

trans_ancestry[trans_ancestry$trait=="LDL",]$PoPS_Top_2_Prioritized_Genes <- as.character(apply(trans_ancestry[trans_ancestry$trait=="LDL",], 1, function(x){
  pops <- c(unlist(strsplit(unlist(strsplit(x[10],":")),";"))[1])
  result <-  intersect(pops,ldl_genes_pops$Genes)
  paste2(unlist(result),collapse=";")
}))

trans_ancestry[trans_ancestry$trait=="nonHDL",]$PoPS_Top_2_Prioritized_Genes <- as.character(apply(trans_ancestry[trans_ancestry$trait=="nonHDL",], 1, function(x){
  pops <- c(unlist(strsplit(unlist(strsplit(x[10],":")),";"))[1])
  result <-  intersect(pops,nonHDL_genes_pops$Genes)
  paste2(unlist(result),collapse=";")
}))

trans_ancestry[trans_ancestry$trait=="logTG",]$PoPS_Top_2_Prioritized_Genes <- as.character(apply(trans_ancestry[trans_ancestry$trait=="logTG",], 1, function(x){
  pops <- c(unlist(strsplit(unlist(strsplit(x[10],":")),";"))[1])
  result <-  intersect(pops,logTG_genes_pops$Genes)
  paste2(unlist(result),collapse=";")
}))

trans_ancestry[trans_ancestry$trait=="TC",]$PoPS_Top_2_Prioritized_Genes <- as.character(apply(trans_ancestry[trans_ancestry$trait=="TC",], 1, function(x){
  pops <- c(unlist(strsplit(unlist(strsplit(x[10],":")),";"))[1])
  result <-  intersect(pops,tc_genes_pops$Genes)
  paste2(unlist(result),collapse=";")
}))

# Closest Genes top 1 #11

trans_ancestry$`Nearest_2_protein_coding_genes_(Ensembl)` <-apply(trans_ancestry, 1, function(x){
  closest <- c(unlist(strsplit(unlist(strsplit(x[11],":")),";"))[1])
  result <- intersect(closest,top2_closest_protein_coding$closest_gene)
  paste2(unlist(result),collapse=";")
})

# Correct several gene names in Closest gene column
trans_ancestry$`Nearest_2_protein_coding_genes_(Ensembl)`[trans_ancestry$`Nearest_2_protein_coding_genes_(Ensembl)` == 43891] <- "MARCHF1"
trans_ancestry$`Nearest_2_protein_coding_genes_(Ensembl)`[trans_ancestry$`Nearest_2_protein_coding_genes_(Ensembl)` == 43898] <- "MARCHF8"
trans_ancestry$`DEPICT (FDR < 0.05)`[trans_ancestry$`DEPICT (FDR < 0.05)` == 43898] <- "MARCHF8"

# Lipid eQTL. # 13
trans_ancestry$Lipid_eQTL <- ifelse(grepl("Liver|Adipose_Subcutaneous|Adipose_Visceral|Whole_Blood|Small_Intestine", trans_ancestry$eQTL),trans_ancestry$eQTL,NA)

lipid_eQTL <- function(x){
  temp <- unlist(strsplit(x[13],","))
  ind <- grepl("Liver|Adipose_Subcutaneous|Adipose_Visceral|Whole_Blood|Small_Intestine",temp)
  result <- temp[ind]
  result2 <- paste2(result,collapse=";")
  paste2(unique(lapply(strsplit(unlist(strsplit(result2,";")),"_"), `[[`, 1)),collapse=";")
}
trans_ancestry$Lipid_eQTL <-apply(trans_ancestry, 1, lipid_eQTL)

# Lipid TWAS. # 14
trans_ancestry$Lipid_TWAS <- ifelse(grepl("Liver|Adipose_Subcutaneous|Adipose_Visceral|Whole_Blood|Small_Intestine", trans_ancestry$TWAS),trans_ancestry$TWAS,NA)
lipid_TWAS <- function(x){
  temp <- unlist(strsplit(x[14],","))
  ind <- grepl("Liver|Adipose_Subcutaneous|Adipose_Visceral|Whole_Blood|Small_Intestine",temp)
  result <- temp[ind]
  result2 <- paste2(result,collapse=";")
  if(is.na(result2)) {
    NA
  } else {
    paste2(unique(lapply(strsplit(unlist(strsplit(result2,";")),":"), `[[`, 2)),collapse=";")
  }
}

trans_ancestry$Lipid_TWAS <-apply(trans_ancestry, 1, lipid_TWAS)

all_eQTL <-function(x){
  temp <- unlist(strsplit(x[8],","))
  result <- temp
  result2 <- paste2(result,collapse=";")
  paste2(unique(lapply(strsplit(unlist(strsplit(result2,";")),"_"), `[[`, 1)),collapse=";")
}

trans_ancestry$eQTL <-apply(trans_ancestry, 1, all_eQTL)

all_TWAS <- function(x){
  if(is.na(x[9])) {
    NA
  } else {
    paste2(unique(lapply(strsplit(unlist(strsplit(x[9],",")),":"), `[[`, 2)),collapse=";")
  }
}
trans_ancestry$TWAS <-apply(trans_ancestry, 1, all_TWAS)

# Star genes
gene_6_star <- function(x){
  pops <- unlist(strsplit(x[10],";"))
  closest <- unlist(strsplit(x[11],";"))
  depict <-unlist(strsplit(x[12],";"))
  eqtl <-unlist(strsplit(x[8],";"))
  twas <-unlist(strsplit(x[9],";"))
  coding <- unlist(strsplit(x[15],";"))
  all <- table(sort(na.omit(c(pops,closest,depict,eqtl,twas,coding))))
  paste2(names(all[all==6]),collapse=";")
}
trans_ancestry$gene_6_star <-as.character(apply(trans_ancestry, 1, gene_6_star))

gene_5_star <- function(x){
  pops <- unlist(strsplit(x[10],";"))
  closest <- unlist(strsplit(x[11],";"))
  depict <-unlist(strsplit(x[12],";"))
  eqtl <-unlist(strsplit(x[8],";"))
  twas <-unlist(strsplit(x[9],";"))
  coding <- unlist(strsplit(x[15],";"))
  all <- table(sort(na.omit(c(pops,closest,depict,eqtl,twas,coding))))
  paste2(names(all[all==5]),collapse=";")
}

trans_ancestry$gene_5_star <-as.character(apply(trans_ancestry, 1, gene_5_star))

gene_4_star <- function(x){
  pops <- unlist(strsplit(x[10],";"))
  closest <- unlist(strsplit(x[11],";"))
  depict <-unlist(strsplit(x[12],";"))
  eqtl <-unlist(strsplit(x[8],";"))
  twas <-unlist(strsplit(x[9],";"))
  coding <- unlist(strsplit(x[15],";"))
  all <- table(sort(na.omit(c(pops,closest,depict,eqtl,twas,coding))))
  paste2(names(all[all==4]),collapse=";")
}

trans_ancestry$gene_4_star <-as.character(apply(trans_ancestry, 1, gene_4_star))

gene_3_star <- function(x){
  pops <- unlist(strsplit(x[10],";"))
  closest <- unlist(strsplit(x[11],";"))
  depict <-unlist(strsplit(x[12],";"))
  eqtl <-unlist(strsplit(x[8],";"))
  twas <-unlist(strsplit(x[9],";"))
  coding <- unlist(strsplit(x[15],";"))
  all <- table(sort(na.omit(c(pops,closest,depict,eqtl,twas,coding))))
  paste2(names(all[all==3]),collapse=";")
}

trans_ancestry$gene_3_star <-as.character(apply(trans_ancestry, 1, gene_3_star))

gene_2_star <- function(x){
  pops <- unlist(strsplit(x[10],";"))
  closest <- unlist(strsplit(x[11],";"))
  depict <-unlist(strsplit(x[12],";"))
  eqtl <-unlist(strsplit(x[8],";"))
  twas <-unlist(strsplit(x[9],";"))
  coding <- unlist(strsplit(x[15],";"))
  all <- table(sort(na.omit(c(pops,closest,depict,eqtl,twas,coding))))
  paste2(names(all[all==2]),collapse=";")
}

trans_ancestry$gene_2_star <-as.character(apply(trans_ancestry, 1, gene_2_star))

gene_1_star <- function(x){
  pops <- unlist(strsplit(x[10],";"))
  closest <- unlist(strsplit(x[11],";"))
  depict <-unlist(strsplit(x[12],";"))
  eqtl <-unlist(strsplit(x[8],";"))
  twas <-unlist(strsplit(x[9],";"))
  coding <- unlist(strsplit(x[15],";"))
  all <- table(sort(na.omit(c(pops,closest,depict,eqtl,twas,coding))))
  paste2(names(all[all==1]),collapse=";")
}

trans_ancestry$gene_1_star <-as.character(apply(trans_ancestry, 1, gene_1_star))

# Coloc
coloc_hdl_cad <- read.delim("./coloc/coloc_results_v2/coloc_HDL_CAD.txt",header=F)
colnames(coloc_hdl_cad) <- c("snp", "chr", "start", "end", "nsnps", "PP.H0.abf",
                             "PP.H1.abf" , "PP.H2.abf" , "PP.H3.abf", "PP.H4.abf")
coloc_hdl_cad$trait <- "HDL"

coloc_ldl_cad <- read.delim("./coloc/coloc_results_v2/coloc_LDL_CAD.txt",header=F)
colnames(coloc_ldl_cad) <- c("snp", "chr", "start", "end", "nsnps", "PP.H0.abf",
                             "PP.H1.abf" , "PP.H2.abf" , "PP.H3.abf", "PP.H4.abf")
coloc_ldl_cad$trait <- "LDL"

coloc_nonhdl_cad <- read.delim("./coloc/coloc_results_v2/coloc_nonHDL_CAD.txt",header=F)
colnames(coloc_nonhdl_cad) <- c("snp", "chr", "start", "end", "nsnps", "PP.H0.abf",
                                "PP.H1.abf" , "PP.H2.abf" , "PP.H3.abf", "PP.H4.abf")
coloc_nonhdl_cad$trait <- "nonHDL"

coloc_logtg_cad <- read.delim("./coloc/coloc_results_v2/coloc_logTG_CAD.txt",header=F)
colnames(coloc_logtg_cad) <- c("snp", "chr", "start", "end", "nsnps", "PP.H0.abf",
                               "PP.H1.abf" , "PP.H2.abf" , "PP.H3.abf", "PP.H4.abf")
coloc_logtg_cad$trait <- "logTG"

coloc_tc_cad <- read.delim("./coloc/coloc_results_v2/coloc_TC_CAD.txt",header=F)
colnames(coloc_tc_cad) <- c("snp", "chr", "start", "end", "nsnps", "PP.H0.abf",
                            "PP.H1.abf" , "PP.H2.abf" , "PP.H3.abf", "PP.H4.abf")
coloc_tc_cad$trait <- "TC"
coloc_lipids_cad <- rbind(coloc_hdl_cad,coloc_ldl_cad,coloc_logtg_cad,coloc_nonhdl_cad,coloc_tc_cad)
coloc_lipids_cad$`pos (build 37)` = coloc_lipids_cad$start+100000
##########
coloc_hdl_t2d <- read.delim("./coloc/coloc_results_v2/coloc_HDL_T2D.txt",header=F)
colnames(coloc_hdl_t2d) <- c("snp", "chr", "start", "end", "nsnps", "PP.H0.abf",
                             "PP.H1.abf" , "PP.H2.abf" , "PP.H3.abf", "PP.H4.abf")
coloc_hdl_t2d$trait <- "HDL"

coloc_ldl_t2d <- read.delim("./coloc/coloc_results_v2/coloc_LDL_T2D.txt",header=F)
colnames(coloc_ldl_t2d) <- c("snp", "chr", "start", "end", "nsnps", "PP.H0.abf",
                             "PP.H1.abf" , "PP.H2.abf" , "PP.H3.abf", "PP.H4.abf")
coloc_ldl_t2d$trait <- "LDL"

coloc_nonhdl_t2d <- read.delim("./coloc/coloc_results_v2/coloc_nonHDL_T2D.txt",header=F)
colnames(coloc_nonhdl_t2d) <- c("snp", "chr", "start", "end", "nsnps", "PP.H0.abf",
                                "PP.H1.abf" , "PP.H2.abf" , "PP.H3.abf", "PP.H4.abf")
coloc_nonhdl_t2d$trait <- "nonHDL"

coloc_logtg_t2d <- read.delim("./coloc/coloc_results_v2/coloc_logTG_T2D.txt",header=F)
colnames(coloc_logtg_t2d) <- c("snp", "chr", "start", "end", "nsnps", "PP.H0.abf",
                               "PP.H1.abf" , "PP.H2.abf" , "PP.H3.abf", "PP.H4.abf")
coloc_logtg_t2d$trait <- "logTG"

coloc_tc_t2d <- read.delim("./coloc/coloc_results_v2/coloc_TC_T2D.txt",header=F)
colnames(coloc_tc_t2d) <- c("snp", "chr", "start", "end", "nsnps", "PP.H0.abf",
                            "PP.H1.abf" , "PP.H2.abf" , "PP.H3.abf", "PP.H4.abf")
coloc_tc_t2d$trait <- "TC"
coloc_lipids_t2d <- rbind(coloc_hdl_t2d,coloc_ldl_t2d,coloc_logtg_t2d,coloc_nonhdl_t2d,coloc_tc_t2d)
coloc_lipids_t2d$`pos (build 37)` = coloc_lipids_t2d$start+100000

# Coding variants (credible set + ld 0.5) damage missense #15
# damage_variants <- read.table("/rprojectnb2/pelosolab/yuxuan/glgc/Credible_set/coding.hg19_ljb23_metasvm_dropped")
# damage_variants <- separate(damage_variants,V2,c("score1", "score2","ind"),sep=",")
# damage_variants <- damage_variants[damage_variants$ind=="D",]
# damage_variants <- damage_variants[,c(10:13)]
# colnames(damage_variants) <- c("trait", "chr","pos (build 37)","damage_missense")
# trans_ancestry <- merge(trans_ancestry,damage_variants,by=c("trait", "chr","pos (build 37)"),sort = F,all.x = T)

coding_variants <- read_csv("credible_set.csv")
trans_ancestry <- merge(trans_ancestry,coding_variants,by=c("trait", "chr","pos (build 37)"),sort = F,all.x = T)
trans_ancestry <- arrange(trans_ancestry,chr,`pos (build 37)`,trait)
ld_0.8 <- read_csv("/rprojectnb2/pelosolab/yuxuan/glgc/Credible_set/ld08.csv")
trans_ancestry <- merge(trans_ancestry,ld_0.8,by=c("trait", "chr","pos (build 37)"),sort = F,all.x = T)
trans_ancestry <- arrange(trans_ancestry,chr,`pos (build 37)`,trait)
trans_ancestry <- trans_ancestry %>% tidyr::unite("coding_variants",credible_set:LD_coding_variants_0.8,sep=";",remove=T,na.rm=T)
trans_ancestry[trans_ancestry$coding_variants=="",]$coding_variants <- NA 
trans_ancestry$coding_variants <- sapply(strsplit(trans_ancestry$coding_variants, ";"), function(x) paste2(unique(x), collapse = ";"))


# Mouse model level 3 results #16
trans_ethnic_mouse_prioritised_genes <- read_excel("trans_ethnic_mouse_prioritised_genes_18Jan2021.xlsx")
trans_ethnic_mouse_prioritised_genes <- trans_ethnic_mouse_prioritised_genes[trans_ethnic_mouse_prioritised_genes$priority_level == 3,]
trans_ethnic_mouse_prioritised_genes <- trans_ethnic_mouse_prioritised_genes[,c("rsid","trait","Symbol")]
trans_ethnic_mouse_prioritised_genes$Symbol <- gsub("'","",trans_ethnic_mouse_prioritised_genes$Symbol)
colnames(trans_ethnic_mouse_prioritised_genes)[3] <- "mouse_model_genes"
trans_ethnic_mouse_prioritised_genes <- separate_rows(trans_ethnic_mouse_prioritised_genes,rsid,sep = ";")
trans_ethnic_mouse_prioritised_genes <-aggregate(mouse_model_genes ~ rsid + trait,data = trans_ethnic_mouse_prioritised_genes, function(x) paste2(x, collapse=";"))
trans_ancestry <- merge(trans_ancestry,trans_ethnic_mouse_prioritised_genes,by.x = c("trait","rsid_dbSNP150"),by.y = c("trait","rsid"),all.x = T)
trans_ancestry <- arrange(trans_ancestry,chr,`pos (build 37)`,trait)


# PoPS and a local. # 17
pops_and_a_local <- function(x){
  pops <- unlist(strsplit(x[10],";"))
  closest <- unlist(strsplit(x[11],";"))
  eqtl <-unlist(strsplit(x[8],";"))
  twas <-unlist(strsplit(x[9],";"))
  credible_set <-unlist(strsplit(x[15],";"))
  pops_closest <- intersect(pops,closest)
  pops_eqtl <-intersect(pops,eqtl)
  pops_twas <- intersect(pops,twas)
  pops_credible_set <- intersect(pops,credible_set)
  result <-unique(c(pops_twas,pops_eqtl,pops_closest,pops_credible_set))
  paste2(unlist(result),collapse=";")
}
trans_ancestry$pops_match_a_local_method <-as.character(apply(trans_ancestry, 1, pops_and_a_local))

# DEPICT and a local. # 18
depict_and_a_local <- function(x){
  depict <-unlist(strsplit(x[12],";"))
  closest <- unlist(strsplit(x[11],";"))
  eqtl <-unlist(strsplit(x[8],";"))
  twas <-unlist(strsplit(x[9],";"))
  credible_set <-unlist(strsplit(x[15],";"))
  depict_closest <- intersect(depict,closest)
  depict_eqtl <-intersect(depict,eqtl)
  depict_twas <- intersect(depict,twas)
  depict_credible_set <- intersect(depict,credible_set)
  result <-unique(c(depict_twas,depict_eqtl,depict_closest,depict_credible_set))
  paste2(unlist(result),collapse=";")
}
trans_ancestry$depict_match_a_local_method <-as.character(apply(trans_ancestry, 1, depict_and_a_local))


# Gene based results 500K # 19
closest_genes <- read.table("/rprojectnb2/pelosolab/yuxuan/glgc/glgc_closest_genes/trans_index_protein_coding_genes_all.txt",header = T)
trans_ancestry <- merge(trans_ancestry,closest_genes,by="rsid_dbSNP150",all.x=T)
trans_ancestry$gene_based_500k <-NA


trans_ancestry[trans_ancestry$trait=="HDL",]$gene_based_500k <- as.character(apply(trans_ancestry[trans_ancestry$trait=="HDL",], 1, function(x){
  closest <- unlist(strsplit(x[19],";"))
  result <-  intersect(closest,hdl_gene_based$gene)
  paste2(unlist(result),collapse=";")
}))

trans_ancestry[trans_ancestry$trait=="LDL",]$gene_based_500k <- as.character(apply(trans_ancestry[trans_ancestry$trait=="LDL",], 1, function(x){
  closest <- unlist(strsplit(x[19],";"))
  result <-  intersect(closest,ldl_gene_based$gene)
  paste2(unlist(result),collapse=";")
}))

trans_ancestry[trans_ancestry$trait=="nonHDL",]$gene_based_500k <- as.character(apply(trans_ancestry[trans_ancestry$trait=="nonHDL",], 1, function(x){
  closest <- unlist(strsplit(x[19],";"))
  result <-  intersect(closest,nonHDL_gene_based$gene)
  paste2(unlist(result),collapse=";")
}))

trans_ancestry[trans_ancestry$trait=="logTG",]$gene_based_500k <- as.character(apply(trans_ancestry[trans_ancestry$trait=="logTG",], 1, function(x){
  closest <- unlist(strsplit(x[19],";"))
  result <-  intersect(closest,logTG_gene_based$gene)
  paste2(unlist(result),collapse=";")
}))

trans_ancestry[trans_ancestry$trait=="TC",]$gene_based_500k <- as.character(apply(trans_ancestry[trans_ancestry$trait=="TC",], 1, function(x){
  closest <- unlist(strsplit(x[19],";"))
  result <-  intersect(closest,tc_gene_based$gene)
  paste2(unlist(result),collapse=";")
}))
trans_ancestry <- trans_ancestry[,-c(19)]

# Update Coding variants (credible set + ld 0.8 + gene based) #15
trans_ancestry <- trans_ancestry %>% tidyr::unite("coding_variants",c(coding_variants,gene_based_500k),sep=";",remove=F,na.rm=T)
trans_ancestry[trans_ancestry$coding_variants=="",]$coding_variants <- NA 
trans_ancestry$coding_variants <- sapply(strsplit(trans_ancestry$coding_variants, ";"), function(x) paste2(unique(x), collapse = ";"))
trans_ancestry <- arrange(trans_ancestry,chr,`pos (build 37)`,trait)

# PoPS ALL
trans_ancestry$PoPS_All <- full_data$PoPS_Top_2_Prioritized_Genes
trans_ancestry$PoPS_All <- as.character(apply(trans_ancestry, 1, function(x){
  pops_all <- c(unlist(strsplit(unlist(strsplit(x[20],":")),";"))[1] , unlist(strsplit(unlist(strsplit(x[20],":")),";"))[3])
  paste2(unlist(pops_all),collapse=";")
}))
trans_ancestry[trans_ancestry$PoPS_All==";",]$PoPS_All <- NA

# T2D
t2d <- read_csv("glgc_variants_diagram_t2d_lookup.csv", 
                col_types = cols(t2d_eur_SNP = col_character()))

colnames(t2d)[c(6,7)] <- c("t2d_eur_NEA","t2d_eur_EA")
t2d$t2d_eur_Beta <- ifelse(t2d$t2d_eur_NEA == trans_ancestry_assign_final$`METAL Ref`,-t2d$t2d_eur_Beta,t2d$t2d_eur_Beta)
t2d <- arrange(t2d,chr,`pos (build 37)`,trait)
# CAD
cad <- read_csv("glgc_variants_mvp_cad_lookup.csv", 
                col_types = cols(cad_mvp_SNP_ID_b37 = col_character()))


colnames(cad)[c(6,7)] <- c("cad_mvp_NEA","cad_mvp_EA")
cad$cad_mvp_Beta <- ifelse(cad$cad_mvp_NEA == trans_ancestry_assign_final$`METAL Ref`,-cad$cad_mvp_Beta,cad$cad_mvp_Beta)
cad <- arrange(cad,chr,`pos (build 37)`,trait)
# NAFLD
nafld_tmp <- read.table(file = "/rprojectnb2/pelosolab/yuxuan/glgc/index_variants_with_HRC_UKB_MGI_meta.txt",header = T)
nafld_tmp$NAFLD_Ref <- toupper(nafld_tmp$NAFLD_Ref)
nafld_tmp$NAFLD_Alt <- toupper(nafld_tmp$NAFLD_Alt)
colnames(nafld_tmp)[c(2:3,57,58,60,62)] <- c("chr","pos (build 37)","nafld_tem_NEA","nafld_tem_EA","nafld_tem_Effect","nafld_tem_Pvalue")
nafld <- merge(trans_ancestry,nafld_tmp,by=c("trait", "chr","pos (build 37)"),sort = F,all.x = T)
nafld <- nafld[,c(4,1:3,80:85)]
nafld <- arrange(nafld,chr,`pos (build 37)`,trait)

identical(t2d$trait,full_data$trait)
identical(cad$trait,full_data$trait)
identical(nafld$trait,full_data$trait)

# GLGC
meta_summary <- fread("meta_summary_stat_all_index_variants.tbl")
meta_summary$Allele1 <- toupper(meta_summary$Allele1)
meta_summary$Allele2 <- toupper(meta_summary$Allele2)

table_s2 <- data.frame(trans_ancestry_assign[,c(1:4)],tmp[,c(33:40)],trans_ancestry_assign[,c(31,7:12,15:16,21:26)],check.names = F)
table_s2$`METAL Ref` <- toupper(table_s2$`METAL Ref`)
table_s2$`METAL Alt` <- toupper(table_s2$`METAL Alt`)

identical(table_s2$trait,full_data$trait)
identical(table_s2$rsid_dbSNP150,full_data$rsid_dbSNP150)
table_s2_ind <- inner_join(table_s2,independent_index_variant)

identical(trans_ancestry_assign_final$trait,t2d$trait)
identical(trans_ancestry_assign_final$rsid_dbSNP150,t2d$rsid_dbSNP150)
table_s3 <- data.frame(trans_ancestry_assign_final,t2d[,c(6,7,9,10,11)],cad[,c(6,7,9,10,11)],nafld[,c(5:6,8,9,10)],check.names = F)
table_s3 <- table_s3[!duplicated(table_s3), ]

table_s3$`Gene Prioritization` <- sapply(strsplit(table_s3$`Gene Prioritization`,split='|', fixed=TRUE), function(x) unlist(strsplit(unlist(x)[1],":"))[2])
table_s3$locus <- full_data$`Locus number accounting for LD between index variants`
table_s3 <- inner_join(table_s3,coloc_lipids_cad[,c(2,10:12)])
colnames(table_s3)[30] <- "pp4_lipids_cad"
table_s3 <- inner_join(table_s3,coloc_lipids_t2d[,c(2,10:12)])
colnames(table_s3)[31] <- "pp4_lipids_t2d"
table_s3 <- arrange(table_s3,chr,`pos (build 37)`,trait)
########################################################################################################
table_s3_unique_index_variants <- table_s3
table_s3_unique_index_variants$trait <- paste0("(", table_s3_unique_index_variants$trait, ")")
table_s3_unique_index_variants$`Gene Prioritization` <- paste0(table_s3_unique_index_variants$`Gene Prioritization`,table_s3_unique_index_variants$trait,sep = "")
table_s3_unique_index_variants_part1 <- table_s3_unique_index_variants[,c(1,3:5)]
table_s3_unique_index_variants_part1$rsid_dbSNP150[is.na(table_s3_unique_index_variants_part1$rsid_dbSNP150)] <- "NA"
table_s3_unique_index_variants_part1 <- aggregate( `Gene Prioritization` ~ .,
                                                   data=table_s3_unique_index_variants_part1,paste,collapse=' | ')
table_s3_unique_index_variants_part1$rsid_dbSNP150[table_s3_unique_index_variants_part1$rsid_dbSNP150 == "NA"] <- NA
table_s3_unique_index_variants_part1 <- arrange(table_s3_unique_index_variants_part1,chr,`pos (build 37)`)

table_s3_unique_index_variants_part2 <- table_s3_unique_index_variants[,c(1,3:4,14:28)]
table_s3_unique_index_variants_part2 <- table_s3_unique_index_variants_part2[!duplicated(table_s3_unique_index_variants_part2),]
table_s3_unique_index_variants_genes <- data.frame(table_s3_unique_index_variants_part1,table_s3_unique_index_variants_part2[,c(4:18)],check.names = F)


table_s3_unique_index_variants_genes$chr <- as.numeric(table_s3_unique_index_variants_genes$chr)
table_s3_unique_index_variants_genes$`pos (build 37)` <- as.numeric(table_s3_unique_index_variants_genes$`pos (build 37)`)

meta_summary <- meta_summary %>% select(chr:MarkerName,LDL_Effect,LDL_StdErr,`LDL_P-value`,
                                        HDL_Effect,HDL_StdErr,`HDL_P-value`,
                                        logTG_Effect,logTG_StdErr,`logTG_P-value`,
                                        nonHDL_Effect,nonHDL_StdErr,`nonHDL_P-value`,
                                        TC_Effect,TC_StdErr,`TC_P-value`)
identical(table_s3_unique_index_variants_genes$rsid_dbSNP150,meta_summary$rsid_dbSNP150)

meta_summary_ldl_lowering <- data.frame(table_s3_unique_index_variants_genes[,c(1:4)],meta_summary[,c(4:5,7:21)],
                                        table_s3_unique_index_variants_genes[,c(5:19)],check.names = F)
meta_summary_ldl_lowering <-meta_summary_ldl_lowering[,-c(22,23,27,28,32,33)]
colnames(meta_summary_ldl_lowering)[c(22:30)] <- c("T2D_Effect","T2D_StdErr","T2D_P-value","CAD_Effect","CAD_StdErr","CAD_P-value","NAFLD_Effect","NAFLD_StdErr","NAFLD_P-value")
meta_summary_ldl_lowering_ind <- inner_join(meta_summary_ldl_lowering,independent_index_variant_unique)

##############################################################

table_s3 <- table_s3 %>% select(rsid_dbSNP150:`METAL Alt`,`METAL Effect`,StdErr, pvalue_neg_log10_METAL,
                                t2d_eur_Beta,t2d_eur_SE,t2d_eur_Pvalue,cad_mvp_Beta,cad_mvp_SE,cad_mvp_pvalue,nafld_tem_Effect,NAFLD_SE,
                                nafld_tem_Pvalue,locus,pp4_lipids_cad,pp4_lipids_t2d)
identical(table_s3$trait,full_data$trait)
identical(table_s3$rsid_dbSNP150,full_data$rsid_dbSNP150)
# Allele 2 is the effect allele
colnames(table_s3)[c(6,7,8,9,10)] <-c("Allele1","Allele2","Lipid_Effect","Lipid_StdErr","Lipid_Neg_Log10_P-value")

table1_a <-table_s3 %>% filter(trait != "HDL") %>% 
  mutate(newAllele1=ifelse(Lipid_Effect > 0, Allele2,Allele1),
         newAllele2=ifelse(Lipid_Effect > 0,Allele1,Allele2),
         Allele1=newAllele1,
         Allele2=newAllele2,
         t2d_eur_Beta=ifelse(Lipid_Effect > 0,-t2d_eur_Beta,t2d_eur_Beta),
         cad_mvp_Beta=ifelse(Lipid_Effect > 0,-cad_mvp_Beta,cad_mvp_Beta),
         nafld_tem_Effect=ifelse(Lipid_Effect > 0,-nafld_tem_Effect,nafld_tem_Effect),
         Lipid_Effect=ifelse(Lipid_Effect > 0,-Lipid_Effect,Lipid_Effect)) %>% select(-newAllele1, -newAllele2)
table1_b <-table_s3 %>% filter(trait == "HDL") %>% 
  mutate(newAllele1=ifelse(Lipid_Effect < 0, Allele2,Allele1),
         newAllele2=ifelse(Lipid_Effect < 0, Allele1,Allele2),
         Allele1=newAllele1,
         Allele2=newAllele2,
         t2d_eur_Beta=ifelse(Lipid_Effect < 0,-t2d_eur_Beta,t2d_eur_Beta),
         cad_mvp_Beta=ifelse(Lipid_Effect < 0,-cad_mvp_Beta,cad_mvp_Beta),
         nafld_tem_Effect=ifelse(Lipid_Effect < 0,-nafld_tem_Effect,nafld_tem_Effect),
         Lipid_Effect=ifelse(Lipid_Effect < 0,-Lipid_Effect,Lipid_Effect)) %>% select(-newAllele1, -newAllele2)
table1 <- rbind(table1_a,table1_b)
colnames(table1)[c(11:19)] <- c("T2D_Effect","T2D_StdErr","T2D_P-value","CAD_Effect","CAD_StdErr","CAD_P-value","NAFLD_Effect","NAFLD_StdErr","NAFLD_P-value")
# Filtering
table1_filtering <- table1 %>% filter(`T2D_P-value` < 1e-06 & `CAD_P-value` < 1e-06)
table1_filtering <- table1_filtering %>% group_by(locus) %>% top_n(1, `Lipid_Neg_Log10_P-value`) %>% ungroup() %>% select(-locus)
table1_filtering <- table1_filtering[table1_filtering$T2D_Effect < 0 & table1_filtering$CAD_Effect < 0, ]
#table1_filtering <- table1_filtering[table1_filtering$pp4_lipids_cad >= 0.8 & table1_filtering$pp4_lipids_t2d >= 0.8, ]
table1_filtering <- data.frame(arrange(table1_filtering,trait,chr,`pos (build 37)`),check.names = F)
table1_ind <- inner_join(table1,independent_index_variant)

# library("xlsx")
# write.xlsx(table1_ind, file = "glgc_tables_1021.xlsx",row.names = F,
#       sheetName = "Table 1 (No Filter)", append = FALSE)

# write.xlsx(table_s2_ind, file = "glgc_tables_1021.xlsx", row.names = F,
#            sheetName="Table S2", append=T)

# write.xlsx(meta_summary_ldl_lowering_ind, file = "glgc_tables_1021.xlsx",row.names = F,
#            sheetName="Table S3", append=T)
#write.xlsx(mouse_genes_list, file = "mouse_model_gene.xlsx",row.names = F,
#      sheetName = "Mouse knockout genes", append = FALSE)