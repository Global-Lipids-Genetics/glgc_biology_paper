# library(data.table)
# 
# LDL <- fread("lipid_results/glgc_2021/glgc_gwas_meta-analysis/trans-ancestry/meta-analysis_AFR_EAS_EUR_HIS_SAS_LDL_INV_ALL_with_N_1.gz")
# LDL$POOLED_REF_AF <- 1- LDL$POOLED_ALT_AF
# 
# fwrite(LDL,file=paste0("glgc/meta-analysis_AFR_EAS_EUR_HIS_SAS_LDL_INV_ALL_with_N_1.gz"),sep = " ",na = "NA",quote = F)
# 
# HDL <- fread("lipid_results/glgc_2021/glgc_gwas_meta-analysis/trans-ancestry/meta-analysis_AFR_EAS_EUR_HIS_SAS_HDL_INV_ALL_with_N_1.gz")
# HDL$POOLED_REF_AF <- 1- HDL$POOLED_ALT_AF
# 
# fwrite(HDL,file=paste0("glgc/meta-analysis_AFR_EAS_EUR_HIS_SAS_HDL_INV_ALL_with_N_1.gz"),sep = " ",na = "NA",quote = F)
# 
# logTG <- fread("lipid_results/glgc_2021/glgc_gwas_meta-analysis/trans-ancestry/meta-analysis_AFR_EAS_EUR_HIS_SAS_logTG_INV_ALL_with_N_1.gz")
# logTG$POOLED_REF_AF <- 1- logTG$POOLED_ALT_AF
# 
# fwrite(logTG,file=paste0("glgc/meta-analysis_AFR_EAS_EUR_HIS_SAS_logTG_INV_ALL_with_N_1.gz"),sep = " ",na = "NA",quote = F)
# 
# nonHDL <- fread("lipid_results/glgc_2021/glgc_gwas_meta-analysis/trans-ancestry/meta-analysis_AFR_EAS_EUR_HIS_SAS_nonHDL_INV_ALL_with_N_1.gz")
# nonHDL$POOLED_REF_AF <- 1- nonHDL$POOLED_ALT_AF
# 
# fwrite(nonHDL,file=paste0("glgc/meta-analysis_AFR_EAS_EUR_HIS_SAS_nonHDL_INV_ALL_with_N_1.gz"),sep = " ",na = "NA",quote = F)
# 
# TC <- fread("lipid_results/glgc_2021/glgc_gwas_meta-analysis/trans-ancestry/meta-analysis_AFR_EAS_EUR_HIS_SAS_TC_INV_ALL_with_N_1.gz")
# TC$POOLED_REF_AF <- 1- TC$POOLED_ALT_AF
# 
# fwrite(TC,file=paste0("glgc/meta-analysis_AFR_EAS_EUR_HIS_SAS_TC_INV_ALL_with_N_1.gz"),sep = " ",na = "NA",quote = F)

library(ggplot2)
library(ggpubr)

lipid <- c("HDL","logTG","LDL","TC","nonHDL")

for (i in lipid){
  
  dat <- read.delim(paste0("geneAssoc_",i,".celltype.txt"))
  x <- stringr::str_split(dat$TissueName,pattern = "-")
  dat$Categories <- unlist(lapply(x, `[[`, 1))
  p1 <- ggplot(dat, aes(x = reorder(TissueName, -Log.p.), y=Log.p., fill=Categories)) + 
    geom_bar(stat = "identity") +
    geom_hline(yintercept=-log10(0.05/54), linetype="dashed", color = "red") +
    guides(fill=guide_legend(ncol=1)) +
    theme_classic() + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1,hjust=1)) + 
    labs(title = paste0(i,"_DESE_GTEx_Gene"),x = "Tissue/Cell Types",y="-log10(P-value)" )
  
  
  dat <- read.delim(paste0("transcriptAssoc_",i,".celltype.txt"))
  x <- stringr::str_split(dat$TissueName,pattern = "-")
  dat$Categories <- unlist(lapply(x, `[[`, 1))
  p2 <- ggplot(dat, aes(x = reorder(TissueName, -Log.p.), y=Log.p., fill=Categories)) + 
    geom_bar(stat = "identity") + 
    geom_hline(yintercept=-log10(0.05/54), linetype="dashed", color = "red") +
    guides(fill=guide_legend(ncol=1)) +
    theme_classic() + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1,hjust=1)) + 
    labs(title = paste0(i,"_DESE_GTEx_Transcript"),x = "Tissue/Cell Types",y="-log10(P-value)" )
  
  p <- ggarrange(p1, p2, ncol=1, nrow=2, common.legend = TRUE, legend="right")
  
  ggsave(filename = paste0(i,".png"),plot = p, width = 18, height = 12)
  
}

library(tidyverse)
lipid <- c("HDL","logTG","LDL","TC","nonHDL")

HDL <- read.delim(paste0("geneAssoc_HDL.celltype.txt"))
HDL <- HDL[,1:2]
HDL <- HDL %>% arrange(TissueName)

logTG <- read.delim(paste0("geneAssoc_logTG.celltype.txt"))
logTG <- logTG[,1:2]
logTG <- logTG %>% arrange(TissueName)


LDL <- read.delim(paste0("geneAssoc_LDL.celltype.txt"))
LDL <- LDL[,1:2]
LDL <- LDL %>% arrange(TissueName)

TC <- read.delim(paste0("geneAssoc_TC.celltype.txt"))
TC <- TC[,1:2]
TC <- TC %>% arrange(TissueName)

nonHDL <- read.delim(paste0("geneAssoc_nonHDL.celltype.txt"))
nonHDL <- nonHDL[,1:2]
nonHDL <- nonHDL %>% arrange(TissueName)

dat <- cbind(HDL,logTG[,2],LDL[,2],TC[,2],nonHDL[,2])
colnames(dat) <- c("TissueName","Pvalue_HDL","Pvalue_logTG","Pvalue_LDL","Pvalue_TC","Pvalue_nonHDL")

write.csv(dat,file = "geneAssoc.celltype.csv")

HDL <- read.delim(paste0("transcriptAssoc_HDL.celltype.txt"))
HDL <- HDL[,1:2]
HDL <- HDL %>% arrange(TissueName)

logTG <- read.delim(paste0("transcriptAssoc_logTG.celltype.txt"))
logTG <- logTG[,1:2]
logTG <- logTG %>% arrange(TissueName)


LDL <- read.delim(paste0("transcriptAssoc_LDL.celltype.txt"))
LDL <- LDL[,1:2]
LDL <- LDL %>% arrange(TissueName)

TC <- read.delim(paste0("transcriptAssoc_TC.celltype.txt"))
TC <- TC[,1:2]
TC <- TC %>% arrange(TissueName)

nonHDL <- read.delim(paste0("transcriptAssoc_nonHDL.celltype.txt"))
nonHDL <- nonHDL[,1:2]
nonHDL <- nonHDL %>% arrange(TissueName)

dat <- cbind(HDL,logTG[,2],LDL[,2],TC[,2],nonHDL[,2])
colnames(dat) <- c("TissueName","Pvalue_HDL","Pvalue_logTG","Pvalue_LDL","Pvalue_TC","Pvalue_nonHDL")

write.csv(dat,file = "transcriptAssoc.celltype.csv")
