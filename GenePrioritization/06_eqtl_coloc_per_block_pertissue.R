library(coloc)
library(qvalue)
library(dplyr)
library(readr)
library(tidyr)


#Functions from Lars to work with very small values
log10toP <- function(log10P){
  log10P <- abs(as.numeric(log10P))
  if(is.na(log10P)) return(NA)
  if(log10P > 300){
    part1 <- log10P%/%100*100
    part2 <- log10P-part1
    if(part2 != 0){
      P <- format(signif(10^-part2,3), scientific = T)
      P <- paste(as.numeric(gsub("e-.+","",P)),"e-",as.numeric(gsub(".+-","",P),sep="")+part1,sep="")
    } else {
      P <- paste("1e-",part1,sep="")
    }
  } else {
    P <- signif(10^-log10P,3)
  }
  return(as.character(P))
}


pheno = "LDL"

args = commandArgs(trailingOnly=TRUE)
print(args)
chr = args[1]
pos = args[2]
sentinel = args[3]

p1 = as.numeric(args[4])
p2 = as.numeric(args[5])
tissue = args[6]
ldblock = c(args[7], as.numeric(args[8]), as.numeric(args[9]))
lower = as.numeric(ldblock[2])
upper = as.numeric(ldblock[3])

print(ldblock)
print(chr)
pos = as.numeric(pos)

#gene coordinates: from gencode v24
#file columns: chr, start, end, gencode IDs
geneids = read_delim("gencode_geneids.txt",comment="#",delim="\t",col_names=FALSE)

#file from GWAS summary stats containing columns: 
#SNP     Chr     Pos     EA      NEA     EAF     Beta    StdErr  Pvalue  N       chisq_association       ndf_association
gwas = read_delim(paste0("forcoloc_",pheno,"_final_",chr,".txt"),delim="\t")
gwas$SNP = gsub(":","_",gwas$SNP)
gwas$EA = toupper(gwas$EA)
gwas$NEA = toupper(gwas$NEA)
gwas$SE = gwas$StdErr

gwas$neg_log10_pvalue <- -pchisq(gwas$chisq_association, df=gwas$ndf_association, lower.tail=F, log.p=T)/log(10)
gwas$pvalue <- as.numeric(sapply(gwas$neg_log10_pvalue,log10toP))
gwas$Pvalue = gwas$pvalue

#find genes within 1Mb
cisgenes = geneids %>% filter(X1 == chr) %>% filter(((X2 > pos - 1000000) & (X2 < pos + 1000000)) | ((X3 > pos - 1000000) & (X3 < pos + 1000000)) | ((X2 <  pos - 1000000) & (X3 > pos + 1000000))) %>% select(X4) %>% distinct()

print(cisgenes)
cisgenes$X4 = gsub("\\..*","",cisgenes$X4)
egenes = read.table(gzfile(paste0("/project/chrbrolab/gtex/results/GTEx_Analysis_v8_eQTL/GTEx_Analysis_v8_eQTL/",tissue,".v8.egenes.txt.gz")), sep="\t",header=T)
egenes[,1] = as.character(egenes[,1])
egenes = as.tbl(egenes)
egenes = egenes %>% filter(qval < 0.05) %>% select(gene_id) %>% distinct()
egenes$gene_id = gsub("\\..*","", egenes$gene_id)
egenes = inner_join(egenes, cisgenes, by=c('gene_id'='X4'))

#get sample size for eqtl in that tissue
command = paste0("cat GTEx_Analysis_v8_eQTL/GTEx_Analysis_v8_eQTL_covariates/",tissue,".v8.covariates.txt | awk '{print NF-1}' | head -1")
num_samples_expn = as.numeric(system(command, intern=TRUE))

coloc_results = tibble(chr="", start="", end="", snp = "", gene="", pval1 = "", pval2 = "", pval3="",pval4="",numsnp="")
for(gene in egenes$gene_id){
  print(gene)
  #read in eqtl data for egenes
  files = list.files(paste0("GTEx_Analysis_v8_eQTL/GTEx_Analysis_v8_eQTL_all_associations_per_gene/",tissue))
  infile = paste0("GTEx_Analysis_v8_eQTL/GTEx_Analysis_v8_eQTL_all_associations_per_gene/",tissue,"/" ,grep(gene, files,value=TRUE)[2])
  eqtl = read_delim(infile,delim="\t")
  eqtl$variant = gsub("_[A-Z]*_[A-Z]*_b38","",eqtl$variant_id)
  
  #merge with gwas data
  m = inner_join(gwas, eqtl, by=c("SNP" = "variant"))
  r1 = which(m$maf == 0)
  r2 = which(m$EAF == 0)
  
  if(length(union(r1, r2)) > 0){
    m = m[-union(r1,r2),]
  }
  
  m = m %>% filter(Pos >= lower) %>% filter(Pos <= upper) %>% filter(!is.na(Pvalue))
  
  if(min(m$Pvalue) == 0){
    m$Pvalue[which(m$Pvalue == 0)] = 1e-320
  }	
  if(length(which(m$Pvalue < 1e-320)) > 0){
    m$Pvalue[which(m$Pvalue < 1e-320)] = 1e-320
  }
  
  print(nrow(m))
  dataset1 = list()
  dataset2 = list()
  
  dataset1$pvalues = m$Pvalue
  dataset1$N = mean(m$N,na.rm=TRUE)
  
  dataset1$MAF = pmin(m$EAF, 1 - m$EAF)
  
  dataset1$type = "quant"
  
  dataset2$pvalues = m$pval_nominal
  dataset2$N = num_samples_expn
  dataset2$MAF = m$maf
  dataset2$beta = m$slope
  dataset2$varbeta = m$slope_se * m$slope_se
  dataset2$type = "quant"
  
  #run coloc
  
  if(nrow(m) == 0){
    next;
  }
  
  COLOC = coloc.abf(dataset1, dataset2)
  coloc_results = add_row(coloc_results, chr=chr, start=as.numeric(ldblock[2]), end = as.numeric(ldblock[3]), snp = sentinel, gene = gene, pval1 = COLOC$summary[3], pval2 = COLOC$summary[4], pval3 = COLOC$summary[5],pval4 = COLOC$summary[6], numsnp = nrow(m))
}

write.table(coloc_results, file=paste0(pheno,"/",tissue,"/egenes_noparameters_coloc_",chr,"_",pos,".txt"),row.names=F,col.names=F,sep="\t",quote=F)
