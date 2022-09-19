library(coloc)
library(qvalue)
library(dplyr)
library(readr)
library(tidyr)

rm(list=ls())
#REQUIRED ARGUMENTS: chr and pos for Sentinel of GWAS 1, name of GWAS1, name of GWAS2

args = commandArgs(trailingOnly=TRUE)
print(args)

chr = args[1]
pos = as.numeric(args[2])

gwas1 = args[3]   #example: HDL
gwas2 = args[4]	  #example: T2D

gwas1file = args[5] 
gwas2file = args[6]

output_path <- args[7]

lower = pos - 100000
upper = pos + 100000

#Read in GWAS 1 file for that chromosome
gwas_1 = read_delim(gwas1file,delim="\t",col_types = cols(
  SNP = col_character(),
  Chr = col_double(),
  Pos = col_double(),
  EA = col_character(),
  NEA = col_character(),
  EAF = col_double(),
  Beta = col_double(),
  SE = col_double(),
  Pvalue = col_double(),
  Neff = col_double()))



#GWAS 2
gwas_2 = read_delim(gwas2file,delim="\t",col_types = cols(
  SNP = col_character(),
  Chr = col_double(),
  Pos = col_double(),
  EA = col_character(),
  NEA = col_character(),
  EAF = col_double(),
  Beta = col_double(),
  SE = col_double(),
  Pvalue = col_double(),
  Neff = col_double()))
gwas_2 = gwas_2  %>% filter(Pos >= lower) %>% filter(Pos <= upper) 

#Now prep for coloc
#merge two gwas datasets
m = merge(gwas_1, gwas_2, by="SNP", all=FALSE)
r1 = which(m$EAF.x == 0)
r2 = which(m$EAF.y == 0)
r3 = which(m$EAF.x == 1)
r4 = which(m$EAF.y == 1)

if(length(union(union(union(r1, r2),r3),r4)) > 0){
  m = m[-union(union(union(r1, r2),r3),r4),]
}

m = m %>% filter(!is.na(Pvalue.x)) %>% filter(!is.na(Pvalue.y))

# Align effect to the same EA and NEA
m$Beta.y <- ifelse(m$EA.x==m$EA.y, m$Beta.y, -m$Beta.y)

print(nrow(m))
dataset1 = list()  #GWAS1
dataset2 = list()  #GWAS2

dataset1$SNP = m$SNP
dataset1$pvalues = m$Pvalue.x
dataset1$N = m$Neff.x
dataset1$MAF = pmin(m$EAF.x, 1 - m$EAF.x)
dataset1$beta = m$Beta.x
dataset1$varbeta = (m$SE.x)^2
dataset1$type = "quant"

dataset2$SNP = m$SNP
dataset2$pvalues = m$Pvalue.y
dataset2$N = m$Neff.y
dataset2$MAF = pmin(m$EAF.y, 1 - m$EAF.y)
dataset2$type = "cc"
#dataset2$s = 74124/824006 # T2D
dataset2$s = 243392/849686 # CAD

#run coloc

if(nrow(m) == 0){
  next;
}
#       COLOC = coloc.abf(dataset1, dataset2, p1=p1, p2=p2,p12=p12)
COLOC = coloc.abf(dataset1, dataset2)
coloc_results = tibble(snp = paste0(chr,":",pos), chr=chr, start=lower, end = upper, nsnps =COLOC$summary[1], PP.H0.abf= COLOC$summary[2],
                       PP.H1.abf = COLOC$summary[3], PP.H2.abf = COLOC$summary[4], PP.H3.abf = COLOC$summary[5],
                       PP.H4.abf = COLOC$summary[6])


write.table(coloc_results, file=paste0(output_path,"coloc_",chr,":",pos,"_", gwas1,"_",gwas2,".txt"),row.names=F,col.names=F,sep="\t",quote=F)

