library(tidyr)
library(readr)
library(dplyr)
library(coloc)
library(qvalue)


overlap = function(chr, start, end){
        toreturn = sentinel %>% filter(X2 == chr) %>% filter((X3 >= start) & (X3 <= end))
        return(nrow(toreturn))
}


args = commandArgs(trailing=TRUE)
tissue = args[1]
pheno = args[2]

#file containing sentinel variants per phenotype
sentinel = read_delim(paste0(pheno,"_sentinelsnps.txt"),delim="\t",col_names=FALSE)

#now submit coloc analyses for each ld block
j = 1

#dummy p1 and p2, because of legacy issues with coloc_per_block_per_tissue script
p1 = 0.0001
p2 = 0.0001

for(i in 1:nrow(sentinel)){
	print(paste0("Starting analysis for LD block ", i))
	lower = as.numeric(sentinel[i,17])
	upper = as.numeric(sentinel[i,16])
	sentinelsnp = paste0(sentinel[i,2], "_",sentinel[i,3], "_")
	chr = strsplit(sentinelsnp, "_")[[1]][1]
	pos = as.numeric(strsplit(sentinelsnp, "_")[[1]][2])
	print(sentinelsnp)
	print(lower)
	print(upper)
	command = paste0('bsub -M 15000 -R "span[hosts=1] rusage [mem=15960]"', ' -J ', chr, ":", pos, " -o /home/shwetar/glgc/analysis_newdata_jan2020_0.25cmor500kboneitherside_2/redone_june2020/",pheno,"/",tissue,"/error",chr,"_",pos, ".o -e /home/shwetar/glgc/analysis_newdata_jan2020_0.25cmor500kboneitherside_2/redone_june2020/",pheno,"/",tissue,"/error",chr,"_",pos, ".e Rscript /home/shwetar/glgc/analysis_newdata_jan2020_0.25cmor500kboneitherside_2/redone_june2020/",pheno,"/",tissue,"/coloc_per_block_pertissue.R ", chr, " ", pos, " ", sentinelsnp, " ", p1, " ", p2, " ", tissue, " ", as.character(chr), " ", as.numeric(lower), " ", as.numeric(upper))
	print(command)
	system(command)
	j = j + 1
#	break
}

