# rm(list = ls())
# 
# genes <- read.delim(file = "high_confidence_gene_wo_mendelian.txt",header = F)
# genes <- genes$V1
# 
# disease <- read.delim(file = "./inputfiles/searchterms.txt",header = F)[1:22,]
# 
# header <- NULL  
# 
# for (i in 1:length(genes)){
#   searchterm <- c(genes[i],disease)
#   write.table(searchterm,file=paste("/rprojectnb2/pelosolab/yuxuan/pubmed/textmine/0508/",genes[i],".txt",sep=""),col.names=F,row.names=F,quote=F)
# }
# 
# 
# for (i in 1:length(genes)){
#   line1 <- paste("python parsesearch.py -d 20210508 -r true -n ",genes[i]," -t ",genes[i],
#                 ".txt -wd /rprojectnb2/pelosolab/yuxuan/pubmed/textmine/0508/ -xd /rprojectnb2/pelosolab/yuxuan/pubmed/",sep = "")
#   line2 <- paste("rm -f /restricted/projectnb/pelosolab/yuxuan/pubmed/*regex.txt")
#   line <- paste(line1,line2,sep = "\n")
#   header <- c(header, line)}
# 
# fileConn<-file("output.txt")
# writeLines(header, fileConn)
# close(fileConn)

library(tidyverse)
dd  <- read_delim("pops_plus_puberrator_regex_output_20220306.txt",delim = "\t",quote = "",col_names = F)
colnames(dd) <- c("Input List","PMID","PublicationDate","Match #","Matching terms","Article Title")

dd <- dd[dd$`Match #` !=1,]
dd <- dd[-1,]

genes <- read.table("../pops_plus.txt")$V1
genes <- tolower(genes)

lipids <- read.delim("lipids.txt",header = F)$V1
lipids <- tolower(lipids)

keep <- function(x){
  any(unlist(strsplit(x[5],split = ", ")) %in% genes) & 
    any(unlist(strsplit(x[5],split = ", ")) %in% lipids) 
}
dd$keep <-as.character(apply(dd, 1, keep))
final <- dd[dd$keep == T,]
final <- unique(final)

terms <- unlist(strsplit(final$`Matching terms`,", "))
count <- as.data.frame(table(terms)) 

genes_count <- count[count$terms %in% genes,]
gene_count_zero <- setdiff(genes,genes_count$terms)
tmp <- data.frame(gene_count_zero,rep(0,23))
colnames(tmp) <- c("terms","Freq")

genes_count <- rbind(genes_count,tmp)
genes_count <- genes_count[order(genes_count$Freq,decreasing = T),]
genes_count$terms <- toupper(genes_count$terms)


final$`Match #` <- as.numeric(final$`Match #`)
final <- final[order(final$`Match #`,decreasing = T),]

write.table(final, "textmine_pops_plus.txt", quote = F, sep = "\t",row.names = F)
write.csv(genes_count,"textmine_pops_plus_count.csv",row.names = F, quote = F)
