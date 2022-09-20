#Step 1
magma --bfile /1KG/ALL.phase3 --gene-annot /POPS_resources/magma_0kb.genes.annot --pval /chrALL_LDL_MR-MEGA.POPS.input.txt ncol=N --gene-model snp-wise=mean --out /LDL_MR-MEGA_magma

#Step 2

python /POPS/pops.feature_selection.py --features /POPS_resources/PoPS.features.txt.gz --gene_results /LDL_MR-MEGA_magma --out /LDL_MR-MEGA_pops_step2

#Step 3

python /POPS/pops.predict_scores.py --gene_loc /POPS_resources/gene_loc.txt --gene_results /LDL_MR-MEGA_magma --features /POPS_resources/PoPS.features.txt.gz --selected_features /LDL_MR-MEGA_pops_step2.features --control_features /POPS_resources/control.features --chromosome ${SGE_TASK_ID} --out /LDL_MR-MEGA_POPS_chr${SGE_TASK_ID}

#################################################
#Need to get all genes within 500kb window of the index variants
R

library(IRanges)
index_variants<-read.table("/trans_ethnic_index_variants.txt",as.is=TRUE,header=TRUE, stringsAsFactors=FALSE)
index1<-subset(index_variants,!(duplicated(rsid)))
biomart<-read.table("/mart_export.txt",as.is=TRUE,header=TRUE, stringsAsFactors=FALSE)
CHR<-1:22

for (i in 1:length(CHR)) {

index2<-subset(index1,chr==CHR[i])
biomart2<-subset(biomart,Chromosome==CHR[i])
snp.gr<-IRanges(index2$pos)
mcols(snp.gr)$rsid <- index2$rsid

annot.gr<-IRanges(biomart2$start,biomart2$end)
mcols(annot.gr)$GeneID <- biomart2$GeneID

ov<-findOverlaps(snp.gr, annot.gr, maxgap=500000)
hits1<-as.data.frame(annot.gr[subjectHits(ov)])
hits2<-as.data.frame(snp.gr[queryHits(ov)])

hits<-cbind(hits2,hits1)
hits<-hits[,c(4,8)]

assign(paste0("chr",CHR[i]),hits)

}

all_hits<-rbind(chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22)
all_hits<-merge(all_hits,biomart[,c(1,2)],by="GeneID")
all_hits$Symbol<-sQuote(all_hits$Symbol)
write.table(all_hits,file="/trans_ethnic_index_variants_biomart_genes.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote = FALSE)


####Need to get the best gene in the locus from the POPS results

chr<-1:22
trait<-c("LDL","HDL","TC","logTG","nonHDL")

for (i in 1:length(chr)) {
for (j in 1:length(trait)) {

pops<-read.table(paste0("/",trait[j],"_MR-MEGA_POPS_chr",chr[i],".",chr[i],".results"),as.is=TRUE,header=TRUE, stringsAsFactors=FALSE)

assign(paste0("chr",chr[i],".",trait[j]),pops)

}
}

ldl_pops<-rbind(chr1.LDL,chr2.LDL,chr3.LDL,chr4.LDL,chr5.LDL,chr6.LDL,chr7.LDL,chr8.LDL,chr9.LDL,chr10.LDL,chr11.LDL,chr12.LDL,chr13.LDL,chr14.LDL,chr15.LDL,chr16.LDL,chr17.LDL,chr18.LDL,chr19.LDL,chr20.LDL,chr21.LDL,chr22.LDL)
hdl_pops<-rbind(chr1.HDL,chr2.HDL,chr3.HDL,chr4.HDL,chr5.HDL,chr6.HDL,chr7.HDL,chr8.HDL,chr9.HDL,chr10.HDL,chr11.HDL,chr12.HDL,chr13.HDL,chr14.HDL,chr15.HDL,chr16.HDL,chr17.HDL,chr18.HDL,chr19.HDL,chr20.HDL,chr21.HDL,chr22.HDL)
tc_pops<-rbind(chr1.TC,chr2.TC,chr3.TC,chr4.TC,chr5.TC,chr6.TC,chr7.TC,chr8.TC,chr9.TC,chr10.TC,chr11.TC,chr12.TC,chr13.TC,chr14.TC,chr15.TC,chr16.TC,chr17.TC,chr18.TC,chr19.TC,chr20.TC,chr21.TC,chr22.TC)
logtg_pops<-rbind(chr1.logTG,chr2.logTG,chr3.logTG,chr4.logTG,chr5.logTG,chr6.logTG,chr7.logTG,chr8.logTG,chr9.logTG,chr10.logTG,chr11.logTG,chr12.logTG,chr13.logTG,chr14.logTG,chr15.logTG,chr16.logTG,chr17.logTG,chr18.logTG,chr19.logTG,chr20.logTG,chr21.logTG,chr22.logTG)
nonhdl_pops<-rbind(chr1.nonHDL,chr2.nonHDL,chr3.nonHDL,chr4.nonHDL,chr5.nonHDL,chr6.nonHDL,chr7.nonHDL,chr8.nonHDL,chr9.nonHDL,chr10.nonHDL,chr11.nonHDL,chr12.nonHDL,chr13.nonHDL,chr14.nonHDL,chr15.nonHDL,chr16.nonHDL,chr17.nonHDL,chr18.nonHDL,chr19.nonHDL,chr20.nonHDL,chr21.nonHDL,chr22.nonHDL)

write.table(ldl_pops,file="/LDL_trans_ethnic_POPS_all_scores_genes.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote = FALSE)
write.table(hdl_pops,file="/HDL_trans_ethnic_POPS_all_scores_genes.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote = FALSE)
write.table(tc_pops,file="/TC_trans_ethnic_POPS_all_scores_genes.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote = FALSE)
write.table(logtg_pops,file="/logTG_trans_ethnic_POPS_all_scores_genes.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote = FALSE)
write.table(nonhdl_pops,file="/nonHDL_trans_ethnic_POPS_all_scores_genes.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote = FALSE)

ldl_pops$trait<-"LDL"
hdl_pops$trait<-"HDL"
tc_pops$trait<-"TC"
logtg_pops$trait<-"logTG"
nonhdl_pops$trait<-"nonHDL"

all_pops<-rbind(ldl_pops,hdl_pops,tc_pops,logtg_pops,nonhdl_pops)

index_variants<-read.table("/trans_ethnic_index_variants.txt",as.is=TRUE,header=TRUE, stringsAsFactors=FALSE)

index_variants_biomart<-read.table("/trans_ethnic_index_variants_biomart_genes.txt",as.is=TRUE,header=TRUE, stringsAsFactors=FALSE)

index_variants2<-merge(index_variants,index_variants_biomart,by="rsid", all.x=TRUE)
index_variants2[is.na(index_variants2)] <- "NA"

index_variants_pops<-merge(index_variants2,all_pops,by.x=c("GeneID","trait"),by.y=c("ENSGID","trait"),all.x=TRUE)
library(dplyr)
lipid<-c("LDL","HDL","TC","logTG","nonHDL")

for (i in 1:length(lipid)) {

select1<-subset(index_variants_pops,trait==lipid[i])
select2<-select1[order( -select1$Score),]
select3<-select2 %>% group_by(rsid) %>% top_n(2, Score)
select3<-as.data.frame(select3)
second<-subset(select3,duplicated(rsid))
select4<-subset(select3,!duplicated(rsid))
select<-merge(select4,second,by="rsid",all.x=TRUE)
select[is.na(select)] <- "NA"
select<-select[,c(3,1,4,5,2,7,8,9,14,15)]
colnames(select)<-c("trait","rsid","chr","pos","GeneID_first","SYMBOL_first","PoPS_Score_first","GeneID_second","SYMBOL_second","PoPS_Score_second")
select$SYMBOL_first<-sQuote(select$SYMBOL_first)
select$SYMBOL_second<-sQuote(select$SYMBOL_second)

assign(paste0(lipid[i]),select)

}

pops_genes_index_variants<-rbind(LDL,HDL,TC,logTG,nonHDL)

write.table(pops_genes_index_variants,file="/trans_ethnic_index_variants_POPS_genes_TopTwo.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote = FALSE)


rm(list=ls(all=TRUE))




