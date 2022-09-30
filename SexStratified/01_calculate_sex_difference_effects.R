#!/usr/bin/Rscript

library(data.table)
library(dplyr)
library(argparse)

parser <- ArgumentParser()
parser$add_argument("-f", "--female_file")
parser$add_argument("-m", "--male_file")
parser$add_argument("-o", "--output_dir")
args <- parser$parse_args()

men <- fread(cmd=paste("zgrep -v '^##'", args$male_file))
women <- fread(cmd=paste("zgrep -v '^##'", args$female_file))

men$SNP <- paste(men$CHROM, men$POS, pmin(men$REF, men$ALT), pmax(men$REF, men$ALT), sep=":")
men$SE <- 1/men$SQRT_V_STAT
men <- men[,c("SNP", "REF", "ALT", "N_INFORMATIVE", "ALT_EFFSIZE", "SE")]
names(men) <- c("SNP", "men_REF", "men_ALT", "men_N_INFORMATIVE", "men_ALT_EFFSIZE", "men_SE")

women$SNP <- paste(women$CHROM, women$POS, pmin(women$REF, women$ALT), pmax(women$REF, women$ALT), sep=":")
women$SE <- 1/women$SQRT_V_STAT
women <- women[,c("SNP", "REF", "ALT", "N_INFORMATIVE", "ALT_EFFSIZE", "SE")]
names(women) <- c("SNP", "women_REF", "women_ALT", "women_N_INFORMATIVE", "women_ALT_EFFSIZE", "women_SE")

sex_combined <- left_join(men, women, by=c("SNP"="SNP"))

#Make sure effect alleles match (it looks like they do, but we'll check beta's just to be sure)
sex_combined$women_ALT_EFFSIZE <- ifelse(sex_combined$men_REF == sex_combined$women_REF & sex_combined$men_ALT == sex_combined$women_ALT, sex_combined$women_ALT_EFFSIZE, -sex_combined$women_ALT_EFFSIZE)

#Calculate overall correlation in effect sizes
r <- cor(sex_combined$women_ALT_EFFSIZE, sex_combined$men_ALT_EFFSIZE, method="pearson", use="complete.obs")

sex_combined$Z_diff = (sex_combined$"men_ALT_EFFSIZE" - sex_combined$"women_ALT_EFFSIZE")/(sex_combined$"men_SE"^2 + sex_combined$"women_SE"^2 - 2*r*sex_combined$"men_SE"*sex_combined$"women_SE")^0.5
sex_combined$p_diff <- 2*pnorm(-abs(sex_combined$Z_diff))
sex_combined$N <- sex_combined$women_N_INFORMATIVE + sex_combined$men_N_INFORMATIVE
sex_combined <- na.omit(sex_combined[,c("SNP", "men_REF", "men_ALT", "N", "Z_diff", "p_diff")])


write.table(sex_combined, paste0(args$output_dir, "/", sub(".gz", "", basename(args$male_file))), col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)


