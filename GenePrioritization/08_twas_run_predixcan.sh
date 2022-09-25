#The following are the two steps in running PrediXcan using mashr models for LDL in adipose. The same code is used for all five lipid phenotypes for all GTEx tissues

#PrediXcan for LDL for Adipose_Subcutaneous

python gwas_parsing.py -gwas_file LDL.txt.gz -liftover hg19ToHg38.over.chain.gz -snp_reference_metadata reference_panel_1000G/variant_metadata.txt.gz METADATA -output_column_map markername variant_id -output_column_map noneffect_allele non_effect_allele -output_column_map effect_allele effect_allele -output_column_map beta effect_size -output_column_map chr chromosome --chromosome_format -output_column_map bp_hg19 position -output_column_map effect_allele_freq frequency -output_column_map se_dgc standard_error -output_order variant_id panel_variant_id chromosome position effect_allele non_effect_allele frequency pvalue zscore effect_size standard_error sample_size n_cases -output LDLparsed.txt.gz

python SPrediXcan.py --gwas_file LDLparsed.txt.gz --snp_column panel_variant_id --effect_allele_column effect_allele --non_effect_allele_column non_effect_allele --zscore_column zscore --model_db_path data/models/eqtl/mashr/mashr_Adipose_Subcutaneous.db --covariance data/models/eqtl/mashr/mashr_Adipose_Subcutaneous.txt.gz --keep_non_rsid --additional_output --model_db_snp_key varID --throw --output_file .//spredix
can/eqtl/LDL_Adipose_Subcutaneous.csv


