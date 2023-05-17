remove.packages('rlang')
remove.packages('Rcpp')
install.packages('Seurat',dependencies=T,source=T)

library(vcfR)
library(dplyr)
library(bio3d) 
#library(Bios2cor) 
library(openxlsx)

install.packages("BiocManager")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("trackViewer")
library(trackViewer)

#clear memory
rm()
rm(list=ls()) 

#set working directory 
setwd("ipt_data")

#=================================CATR===========================================
# Input the SNP vcf files.
vcf_file<-"CATR_HYHYVDSX2_raw_snp.vcf"
vcf_snps <- read.vcfR(vcf_file, verbose = FALSE)

# Input the Indel vcf files.
vcf_file<-"CATR_HYHYVDSX2_raw_indel.vcf"
vcf_indel <- read.vcfR(vcf_file, verbose = FALSE)

#transfer vcf into list
table_snps <- vcfR2tidy(vcf_snps,single_frame = FALSE)
table_indel <- vcfR2tidy(vcf_indel,single_frame = FALSE)

#transfer list into frame
tableFrame_snps_CATR <- as.data.frame(table_snps$fix, row.names = TRUE)
tableFrame_indel_CATR <- as.data.frame(table_indel$fix, row.names = TRUE)
#=================================CATY===========================================
# Input the SNP vcf files.
vcf_file<-"CATY_HYMGMDSX2_raw_snp.vcf"
vcf_snps <- read.vcfR(vcf_file, verbose = FALSE)
``
# Input the Indel vcf files.
vcf_file<-"CATY_HYMGMDSX2_raw_indel.vcf"
vcf_indel <- read.vcfR(vcf_file, verbose = FALSE)

#transfer vcf into list
table_snps <- vcfR2tidy(vcf_snps,single_frame = FALSE)
table_indel <- vcfR2tidy(vcf_indel,single_frame = FALSE)

#transfer list into frame
tableFrame_snps_CATY <- as.data.frame(table_snps$fix, row.names = TRUE)
tableFrame_indel_CATY <- as.data.frame(table_indel$fix, row.names = TRUE)
#=================================COFC===========================================
# Input the SNP vcf files.
vcf_file<-"COFC5_H3YFFDSX2_raw_snp.vcf"
vcf_snps <- read.vcfR(vcf_file, verbose = FALSE)

# Input the Indel vcf files.
vcf_file<-"COFC5_H3YFFDSX2_raw_indel.vcf"
vcf_indel <- read.vcfR(vcf_file, verbose = FALSE)

#transfer vcf into list
table_snps <- vcfR2tidy(vcf_snps,single_frame = FALSE)
table_indel <- vcfR2tidy(vcf_indel,single_frame = FALSE)

#transfer list into frame
tableFrame_snps_COFC <- as.data.frame(table_snps$fix, row.names = TRUE)
tableFrame_indel_COFC <- as.data.frame(table_indel$fix, row.names = TRUE)
#=================================COFT===========================================
# Input the SNP vcf files.
vcf_file<-"COFT2_H3YFFDSX2_raw_snp.vcf"
vcf_snps <- read.vcfR(vcf_file, verbose = FALSE)

# Input the Indel vcf files.
vcf_file<-"COFT2_H3YFFDSX2_raw_indel.vcf"
vcf_indel <- read.vcfR(vcf_file, verbose = FALSE)

#transfer vcf into list
table_snps <- vcfR2tidy(vcf_snps,single_frame = FALSE)
table_indel <- vcfR2tidy(vcf_indel,single_frame = FALSE)

#transfer list into frame
tableFrame_snps_COFT <- as.data.frame(table_snps$fix, row.names = TRUE)
tableFrame_indel_COFT <- as.data.frame(table_indel$fix, row.names = TRUE)
#=================================G351===========================================
# Input the SNP vcf files.
vcf_file<-"G351_HYHYVDSX2_raw_snp.vcf"
vcf_snps <- read.vcfR(vcf_file, verbose = FALSE)

# Input the Indel vcf files.
vcf_file<-"G351_HYHYVDSX2_raw_indel.vcf"
vcf_indel <- read.vcfR(vcf_file, verbose = FALSE)

#transfer vcf into list
table_snps <- vcfR2tidy(vcf_snps,single_frame = FALSE)
table_indel <- vcfR2tidy(vcf_indel,single_frame = FALSE)

#transfer list into frame
tableFrame_snps_G351 <- as.data.frame(table_snps$fix, row.names = TRUE)
tableFrame_indel_G351 <- as.data.frame(table_indel$fix, row.names = TRUE)
#=================================HDT===========================================
# Input the SNP vcf files.
vcf_file<-"HDT_HYMGMDSX2_raw_snp.vcf"
vcf_snps <- read.vcfR(vcf_file, verbose = FALSE)

# Input the Indel vcf files.
vcf_file<-"HDT_HYMGMDSX2_raw_indel.vcf"
vcf_indel <- read.vcfR(vcf_file, verbose = FALSE)

#transfer vcf into list
table_snps <- vcfR2tidy(vcf_snps,single_frame = FALSE)
table_indel <- vcfR2tidy(vcf_indel,single_frame = FALSE)

#transfer list into frame
tableFrame_snps_HDT <- as.data.frame(table_snps$fix, row.names = TRUE)
tableFrame_indel_HDT <- as.data.frame(table_indel$fix, row.names = TRUE)
#=============================================TPS1========================================================

source("K:/SNP_process/compare_with_ref/code/variants_statistic_activity/TPS1_GENE_VARIANTS_ANALYSIS.R")
TPS1_GENE_VARIANTS_ANALYSIS(0,0,tableFrame_snps_CATR,tableFrame_indel_CATR)

source("L:/SNP_process/compare_with_ref/code/variants_statistic_activity/TPS1_GENE_VARIANTS_ANALYSIS.R")
TPS1_GENE_VARIANTS_ANALYSIS(0,0,tableFrame_snps_CATY,tableFrame_indel_CATY)

source("L:/SNP_process/compare_with_ref/code/variants_statistic_activity/TPS1_GENE_VARIANTS_ANALYSIS.R")
TPS1_GENE_VARIANTS_ANALYSIS(0,0,tableFrame_snps_COFC,tableFrame_indel_COFC)

source("L:/SNP_process/compare_with_ref/code/variants_statistic_activity/TPS1_GENE_VARIANTS_ANALYSIS.R")
TPS1_GENE_VARIANTS_ANALYSIS(0,0,tableFrame_snps_COFT,tableFrame_indel_COFT)

source("L:/SNP_process/compare_with_ref/code/variants_statistic_activity/TPS1_GENE_VARIANTS_ANALYSIS.R")
TPS1_GENE_VARIANTS_ANALYSIS(0,0,tableFrame_snps_G351,tableFrame_indel_G351)

source("L:/SNP_process/compare_with_ref/code/variants_statistic_activity/TPS1_GENE_VARIANTS_ANALYSIS.R")
TPS1_GENE_VARIANTS_ANALYSIS(0,0,tableFrame_snps_HDT,tableFrame_indel_HDT)
#=============================================TPS2========================================================

source("L:/SNP_process/compare_with_ref/code/variants_statistic_activity/TPS2_GENE_VARIANTS_ANALYSIS.R")
TPS2_GENE_VARIANTS_ANALYSIS(1,0,tableFrame_snps_CATR,tableFrame_indel_CATR)

source("L:/SNP_process/compare_with_ref/code/variants_statistic_activity/TPS2_GENE_VARIANTS_ANALYSIS.R")
TPS2_GENE_VARIANTS_ANALYSIS(1,0,tableFrame_snps_CATY,tableFrame_indel_CATY)

source("L:/SNP_process/compare_with_ref/code/variants_statistic_activity/TPS2_GENE_VARIANTS_ANALYSIS.R")
TPS2_GENE_VARIANTS_ANALYSIS(1,0,tableFrame_snps_COFC,tableFrame_indel_COFC)

source("L:/SNP_process/compare_with_ref/code/variants_statistic_activity/TPS2_GENE_VARIANTS_ANALYSIS.R")
TPS2_GENE_VARIANTS_ANALYSIS(1,0,tableFrame_snps_COFT,tableFrame_indel_COFT)

source("L:/SNP_process/compare_with_ref/code/variants_statistic_activity/TPS2_GENE_VARIANTS_ANALYSIS.R")
TPS2_GENE_VARIANTS_ANALYSIS(1,0,tableFrame_snps_G351,tableFrame_indel_G351)

source("L:/SNP_process/compare_with_ref/code/variants_statistic_activity/TPS2_GENE_VARIANTS_ANALYSIS.R")
TPS2_GENE_VARIANTS_ANALYSIS(1,0,tableFrame_snps_HDT,tableFrame_indel_HDT)



