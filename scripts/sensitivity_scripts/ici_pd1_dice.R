# Load required packages
library(tidyverse)
library(data.table)
library(TwoSampleMR)
library(here)
readRenviron(here("config.env"))

datadir <- Sys.getenv("datadir")
setwd(datadir)
setwd("ici_analyses/UKBB/eQTL/eQTL_DICE_")

protein <- c("PD1")
window <- 500
pval <- 5e-06
r <- 0.3
measure <- c("survival")
survivals <- c("breast", "colorectal", "lung", "melanoma", "ovarian", "prostate")

## PDCD1
DICE_PDCD1 <- data.table::fread("PDCD1_SNPs_eQTL_DICE_130625.csv", sep=",", header=T)
colnames(DICE_PDCD1) <- c("rsid", "Cell_type", "Pos", "Padj_eQTL", "Effect_size", "GWAS", "pieQTL")
DICE_PDCD1$CHR <- "2"
path <- c("PDCD1_Q15116_OID21396_v1_Oncology_COMBINED")

p_dis <- gsub('-', '', pval)

exp_dat_DICE_PDCD1_all <- data.frame()

for (survival in survivals) {

### Generate directory to store output files for combination of instrument selection thresholds in
dir <- paste0(path, "/", window, "kb_p", p_dis, "_r", r, "_240925")

exp_dat <- TwoSampleMR::read_exposure_data(
  filename=paste0("../../", protein, "/",dir,"/",protein,survival, measure, window,"kb_p",p_dis,"r",r,"inst_100525.csv"),
  sep = ",",
  snp_col = "rsid",
  beta_col = "BETA",
  se_col = "SE",
  eaf_col = "A1FREQ",
  effect_allele_col = "ALLELE1",
  other_allele_col = "ALLELE0",
  pval_col = "exp_pval",
  ncase_col = "N")

exp_dat_DICE_PDCD1 <- dplyr::inner_join(DICE_PDCD1, exp_dat, by = c("rsid" = "SNP")) %>%
dplyr::select(all_of(colnames(DICE_PDCD1))) %>%
dplyr::mutate(inst_set = survival)

exp_dat_DICE_PDCD1_all <- rbind(exp_dat_DICE_PDCD1_all, exp_dat_DICE_PDCD1)

}

resultsdir <- Sys.getenv("resultsdir")
setwd(resultsdir)

write.table(exp_dat_DICE_PDCD1_all, "PD1_dice_assoc_240925.tsv", row.names=F, quote=F, sep = "\t")


