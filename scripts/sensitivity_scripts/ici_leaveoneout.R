# Load required packages
library(tidyverse)
library(data.table)

library(here)
readRenviron(here("config.env"))

datadir <- Sys.getenv("datadir")
setwd(datadir)
setwd("ici_analyses/UKBB")

resultsdir <- Sys.getenv("resultsdir")

# Load required packages
library(TwoSampleMR)
library(tidyverse)
library(MendelianRandomization)

# Set up variables which loop will be performed over (including instrument selection window (kb), p-value and LD r2 criteria)
proteins <- c("PD1", "PDL1")
survivals <- c("breast", "colorectal", "lung", "melanoma", "ovarian", "prostate")
window <- c(500)
pval <- c(5e-6)
r <- c(0.3)
measure <- c("survival")

# Loop over proteins whose expression will be genetically proxied by the instruments
for (protein in proteins) {
## 1. Instrument construction:
if(protein=="PD1") {
path <- c("PD1/PDCD1_Q15116_OID21396_v1_Oncology_COMBINED")
} else if (protein=="PDL1") {
path <- c("PDL1/CD274_Q9NZQ7_OID20966_v1_Neurology_COMBINED")
}

### 1.1: Load outcome GWAS summary statistics depending on which cancer site and outcome measure required (all formatted to have same columns)
# Loop over cancer sites 
for (survival in survivals) {

if (survival == "colorectal") {
outcome_sumstats_path <- c("colorectalcancer_sumstats")

} else if (survival == "prostate") {
outcome_sumstats_path <- c("prostatecancer_sumstats") 

### Load breast cancer survival GWAS summary statistics
} else if (survival=="breast") {
outcome_sumstats_path <- c("survival_2021_gwas_results_files_BCAC") 

### Load ovarian cancer survival GWAS summary statistics
} else if (survival== "ovarian") {
outcome_sumstats_path <- c("OCAC_sumstats")

} else if (survival == "lung") {
outcome_sumstats_path <- c("lungcancer_sumstats")

} else if (survival == "melanoma") {
outcome_sumstats_path <- c("melanoma_sumstats")
}

### Remove - from p-value threshold input so can be used in file/directory names
p_dis <- gsub('-', '', pval)

### Generate directory to store output files for combination of instrument selection thresholds in
dir <- paste0(path, "/", window, "kb_p", p_dis, "_r", r, "_240925")

### 1.6: Read exposure data for genetic instruments selected using the above p-value, window and LD r2 thresholds
exp_dat <- TwoSampleMR::read_exposure_data(
  filename=paste0(dir,"/",protein, survival, measure, window,"kb_p",p_dis,"r",r,"inst_100525.csv"),
  sep = ",",
  snp_col = "rsid",
  beta_col = "BETA",
  se_col = "SE",
  eaf_col = "A1FREQ",
  effect_allele_col = "ALLELE1",
  other_allele_col = "ALLELE0",
  pval_col = "exp_pval",
  ncase_col = "N")

allele_freqs <- read.table(paste0(dir, "/" ,protein, window, "kb_p",p_dis,"_r",r, survival, measure, "_minor_100525.frq"), header=T)

### List of SNPs where A1 in .frq file (assume minor allele) is not the effect allele
mismatch_minoralleles <- full_join(exp_dat, allele_freqs, by = c("SNP")) %>%
dplyr::filter(!((A1==effect_allele.exposure)&(A2==other_allele.exposure)))

exp_dat_minorallele <- dplyr::mutate(exp_dat, effectallele.exposure_min = case_when(exp_dat$SNP %in% mismatch_minoralleles$SNP ~ other_allele.exposure, T ~ effect_allele.exposure),
otherallele.exposure_min = case_when(exp_dat$SNP %in% mismatch_minoralleles$SNP ~ effect_allele.exposure, T ~ other_allele.exposure),
beta.exposure_min = case_when(exp_dat$SNP %in% mismatch_minoralleles$SNP ~ (beta.exposure*-1), T ~ beta.exposure),
eaf.exposure_min = case_when(exp_dat$SNP %in% mismatch_minoralleles$SNP ~ (1-eaf.exposure), T ~ eaf.exposure)) %>%
dplyr::select(-c(eaf.exposure, effect_allele.exposure, other_allele.exposure, beta.exposure)) %>%
dplyr::rename(eaf.exposure = eaf.exposure_min, effect_allele.exposure = effectallele.exposure_min, other_allele.exposure = otherallele.exposure_min, beta.exposure = beta.exposure_min)

## 2: Load outcome data for genetic instruments using TwoSampleMR package
survival_data <- TwoSampleMR::read_outcome_data(
  snps = exp_dat$SNP,
  filename = paste0("../../cancer_sumstats/", outcome_sumstats_path, "/", survival, measure, "_sumstats_", protein, ".csv"),
  sep = ",",
  snp_col = "rsid",
  beta_col = "Beta",
  se_col = "SE",
  effect_allele_col = "a1",
  other_allele_col = "a0",
  eaf_col = "a1_freq",
  pval_col = "P")

## 3: Harmonise exposure and outcome data using TwoSampleMR package
package_harm_all <- TwoSampleMR::harmonise_data(exp_dat_minorallele, survival_data, action=2)
package_harm_all <- package_harm_all[package_harm_all$mr_keep==T, ] 

rsids <- package_harm_all$SNP
twoSMR_UKBld_tbl_LOOres_all <- data.frame()

for (rsid in rsids) {

package_harm <- dplyr::filter(package_harm_all, !SNP == rsid)

### Write table of harmonised SNPs rsids to use as PLINK input
write.table(package_harm$SNP, paste0(dir, "/", protein,  window, "kb_p", p_dis, "_r", r, survival, measure, "harm_minorallele_excl_", rsid, "_010725.txt"), sep="\t", row.names=F, quote=F)

### Find SNPs which are present in the UK biobank reference panel
system(paste0("plink --bfile ../../ukb_bed/UKBB_10K --extract ", dir, "/", protein, window, "kb_p", p_dis, "_r", r, survival, measure, "harm_minorallele_excl_", rsid, "_010725.txt --write-snplist --out ", dir, "/", protein, window, "kb_p", p_dis, "_r", r, survival, measure, "harm_minorallele_excl_", rsid, "_010725"))

### Generate LD matrix for harmonised SNPs using PLINK (signed r LD measure)
system(paste0("plink --bfile ../../ukb_bed/UKBB_10K --extract ", dir, "/", protein, window, "kb_p", p_dis, "_r", r, survival, measure, "harm_minorallele_excl_", rsid, "_010725.txt --out ", dir, "/", protein, window, "kb_p", p_dis, "_r", r, survival, measure, "harm_minorallele_excl_", rsid, "_010725 --r square"))

## 6. Perform two-sample MR using TwoSampleMR package, accounting for LD between SNPs using UK Biobank reference panel
### Read list of SNPs which are present in the UK biobank reference panel
SNP_list_UKB <- read.table(paste0(dir, "/", protein, window, "kb_p", p_dis, "_r", r, survival, measure, "harm_minorallele_excl_", rsid, "_010725.snplist"), sep = "\t")
colnames(SNP_list_UKB) <- c("SNP")

### Read LD matrix and add SNPs as column/row names
UKBB_LD <- read.table(paste0(dir, "/", protein, window, "kb_p", p_dis, "_r", r, survival, measure, "harm_minorallele_excl_", rsid, "_010725.ld"))
colnames(UKBB_LD) <- SNP_list_UKB$SNP
row.names(UKBB_LD) <- SNP_list_UKB$SNP
UKBB_LD_matrix <- data.matrix(UKBB_LD)

### Restrict harmonised data to SNPs present in the UK biobank reference panel (should be in same order as the LD matrix)
harm_ld_order_UKB <- inner_join(SNP_list_UKB, package_harm, by = c("SNP"))

### Format harmonised data and LD matrix as MendelianRandomization input
MRcor_LOO_cor <- MendelianRandomization::mr_input(
  bx = harm_ld_order_UKB$beta.exposure,
  bxse = harm_ld_order_UKB$se.exposure,
  by = harm_ld_order_UKB$beta.outcome,
  byse = harm_ld_order_UKB$se.outcome,
  correlation = UKBB_LD_matrix,
  exposure = harm_ld_order_UKB$id.exposure,
  outcome = harm_ld_order_UKB$id.outcome,
  snps = harm_ld_order_UKB$SNP,
  effect_allele = harm_ld_order_UKB$effect_allele.exposure,
  other_allele = harm_ld_order_UKB$other_allele.exposure,
  eaf = harm_ld_order_UKB$eaf.exposure)

### Perform random-effects IVW two-sample MR adjusting for LD between SNPs
twoSMR_ld_LOO_cor <- MendelianRandomization::mr_ivw(MRcor_LOO_cor, correl = MRcor_LOO_cor$correlation, model = "random")

twoSMR_UKBld_tbl_LOOcor <- data.frame("SNP excluded"=rsid, "Method" = "IVW accounting for LD","nsnp"= twoSMR_ld_LOO_cor@SNPs,"b"=(twoSMR_ld_LOO_cor@Estimate*-1), "CI lower"=(twoSMR_ld_LOO_cor@CIUpper*-1), "CI upper"=(twoSMR_ld_LOO_cor@CILower*-1),"se"= twoSMR_ld_LOO_cor@StdError, "pval"= twoSMR_ld_LOO_cor@Pvalue, "Protein"=protein, "Survival"=survival,"Residual SE"= twoSMR_ld_LOO_cor@RSE, "Cochran's Q"= twoSMR_ld_LOO_cor@Heter.Stat[1], "Heterogeneity p-val"= twoSMR_ld_LOO_cor@Heter.Stat[2], "p_thresh"=p_dis, "window_thresh"=window, "r2_thresh"=0.3, "GWASoutcome"=measure)

twoSMR_UKBld_tbl_LOOres_all <- rbind(twoSMR_UKBld_tbl_LOOres_all, twoSMR_UKBld_tbl_LOOcor)

}

setwd(resultsdir)
write.table(twoSMR_UKBld_tbl_LOOres_all, paste0("leaveoneout_290925/", protein, "_maininst_", survival, "_survival_loo_estimates_290925.csv"), sep = ",", quote=F, row.names = F)

setwd(datadir)
setwd("ici_analyses/UKBB")

}}


setwd(resultsdir)
leaveoneout_files <- list.files("leaveoneout_290925")

for (file in leaveoneout_files) {

twoSMR_UKBld_tbl_LOOres_all <- read.table(paste0("leaveoneout_290925/", file), sep = ",", header = T) %>%
dplyr::arrange(b)

protein <- strsplit(file, "_")[[1]][1]
site <- strsplit(file, "_")[[1]][3]

twoSMR_UKBld_tbl <- read.table("PD1_PDL1_maininst_surv_res_240925.csv", sep = ",", header = T) %>%
dplyr::filter(Survival == site) %>%
dplyr::filter(Protein == protein) %>%
dplyr::filter(method == "IVW accounting for LD using UKB") %>%
dplyr::mutate(SNP.excluded = "none")

twoSMR_UKBld_tbl_LOOres_all <- dplyr::full_join(twoSMR_UKBld_tbl_LOOres_all, twoSMR_UKBld_tbl)
twoSMR_UKBld_tbl_LOOres_all$SNP.excluded <- factor(twoSMR_UKBld_tbl_LOOres_all$SNP.excluded, levels = twoSMR_UKBld_tbl_LOOres_all$SNP.excluded)

## set constant for arms of forest plot
constant_for95CI <- stats::qnorm(1 - (1 - 0.95) / 2)

pdf(paste0("leaveoneout_290925_forplots/", protein, "_", site, "_looforplot_171225.pdf"), width = 12, height = 18)
surv_plot_loo <- ggplot(data = twoSMR_UKBld_tbl_LOOres_all) +
  geom_point(aes(x = exp(b), y = SNP.excluded), size = 1.75, show.legend = F) +
  geom_vline(xintercept = exp(0))+
  geom_errorbarh(aes(y = SNP.excluded, xmin = exp(b-(constant_for95CI*se)), 
                     xmax = exp(b+(constant_for95CI*se)),
                     height = 0.2), linewidth = 1, show.legend = F) +
  theme_bw(base_size = 15) +
  labs(x = "HR (95% CI)", y = "Instrument excluded from set") +
ggtitle(paste0("Leave one out: ", protein, " instruments on risk of mortality for ", site, " cancer patients"))
print(surv_plot_loo)
dev.off()

}






