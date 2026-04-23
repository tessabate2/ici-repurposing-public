# Load required packages
library(TwoSampleMR)
library(tidyverse)
library(MendelianRandomization)

library(here)
readRenviron(here("config.env"))

datadir <- Sys.getenv("datadir")
resultsdir <- Sys.getenv("resultsdir")

setwd(datadir)
setwd("ici_analyses/UKBB")

# Set up variables which loop will be performed over (including instrument selection window (kb), p-value and LD r2 criteria)
proteins <- c("PD1", "PDL1")
survival <- c("lung")
subtypes <- c("allstages", "allstages_stratified1234", "stage12", "stage34")
windows <- c(500)
pvalues <- c(5e-6)
r_lds <- c(0.3)
measure <- c("survival")

# Create data frame for all results to be stored in
all_MR_res <- data.frame()

# Loop over proteins whose expression will be genetically proxied by the instruments
for (protein in proteins) {
## 1. Instrument construction:
if(protein=="PD1") {
path <- c("PDCD1_Q15116_OID21396_v1_Oncology_COMBINED")
rsids <- read.table("UKBB_pQTL_rsids/olink_rsids/olink_rsid_map_mac5_info03_b0_7_chr2_patched_v2.tsv", header = T)
chrom <- 2
} else if (protein=="PDL1") {
path <- c("CD274_Q9NZQ7_OID20966_v1_Neurology_COMBINED")
rsids <- read.table("UKBB_pQTL_rsids/olink_rsids/olink_rsid_map_mac5_info03_b0_7_chr9_patched_v2.tsv", header = T)
chrom <- 9
}

### 1.1: Load outcome GWAS summary statistics depending on which cancer site and outcome measure required (all formatted to have same columns)

outcome_sumstats_path <- c("lungcancer_sumstats")

# Loop over instrument selection p-value thresholds
for (pval in pvalues) {

# Loop over instrument selection window criteria
for (window in windows) {

### Remove - from p-value threshold input so can be used in file/directory names
p_dis <- gsub('-', '', pval)

# Loop over LD r2 thresholds
for (r in r_lds) {

### Generate directory to store output files for combination of instrument selection thresholds in
dir <- paste0(protein, "/",  path, "/", window, "kb_p", p_dis, "_r", r, "_240925")

### 1.6: Read exposure data for genetic instruments selected using the above p-value, window and LD r2 thresholds

for (subtype in subtypes) {

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

write.table(exp_dat$SNP, paste0(dir, "/", protein,  window, "kb_p", p_dis, "_r", r, survival, measure, "harm_minorallele_100525.txt"), quote=F, row.names=F)

### Load .frq file
system(paste0("plink --bfile ../../ukb_bed/UKBB_10K --extract ", dir, "/", protein,  window, "kb_p", p_dis, "_r", r, survival, measure, "harm_minorallele_100525.txt --freq --out ", dir, "/", protein, window, "kb_p", p_dis, "_r", r, survival, measure, "_minor_100525"))

allele_freqs <- read.table(paste0(dir, "/" ,protein, window, "kb_p",p_dis,"_r",r, survival, measure, "_minor_100525.frq"), header=T)

### List of SNPs where A1 in .frq file (assume minor allele) is not the effect allele
mismatch_minoralleles <- full_join(exp_dat, allele_freqs, by = c("SNP")) %>%
dplyr::filter(!((A1==effect_allele.exposure)&(A2==other_allele.exposure)))

exp_dat_minorallele <- dplyr::mutate(exp_dat, effectallele.exposure_min = case_when(exp_dat$SNP %in% mismatch_minoralleles$SNP ~ other_allele.exposure, T ~ effect_allele.exposure),
otherallele.exposure_min = case_when(exp_dat$SNP %in% mismatch_minoralleles$SNP ~ effect_allele.exposure, T ~ other_allele.exposure),
beta.exposure_min = case_when(exp_dat$SNP %in% mismatch_minoralleles$SNP ~ (beta.exposure*-1), T ~ beta.exposure),
eaf.exposure_min = case_when(exp_dat$SNP %in% mismatch_minoralleles$SNP ~ (1-eaf.exposure), T ~ eaf.exposure)) %>%
dplyr::select(-c(eaf.exposure, effect_allele.exposure, other_allele.exposure, beta.exposure)) %>%
dplyr::rename(eaf.exposure = eaf.exposure_min, effect_allele.exposure = effectallele.exposure_min, other_allele.exposure = otherallele.exposure_min, beta.exposure = beta.exposure_min) %>%
dplyr::select(!z.exposure)

## 2: Load outcome data for genetic instruments using TwoSampleMR package
out_data <- read.table(paste0("../../cancer_sumstats/lungcancer_sumstats/meta_analysis/METAANALYSIS_survival_", subtype, "_ilcco_dfci_ge_mgi1_formatted.txt"), sep = "\t", header = T) %>%
dplyr::filter(chr==chrom) %>%
inner_join(rsids, by = c("bp"="POS38"), relationship = "many-to-many") %>%
dplyr::mutate(Allele1=toupper(Allele1)) %>%
dplyr::mutate(Allele2=toupper(Allele2)) %>%
dplyr::mutate(allele_match=ifelse(REF==Allele1&ALT==Allele2 | REF==Allele2&ALT==Allele1, T, F)) %>%
dplyr::filter(allele_match==T) %>%
dplyr::select(!c(ID, REF, ALT, POS19, allele_match))

survival_data <- TwoSampleMR::format_data(out_data, type = "outcome", beta_col = "Effect", se_col = "StdErr", eaf_col = "Freq1", pval_col = "P.value", effect_allele_col = "Allele1", other_allele_col = "Allele2",
snp_col = "rsid")

## 3: Harmonise exposure and outcome data using TwoSampleMR package
package_harm <- TwoSampleMR::harmonise_data(exp_dat_minorallele, survival_data, action=2)
package_harm <- package_harm[package_harm$mr_keep==T, ] 

## 4. Perform two-sample MR using TwoSampleMR package (not accounting for LD between SNPs)
### If there were SNPs present in the harmonised dataset, perform two-sample MR
if (!dim(package_harm)[1] ==0) {
twoSMR <- TwoSampleMR::mr(package_harm)

### If there was only 1 SNP present in the harmonised dataset, use the Wald ratio result
if (dim(package_harm)[1] == 1) {
twoSMR_res <- twoSMR[twoSMR$method=="Wald ratio", ]

### If there were more than 1 SNPs present in the harmonised dataset, use the inverse-variance weighted result
} else if (dim(package_harm)[1] > 1) {
twoSMR_res <- twoSMR[twoSMR$method=="Inverse variance weighted", ]
}

### If there were no SNPs present in the harmonised dataset, use NA as results
} else if (dim(package_harm)[1] ==0) {
twoSMR_res$method <- NA
twoSMR_res$nsnp <- 0
twoSMR_res$b <- NA
twoSMR_res$se <- NA
twoSMR_res$pval <- NA
twoSMR_res$id.exposure <- NA
twoSMR_res$id.outcome <- NA
twoSMR_res$outcome <- NA
twoSMR_res$exposure <- NA
}

### Store parameter information (protein expression proxied, cancer site, criteria used in instrument selection, and outcome measure)
twoSMR_res$b <- twoSMR_res$b*-1
twoSMR_res$Protein <- protein
twoSMR_res$Survival <- survival
twoSMR_res$p_thresh <- p_dis
twoSMR_res$window_thresh <- window
twoSMR_res$r2_thresh <- r
twoSMR_res$GWASoutcome <- measure
twoSMR_res$subgroup <- subtype
twoSMR_res <- dplyr::select(twoSMR_res, c("Protein", "Survival","subgroup", "GWASoutcome","p_thresh",  "window_thresh", "r2_thresh", "method", "nsnp", "b", "se", "pval"))

## 6. Perform two-sample MR using TwoSampleMR package, accounting for LD between SNPs using UK Biobank reference panel
### Write table of harmonised SNPs ^`^y rsids to use as PLINK input
write.table(package_harm$SNP, paste0(dir, "/", protein,  window, "kb_p", p_dis, "_r", r, survival,"_",subtype, measure, "harm_minorallele_100525.txt"), sep="\t", row.names=F, quote=F)

### Find SNPs which are present in the UK biobank reference panel
system(paste0("plink --bfile ../../ukb_bed/UKBB_10K --extract ", dir, "/", protein, window, "kb_p", p_dis, "_r", r, survival, "_", subtype, measure, "harm_minorallele_100525.txt --write-snplist --out ", dir, "/", protein, window, "kb_p", p_dis, "_r", r, survival, "_", subtype, measure, "harm_minorallele_100525"))

### Generate LD matrix for harmonised SNPs using PLINK (signed r LD measure)
system(paste0("plink --bfile ../../ukb_bed/UKBB_10K --extract ", dir, "/", protein, window, "kb_p", p_dis, "_r", r, survival,"_", subtype, measure, "harm_minorallele_100525.txt --out ", dir, "/", protein, window, "kb_p", p_dis, "_r", r, survival,"_", subtype, measure, "harm_minorallele_100525 --r square"))

### Read list of SNPs which are present in the UK biobank reference panel
SNP_list_UKB <- read.table(paste0(dir, "/", protein, window, "kb_p", p_dis, "_r", r, survival, "_", subtype, measure, "harm_minorallele_100525.snplist"), sep = "\t")
colnames(SNP_list_UKB) <- c("SNP")

### Read LD matrix and add SNPs as column/row names
UKBB_LD <- read.table(paste0(dir, "/", protein, window, "kb_p", p_dis, "_r", r, survival,"_", subtype, measure, "harm_minorallele_100525.ld"))
colnames(UKBB_LD) <- SNP_list_UKB$SNP
row.names(UKBB_LD) <- SNP_list_UKB$SNP
UKBB_LD_matrix <- data.matrix(UKBB_LD)

### Restrict harmonised data to SNPs present in the UK biobank reference panel (should be in same order as the LD matrix)
harm_ld_order_UKB <- inner_join(SNP_list_UKB, package_harm, by = c("SNP"))

### Format harmonised data and LD matrix as MendelianRandomization input
MRcor <- mr_input(
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
twoSMR_ld <- MendelianRandomization::mr_ivw(MRcor, correl = MRcor$correlation, model = "random")

twoSMR_UKBld_tbl <- data.frame("method"="IVW accounting for LD using UKB","nsnp"=twoSMR_ld@SNPs,"b"=(twoSMR_ld@Estimate*-1), "CI lower"=(twoSMR_ld@CIUpper*-1), "CI upper"=(twoSMR_ld@CILower*-1),"se"=twoSMR_ld@StdError, "pval"=twoSMR_ld@Pvalue, "Protein"=protein, "Survival"=survival,"Residual SE"=twoSMR_ld@RSE, "Cochran's Q"=twoSMR_ld@Heter.Stat[1], "Heterogeneity p-val"=twoSMR_ld@Heter.Stat[2], "p_thresh"=p_dis, "window_thresh"=window, "r2_thresh"=r, "GWASoutcome"=measure, "subgroup" = subtype)

all_MR <- full_join(twoSMR_res, twoSMR_UKBld_tbl, by=colnames(twoSMR_res))

### Join MR results to table of all MR results
all_MR_res <- rbind(all_MR_res, all_MR)


}}}}}

setwd(resultsdir)

write.table(all_MR_res, "PD1_PDL1_maininst_surv_res_lungsubgroups_240925.csv", sep = ",", quote = F, row.names = F)

all_MR_res_lungmain <- read.table("PD1_PDL1_maininst_surv_res_lungsubgroups_240925.csv", sep = ",", header = T)

constant_for95CI <- stats::qnorm(1 - (1 - 0.95) / 2)

all_MR_res_lungmain <- dplyr::mutate(all_MR_res_lungmain, sub_form = 
                                      case_when(subgroup == "allstages" ~ "all stages",
                                                subgroup == "allstages_stratified1234" ~ "all stages stratified",
                                                subgroup == "stage12" ~ "stages 1 and 2",
                                                subgroup == "stage34" ~ "stages 3 and 4",
                                                T ~ subgroup))

all_MR_res_lungmain$sub_form <- factor(all_MR_res_lungmain$sub_form,
                                      levels = c("stages 3 and 4", "stages 1 and 2", 
                                                 "all stages stratified",
                                                 "all stages"))

pdf("PD1_PDL1_lungsubgroups_surv_forplot_240925.pdf", height = 12, width = 12)
surv_plot <- ggplot(data = dplyr::filter(all_MR_res_lungmain, method == "IVW accounting for LD using UKB")) +
  geom_point(aes(x = exp(b), y = sub_form, colour = sub_form), size = 1.75, show.legend = F) +
  geom_vline(xintercept = exp(0))+
  geom_errorbarh(aes(y = sub_form, xmin = exp(b-(constant_for95CI*se)),
                     xmax = exp(b+(constant_for95CI*se)),
                     height = 0.2, colour = sub_form), linewidth = 1, show.legend = F) +
  facet_grid(rows = vars(Protein), space = "free", scales = "free") +
  theme_bw(base_size = 30) +
  labs(x = "Hazard ratio (95% CI) for mortality per \ngenetically proxied protein NPX unit \ndecrease", 
       y = "Lung cancer stage")
print(surv_plot)
dev.off()
