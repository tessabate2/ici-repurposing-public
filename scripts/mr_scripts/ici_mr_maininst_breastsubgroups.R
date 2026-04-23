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
survival <- c("breast")
window <- c(500)
pval <- c(5e-6)
r <- c(0.3)
measure <- c("survival")

# breast cancer survival
bcac_surv_path <- c("../../cancer_sumstats/survival_2021_gwas_results_files_BCAC")
BCAC_subgroups <- list.files(bcac_surv_path)
BCAC_subgroups_df <- data.frame(file = BCAC_subgroups) %>%
  dplyr::mutate(file_description = case_when(grepl("_subgroup_a_", file) ~ "younger than age 40 years at diagnosis",
                                             grepl("_subgroup_b_", file) ~ "diagnosed with grade 3 tumours",
                                             grepl("_subgroup_c_", file) ~ "diagnosed with ER+ tumours, who received endocrine therapy any kind",
                                             grepl("_subgroup_d_", file) ~ "diagnosed with ER- tumours, who received chemotherapy any kind",
                                             grepl("_subgroup_e_", file) ~ "diagnosed with tumours that were ER+ or PR+, and HER2-",
                                             grepl("_subgroup_f_", file) ~ "diagnosed with ER+ or PR+, and HER2- tumours, who received chemotherapy any kind",
                                             grepl("_subgroup_g_", file) ~ "diagnosed with ER+ or PR+, and HER2- tumours, who did not receive chemotherapy",
                                             grepl("_subgroup_h_", file) ~ "diagnosed with ER+ or PR+, and HER2+ tumours",
                                             grepl("_subgroup_i_", file) ~ "diagnosed with ER- and PR- and HER2+ tumours",
                                             grepl("_subgroup_j_", file) ~ "diagnosed with ER- and PR- and HER2- tumours",
                                             grepl("_subgroup_k_", file) ~ "who received Tamoxifen",
                                             grepl("_subgroup_l_", file) ~ "who received an aromatase inhibitor",
                                             grepl("_subgroup_m_", file) ~ "who received a CMF-like chemotherapy regimen",
                                             grepl("_subgroup_n_", file) ~ "who received Taxanes",
                                             grepl("_subgroup_o_", file) ~ "who received anthracyclines",
                                             T ~ NA)) %>%
  dplyr::mutate(follow_up = case_when(grepl("_15_years.txt", file) ~ "15",
                                      grepl("_15 years.txt", file) ~ "15",
                                      grepl("_5_years.txt", file) ~ "5",
                                      T ~ NA)) %>%
  drop_na() %>%
  dplyr::filter(follow_up == 15)

subgroup_files <- BCAC_subgroups_df$file

# Create data frame for all results to be stored in
all_MR_res <- data.frame()

rsids <- data.table::fread("UKBB_pQTL_rsids/olink_rsids/all_rsids.csv", sep = ",", header = T)

# Loop over proteins whose expression will be genetically proxied by the instruments
for (protein in proteins) {
## 1. Instrument construction:
if(protein=="PD1") {
path <- c("PDCD1_Q15116_OID21396_v1_Oncology_COMBINED")
chrom <- 2
} else if (protein=="PDL1") {
path <- c("CD274_Q9NZQ7_OID20966_v1_Neurology_COMBINED")
chrom <- 9
}

### Remove - from p-value threshold input so can be used in file/directory names
p_dis <- gsub('-', '', pval)

### Generate directory to store output files for combination of instrument selection thresholds in
dir <- paste0(protein, "/",  path, "/", window, "kb_p", p_dis, "_r", r, "_240925")

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

for (file in subgroup_files) {
survival_df <- data.table::fread(paste0("../../cancer_sumstats/survival_2021_gwas_results_files_BCAC/", file), sep = " ", header = T) %>%
dplyr::filter(Chromosome==chrom) %>%
dplyr::mutate(prioritised_EAF = coalesce(exp_freq_a1_OncoArray, exp_freq_a1_iCOGS)) %>%
inner_join(rsids, by = c("Position"="POS19"), relationship = "many-to-many") %>%
dplyr::mutate(allele_match=ifelse(REF==a1&ALT==a0 | REF==a0&ALT==a1, T, F)) %>%
dplyr::filter(allele_match==T) %>%
dplyr::select(!c(ID, REF, ALT, POS38, allele_match))

survival_data <- TwoSampleMR::format_data(dat = as.data.frame(survival_df), 
type = "outcome", beta_col = "Beta", 
se_col = "SE", pval_col = "P", effect_allele_col = "a1", 
other_allele_col = "a0", eaf_col = "prioritised_EAF", snp_col = "rsid")

subgroup <- BCAC_subgroups_df[BCAC_subgroups_df$file == file, ]$file_description
subgroup <- gsub(",", "", subgroup)
subgroup <- gsub("[+]", "pos", subgroup)
subgroup <- gsub("[-]", "neg", subgroup)
subgroup <- gsub(" ", "_", subgroup)

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
twoSMR_res$subgroup <- subgroup
twoSMR_res$file <- file
twoSMR_res <- dplyr::select(twoSMR_res, c("Protein", "Survival","subgroup", "file", "GWASoutcome","p_thresh",  "window_thresh", "r2_thresh", "method", "nsnp", "b", "se", "pval"))


## 6. Perform two-sample MR using TwoSampleMR package, accounting for LD between SNPs using UK Biobank reference panel
### Write table of harmonised SNPs ^`^y rsids to use as PLINK input
write.table(package_harm$SNP, paste0(dir, "/", protein,  window, "kb_p", p_dis, "_r", r, survival,"_",subgroup, measure, "harm_minorallele_100525.txt"), sep="\t", row.names=F, quote=F)

### Find SNPs which are present in the UK biobank reference panel
system(paste0("plink --bfile ../../ukb_bed/UKBB_10K --extract ", dir, "/", protein, window, "kb_p", p_dis, "_r", r, survival,"_", subgroup, measure, "harm_minorallele_100525.txt --write-snplist --out ", dir, "/", protein, window, "kb_p", p_dis, "_r", r, survival, "_", subgroup, measure, "harm_minorallele_100525"))

### Generate LD matrix for harmonised SNPs using PLINK (signed r LD measure)
system(paste0("plink --bfile ../../ukb_bed/UKBB_10K --extract ", dir, "/", protein, window, "kb_p", p_dis, "_r", r, survival,"_", subgroup, measure, "harm_minorallele_100525.txt --out ", dir, "/", protein, window, "kb_p", p_dis, "_r", r, survival,"_", subgroup, measure, "harm_minorallele_100525 --r square"))

### Read list of SNPs which are present in the UK biobank reference panel
SNP_list_UKB <- read.table(paste0(dir, "/", protein, window, "kb_p", p_dis, "_r", r, survival, "_", subgroup,measure, "harm_minorallele_100525.snplist"), sep = "\t")
colnames(SNP_list_UKB) <- c("SNP")

### Read LD matrix and add SNPs as column/row names
UKBB_LD <- read.table(paste0(dir, "/", protein, window, "kb_p", p_dis, "_r", r, survival,"_", subgroup, measure, "harm_minorallele_100525.ld"))
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

### Format results of two-sample MR
twoSMR_UKBld_tbl <- data.frame("method"="IVW accounting for LD using UKB","nsnp"=twoSMR_ld@SNPs,"b"=(twoSMR_ld@Estimate*-1), "CI lower"=(twoSMR_ld@CIUpper*-1), "CI upper"=(twoSMR_ld@CILower*-1),"se"=twoSMR_ld@StdError, "pval"=twoSMR_ld@Pvalue, "Protein"=protein, "Survival"=survival,"Residual SE"=twoSMR_ld@RSE, "Cochran's Q"=twoSMR_ld@Heter.Stat[1], "Heterogeneity p-val"=twoSMR_ld@Heter.Stat[2], "p_thresh"=p_dis, "window_thresh"=window, "r2_thresh"=r, "GWASoutcome"=measure, "subgroup" = subgroup, "file" = file)

all_MR <- full_join(twoSMR_res, twoSMR_UKBld_tbl, by=colnames(twoSMR_res))

### Join MR results to table of all MR results
all_MR_res <- rbind(all_MR_res, all_MR)

}}

setwd(resultsdir)

write.table(all_MR_res, "PD1_PDL1_maininst_surv_res_breastsubgroups_260925.csv", sep = ",", quote = F, row.names = F)

all_MR_res_breastmain <- read.table("PD1_PDL1_maininst_surv_res_breastsubgroups_260925.csv",sep = ",", header = T) %>%
  dplyr::mutate(subgroup_desc = case_when(grepl("_subgroup_a_", file) ~ "younger than age 40 years at diagnosis",
                                             grepl("_subgroup_b_", file) ~ "diagnosed with grade 3 tumours",
                                             grepl("_subgroup_c_", file) ~ "diagnosed with ER+ tumours, received endocrine therapy any kind",
                                             grepl("_subgroup_d_", file) ~ "diagnosed with ER- tumours, received chemotherapy any kind",
                                             grepl("_subgroup_e_", file) ~ "diagnosed with tumours that were ER+ or PR+, and HER2-",
                                             grepl("_subgroup_f_", file) ~ "diagnosed with ER+ or PR+, and HER2- tumours, received chemotherapy any kind",
                                             grepl("_subgroup_g_", file) ~ "diagnosed with ER+ or PR+, and HER2- tumours, did not receive chemotherapy",
                                             grepl("_subgroup_h_", file) ~ "diagnosed with ER+ or PR+, and HER2+ tumours",
                                             grepl("_subgroup_i_", file) ~ "diagnosed with ER- and PR- and HER2+ tumours",
                                             grepl("_subgroup_j_", file) ~ "diagnosed with ER- and PR- and HER2- tumours",
                                             grepl("_subgroup_k_", file) ~ "received Tamoxifen",
                                             grepl("_subgroup_l_", file) ~ "received an aromatase inhibitor",
                                             grepl("_subgroup_m_", file) ~ "received a CMF-like chemotherapy regimen",
                                             grepl("_subgroup_n_", file) ~ "received Taxanes",
                                             grepl("_subgroup_o_", file) ~ "received anthracyclines",
                                             T ~ NA))

constant_for95CI <- stats::qnorm(1 - (1 - 0.95) / 2)

pdf("PD1_PDL1_breastsubtypes_surv_forplot_260925.pdf", height = 12, width = 12)
surv_plot <- ggplot(data = dplyr::filter(all_MR_res_breastmain, method == "IVW accounting for LD using UKB")) +
  geom_point(aes(x = exp(b), y = subgroup_desc, colour = subgroup_desc), size = 1.75, show.legend = F) +
  geom_vline(xintercept = exp(0))+
  geom_errorbarh(aes(y = subgroup_desc, xmin = exp(b-(constant_for95CI*se)),
                     xmax = exp(b+(constant_for95CI*se)),
                     height = 0.2, colour = subgroup_desc), linewidth = 1, show.legend = F) +
  facet_grid(rows = vars(Protein), space = "free", scales = "free") +
  theme_bw(base_size = 10) +
  labs(x = "Hazard ratio (95% CI) for mortality per genetically proxied protein NPX unit decrease", 
       y = "Breast cancer subgroup")
print(surv_plot)
dev.off()
