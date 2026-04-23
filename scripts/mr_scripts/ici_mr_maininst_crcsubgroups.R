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
survival <- c("colorectal")
window <- c(500)
pval <- c(5e-6)
r <- c(0.3)
measure <- c("survival")
subgroups <- c("stage_1", "stage_23", "stage_4", "proximal_colon", "distal_colon", "rectal")

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

for (subgroup in subgroups) {
## 2: Load outcome data for genetic instruments using TwoSampleMR package

if (subgroup == "stage_1" | subgroup == "stage_23") {
survival_df <- data.table::fread(paste0("../../cancer_sumstats/crc_PD_L1_PD1SumStats/GWAS_CRC_", subgroup, ".csv"), sep=",", header=T) %>%
dplyr::filter(chromosome == chrom) %>%
dplyr::inner_join(rsids, by = c("position"="POS19"), relationship = "many-to-many") %>%
dplyr::filter((REF==ref_allele&ALT==alt_allele | REF==alt_allele&ALT==ref_allele)) %>%
dplyr::select(!c(ID, REF, ALT, POS38))

} else if (subgroup == "stage_4") {
survival_df <- data.table::fread(paste0("../../cancer_sumstats/crc_PD_L1_PD1SumStats/GWAS_CRC_", subgroup), sep=",", header=T) %>%
dplyr::filter(chromosome == chrom) %>%
dplyr::inner_join(rsids, by = c("position"="POS19"), relationship = "many-to-many") %>%
dplyr::filter((REF==ref_allele&ALT==alt_allele | REF==alt_allele&ALT==ref_allele)) %>%
dplyr::select(!c(ID, REF, ALT, POS38))

} else if (subgroup == "distal_colon" | subgroup == "proximal_colon" | subgroup == "rectal") {
survival_df <- data.table::fread(paste0("../../cancer_sumstats/crc_PD_L1_PD1SumStats/GWAS_CRC_survival_", subgroup, ".csv"), sep=",", header=T) %>%
dplyr::filter(chromosome == chrom) %>%
dplyr::inner_join(rsids, by = c("position"="POS19"), relationship = "many-to-many") %>%
dplyr::filter((REF==ref_allele&ALT==alt_allele | REF==alt_allele&ALT==ref_allele)) %>%
dplyr::select(!c(ID, REF, ALT, POS38))
}

survival_data <- TwoSampleMR::format_data(dat = as.data.frame(survival_df), type = "outcome", beta_col = "crc_gbeta", se_col = "crc_gse", pval_col = "crc_gp", effect_allele_col = "alt_allele", other_allele_col = "ref_allele", eaf_col = "CAF", snp_col = "rsid")

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
twoSMR_res <- dplyr::select(twoSMR_res, c("Protein", "Survival","subgroup", "GWASoutcome","p_thresh",  "window_thresh", "r2_thresh", "method", "nsnp", "b", "se", "pval"))

## 6. Perform two-sample MR using TwoSampleMR package, accounting for LD between SNPs using UK Biobank reference panel
### Write table of harmonised SNPs’ rsids to use as PLINK input
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
twoSMR_UKBld_tbl <- data.frame("method"="IVW accounting for LD using UKB","nsnp"=twoSMR_ld@SNPs,"b"=(twoSMR_ld@Estimate*-1), "CI lower"=(twoSMR_ld@CIUpper*-1), "CI upper"=(twoSMR_ld@CILower*-1),"se"=twoSMR_ld@StdError, "pval"=twoSMR_ld@Pvalue, "Protein"=protein, "Survival"=survival,"Residual SE"=twoSMR_ld@RSE, "Cochran's Q"=twoSMR_ld@Heter.Stat[1], "Heterogeneity p-val"=twoSMR_ld@Heter.Stat[2], "p_thresh"=p_dis, "window_thresh"=window, "r2_thresh"=r, "GWASoutcome"=measure, "subgroup" = subgroup)

all_MR <- full_join(twoSMR_res, twoSMR_UKBld_tbl, by=colnames(twoSMR_res)) 

### Join MR results to table of all MR results
all_MR_res <- rbind(all_MR_res, all_MR) 

}}

setwd(resultsdir)

write.table(all_MR_res, "PD1_PDL1_maininst_surv_res_crcsubgroups_260925.csv", sep = ",", quote = F, row.names = F)

all_MR_res <- read.table("PD1_PDL1_maininst_surv_res_crcsubgroups_260925.csv", sep = ",", header = T) %>%
dplyr::filter(method ==	"IVW accounting for LD using UKB")

crc_main_MR_res <- read.table("PD1_PDL1_maininst_surv_res_240925.csv", sep = ",", header = T) %>%
dplyr::mutate(subgroup = "all")  %>%
dplyr::filter(Survival == "colorectal") %>%
dplyr::filter(method == "IVW accounting for LD using UKB")

all_MR_res_crcmain <- dplyr::full_join(all_MR_res, crc_main_MR_res)

## forest plot

constant_for95CI <- stats::qnorm(1 - (1 - 0.95) / 2)

all_MR_res_crcmain <- dplyr::mutate(all_MR_res_crcmain, sub_form = 
                                      case_when(subgroup == "proximal_colon" ~ "proximal colon",
                                                subgroup == "distal_colon" ~ "distal colon",
                                                subgroup == "stage_4" ~ "stage 4",
                                                subgroup == "stage_23" ~ "stages 2 and 3",
                                                subgroup == "stage_1" ~ "stage 1",
                                                T ~ subgroup))

all_MR_res_crcmain$sub_form <- factor(all_MR_res_crcmain$sub_form,
                                      levels = c("rectal", "proximal colon", 
                                                 "distal colon", "stage 4", "stages 2 and 3", 
                                                 "stage 1", "all"))

pdf("PD1_PDL1_crcsubtypes_surv_forplot_260925.pdf", height = 12, width = 12)
surv_plot <- ggplot(data = all_MR_res_crcmain) +
  geom_point(aes(x = exp(b), y = sub_form, colour = sub_form), size = 1.75, show.legend = F) +
  geom_vline(xintercept = exp(0))+
  geom_errorbarh(aes(y = sub_form, xmin = exp(b-(constant_for95CI*se)),
                     xmax = exp(b+(constant_for95CI*se)),
                     height = 0.2, colour = sub_form), linewidth = 1, show.legend = F) +
  facet_grid(rows = vars(Protein), space = "free", scales = "free") +
  theme_bw(base_size = 30) +
labs(x = "Hazard ratio (95% CI) for mortality per \ngenetically proxied protein NPX unit \ndecrease", 
       y = "Colorectal cancer subgroup")
print(surv_plot)
dev.off()
