# Load required packages
library(tidyverse)
library(data.table)

library(here)
readRenviron(here("config.env"))

datadir <- Sys.getenv("datadir")

setwd(datadir)

setwd("ici_analyses/UKBB")

# Set up variables which loop will be performed over (including instrument selection window (kb), p-value and LD r2 criteria)
proteins <- c("PD1", "PDL1")
survivals <- c("breast", "colorectal", "lung", "melanoma", "ovarian", "prostate")
windows <- c(500, 250, 100)
pvalues <- c(5e-6, 5e-7, 5e-8)
r_lds <- c(0.3, 0.2, 0.1, 0.001)
measures <- c("survival")

N_SNPs_tbl_all <- data.frame()

# Load breast cancer summary statistics
breast_sumstats <- data.table::fread("../../cancer_sumstats/breastcancer_sumstats/gwas_summary_estimates_all_patients_15_years.txt", header = T)

# Loop over proteins whose expression will be genetically proxied by the instruments
for (protein in proteins) {
## 1. Instrument construction:
### Set up values for variables soft coded across the loop 
### path = name of directory containing UK Biobank protein expression results
### rsids = file of rsids for a specific chromosome (obtained from the UK Biobank protein expression GWAS)
### sumstats_path = path to folder of summary statistics for the specific protein expression outcome and the SNPs on the same chromosome as the gene encoding the protein
### chrom = chromosome which the gene encoding the protein is found on
### lower_pos and upper_pos = genomic co-ordinates of the gene encoding the protein in hg19/GRCh37

if (protein=="PD1") {
path <- c("PDCD1_Q15116_OID21396_v1_Oncology_COMBINED")
rsids <- data.table::fread("UKBB_pQTL_rsids/olink_rsids/olink_rsid_map_mac5_info03_b0_7_chr2_patched_v2.tsv", header = T)
sumstats_path <- c("PDCD1_sumstats/combined_chr2_PDCD1:Q15116:OID21396:v1:Oncology")
chrom <- 2
lower_pos <- 242792033
upper_pos <- 242801060

} else if (protein=="PDL1") {
path <- c("CD274_Q9NZQ7_OID20966_v1_Neurology_COMBINED")
rsids <- read.table("UKBB_pQTL_rsids/olink_rsids/olink_rsid_map_mac5_info03_b0_7_chr9_patched_v2.tsv", header = T)
sumstats_path <- c("CD274_sumstats/combined_chr9_CD274:Q9NZQ7:OID20966:v1:Neurology")
chrom <- 9
lower_pos <- 5450503
upper_pos <- 5470566
}

### 1.1: Load outcome GWAS summary statistics depending on which cancer site and outcome measure required (all formatted to have same columns)
# Loop over cancer sites 
for (survival in survivals) {

# Loop over outcome measure to be analysed (risk or survival)
for (measure in measures) {
if (measure=="survival") {

if (survival == "colorectal") {

outcome_sumstats_path <- c("../../cancer_sumstats/colorectalcancer_sumstats") 

outcome_sumstats <- read.table("../../cancer_sumstats/crc_PD_L1_PD1SumStats/GWAS_CRC_survival.csv", sep=",", header=T)
colnames(outcome_sumstats)[colnames(outcome_sumstats) == "crc_gbeta"] <- "Beta"
colnames(outcome_sumstats)[colnames(outcome_sumstats) == "crc_gse"] <- "SE"
colnames(outcome_sumstats)[colnames(outcome_sumstats) == "crc_gp"] <- "P"
colnames(outcome_sumstats)[colnames(outcome_sumstats) == "ref_allele"] <- "a0"
colnames(outcome_sumstats)[colnames(outcome_sumstats) == "alt_allele"] <- "a1"
colnames(outcome_sumstats)[colnames(outcome_sumstats) == "CAF"] <- "a1_freq"

### add rsids
outcome_sumstats <- dplyr::filter(outcome_sumstats, chromosome == chrom) %>%
inner_join(rsids, by = c("position"="POS19"), relationship = "many-to-many") %>%
dplyr::filter((REF==a1&ALT==a0 | REF==a0&ALT==a1)) %>%
dplyr::select(!c(ID, REF, ALT, POS38))

} else if (survival == "prostate"){

outcome_sumstats_path <- c("../../cancer_sumstats/prostatecancer_sumstats") 

outcome_sumstats <- data.table::fread("survival_meta_iCOGS_Onco_BPC3_CAPS_INFO.4_2019-08-22.tsv",  header=T)

colnames(outcome_sumstats)[colnames(outcome_sumstats) == "B.meta.fixed"] <- "Beta"
colnames(outcome_sumstats)[colnames(outcome_sumstats) == "SE.meta.fixed"] <- "SE"
colnames(outcome_sumstats)[colnames(outcome_sumstats) == "P.meta.fixed"] <- "P"
colnames(outcome_sumstats)[colnames(outcome_sumstats) == "AlleleA"] <- "a0"
colnames(outcome_sumstats)[colnames(outcome_sumstats) == "AlleleB"] <- "a1"
colnames(outcome_sumstats)[colnames(outcome_sumstats) == "F_B"] <- "a1_freq"

### add rsids
outcome_sumstats <- dplyr::filter(outcome_sumstats, CHR == chrom) %>%
inner_join(rsids, by = c("BP"="POS19"), relationship = "many-to-many") %>%
dplyr::filter((REF==a1&ALT==a0 | REF==a0&ALT==a1)) %>%
dplyr::select(!c(ID, REF, ALT, POS38))

### Load breast cancer survival GWAS summary statistics
} else if (survival=="breast") {
outcome_sumstats_path <- c("../../cancer_sumstats/breastcancer_sumstats") 

outcome_sumstats <- dplyr::filter(breast_sumstats, Chromosome==chrom) %>%
dplyr::mutate(prioritised_EAF = coalesce(exp_freq_a1_OncoArray, exp_freq_a1_iCOGS)) %>%
inner_join(rsids, by = c("Position"="POS19"), relationship = "many-to-many") %>%
dplyr::mutate(allele_match=ifelse(REF==a1&ALT==a0 | REF==a0&ALT==a1, T, F)) %>%
dplyr::filter(allele_match==T) %>%
dplyr::select(!c(ID, REF, ALT, POS38, allele_match))
colnames(outcome_sumstats)[colnames(outcome_sumstats) == "prioritised_EAF"] <- "a1_freq"

### Load ovarian cancer survival GWAS summary statistics
} else if (survival== "ovarian") {
outcome_sumstats_path <- c("../../cancer_sumstats/OCAC_sumstats")

outcome_sumstats <- read.table(paste0("../../cancer_sumstats/OCAC_sumstats/ocac_imputed_results_",protein, ".csv"), sep=",", header=T)
colnames(outcome_sumstats)[colnames(outcome_sumstats) == "X1KG_Phase1_ID"] <- "rsid"
colnames(outcome_sumstats)[colnames(outcome_sumstats) == "ocac_OverallSurvival_OR"] <- "Beta"
colnames(outcome_sumstats)[colnames(outcome_sumstats) == "ocac_OverallSurvival_se"] <- "SE"
colnames(outcome_sumstats)[colnames(outcome_sumstats) == "ocac_eaf"] <- "a1_freq"
colnames(outcome_sumstats)[colnames(outcome_sumstats) == "ocac_OverallSurvival_pvalue"] <- "P"

} else if (survival == "lung") {
outcome_sumstats_path <- c("../../cancer_sumstats/lungcancer_sumstats")

outcome_sumstats <- read.table("../../cancer_sumstats/lungcancer_sumstats/meta_analysis/METAANALYSIS_survival_allstages_ilcco_dfci_ge_mgi1_formatted.txt", header=T) %>%
dplyr::filter(chr==chrom) %>%
inner_join(rsids, by = c("bp"="POS38"), relationship = "many-to-many") %>%
dplyr::mutate(Allele1=toupper(Allele1)) %>%
dplyr::mutate(Allele2=toupper(Allele2)) %>%
dplyr::mutate(allele_match=ifelse(REF==Allele1&ALT==Allele2 | REF==Allele2&ALT==Allele1, T, F)) %>%
dplyr::filter(allele_match==T) %>%
dplyr::select(!c(ID, REF, ALT, POS19, allele_match))

colnames(outcome_sumstats)[colnames(outcome_sumstats) == "Effect"] <- "Beta"
colnames(outcome_sumstats)[colnames(outcome_sumstats) == "StdErr"] <- "SE"
colnames(outcome_sumstats)[colnames(outcome_sumstats) == "Freq1"] <- "a1_freq"
colnames(outcome_sumstats)[colnames(outcome_sumstats) == "P-value"] <- "P"
colnames(outcome_sumstats)[colnames(outcome_sumstats) == "Allele1"] <- "a1"
colnames(outcome_sumstats)[colnames(outcome_sumstats) == "Allele2"] <- "a0"

} else if (survival == "melanoma") {
outcome_sumstats_path <- c("../../cancer_sumstats/melanoma_sumstats")

outcome_sumstats <- read.table("../../cancer_sumstats/melanoma_sumstats/PD1_PDL1_extraction_06072023", header=T) %>%
dplyr::filter(CHR==chrom)

colnames(outcome_sumstats)[colnames(outcome_sumstats) == "BETA"] <- "Beta"
colnames(outcome_sumstats)[colnames(outcome_sumstats) == "FRQ1"] <- "a1_freq"
colnames(outcome_sumstats)[colnames(outcome_sumstats) == "A1"] <- "a1"
colnames(outcome_sumstats)[colnames(outcome_sumstats) == "A2"] <- "a0"
colnames(outcome_sumstats)[colnames(outcome_sumstats) == "rsID"] <- "rsid"

}
}

if("a1_freq" %in% colnames(outcome_sumstats)){
### Filter 0.01 < EAF < 0.99 for outcome summary statistics and remove SNPs with multiple alleles
outcome_sumstats <- outcome_sumstats[outcome_sumstats$a1_freq > 0.01, ]
outcome_sumstats <- outcome_sumstats[outcome_sumstats$a1_freq < 0.99, ]
}

alleles <- c("A", "C", "G", "T")
outcome_sumstats <- dplyr::filter(outcome_sumstats, a0 %in% alleles) %>%
dplyr::filter(a1 %in% alleles)

### Write table of summary statistics for the outcome so can be read by TwoSampleMR package
write.table(outcome_sumstats, paste0(outcome_sumstats_path, "/",  survival, measure, "_sumstats_", protein, "_100525.csv"), sep=",", row.names=F, quote=F)

### 1.2: Load UK biobank protein expression summary statistics
sumstats <- read.table(paste0(protein, "/", path, "/", sumstats_path), header=T) %>%
### Add corresponding rsids to summary statistics
inner_join(rsids, by = "ID") %>%
dplyr::mutate(exp_pval = (10^-(LOG10P)))

input_sumstats <- dim(sumstats)[1]

### Restrict to SNPs present in the outcome dataset
sumstats <- sumstats[sumstats$rsid %in% outcome_sumstats$rsid, ]

input_surv_sumstats <- dim(sumstats)[1]

### Remove SNPs with multiple alleles
sumstats <- dplyr::filter(sumstats, ALLELE0 %in% alleles) %>%
dplyr::filter(ALLELE1 %in% alleles)

input_surv_sumstats_singleall <- dim(sumstats)[1]

### Find the min and max EAFs
#min(sumstats$A1FREQ)
#max(sumstats$A1FREQ)

### (if necessary) Restrict to SNPs with 0.01 < EAF < 0.99
sumstats <- sumstats[!((sumstats$A1FREQ < 0.01)|(sumstats$A1FREQ > 0.99)), ]
#min(sumstats$A1FREQ)
#max(sumstats$A1FREQ)

input_surv_sumstats_singleall_maf <- dim(sumstats)[1]

### Calculate inverse z-scores for each SNP:
#### Calculate z-score = beta/SE
sumstats$z <- (sumstats$BETA)/(sumstats$SE)
#### Find absolute z-score by squaring, then square rooting z-scores
sumstats$z2 <- (sumstats$z)^2
sumstats$sqrtz2 <- sqrt(sumstats$z2)
#### Calculate inverse z-score = 1/z-score
sumstats$inv_sqrtz2 <- 1/(sumstats$sqrtz2)

### 1.3: Filter by instrument p-value threshold using corresponding z-scores found using https://www.gigacalculator.com/calculators/p-value-to-z-score-calculator.php (two-tailed)
# Loop over instrument selection p-value thresholds
for (pval in pvalues) {
if (pval == "5e-06"){
p_val_zthresh <- 1/4.564798
} else if (pval == "5e-07") {
p_val_zthresh <- 1/5.026337
} else if (pval == "5e-08") {
p_val_zthresh <- 1/5.451355
}

### Restrict to SNPs with inverse z-scores below the inverse z-score threshold
sumstats_pz <- sumstats[sumstats$inv_sqrtz2 < p_val_zthresh, ]

sumstats_qc_p <- dim(sumstats_pz)[1]

### 1.4: Filter by window size (kb)
# Loop over instrument selection window criteria
for (window in windows) {
### Restrict to SNPs with genomic position (hg19) within 500kb below or above the gene encoding the protein of interest
sumstats_pz_window <- sumstats_pz[sumstats_pz$POS19 > (lower_pos - (window*1000)), ]
sumstats_pz_window <- sumstats_pz_window[sumstats_pz_window$POS19 < (upper_pos+ (window*1000)), ]

sumstats_qc_p_window <- dim(sumstats_pz_window)[1]

### Sort inverse z-scores for retained SNPs
sumstats_pz_window_order <- sumstats_pz_window[order(sumstats_pz_window$inv_sqrtz2), ] %>%
### Keep rsid and inverse z-score columns for PLINK input
dplyr::select(c("rsid", "inv_sqrtz2"))

### Remove - from p-value threshold input so can be used in file/directory names
p_dis <- gsub('-', '', pval)

# Loop over LD r2 thresholds
for (r in r_lds) {

### Generate directory to store output files for combination of instrument selection thresholds in
dir <- paste0(protein, "/", path, "/", window, "kb_p", p_dis, "_r", r, "_240925")
system(paste0("mkdir ", dir))

### Save restricted SNPs so can be clumped using PLINK
write.table(sumstats_pz_window_order, paste0(dir, "/", protein, survival,  measure, window,"kb_p", p_dis, "inst_100525.tsv"), row.names = FALSE, sep = "\t", quote=F)

### 1.5: Filter by LD r2 threshold using PLINK clumping based on inverse z-scores
#system(paste0("plink --bfile ../../ukb_bed/UKBB_10K --clump ", dir, "/", protein, survival, measure, window,"kb_p", p_dis, "inst_100525.tsv --clump-p1 ", p_val_zthresh,  " --clump-field inv_sqrtz2 --clump-snp-field rsid --clump-r2 ",r," --clump-kb 1050  --out ", dir, "/", protein, survival, measure, window,"kb_p", p_dis,"r",r, "_100525"))

### Read table of clumped SNPs generated using PLINK
inst_window_p_clump <- read.table(paste0(dir,"/", protein,survival,measure, window,"kb_p",p_dis,"r",r, "_100525.clumped"), header=T) 
### Retain summary statistics for clumped SNPs only
inst_window_p_clump_sumstats <- sumstats[sumstats$rsid %in% inst_window_p_clump$SNP, ]

sumstats_qc_p_window_r2 <- dim(inst_window_p_clump_sumstats)[1]

N_SNPs_tbl <- data.frame(protein = protein, survival = survival, measure = measure, window_thresh = window, pval_thresh = pval, r2_thresh = r,
input_N = input_sumstats, survival_dat_N = input_surv_sumstats, single_allele_N = input_surv_sumstats_singleall, maf_N = input_surv_sumstats_singleall_maf, pval = sumstats_qc_p, window = sumstats_qc_p_window, r2 = sumstats_qc_p_window_r2)

N_SNPs_tbl_all <- rbind(N_SNPs_tbl, N_SNPs_tbl_all)

### Save table of genetic instruments so can be read in as exposure data using TwoSampleMR package
write.csv(inst_window_p_clump_sumstats, paste0(dir,"/", protein, survival, measure, window,"kb_p",p_dis,"r",r,"inst_100525.csv"), row.names=F)

}}}}}}

resultsdir <- Sys.getenv("resultsdir")
setwd(resultsdir)
write.csv(N_SNPs_tbl_all, "N_SNPs_PD1_PDL1_combined_240925.csv", row.names=F, quote=F)

