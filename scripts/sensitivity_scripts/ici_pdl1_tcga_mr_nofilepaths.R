library(data.table)

#1. Look up whether lead PD-L1 instrument is a CD274 eQTL in pancanQTL any cancer site (FDR p<0.05)
cis_pancanqtl <- data.table::fread("cis_eQTLs_all_re", sep = "\t", header = F)
colnames(cis_pancanqtl) <- c("cancer_type", "rsid", "chr", "position", "alleles", "gene", "gene_position", "beta", "t-stat", "p-value")

cis_pancanqtl_cd274_leadpdl1inst <- dplyr::filter(cis_pancanqtl, gene == "CD274" & rsid == "rs822341")

#2. Look up whether any proxies (R2>0.6) of the lead PD-L1 instrument are a CD274 eQTL in pancanQTL any cancer site (FDR p<0.05)
library(gwasvcf)
gwasvcf::set_plink(path = "")
gwasvcf::check_plink()
pdl1inst_proxy <- gwasvcf::get_ld_proxies(
  rsid = "rs822341",
  bfile = "/data/ukb_bed/UKBB_10K",
  searchspace = NULL,
  tag_kb = 5000,
  tag_nsnp = 5000,
  tag_r2 = 0.6,
  threads = 1,
  out = tempfile())

cis_pancanqtl_cd274_leadpdl1inst_ldproxies <- dplyr::filter(cis_pancanqtl, gene == "CD274" & rsid %in% pdl1inst_proxy$SNP_B)


#3. Look up whether any pancanQTL CD274 eQTLs (at any site) are in LD with the lead PD-L1 instrument
cis_pancanqtl_cd274 <- dplyr::filter(cis_pancanqtl, gene == "CD274")
pancanqtl_cis_vars <- unique(cis_pancanqtl_cd274$rsid)

library(ieugwasr)
library(tidyverse)
library(tidyr)

ld_ukb <- ieugwasr::ld_matrix(
  variants = c(pancanqtl_cis_vars, "rs822341"),
  with_alleles = TRUE,
  bfile = "/data/ukb_bed/UKBB_10K",
  plink_bin = genetics.binaRies::get_plink_binary()
) %>%
### output = signed, non-squared R values
as.data.frame() %>%
tibble::rownames_to_column("rsids") %>%
pivot_longer(cols = -c(rsids)) %>%
dplyr::filter(!rsids == name) %>%
dplyr::mutate(r2 = value^2) %>%
dplyr::filter(rsids == "rs822341_T_C") %>%
dplyr::rename("lead_inst_rsid" = "rsids", "pancaneqtl_rsid" = "name", "r_ukb" = "value", "r2_ukb" = "r2")

ld_1kgp <- ld_matrix(
  variants = c(pancanqtl_cis_vars, "rs822341"),
  with_alleles = TRUE,
  bfile = NULL,
  pop = "EUR",
  opengwas_jwt = get_opengwas_jwt(),
  plink_bin = genetics.binaRies::get_plink_binary()
) %>%
as.data.frame() %>%
tibble::rownames_to_column("rsids") %>%
pivot_longer(cols = -c(rsids)) %>%
dplyr::filter(!rsids == name) %>%
dplyr::mutate(r2 = value^2) %>%
dplyr::filter(rsids == "rs822341_T_C") %>%
dplyr::rename("lead_inst_rsid" = "rsids", "pancaneqtl_rsid" = "name", "r_1kgp_eur" = "value", "r2_1kgp_eur" = "r2")

ld <- dplyr::full_join(ld_ukb, ld_1kgp)

## range ld w lead PD-l1 instrument
ld[which.min(ld$r2_ukb), ]
ld[which.max(ld$r2_ukb), ]

write.table(ld, "pancanqtl_variants_ldwleadpdl1inst_041225.csv", sep = ",", quote = F, row.names = F)


#4. MR: tumour CD274 gene expression-cancer survival (for same respective sites)
sites <- c("BRCA", "COAD", "READ", "LUAD", "LUSC", "SKCM", "OV", "PRAD")
cd274_cispancanqtl_sites <- dplyr::filter(cis_pancanqtl, gene == "CD274" & cancer_type %in% sites)

## find lead cis CD274 eQTL for each included site (which has at least one cis CD274 eQTL)
cd274_cispancanqtl_site_lead <- cis_pancanqtl_cd274 %>% 
             group_by(cancer_type) %>%
             filter(`p-value` == min(`p-value`)) %>%
             arrange()
cd274_cispancanqtl_site_lead$rsid

## try to extract lead cis CD274 eQTL for each included site from respective cancer survival GWAS

## check whether cis eqtls present in cancer survival summary statistics
breast_sumstats <- data.table::fread("cancer_sumstats/breastcancer_sumstats/gwas_summary_estimates_all_patients_15_years_RSIDs.txt",
sep = ",", header = T) %>%
dplyr::filter(rsid %in% dplyr::filter(cd274_cispancanqtl_sites, cancer_type == "BRCA")$rsid)

lung_sumstats <- data.table::fread("cancer_sumstats/lungcancer_sumstats/meta_analysis/METAANALYSIS_survival_allstages_ilcco_dfci_ge_mgi1_formatted.txt",
sep = "\t", header = T)
rsids <- data.table::fread("UKBB/UKBB_pQTL_rsids/olink_rsids/olink_rsid_map_mac5_info03_b0_7_chr9_patched_v2.tsv", header = T)
lung_sumstats <- dplyr::filter(lung_sumstats, chr == 9) %>%
dplyr::inner_join(rsids, by = c("bp" = "POS38")) %>%
dplyr::filter(rsid %in% dplyr::filter(cd274_cispancanqtl_sites, cancer_type == "LUSC")$rsid)

ovarian_sumstats <- data.table::fread("cancer_sumstats/OCAC_sumstats/ovariansurvival_sumstats_PDL1_190925.csv", sep = ",", header= T) %>%
dplyr::filter(rsid %in% dplyr::filter(cd274_cispancanqtl_sites, cancer_type == "OV")$rsid)
#dplyr::filter(chr == 9) %>%
#dplyr::filter(position_b37 == dplyr::filter(cd274_cispancanqtl_sites, cancer_type == "OV")$position)

system("mkdir ici_analyses/TCGA")
write.table(cd274_cispancanqtl_sites$rsid, "ici_analyses/TCGA/TCGA_cis_CD274_eqtls_130126.txt", quote=F, row.names=F)
system("plink --bfile ukb_bed/UKBB_10K --extract ici_analyses/TCGA/TCGA_cis_CD274_eqtls_130126.txt --freq --out ici_analyses/TCGA/TCGA_cis_CD274_eqtls_maf_130126")

cd274_cispancanqtl_sites_maf <- read.table("ici_analyses/TCGA/TCGA_cis_CD274_eqtls_maf_130126.frq", header = T)


## ld proxies of (2nd) lead LUSC eqtl
rs10975588_proxy <- gwasvcf::get_ld_proxies(
  rsid = "rs10975588",
  bfile = "/data/ukb_bed/UKBB_10K",
  searchspace = NULL,
  tag_kb = 5000,
  tag_nsnp = 5000,
  tag_r2 = 0.6,
  threads = 1,
  out = tempfile())
rs10975588_proxy[which.max(rs10975588_proxy$R), ]$SNP_B
rs10975588_proxy[which.max(rs10975588_proxy$R), ]$R^2

lung_sumstats <- data.table::fread("cancer_sumstats/lungcancer_sumstats/meta_analysis/METAANALYSIS_survival_allstages_ilcco_dfci_ge_mgi1_formatted.txt",
sep = "\t", header = T)
rsids <- data.table::fread("ici_analyses/UKBB/UKBB_pQTL_rsids/olink_rsids/olink_rsid_map_mac5_info03_b0_7_chr9_patched_v2.tsv", header = T)
lung_sumstats <- dplyr::filter(lung_sumstats, chr == 9) %>%
dplyr::inner_join(rsids, by = c("bp" = "POS38")) %>%
dplyr::filter(rsid==rs10975588_proxy[which.max(rs10975588_proxy$R), ]$SNP_B)

system(paste0("plink --bfile ukb_bed/UKBB_10K --ld rs118166192 rs822341"))
system(paste0("plink --bfile ukb_bed/UKBB_10K --ld rs17804441 rs822341"))

library(gwasvcf)
gwasvcf::set_plink(path = "")
gwasvcf::check_plink()
rs822337_proxy <- gwasvcf::get_ld_proxies(
  rsid = "rs822337",
  bfile = "/data/ukb_bed/UKBB_10K",
  searchspace = NULL,
  tag_kb = 5000,
  tag_nsnp = 5000,
  tag_r2 = 0.6,
  threads = 1,
  out = tempfile())

breast_sumstats <- data.table::fread("cancer_sumstats/breastcancer_sumstats/gwas_summary_estimates_all_patients_15_years_RSIDs.txt",
sep = ",", header = T) %>%
dplyr::filter(rsid %in% dplyr::filter(cd274_cispancanqtl_sites, cancer_type == "BRCA")$rsid | rsid %in% rs822337_proxy$SNP_B) %>%
dplyr::mutate(outcome = "breast_cancer_survival")

brca_leadinst <- dplyr::filter(cd274_cispancanqtl_sites, cancer_type == "BRCA") %>%
dplyr::mutate(effect_allele = stringr::str_split(alleles, "/", simplify = TRUE)[ , 2]) %>%
dplyr::mutate(other_allele = stringr::str_split(alleles, "/", simplify = TRUE)[ , 1]) %>%
dplyr::mutate(se = beta/(abs(`t-stat`))) %>%
as.data.frame() %>%
TwoSampleMR::format_data(type = "exposure", snp_col = "rsid", pval_col = "p-value", phenotype = "BRCA CD274 gene expression")

breast_outdat <- TwoSampleMR::format_data(as.data.frame(breast_sumstats), type = "outcome", 
snp_col = "rsid",
beta_col = "Beta",
se_col = "SE",
pval_col = "P",
effect_allele_col = "a1",
other_allele_col = "a0",
eaf_col = "exp_freq_a1_iCOGS",
phenotype = "breast cancer survival")

eqtl_breastsurv_harm <- TwoSampleMR::harmonise_data(brca_leadinst, breast_outdat) %>%
dplyr::filter(mr_keep == TRUE)
eqtl_breastsurv_harm <- eqtl_breastsurv_harm[which.min(eqtl_breastsurv_harm$`pval.exposure`), ]  %>%
dplyr::mutate(exposure = "BRCA CD274 gene expression", outcome = "breast cancer survival")
eqtl_breastsurv_mrres <- TwoSampleMR::mr(eqtl_breastsurv_harm)

## lead BRCA eQTL LD proxy:
brca_leadinst_proxy <- dplyr::filter(brca_leadinst, SNP == "rs822337")
breast_outdat_proxy <- dplyr::filter(breast_outdat, SNP == rs822337_proxy$SNP_B)
brca_leadinst_proxy$SNP <- breast_outdat_proxy$SNP
brca_leadinst_proxy$effect_allele.exposure <- breast_outdat_proxy$effect_allele.outcome
brca_leadinst_proxy$other_allele.exposure <- breast_outdat_proxy$other_allele.outcome
eqtl_breastsurv_harm_proxy <- TwoSampleMR::harmonise_data(brca_leadinst_proxy, breast_outdat_proxy)
eqtl_breastsurv_harm_proxy_mrres <- TwoSampleMR::mr(eqtl_breastsurv_harm_proxy) %>%
dplyr::mutate(exposure = "BRCA CD274 gene expression", outcome = "breast cancer survival")

## usable LUSC eQTL
lusc_leadinst <- dplyr::filter(cd274_cispancanqtl_sites, cancer_type == "LUSC") %>% 
dplyr::mutate(effect_allele = stringr::str_split(alleles, "/", simplify = TRUE)[ , 2]) %>%
dplyr::mutate(other_allele = stringr::str_split(alleles, "/", simplify = TRUE)[ , 1]) %>%
dplyr::mutate(se=beta/(abs(`t-stat`))) %>%
as.data.frame() %>%
TwoSampleMR::format_data(type = "exposure", snp_col = "rsid", pval_col = "p-value")
lusc_leadinst_se <- dplyr::mutate(cd274_cispancanqtl_sites, `se.exposure`=beta/(abs(`t-stat`))) %>%
dplyr::rename("SNP" = "rsid") %>%
dplyr::select(c(SNP, `se.exposure`))
lusc_leadinst <- dplyr::select(lusc_leadinst, -c(`se.exposure`, `mr_keep.exposure`)) %>%
dplyr::inner_join(lusc_leadinst_se)
lusc_leadinst$`mr_keep.exposure` <- TRUE

lung_sumstats <- data.table::fread("cancer_sumstats/lungcancer_sumstats/meta_analysis/METAANALYSIS_survival_allstages_ilcco_dfci_ge_mgi1_formatted.txt",
sep = "\t", header = T)
rsids <- data.table::fread("ici_analyses/UKBB/UKBB_pQTL_rsids/olink_rsids/olink_rsid_map_mac5_info03_b0_7_chr9_patched_v2.tsv", header = T)
lung_sumstats <- dplyr::filter(lung_sumstats, chr == 9) %>%
dplyr::inner_join(rsids, by = c("bp" = "POS38")) %>%
dplyr::filter(rsid %in% lusc_leadinst$SNP) %>%
dplyr::mutate(outcome = "lung_cancer_survival")
lung_sumstats$chr <- as.numeric(lung_sumstats$chr)

lung_sumstats <- TwoSampleMR::format_data(as.data.frame(lung_sumstats), type = "outcome", snp_col = "rsid", beta_col = "Effect",
se_col = "StdErr", pval_col = "P-value", effect_allele_col = "Allele1", other_allele_col = "Allele2", eaf_col = "Freq1")

eqtl_lungsurv_harm <- TwoSampleMR::harmonise_data(lusc_leadinst, lung_sumstats)

rs10975588_proxy <- gwasvcf::get_ld_proxies(
  rsid = "rs10975588",
  bfile = "/data/ukb_bed/UKBB_10K",
  searchspace = NULL,
  tag_kb = 5000,
  tag_nsnp = 5000,
  tag_r2 = 0.6,
  threads = 1,
  out = tempfile())

lung_sumstats <- data.table::fread("cancer_sumstats/lungcancer_sumstats/meta_analysis/METAANALYSIS_survival_allstages_ilcco_dfci_ge_mgi1_formatted.txt",
sep = "\t", header = T)
lung_sumstats_rs10975588_proxy <- dplyr::filter(lung_sumstats, chr == 9) %>%
dplyr::inner_join(rsids, by = c("bp" = "POS38")) %>%
dplyr::filter(rsid == rs10975588_proxy[which.max(rs10975588_proxy$R), ]$SNP_B) %>%
as.data.frame() %>%
TwoSampleMR::format_data(type = "outcome", snp_col = "rsid", beta_col = "Effect",
se_col = "StdErr", pval_col = "P-value", effect_allele_col = "Allele1", other_allele_col = "Allele2", eaf_col = "Freq1")

lusc_leadinst_rs10975588_proxy <- dplyr::filter(lusc_leadinst, SNP == "rs10975588")
lusc_leadinst_rs10975588_proxy$SNP <- rs10975588_proxy[which.max(rs10975588_proxy$R), ]$SNP_B
lusc_leadinst_rs10975588_proxy$`effect_allele.exposure` <- rs10975588_proxy[which.max(rs10975588_proxy$R), ]$B1
lusc_leadinst_rs10975588_proxy$`other_allele.exposure`	<- rs10975588_proxy[which.max(rs10975588_proxy$R), ]$B2  

eqtl_lungsurv_harm_rs10975588_proxy <- TwoSampleMR::harmonise_data(lusc_leadinst_rs10975588_proxy, lung_sumstats_rs10975588_proxy) %>%
dplyr::mutate(exposure = "LUSC CD274 gene expression", outcome = "lung cancer survival")

#eqtl_lungsurv_harm_rs10975588_proxy_mrres <- TwoSampleMR::mr(eqtl_lungsurv_harm_rs10975588_proxy)

eqtl_harm <- dplyr::full_join(eqtl_breastsurv_harm, eqtl_lungsurv_harm_rs10975588_proxy)

#TwoSampleMR::mr(eqtl_harm)

rsids <- data.table::fread("ici_analyses/UKBB/UKBB_pQTL_rsids/olink_rsids/olink_rsid_map_mac5_info03_b0_7_chr9_patched_v2.tsv", header = T)
pdl1_chr9_sumstats <- data.table::fread("ici_analyses/UKBB/PDL1/CD274_Q9NZQ7_OID20966_v1_Neurology_COMBINED/CD274_sumstats/combined_chr9_CD274:Q9NZQ7:OID20966:v1:Neurology", header = T) %>%
dplyr::inner_join(rsids) %>%
dplyr::filter(rsid %in% eqtl_harm$SNP) %>%
dplyr::mutate(outcome = "PDL1_protein_level") %>%
dplyr::select(c(outcome, CHROM, POS38, BETA, SE, LOG10P, ALLELE0, ALLELE1, A1FREQ, rsid))

pdl1_chr9_sumstats_exp <- TwoSampleMR::format_data(as.data.frame(pdl1_chr9_sumstats), snp_col = "rsid", beta_col = "BETA", se_col = "SE", pval_col = "LOG10P", log_pval = TRUE, 
effect_allele_col = "ALLELE1", other_allele_col = "ALLELE0")

lung_sumstats_pdl1 <- dplyr::filter(lung_sumstats, chr == 9) %>%
dplyr::inner_join(rsids, by = c("bp" = "POS38")) %>%
dplyr::filter(rsid %in% dplyr::filter(eqtl_harm, outcome == "lung cancer survival")$SNP) %>%
as.data.frame() %>%
TwoSampleMR::format_data(type = "outcome", snp_col = "rsid", beta_col = "Effect",
se_col = "StdErr", pval_col = "P-value", effect_allele_col = "Allele1", other_allele_col = "Allele2", eaf_col = "Freq1")

pdl1_harm_lung <- TwoSampleMR::harmonise_data(pdl1_chr9_sumstats_exp, lung_sumstats_pdl1) %>%
dplyr::mutate(exposure = "plasma PD-L1 level", outcome = "lung cancer survival")

system(paste0("plink --bfile ukb_bed/UKBB_10K --ld rs118166192 rs822341"))

breast_sumstats_pdl1 <- breast_sumstats <- data.table::fread("cancer_sumstats/breastcancer_sumstats/gwas_summary_estimates_all_patients_15_years_RSIDs.txt",
sep = ",", header = T) %>%
dplyr::filter(rsid %in% dplyr::filter(eqtl_harm, outcome == "breast cancer survival")$SNP) %>%
as.data.frame() %>%
TwoSampleMR::format_data(type = "outcome", snp_col = "rsid",
beta_col = "Beta",
se_col = "SE",
pval_col = "P",
effect_allele_col = "a1",
other_allele_col = "a0",
eaf_col = "exp_freq_a1_iCOGS",
phenotype = "breast cancer survival")

pdl1_harm_breast <- TwoSampleMR::harmonise_data(pdl1_chr9_sumstats_exp, breast_sumstats_pdl1) %>%
dplyr::mutate(exposure = "plasma PD-L1 level", outcome = "breast cancer survival")

pqtl_harm <- dplyr::full_join(pdl1_harm_lung, pdl1_harm_breast)

harm <- dplyr::full_join(eqtl_harm, pqtl_harm)
eqtl_pqtl_res <- TwoSampleMR::mr(harm)
eqtl_pqtl_res$b <- eqtl_pqtl_res$b*-1
eqtl_pqtl_res$hr <- exp(eqtl_pqtl_res$b)


write.table(harm, "/tcga_pancanqtl/pancanqtl_cd274_vars_survival_harmdat_061225.csv", sep = ",", quote = F, row.names = F)
write.table(eqtl_pqtl_res, "/tcga_pancanqtl/pancanqtl_cd274_vars_survival_mrres_061225.csv", sep = ",", quote = F, row.names = F)
write.table(eqtl_breastsurv_harm_proxy_mrres, "/tcga_pancanqtl/pancanqtl_cd274_vars_survival_mrres_061225_leadBRCAproxy.csv", sep = ",", 
quote = F, row.names = F)
