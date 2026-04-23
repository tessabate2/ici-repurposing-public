#module load languages/R/4.3.3
#module load plink/1.9-beta6.27-openblas
#module load  apps/plink1.9/1.90-b77
#module load gcc/12.3.0
#module load cmake/3.27.9

library(tidyverse)
library(data.table)
library(geni.plots)
library(here)

readRenviron(here("config.env"))
datadir <- Sys.getenv("datadir")
setwd(datadir)
setwd("ici_analyses/colocalisation_analyses")
resultsdir <- Sys.getenv("resultsdir")

ici_coloc <- function(protein, cancer_site) {

alleles <- c("A", "C", "G", "T")

if (protein == "PD1") {
PD1_sumstats_chr <- data.table::fread("raw_input_files/combined_chr2_PDCD1_Q15116_OID21396_v1_Oncology", header=T)
chr_sumstats <- data.table::fread("raw_input_files/olink_rsid_map_mac5_info03_b0_7_chr2_patched_v2.tsv", header=T)
PD1_sumstats_chr_rsids <- dplyr::inner_join(PD1_sumstats_chr, chr_sumstats, by = c("ID"))

PD1_inst <- read.table(paste0("../UKBB/PD1/PDCD1_Q15116_OID21396_v1_Oncology_COMBINED/500kb_p5e06_r0.3_240925/PD1500kb_p5e06_r0.3", cancer_site, "survivalharm_minorallele_100525.snplist"), header=F)
colnames(PD1_inst) <- c("rsid")

PD1_sumstats_inst <- PD1_sumstats_chr_rsids[PD1_sumstats_chr_rsids$rsid %in% PD1_inst$rsid, ]
PD1_sumstats_leadinst <- PD1_sumstats_inst[which.max(PD1_sumstats_inst$LOG10P), ]

## restrict to lead PD-1 instrument +/- 500kb
chrom <- c("2")
pos <- 242804707

PD1_sumstats_chr_rsids <- dplyr::filter(PD1_sumstats_chr_rsids, POS19 < pos+500000) %>%
dplyr::filter(POS19 > pos- 500000)

PD1_sumstats_QC <- dplyr::filter(PD1_sumstats_chr_rsids, ALLELE0 %in% alleles) %>%
  dplyr::filter(ALLELE1 %in% alleles) %>%
  dplyr::filter(A1FREQ <= 0.99) %>%
  dplyr::filter(A1FREQ >= 0.01)

PD1_sumstats_leadSNPwindow_regplot <- dplyr::mutate(PD1_sumstats_QC, p = 10^-LOG10P) %>%
                                      dplyr::select(c(rsid, CHROM, POS19, p, REF, ALT)) %>%
                                      dplyr::rename(marker = rsid, chr = CHROM,
                                                    pos = POS19, pvalue_1 = p)

protein_sumstats_QC <- PD1_sumstats_QC

exp_regplot_df <- PD1_sumstats_leadSNPwindow_regplot

lead_SNP <- c("rs75960776")
protein_dis <- c("PD-1")

} else if (protein == "PDL1") {

PDL1_sumstats_chr <- data.table::fread("../UKBB/PDL1/CD274_Q9NZQ7_OID20966_v1_Neurology_COMBINED/CD274_sumstats/combined_chr9_CD274:Q9NZQ7:OID20966:v1:Neurology", header=T)
chr_sumstats <- data.table::fread("../UKBB/UKBB_pQTL_rsids/olink_rsids/olink_rsid_map_mac5_info03_b0_7_chr9_patched_v2.tsv", header=T)
PDL1_sumstats_chr_rsids <- dplyr::inner_join(PDL1_sumstats_chr, chr_sumstats, by = c("ID"))

PDL1_inst <- read.table(paste0("../UKBB/PDL1/CD274_Q9NZQ7_OID20966_v1_Neurology_COMBINED/500kb_p5e06_r0.3_240925/PDL1500kb_p5e06_r0.3", cancer_site, "survivalharm_minorallele_100525.snplist"), header=F)
colnames(PDL1_inst) <- c("rsid")

PDL1_sumstats_inst <- PDL1_sumstats_chr_rsids[PDL1_sumstats_chr_rsids$rsid %in% PDL1_inst$rsid, ]
PDL1_sumstats_leadinst <- PDL1_sumstats_inst[which.max(PDL1_sumstats_inst$LOG10P), ]

## restrict to lead PD-1 instrument +/- 200kb
pos <- 5453396
chrom <- c("9")
PDL1_sumstats_chr_rsids <- dplyr::filter(PDL1_sumstats_chr_rsids, POS19 < pos +500000) %>%
dplyr::filter(POS19 > pos - 500000)

PDL1_sumstats_QC <- dplyr::filter(PDL1_sumstats_chr_rsids, ALLELE0 %in% alleles) %>%
  dplyr::filter(ALLELE1 %in% alleles) %>%
  dplyr::filter(A1FREQ <= 0.99) %>%
  dplyr::filter(A1FREQ >= 0.01)

PDL1_sumstats_leadSNPwindow_regplot <- dplyr::mutate(PDL1_sumstats_QC, p = 10^-LOG10P) %>%
                                      dplyr::select(c(rsid, CHROM, POS19, p, REF, ALT)) %>%
                                      dplyr::rename(marker = rsid, chr = CHROM,
                                                    pos = POS19, pvalue_1 = p)

protein_sumstats_QC <- PDL1_sumstats_QC

exp_regplot_df <- PDL1_sumstats_leadSNPwindow_regplot

lead_SNP <- c("rs822341")
protein_dis <- c("PD-L1")

}

leadSNP_window_pwcoco <- dplyr::mutate(protein_sumstats_QC, p = 10^-LOG10P) %>%
  dplyr::select(c(rsid, ALLELE1, ALLELE0, A1FREQ, BETA, SE, p, N)) %>%
  dplyr::rename(SNP = rsid, effect_allele.exposure = ALLELE1,
                other_allele.exposure = ALLELE0, eaf.exposure = A1FREQ,
                beta.exposure = BETA, se.exposure = SE, pval.exposure = p,
                samplesize.exposure = N)

if (cancer_site=="colorectal") {

cancer_dis <- c("Colorectal")
CRC_sumstats <- data.table::fread("../../cancer_sumstats/crc_PD_L1_PD1SumStats/GWAS_CRC_survival.csv", 
                                   header = T) %>%
dplyr::filter(position < (pos+500000)) %>%
dplyr::filter(position > (pos-500000)) %>%
dplyr::filter(chromosome==chrom)

CRC_sumstats_QC <- dplyr::filter(CRC_sumstats, ref_allele %in% alleles) %>%
  dplyr::filter(alt_allele %in% alleles) %>%
  dplyr::filter(CAF <= 0.99) %>%
  dplyr::filter(CAF >= 0.01) %>%
dplyr::inner_join(chr_sumstats, by = c("position"="POS19"), relationship = "many-to-many") %>%
dplyr::filter((REF==ref_allele&ALT==alt_allele | REF==alt_allele&ALT==ref_allele)) %>%
dplyr::rename(SNP = "rsid")

CRC_sumstats_leadSNPwindow_regplot <- dplyr::select(CRC_sumstats_QC,
                                                     c(SNP, chromosome, position, 
                                                       crc_gp, REF, ALT)) %>%
                                      dplyr::rename(marker = SNP, pos = position,
                                                    pvalue_2 = crc_gp, chr = chromosome)

surv_sumstats_QC <- CRC_sumstats_QC
out_regplot_df <- CRC_sumstats_leadSNPwindow_regplot

surv_leadSNPwindow_pwcoco <- dplyr::select(surv_sumstats_QC,
                                                    c(SNP, ref_allele, alt_allele, CAF,
                                                      crc_gbeta,
                                                      crc_gse,
                                                      crc_gp)) %>%
  dplyr::rename(effect_allele.outcome = alt_allele, other_allele.outcome = ref_allele,
                eaf.outcome = CAF, beta.outcome = crc_gbeta,
                se.outcome = crc_gse, pval.outcome = crc_gp)

} else if (cancer_site == "ovarian") {
cancer_dis <- c("Ovarian")
OCAC_sumstats <- data.table::fread("raw_input_files/ocac_imputed_results_2_242292033_243199373_PD1chr2.csv", 
                                   header = T) %>%
dplyr::filter(position_b37 < pos+500000) %>%
dplyr::filter(position_b37 > pos - 500000) %>%
dplyr::filter(chr==chrom)

OCAC_sumstats_QC <- dplyr::filter(OCAC_sumstats, a0 %in% alleles) %>%
  dplyr::filter(a1 %in% alleles) %>%
  dplyr::filter(ocac_eaf <= 0.99) %>%
  dplyr::filter(ocac_eaf >= 0.01) %>%
  dplyr::rename(SNP = "1KG_Phase1_ID")

OCAC_sumstats_leadSNPwindow_regplot <- dplyr::select(OCAC_sumstats_QC,
                                                     c(SNP, chr, position_b37, 
                                                       ocac_OverallSurvival_pvalue, a0, a1)) %>%
                                      dplyr::rename(marker = SNP, pos = position_b37,
                                                    pvalue_2 = ocac_OverallSurvival_pvalue)   


surv_sumstats_QC <- OCAC_sumstats_QC
out_regplot_df <- OCAC_sumstats_leadSNPwindow_regplot

surv_leadSNPwindow_pwcoco <- dplyr::select(OCAC_sumstats_QC,
                                                    c(SNP, a1, a0, ocac_eaf, 
                                                      ocac_OverallSurvival_OR, 
                                                      ocac_OverallSurvival_se, 
                                                      ocac_OverallSurvival_pvalue)) %>%
  dplyr::rename(effect_allele.outcome = a1, other_allele.outcome = a0,
                eaf.outcome = ocac_eaf, beta.outcome = ocac_OverallSurvival_OR,
                se.outcome = ocac_OverallSurvival_se, pval.outcome = ocac_OverallSurvival_pvalue)

}

regplot_input <- dplyr::inner_join(exp_regplot_df, out_regplot_df,  by = c("marker", "chr", "pos"))

if (cancer_site == "colorectal") {

regplot_input <- dplyr::filter(regplot_input, ((REF.x == ALT.y)&(ALT.x == REF.y)) | ((REF.x == REF.y) &(ALT.x==ALT.y)))

} else if (cancer_site == "ovarian") {

regplot_input <- dplyr::filter(regplot_input, ((REF == a0)&(ALT == a1)) | ((REF == a1) &(ALT==a0)))
}


write.table(regplot_input$marker, paste0(protein, "_", cancer_site, "_regplot_input_SNPs_300925.txt"), sep = "\t", quote=F, row.names=F)

system(paste0("plink --bfile ../../ukb_bed/UKBB_10K --extract ", protein, "_", cancer_site, "_regplot_input_SNPs_300925.txt --write-snplist --out ", protein, "_", cancer_site, "_regplot_input_SNPs_300925"))

system(paste0("plink --bfile ../../ukb_bed/UKBB_10K --extract ", protein, "_", cancer_site, "_regplot_input_SNPs_300925.txt --out ", protein, "_", cancer_site, "_regplot_input_SNPs_300925 --r square"))


### Read list of SNPs which are present in the UK biobank reference panel
SNP_list_UKB <- read.table(paste0(protein, "_", cancer_site, "_regplot_input_SNPs_300925.snplist"), sep = "\t")
colnames(SNP_list_UKB) <- c("SNP")

### Read LD matrix and add SNPs as column/row names
UKBB_LD <- read.table(paste0(protein, "_", cancer_site, "_regplot_input_SNPs_300925.ld"))
colnames(UKBB_LD) <- SNP_list_UKB$SNP
row.names(UKBB_LD) <- SNP_list_UKB$SNP
UKBB_LD_matrix <- data.matrix(UKBB_LD)


### Restrict harmonised data to SNPs present in the UK biobank reference panel (should be in same order as the LD matrix)
regplot_input_ukBorder <- dplyr::inner_join(SNP_list_UKB, regplot_input, 
                                              by = c("SNP" = "marker")) %>%
                                   dplyr::rename(marker = SNP)

stack_region <- list()
stack_region$assoc <- regplot_input_ukBorder
stack_region$corr <- UKBB_LD_matrix

## ld bw lead protein snp and lead survival snp

reg_sumstats <- stack_region$assoc
reg_sumstats[which.min(reg_sumstats$pvalue_2), ]

system(paste0("plink --bfile ../../ukb_bed/UKBB_10K --ld ", reg_sumstats[which.min(reg_sumstats$pvalue_1), ]$marker, " ",reg_sumstats[which.min(reg_sumstats$pvalue_2), ]$marker, " > ", resultsdir, "/colocalisation_results/", protein, "_", cancer_site, "_UKB_ld_300925.txt"))

lead_snps_ld <- plinkr::read_plink_log_file(paste0(resultsdir, "/colocalisation_results/", protein, "_", cancer_site, "_UKB_ld_300925.txt"))
lead_snps_ld <- as.data.frame(lead_snps_ld)
lead_snps_ld <- dplyr::filter(lead_snps_ld, grepl("R-sq = ", lead_snps_ld)==TRUE)

tmp <- sub("*\\   R-sq = [0-9]*", "", lead_snps_ld$lead_snps_ld)

lead_snps_ld_df <- data.frame(Protein = protein, Site = cancer_site, "R_sq" = round(as.numeric(sub("*\\    D' = .*", "", tmp)), 4))


## pwcoco

write.table(leadSNP_window_pwcoco, paste0(protein, "sumstats_leadSNPwindow_PWcocoinput_300925.txt"), sep="\t", row.names=F, quote=F)

write.table(surv_leadSNPwindow_pwcoco, paste0(cancer_site, "sumstats_leadSNPwindow_PWCOCOinput_300925.txt"), sep="\t", row.names=F, quote=F)

system(paste0("pwcoco/build/pwcoco --bfile ../../ukb_bed/UKBB_10K --sum_stats1 ", protein, "sumstats_leadSNPwindow_PWcocoinput_300925.txt --sum_stats2 ", cancer_site, "sumstats_leadSNPwindow_PWCOCOinput_300925.txt --maf 0.01 --out ", protein, "_", cancer_site, "_pwcoco_240925"))

pwcoco_res <- read.table(paste0(protein, "_", cancer_site, "_pwcoco_240925.coloc"), header = T)
pwcoco_res <- pwcoco_res[!duplicated(pwcoco_res), ]

pwcoco_res_df <- dplyr::filter(pwcoco_res, SNP1 == "unconditioned" & SNP2 == "unconditioned")
pwcoco_res_df$Protein <- protein
pwcoco_res_df$Site <- cancer_site

colocalisation_res <- dplyr::inner_join(lead_snps_ld_df, pwcoco_res_df, by = c("Protein", "Site")) %>%
dplyr::select(-c(Dataset1, Dataset2))

pdf(paste0(resultsdir, "/colocalisation_results/regassocplot_", protein, "_", cancer_site, "surv_090326_addH4.pdf"), width=10, height = 12)
regassoc_plot <- geni.plots::fig_region_stack(
  data = stack_region$assoc,
  traits = c(paste0("Plasma ", protein_dis, " level"), paste0(cancer_dis,"  cancer survival")),
  corr = stack_region$corr,
  build = 37,
  highlights = paste0(lead_SNP),
  title_center = TRUE)
print(regassoc_plot + labs(caption = paste0("H4=", round(colocalisation_res$H4, 3))) + theme(plot.caption = element_text(colour = "blue", hjust = 1, size = 20)))
dev.off()

return(colocalisation_res)

}

resultsdir <- Sys.getenv("resultsdir")

pd1_crc_coloc <- ici_coloc(protein = "PD1", cancer_site = "colorectal")
pdl1_crc_coloc <- ici_coloc(protein = "PDL1", cancer_site = "colorectal")
pd1_ov_coloc <- ici_coloc(protein = "PD1", cancer_site = "ovarian")

coloc_all_res <- dplyr::full_join(pd1_crc_coloc, pd1_ov_coloc) %>%
dplyr::full_join(pdl1_crc_coloc) %>%
dplyr::rename("LD R^2 between lead protein and cancer survival SNPs" = "R_sq") %>%
dplyr::select(-c(SNP1, SNP2, nsnps, log_abf_all))

coloc_all_res$H0 <- signif(coloc_all_res$H0, digits = 3)
coloc_all_res$H1 <- signif(coloc_all_res$H1, digits = 3)
coloc_all_res$H2 <- signif(coloc_all_res$H2, digits = 3)
coloc_all_res$H3 <- signif(coloc_all_res$H3, digits = 3)
coloc_all_res$H4 <- signif(coloc_all_res$H4, digits = 3)

write.table(coloc_all_res, paste0(resultsdir, "/protein_cancersurv_colocalisationresults_250925.csv"), sep = ",", quote = F, row.names = F)
