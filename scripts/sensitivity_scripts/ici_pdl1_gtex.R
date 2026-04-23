
# Load required packages
library(tidyverse)
library(data.table)

library(here)
readRenviron(here("config.env"))
gtexdatadir <- Sys.getenv("gtexdatadir")

datadir <- Sys.getenv("datadir")
setwd(datadir)
setwd("ici_analyses/UKBB")

protein <- c("PDL1")
survivals <- c("breast", "colorectal", "lung", "melanoma", "ovarian", "prostate")
window <- c(500)
pval <- c(5e-6)
r <- c(0.3)
measure <- c("survival")

# Create data frame for all results to be stored in
gtex_eQTL_protein_leadsnp_all <- data.frame()

## 1. Instrument:
if (protein=="PDL1") {
path <- c(paste0("PDL1/CD274_Q9NZQ7_OID20966_v1_Neurology_COMBINED/", window, "kb_p", pval, "_r",r))
ensg_id <- c("ENSG00000120217")
rsids <- read.table("UKBB_pQTL_rsids/olink_rsids/olink_rsid_map_mac5_info03_b0_7_chr9_patched_v2.tsv", header = T)
lead_snp <- c("rs822341")
}

### 1.1: Load outcome GWAS summary statistics depending on which cancer site and outcome measure required (all formatted to have same columns)
# Loop over cancer sites 
for (survival in survivals) {

if (survival == "prostate") {
gtex_tissues <- c("Prostate")

### Load breast cancer survival GWAS summary statistics
} else if (survival=="breast") {
gtex_tissues <- c("Breast_Mammary_Tissue")

### Load ovarian cancer survival GWAS summary statistics
} else if (survival== "ovarian") {
gtex_tissues <- c("Ovary")

} else if (survival == "lung") {
gtex_tissues <- c("Lung")

} else if (survival == "melanoma") {
gtex_tissues <- c("Skin_Not_Sun_Exposed_Suprapubic", "Skin_Sun_Exposed_Lower_leg")

} else if (survival == "colorectal") {
gtex_tissues <- c("Colon_Sigmoid", "Colon_Transverse")

}

setwd(gtexdatadir)

for (tissue in gtex_tissues) {
gtex_eQTL <- data.table::fread(paste0(tissue, ".allpairs.txt"), sep = "\t", header=T)

gtex_eQTL_protein <- dplyr::filter(gtex_eQTL, grepl(ensg_id, gene_id)) %>%
separate(variant_id, into = c('chr', 'pos', 'ref', 'alt', 'build'), '_', convert = TRUE) %>%
dplyr::full_join(rsids, by = c("pos" = "POS38"), relationship = "many-to-many") %>%
dplyr::filter((ref == REF & alt == ALT) | (ref == ALT & alt ==REF))

gtex_eQTL_protein_leadsnp <- dplyr::filter(gtex_eQTL_protein, rsid == lead_snp)

gtex_eQTL_protein_leadsnp$protein <- protein
gtex_eQTL_protein_leadsnp$tissue <- tissue

gtex_eQTL_protein_leadsnp_all <- rbind(gtex_eQTL_protein_leadsnp_all, gtex_eQTL_protein_leadsnp)

}}

resultsdir <- Sys.getenv("resultsdir")
setwd(resultsdir)

write.csv(gtex_eQTL_protein_leadsnp_all, "PDL1_gtex_assoc_130625.csv", row.names=F, quote=F)

gtex_res <- read.table("PDL1_gtex_assoc_130625.csv", sep = ",", header = T) %>%
dplyr::mutate(gtex_beta = case_when(!(alt == ALT) ~ slope*-1,
                                    (alt == ALT) ~ slope)) %>%
  dplyr::filter(!tissue == "Skin_Not_Sun_Exposed_Suprapubic") %>%
  dplyr::select(c(rsid, ref, alt, maf, gtex_beta, slope_se, pval_nominal, protein,
                  tissue)) %>%
  dplyr::rename(beta = gtex_beta, se = slope_se) %>%
  dplyr::mutate(comparison_cancersite = case_when(tissue == "Breast_Mammary_Tissue" ~ "breast",
                                                  tissue == "Colon_Sigmoid" ~ "colorectal",
                                                  tissue == "Colon_Transverse" ~ "colorectal",
                                                  tissue == "Lung" ~ "lung",
                                                  tissue == "Skin_Sun_Exposed_Lower_leg" ~"melanoma",
                                                  tissue == "Ovary" ~ "ovarian",
                                                  tissue == "Prostate" ~ "prostate"))

gtex_res$tissue <- sub("_", " ", gtex_res$tissue)
gtex_res$tissue <- sub("_", " ", gtex_res$tissue)
gtex_res$tissue <- sub("_", " ", gtex_res$tissue)
gtex_res$tissue <- sub("_", " ", gtex_res$tissue)

myColors <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
names(myColors) <- levels(gtex_res$comparison_cancersite)
colScale <- scale_colour_manual(name = "class",values = myColors)
constant_for95CI <- stats::qnorm(1 - (1 - 0.95) / 2)

gtex_res$tissue <- factor(gtex_res$tissue,
                          levels = c("Prostate", "Ovary", "Skin Sun Exposed Lower leg",
                                     "Lung", "Colon Transverse", "Colon Sigmoid",
                                     "Breast Mammary Tissue"))

pdf("PDL1_gtex_assoc_130625.pdf", width = 10, height = 10)
gtex_forplot <- ggplot(data = gtex_res) +
  geom_point(aes(x = beta, y = tissue, colour = comparison_cancersite),
             show.legend = F) +
  geom_vline(xintercept = 0, linetype = "dotted")+
  geom_errorbarh(aes(y = tissue, xmin = beta-(constant_for95CI*se), 
                     xmax = beta+(constant_for95CI*se),
                     height = 0.2, colour = comparison_cancersite), show.legend = F) +
  theme_bw(base_size = 16) +
  xlab("Beta (95% CI)") +
  ylab("Tissue") +
  ggtitle("Lead PD-L1 instrument effect on CD274 gene expression \nin different tissues")
print(gtex_forplot)
dev.off()
