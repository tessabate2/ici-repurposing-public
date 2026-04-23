# Load required packages
library(tidyverse)
library(data.table)
library(here)
readRenviron(here("config.env"))

datadir <- Sys.getenv("datadir")
setwd(datadir)
setwd("ici_analyses/dice_t_cell_vcf")

tcell_pops <- c("CD4_memory", "CD4_naive", "CD4_naive_activated", "CD4_naive_Treg",
                "CD4_TFH", "CD4_TH1", "CD4_TH1_TH17", "CD4_TH17", "CD4_TH2", "CD8_naive",
                "CD8_naive_activated")

vcf_res <- data.frame()

for (tcell_pop in tcell_pops) {
if (tcell_pop == "CD4_memory") {
vcf_file <- c("TREG_MEM.vcf")
} else if (tcell_pop == "CD4_naive") {
vcf_file <- c("CD4_NAIVE.vcf")
} else if (tcell_pop == "CD4_naive_activated") {
vcf_file <- c("CD4_STIM.vcf")
} else if (tcell_pop == "CD4_naive_Treg") {
vcf_file <- c("TREG_NAIVE.vcf") 
} else if (tcell_pop == "CD4_TFH") {
vcf_file <- c("TFH.vcf") 
} else if (tcell_pop == "CD4_TH1") {
vcf_file <- c("TH1.vcf")
} else if (tcell_pop == "CD4_TH1_TH17") {
vcf_file <- c("THSTAR.vcf")
} else if (tcell_pop == "CD4_TH17") {
vcf_file <- c("TH17.vcf")
} else if (tcell_pop == "CD4_TH2") {
vcf_file <- c("TH2.vcf")
} else if (tcell_pop == "CD8_naive") {
vcf_file <- c("CD8_NAIVE.vcf")
} else if (tcell_pop == "CD8_naive_activated") {
vcf_file <- c("CD8_STIM.vcf")
}

vcf <- data.table::fread(vcf_file, skip = 6, sep = "\t") %>%
  dplyr::filter(ID == "rs75960776") %>%
  dplyr::filter(grepl("ENSG00000188389", INFO))

vcf <-  separate(data = vcf, col = INFO,
         sep = ";", into = c("Gene", "GeneSymbol", "Pvalue", "Beta", "Statistic", "FDR")) %>%
  dplyr::mutate(Gene = gsub(".*\\=", "", Gene)) %>%
  dplyr::mutate(GeneSymbol = gsub(".*\\=", "", GeneSymbol)) %>%
  dplyr::mutate(Pvalue = gsub(".*\\=", "", Pvalue)) %>%
  dplyr::mutate(Beta = gsub(".*\\=", "", Beta)) %>%
  dplyr::mutate(Statistic = gsub(".*\\=", "", Statistic)) %>%
  dplyr::mutate(FDR = gsub(".*\\=", "", FDR))
vcf$cell_population <- tcell_pop

vcf_res <- rbind(vcf_res, vcf)
#print(vcf_res)
}
write.table(vcf_res, "dice_pdcd1_leadpd1inst_cellpops_181125.csv",
            sep = ",", quote = F, row.names = F)


## format output
vcf_res <- read.table("dice_pdcd1_leadpd1inst_cellpops_181125.csv", sep = ",", header = F)
colnames(vcf_res) <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","Gene","GeneSymbol","Pvalue","Beta","Statistic","FDR","cell_population")

vcf_res <- dplyr::mutate(vcf_res, diceTcell_pop = case_when(cell_population == "CD4_memory" ~ "T cell, CD4, memory TREG",
cell_population == "CD4_naive" ~ "T cell, CD4, naive",
cell_population == "CD4_naive_activated" ~ "T cell, CD4, naive [activated]",
cell_population == "CD4_naive_Treg" ~ "T cell, CD4, naive TREG",
cell_population == "CD4_TFH" ~ "T cell, CD4, TFH",
cell_population == "CD4_TH1" ~ "T cell, CD4, TH1",
cell_population == "CD4_TH1_TH17" ~ "T cell, CD4, TH1/17",
cell_population == "CD4_TH17" ~ "T cell, CD4, TH17",
cell_population == "CD4_TH2" ~ "T cell, CD4, TH2",
cell_population == "CD8_naive" ~ "T cell, CD8, naive",
cell_population == "CD8_naive_activated" ~ "T cell, CD8, naive [activated]")) %>%
dplyr::mutate(ref = case_when(REF == "TRUE" ~ "T")) %>%
dplyr::select(c(diceTcell_pop, ID, ref, ALT, Gene, GeneSymbol, Beta, Pvalue, FDR)) %>%
dplyr::rename("DICE T cell population" = "diceTcell_pop", "rsid" = "ID", "REF" = "ref")

resultsdir <- Sys.getenv("resultsdir")
setwd(resultsdir)
write.table(vcf_res, "dice_unfiltered_res_201125.tsv", sep="\t", quote = F, row.names = F)



