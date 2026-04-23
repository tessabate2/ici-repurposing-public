library(ieugwasr)
library(tidyverse)

library(here)
readRenviron(here("config.env"))
resultsdir <- Sys.getenv("resultsdir")
setwd(resultsdir)

#ieugwasr::batches()

##default pval
#pd1_lead_phewas <- ieugwasr::phewas(variants = "rs75960776", pval = 1e-05)
#pdl1_lead_phewas <- ieugwasr::phewas(variants = "rs822341", pval = 1e-05)

##n comparisons = n GWAS on OpenGWAS - outcomes on OpenGWAS already used
##outcomes on OpenGWAS used in previous analyses = ovarian cancer survival (1:ieu-a-1234), colorectal cancer risk (2:ebi-a-GCST90018808),...
##..., lung cancer risk (3:ieu-a-987), ovarian cancer risk (4:ieu-a-1120), prostate cancer risk (5:ieu-b-85) 
ids_exclude <- c("ieu-a-1234", "ebi-a-GCST90018808", "ieu-a-987", "ieu-a-1120", "ieu-b-85")
n_gwas <- ieugwasr::gwasinfo()

pd1_pdl1_proteins <- dplyr::filter(n_gwas, grepl('PD-1|PD-L1|Programmed cell death protein 1|Programmed cell death 1 ligand 1|Programmed Cell Death 1 Protein|Programmed Death Ligand 1', trait))$id

n_gwas <- n_gwas[!n_gwas$id %in% ids_exclude, ]
n_gwas <- n_gwas[!n_gwas$id %in% pd1_pdl1_proteins, ]

##bonferroni correction (overly stringent so not using results)
#pd1_lead_phewas_bfcor <- ieugwasr::phewas(variants = "rs75960776", pval = (5e-08/length(unique(n_gwas$id))))
#pdl1_lead_phewas_bfcor <- ieugwasr::phewas(variants = "rs822341",pval = (5e-08/length(unique(n_gwas$id))))


##FDR correction
library(stats)
### PD-1
pd1_lead_phewas_fdrcor <- ieugwasr::phewas(variants = "rs75960776", pval = 1e-04)
pd1_lead_phewas_fdrcor$fdr_p <- stats::p.adjust(pd1_lead_phewas_fdrcor$p, 
                                         method = "fdr", 
                                         n = length(unique(n_gwas$id))) 
pd1_lead_phewas_fdrcor <- pd1_lead_phewas_fdrcor[pd1_lead_phewas_fdrcor$id %in% n_gwas$id, ]
dec_prot_PD1 <- read.table("../../data/ici_analyses/UKBB/PD1/PDCD1_Q15116_OID21396_v1_Oncology_COMBINED/500kb_p5e06_r0.001/PD1ovariansurvival500kb_p5e06r0.001inst_100525.csv", 
sep = ",",header=T) %>%
dplyr::filter(rsid == "rs75960776") %>%
dplyr::mutate(beta_dec = case_when(BETA > 0 ~ BETA*-1, BETA < 0	~ BETA, T ~ NA)) %>%
dplyr::mutate(allele1_dec = case_when(BETA > 0 ~ ALLELE0, BETA < 0 ~ ALLELE1, T ~ NA)) %>%
dplyr::mutate(allele0_dec = case_when(BETA > 0 ~ ALLELE1, BETA < 0 ~ ALLELE0, T ~ NA)) %>%
dplyr::select(-c(BETA, ALLELE0, ALLELE1)) %>%
dplyr::rename(BETA = beta_dec, ALLELE0 = allele0_dec, ALLELE1 = allele1_dec)

pd1_lead_phewas_fdrcor_decPD1 <- dplyr::mutate(pd1_lead_phewas_fdrcor, beta_decprotein = case_when((dec_prot_PD1$ALLELE1 == ea & dec_prot_PD1$ALLELE0 == nea) ~ beta, 
(dec_prot_PD1$ALLELE1 == nea & dec_prot_PD1$ALLELE0 == ea) ~ beta*-1, T~NA)) %>%
dplyr::mutate(a1_decprotein = case_when((dec_prot_PD1$ALLELE1 == ea) ~ ea, (dec_prot_PD1$ALLELE1 == nea) ~ nea)) %>%
dplyr::mutate(a0_decprotein = case_when((dec_prot_PD1$ALLELE0 == ea) ~ ea, (dec_prot_PD1$ALLELE0 == nea) ~ nea)) %>%
dplyr::mutate(eaf_decprotein = case_when((dec_prot_PD1$ALLELE1 == ea) ~ eaf, (dec_prot_PD1$ALLELE1 == nea) ~ (1-eaf))) %>%
dplyr::select(-c(beta, ea, nea, eaf)) %>%
dplyr::rename(beta = beta_decprotein, ea = a1_decprotein, nea = a0_decprotein, eaf = eaf_decprotein)
pd1_lead_phewas_fdrcor_decPD1$instrument <- c("PD1")
#write.table(pd1_lead_phewas_fdrcor, "leadpd1inst_phewas_090625.csv", sep = ",", quote = F, row.names = F)

### PD-L1
pdl1_lead_phewas_fdrcor <- ieugwasr::phewas(variants = "rs822341", pval = 1e-04)
pdl1_lead_phewas_fdrcor$fdr_p <- stats::p.adjust(pdl1_lead_phewas_fdrcor$p, 
                                         method = "fdr", 
                                         n = length(unique(n_gwas$id)))
pdl1_lead_phewas_fdrcor <- pdl1_lead_phewas_fdrcor[pdl1_lead_phewas_fdrcor$id %in% n_gwas$id, ]

dec_prot_PDL1 <- read.table("../../data/ici_analyses/UKBB/PDL1/CD274_Q9NZQ7_OID20966_v1_Neurology_COMBINED/500kb_p5e06_r0.001/PDL1ovariansurvival500kb_p5e06r0.001inst_100525.csv",
sep = ",",header=T) %>%
dplyr::filter(rsid == "rs822341") %>%
dplyr::mutate(beta_dec = case_when(BETA	> 0 ~ BETA*-1, BETA < 0 ~ BETA, T ~ NA)) %>%
dplyr::mutate(allele1_dec = case_when(BETA > 0 ~ ALLELE0, BETA < 0 ~ ALLELE1, T ~ NA)) %>%
dplyr::mutate(allele0_dec = case_when(BETA > 0 ~ ALLELE1, BETA < 0 ~ ALLELE0, T ~ NA)) %>%
dplyr::select(-c(BETA, ALLELE0, ALLELE1)) %>%
dplyr::rename(BETA = beta_dec, ALLELE0 = allele0_dec, ALLELE1 = allele1_dec)

pdl1_lead_phewas_fdrcor_decPDL1 <- dplyr::mutate(pdl1_lead_phewas_fdrcor, beta_decprotein = case_when((dec_prot_PDL1$ALLELE1 == ea & dec_prot_PDL1$ALLELE0 == nea) ~ beta,
(dec_prot_PDL1$ALLELE1 == nea & dec_prot_PDL1$ALLELE0 == ea) ~ beta*-1, T~NA)) %>%
dplyr::mutate(a1_decprotein = case_when((dec_prot_PDL1$ALLELE1 == ea) ~ ea, (dec_prot_PDL1$ALLELE1 == nea) ~ nea)) %>%
dplyr::mutate(a0_decprotein = case_when((dec_prot_PDL1$ALLELE0 == ea) ~ ea, (dec_prot_PDL1$ALLELE0 == nea) ~ nea)) %>%
dplyr::mutate(eaf_decprotein = case_when((dec_prot_PDL1$ALLELE1 == ea) ~ eaf, (dec_prot_PDL1$ALLELE1 == nea) ~ (1-eaf))) %>%
dplyr::select(-c(beta, ea, nea, eaf)) %>%
dplyr::rename(beta = beta_decprotein, ea = a1_decprotein, nea = a0_decprotein, eaf = eaf_decprotein)
pdl1_lead_phewas_fdrcor_decPDL1$instrument <- c("PDL1")
#write.table(pdl1_lead_phewas_fdrcor, "leadpdl1inst_phewas_090625.csv", sep = ",",	quote =	F, row.names = F)

pd1_pdl1_phewas_fdr <- dplyr::full_join(pd1_lead_phewas_fdrcor_decPD1, pdl1_lead_phewas_fdrcor_decPDL1) %>%
dplyr::filter(fdr_p <0.05) %>%
dplyr::select(c(instrument, rsid, chr, position, id, trait, ea, nea, eaf, n, beta, se, p, fdr_p))

write.table(pd1_pdl1_phewas_fdr, "lead_pd1_pdl1_instphewas_271025.csv", sep = ",", quote = F, row.names = F)
