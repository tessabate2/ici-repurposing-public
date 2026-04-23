# Load required packages
library(tidyverse)

library(here)
readRenviron(here("config.env"))

datadir <- Sys.getenv("datadir")
setwd(datadir)
setwd("ici_analyses/UKBB")

resultsdir <- Sys.getenv("resultsdir")

proteins <- c("PD1", "PDL1")

for (protein in proteins) {

if (protein == "PD1") {
setwd("PD1/PDCD1_Q15116_OID21396_v1_Oncology_COMBINED/500kb_p5e06_r0.3_240925")
combined_files <- c("PD1breastsurvival500kb_p5e06r0.3inst_100525.csv", "PD1colorectalsurvival500kb_p5e06r0.3inst_100525.csv", "PD1lungsurvival500kb_p5e06r0.3inst_100525.csv", 
"PD1melanomasurvival500kb_p5e06r0.3inst_100525.csv", "PD1ovariansurvival500kb_p5e06r0.3inst_100525.csv", "PD1prostatesurvival500kb_p5e06r0.3inst_100525.csv")

} else if (protein == "PDL1") {
setwd("PDL1/CD274_Q9NZQ7_OID20966_v1_Neurology_COMBINED/500kb_p5e06_r0.3_240925")
combined_files <- c("PDL1breastsurvival500kb_p5e06r0.3inst_100525.csv", "PDL1colorectalsurvival500kb_p5e06r0.3inst_100525.csv", "PDL1lungsurvival500kb_p5e06r0.3inst_100525.csv",
"PDL1melanomasurvival500kb_p5e06r0.3inst_100525.csv", "PDL1ovariansurvival500kb_p5e06r0.3inst_100525.csv", "PDL1prostatesurvival500kb_p5e06r0.3inst_100525.csv")
}

combined_allinst <- lapply(combined_files, read.csv)

combined_allinst[[1]]$breast_surv_inst <- c("Yes")
combined_allinst[[2]]$crc_surv_inst <- c("Yes")
combined_allinst[[3]]$lung_surv_inst <- c("Yes")
combined_allinst[[4]]$melanoma_surv_inst <- c("Yes")
combined_allinst[[5]]$ovarian_surv_inst <- c("Yes")
combined_allinst[[6]]$prostate_surv_inst <- c("Yes")

inst <- dplyr::full_join(combined_allinst[[1]], combined_allinst[[2]]) %>%
dplyr::full_join(combined_allinst[[3]]) %>%
dplyr::full_join(combined_allinst[[4]]) %>%
dplyr::full_join(combined_allinst[[5]]) %>%
dplyr::full_join(combined_allinst[[6]]) 

inst <- dplyr::mutate(inst, maf = case_when(A1FREQ>0.5 ~ 1-A1FREQ, A1FREQ<0.5 ~ A1FREQ))
inst$b_sd <- inst$z/sqrt(2*inst$maf*(1-inst$maf)*(inst$N+inst$z^2))

var <- 1
inst$r2 <- 2*inst$b_sd^2*inst$maf*(1-inst$maf)/var

inst$r2_perc <- inst$r2*100

inst$F <- (inst$r2)*(inst$N-2)/(1-inst$r2)

min(inst$F)
max(inst$F)

max(inst$r2_perc)

inst_df <- dplyr::select(inst, c(rsid, breast_surv_inst, crc_surv_inst, lung_surv_inst, melanoma_surv_inst, ovarian_surv_inst, prostate_surv_inst, r2_perc, F)) %>%
replace_na(list(breast_surv_inst = 'No', crc_surv_inst = "No", lung_surv_inst = "No", melanoma_surv_inst = "No", ovarian_surv_inst = "No", prostate_surv_inst = "No")) 

inst_df$protein <- protein
inst_df$r2_perc <- round(inst_df$r2_perc, 2)
inst_df$F	<- round(inst_df$F, 1)

inst_df <- dplyr::select(inst_df,
c(protein, rsid, r2_perc, F, breast_surv_inst, crc_surv_inst, lung_surv_inst, melanoma_surv_inst, ovarian_surv_inst, prostate_surv_inst))

inst_df <- dplyr::rename(inst_df,
"Protein instrumented" = "protein",
"Variant in breast cancer site-specific instrument set?" = "breast_surv_inst",
"Variant in colorectal cancer site-specific instrument set?" = "crc_surv_inst",
"Variant in lung cancer site-specific instrument set?" = "lung_surv_inst",
"Variant in melanoma site-specific instrument set?" = "melanoma_surv_inst",
"Variant in ovarian cancer site-specific instrument set?" = "ovarian_surv_inst",
"Variant in prostate cancer site-specific instrument set?" = "prostate_surv_inst",
"Variance explained in protein level (%)" = "r2_perc",
"F-statistic" = "F")

setwd(resultsdir)
write.table(inst_df, paste0(protein, "_maininst_varexp_F_240925.csv"), sep = ",", quote = F, row.names = F)

setwd(datadir)
setwd("ici_analyses/UKBB")

}


