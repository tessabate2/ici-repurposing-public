#module load languages/R/4.1.2 

# load packages
library(survSNP)
library(tidyverse)

library(here)
readRenviron(here("config.env"))

resultsdir <- Sys.getenv("resultsdir")
setwd(resultsdir)
setwd("power_estimates")

### POWER CALCULATIONS FOR DISCOVERY+REPLICATION (COMBINED) pQTL
### (have excluded GE from power calculations as not including in main manuscript)

### calculate power to detect GHR
# n  = combined number of participants (across consortium GWAS)
# event_rate = death proportion (across consortium GWAS)
# fiveyr_surv = 5-year survival (as proportion of 1) from CRUK
# allele_freq = relative allele frequency for lead cis SNP (UKB pQTL GWAS combined)
# SNP_beta = beta for lead cis SNP (UKB pQTL GWAS combined) (*-1 if is negative)
power_calc <- function(protein, survival, n, event_rate, fiveyr_surv, allele_freq, SNP_beta) {
  power_PD1 <- survSNP.power.table(GHRs = c(seq(1.00, 1.3, 0.00005)), ns = n, 
                                   rafs = allele_freq, erates = event_rate, model = "additive", 
                                   test = "additive", alpha = 0.05, pilm = fiveyr_surv, 
                                   lm = 5) %>%
    ### loge transform GHRs
    mutate(log_HR = log(GHR)) %>%
    ### divide log transformed GHRs by SNP beta to scale per SD change
    mutate(scaled_log_HR = log_HR/SNP_beta) %>%
    ### take exponent scaled HRs to find HR per SD change
    mutate(scaled_HR = exp(scaled_log_HR)) %>%
    ### round calculated power values to 3dp
    mutate(pow0_3dp = round(pow0, 3))
  ### find scaled HR for calculated power = 0.800 (to 3dp)
  power_PD1_p80 <- power_PD1[power_PD1$pow0_3dp == 0.800, ]

if (dim(power_PD1_p80)[1] > 1) {
power_PD1_p80 <- power_PD1_p80[which.max(power_PD1_p80$GHR), ]

} else {
power_PD1_p80 <- power_PD1_p80
}

power_tbl <- data.frame(protein = protein, survival = survival, n=n, event_rate = round(event_rate, 3), fiveyr_survival = fiveyr_surv, power = 0.80, lower_HR = round(1/(power_PD1_p80$scaled_HR), 2), 
upper_HR = round(power_PD1_p80$scaled_HR, 2))

power_tbl <- dplyr::rename(power_tbl, "Protein" = "protein",
"Survival outcome" = "survival",
"Survival GWAS N participants" = "n",
"Survival GWAS mortality event rate" = "event_rate",
"Five year survival" = "fiveyr_survival",
"Estimated lower HR detectable at 80% power" = "lower_HR",
"Estimated upper HR detectable at 80% power" = "upper_HR") %>%
dplyr::select(-power)

return(power_tbl)
}

## PD-1
breast_PD1 <- power_calc(protein = "PD1", survival = "breast", n = 91686, event_rate = (7531/91686), fiveyr_surv = 0.859, allele_freq = 0.2119, SNP_beta = -0.285*-1)
crc_PD1 <- power_calc(protein = "PD1", survival = "colorectal", n = 16964, event_rate = (4010/16964), fiveyr_surv = 0.584, allele_freq = 0.2119, SNP_beta = -0.285*-1)
lung_PD1 <- power_calc(protein = "PD1", survival = "lung", n = 7352, event_rate = (4598/7352), fiveyr_surv = 0.21, allele_freq = 0.2119, SNP_beta = -0.285*-1)
mel_PD1 <- power_calc(protein = "PD1", survival = "melanoma", n = 10982, event_rate = (1041/10982), fiveyr_surv = 0.926, allele_freq = 0.2119, SNP_beta = -0.285*-1)
ov_PD1 <- power_calc(protein = "PD1", survival = "ovarian", n = 2901, event_rate = (1656/2901), fiveyr_surv = 0.45, allele_freq = 0.2119, SNP_beta = -0.285*-1)
pr_PD1 <- power_calc(protein = "PD1", survival = "prostate", n = 67758, event_rate = (7914/67758), fiveyr_surv = 0.885, allele_freq = 0.2119, SNP_beta = -0.285*-1)

PD1_powerest <- rbind(breast_PD1, crc_PD1, lung_PD1, mel_PD1, ov_PD1, pr_PD1)
PD1_powerest <- PD1_powerest[!duplicated(PD1_powerest), ]

## PD-L1
breast_PDL1 <- power_calc(protein = "PDL1", survival = "breast", n = 91686, event_rate = (7531/91686), fiveyr_surv = 0.859, allele_freq = 0.7486, SNP_beta = 0.340)
crc_PDL1 <- power_calc(protein = "PDL1", survival = "colorectal", n = 16964, event_rate = (4010/16964), fiveyr_surv = 0.584, allele_freq = 0.7486, SNP_beta = 0.340)
lung_PDL1 <- power_calc(protein = "PDL1", survival = "lung", n = 7352, event_rate = (4598/7352), fiveyr_surv = 0.21, allele_freq = 0.7486, SNP_beta = 0.340)
mel_PDL1 <- power_calc(protein = "PDL1", survival = "melanoma", n = 10982, event_rate = (1041/10982), fiveyr_surv = 0.926, allele_freq = 0.7486, SNP_beta = 0.340)
ov_PDL1 <- power_calc(protein = "PDL1", survival = "ovarian", n = 2901, event_rate = (1656/2901), fiveyr_surv = 0.45, allele_freq = 0.7486, SNP_beta = 0.340)
pr_PDL1 <- power_calc(protein = "PDL1", survival = "prostate", n = 67758, event_rate = (7914/67758), fiveyr_surv = 0.885, allele_freq = 0.7486, SNP_beta = 0.340)

PDL1_powerest <- rbind(breast_PDL1, crc_PDL1, lung_PDL1, mel_PDL1, ov_PDL1, pr_PDL1)
PDL1_powerest <- PDL1_powerest[!duplicated(PDL1_powerest), ]

PD1_PDL1_powerest_main <- dplyr::full_join(PD1_powerest, PDL1_powerest)

write.table(PD1_PDL1_powerest_main, "main_power_estimates_pd1_pdl1_240925.csv", sep = ",", quote = F, row.names = F)

#### SENSITIVITY method
## Using powerSurvEpi package
#install.packages("powerSurvEpi")
library(powerSurvEpi)

survival <- c("breast", "colorectal", "lung", "melanoma", "ovarian", "prostate")
n <- c(91686, 16964, 7352, 10982, 2901, 67758)
prop_died <- c((7531/91686), (4010/16964), (4598/7352), (1041/10982), (1656/2901), (7914/67758))
pd1_r2 <- 2.48/100
pdl1_r2 <- 4.47/100
# r2<-c(pd1_r2,pdl1_r2)
n2_pd1 <- n*pd1_r2
n2_pdl1 <- n*pdl1_r2
hr_pd1 <- dplyr::filter(PD1_PDL1_powerest_main, Protein == "PD1")$'Estimated upper HR detectable at 80% power'
hr_pdl1 <- dplyr::filter(PD1_PDL1_powerest_main, Protein	== "PDL1")$'Estimated upper HR detectable at 80% power'
sd <- 1
conf <- 0
p <- 0.05

pd1_powest <- powerEpiCont.default(n = n2_pd1, theta = hr_pd1, psi = prop_died, rho2 = conf, alpha = p, sigma2 = sd)*100
pdl1_powest <- powerEpiCont.default(n = n2_pdl1, theta = hr_pdl1, psi = prop_died, rho2 = conf, alpha = p, sigma2 = sd)*100

pd1_powest_df <- data.frame(protein = "PD1", survival = survival, n = n, event_rate = round(prop_died,3), n_r2 = round(n2_pd1, 0), hr = hr_pd1, pow_to_detect_hr = round(pd1_powest, 1))
pdl1_powest_df <- data.frame(protein = "PDL1", survival = survival, n = n,	event_rate = round(prop_died,3),	n_r2 = round(n2_pdl1, 0), hr = hr_pdl1, pow_to_detect_hr = round(pdl1_powest, 1))

PD1_PDL1_powerest_sensitivity <- dplyr::full_join(pd1_powest_df, pdl1_powest_df) %>%
dplyr::rename("Protein" = "protein", "Survival outcome" = "survival", "Survival GWAS N participants" = "n", "Survival GWAS mortality event rate" = "event_rate",
"Scaled N participants" = "n_r2", "Estimated HR detectable at 80% power" = "hr", "Power to detect estimated HR (%)" = "pow_to_detect_hr")

write.table(PD1_PDL1_powerest_sensitivity, "sensitivity_power_estimates_pd1_pdl1_240925.csv",	sep = ",", quote = F, row.names	= F)
