library(tidyverse)
library(ggpubr)
library(data.table)
library(readr)
library(readxl)
library(formattable)
library(OlinkAnalyze)
library(gtsummary)
library(flextable)
library(ggplotify)
library(pheatmap)
library(plotly)
library(gapminder)
library(manipulateWidget)
library(ggfortify)
library(reactable)
library(mice)
library(naniar)
library(ggrepel)
library(writexl)
library(qvalue)

# setwd("/home/meizhen/proteomics/data")


setwd("~/Documents/projects/proteomics_ht/data")

## import data (clean data - remove controls in sample and assay)
data <- fread("intermediate/whole_inverserank_data.txt")

## Proteins with <LOD% > 75%
OlinkID_lod_remove75<- fread("intermediate/protein_lod_remove75.txt")$OlinkID


chemical_name <- names(data %>% 
                         dplyr::select(starts_with("pfas_"), starts_with("pth_") , starts_with("pah_") , starts_with("parab_") , starts_with("bis_") , starts_with("metal_") , -ends_with("_blod") ) )


data_ready<-  data %>% 
              mutate(ptb_type = case_when(ptb == 0 ~ "Not preterm birth",
                                          ptb == 1 & ptb_spont == 1 ~ "Spontaneous preterm birth",
                                          ptb == 1 & ptb_placental == 1 ~ "Placental preterm birth",
                                          ptb == 1 & ptb_protocol == 1 ~ "Protocol preterm birth"),
                     ptb_cat = factor(ptb, levels = c(0,1), labels = c("No", "Yes")),
                     hispanic = factor(hispanic, levels = c(0,1), labels = c("No", "Yes")),
                     race_recat = factor(ifelse(race == "Caucasian", "Caucasian", "Non-Caucasian"), levels = c("Caucasian", "Non-Caucasian")),
                     parity_recat = factor(case_when(parity == 0 ~ "0",
                                                     parity == 1 ~ "1",
                                                     .default = ">1"), levels = c("0", "1", ">1"))) %>% 
                    mutate(across(chemical_name, ~(log2(.x)), .names = '{col}_l'),
                           weight = ifelse(ptb == 1, 1/(43/143), 1/(43/1038)))


################################################
## function for robust linear regression model
################################################
lm_fit_info<- function(outcome, expo_names, covariates, weight_name, data, family){
  
  lm <- data.frame(Outcome = NA_character_, Exposure = NA_character_, Value = NA_real_, Std.Error = NA_real_, Lower.limit = NA_real_, Upper.limit = NA_real_, p.value = NA_real_)
  
  for(i in 1:length(outcome)){
    
    for (j in 1:length(expo_names)) {
      
      s_lm <- (glm(as.formula(paste0(outcome[i], "~", expo_names[j],  paste(covariates, collapse = " +"))), 
                   weights=weight_name, 
                   data = data,
                   family = family))
      
      # robust vcov and inference
      Vrob <- sandwich::vcovHC(s_lm, type = "HC3")
      ct   <- lmtest::coeftest(s_lm, vcov. = Vrob)
      ci   <- lmtest::coefci(s_lm, vcov. = Vrob)
      
      r.est <- cbind(
        Estimate     = coef(s_lm),
        Std.Error  = sqrt(diag(Vrob)),
        lower_limit  = ci[, 1],
        upper_limit  = ci[, 2],
        p.value  = ct[, 4]
      )
      
      lm <- rbind(lm, c(outcome[i], expo_names[j], as.numeric(r.est[2, c(1, 2, 3, 4, 5)])))
    }
    
  }
  
  lm$Value <- as.numeric(lm$Value)
  lm$Std.Error <- as.numeric(lm$Std.Error)
  lm$Lower.limit <- as.numeric(lm$Lower.limit)
  lm$Upper.limit <- as.numeric(lm$Upper.limit)
  lm$p.value <- as.numeric(lm$p.value)
  
  lm = lm[-1,]  
  
}


protein<- setdiff(names(data_ready %>%
                          dplyr::select(starts_with("OID"))), OlinkID_lod_remove75)

# protein<- names(data_ready[, 3:4])

covariates_name<- c("+", "age", "bmi_prepreg", "race_recat")


## pth - without weight, GM
exposure<- names(data_ready %>% 
                  dplyr::select(starts_with("pth") & ends_with("_SG_gm_l")))

model<- lm_fit_info(protein, exposure, covariates_name, weight_name = NULL, data_ready, "gaussian") %>% 
           mutate(fdr_bh = p.adjust(p.value, method = "BH"),
                  fdr_by = p.adjust(p.value, method = "BY"),
                  qvalue = qvalue(p.value)$qvalues) %>% 
            group_by(Exposure) %>% 
            mutate(fdr_bh_single = p.adjust(p.value, method = "BH"),
                   fdr_by_single = p.adjust(p.value, method = "BY"),
                   qvalue_single = qvalue(p.value)$qvalues) %>% 
            ungroup()


write.csv(model, 'model/paper1/pth_gm_ad_noweight_robust.csv')




