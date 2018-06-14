# Workspace ---------------------------------------------------------------
library(data.table)
library(dplyr)
library(tidyr)
library(igraph)
library(visNetwork)
library(ggplot2)
library(ggrepel)
library(DT)

set.seed(12)

# Data prep ---------------------------------------------------------------
load("Data/20171111.RData")
setDT(hitsLD)

hitsType <- merge(hitsType, MAP[ , list(rsid, BP)], by.x = "oncoID", by.y = "rsid")

# update hits' SNP names
metaRegionClean <- merge(metaRegionClean, 
                         unique(hitsType[ , list(SNPid, temp_SNP = SNP)]),
                         by = "SNPid", all.x = TRUE)
metaRegionClean$SNP <- ifelse(!is.na(metaRegionClean$temp_SNP),
                              metaRegionClean$temp_SNP,
                              metaRegionClean$SNP)
metaRegionClean$temp_SNP <- NULL

# for network plot, "hover mouse text"
metaRegionClean$title <- 
  paste0("<p><b>rsid: ", metaRegionClean$rsid, "</b><br></p>",
         "<p><b>SNP: ", metaRegionClean$SNP, "</b><br>",
         "<b>BP: ", metaRegionClean$BP, "</b><br>",
         "<b>A1: ", toupper(metaRegionClean$Allele1), "</b><br>",
         "<b>A2: ", toupper(metaRegionClean$Allele2), "</b><br>",
         "<b>Type: ", if_else(metaRegionClean$TYPED == 2, "typed", "imputed"), "</b><br>",
         "<b>MAF: ", metaRegionClean$Freq1,"</b><br>",
         "<b>Info: ", metaRegionClean$info, "</b><br>",
         "<b>Effect: ", metaRegionClean$Effect, "</b><br>",
         "<b>P: ", metaRegionClean$`P-value`, "</b><br>",
         "<b>MaxPostProb: ", metaRegionClean$MaxPostProb, "</b><br>",
         "<b>MaxBF: ", metaRegionClean$MaxBF, "</b><br>",
         "</p>")

# Clean MAP ---------------------------------------------------------------
MAP$SNP <- sapply(strsplit(MAP$rsid, ":"), function(i){
  if_else(substr(i[1], 1, 2) == "rs", i[1], paste(i, collapse = ":"))
  }) 
  
MAP <- MAP %>% 
  group_by(SNP) %>% 
  mutate(SNPcnt = n()) %>% 
  ungroup() %>% 
  mutate(SNP = if_else(SNPcnt > 1, rsid, SNP),
         SNPid = paste0("SNP_", rn)) %>% 
  as.data.table()

# LD clean with rs numbesr ------------------------------------------------
x <- merge(hitsLD, unique(hitsType[ , list(SNPid, SNP_A = SNP)]),
           by.x = "SNP_hit", by.y = "SNPid", all.x = TRUE)
xx <- merge(x, unique(hitsType[ , list(SNPid, SNP_B = SNP)]),
            by.x = "SNP", by.y = "SNPid", all.x = TRUE)
xxx <- merge(xx, MAP[, list(SNPid, MAPrsid_A = SNP, BP_A = BP)],
             by.x = "SNP_hit", by.y = "SNPid", all.x = TRUE)
xxxx <- merge(xxx, MAP[, list(SNPid, MAPrsid_B = SNP, BP_B = BP)],
              by.x = "SNP", by.y = "SNPid", all.x = TRUE)
#head(xxxx)
hitsLDclean <-
  xxxx %>% 
  transmute(SNP_A = if_else(is.na(SNP_A), MAPrsid_A, SNP_A),
            SNP_B = if_else(is.na(SNP_B), MAPrsid_B, SNP_B),
            SNPid_A = SNP_hit,
            SNPid_B = SNP,
            BP_A, BP_B,
            R2)
rm(x, xx, xxx, xxxx)

# Custom functions --------------------------------------------------------
source("UDF.R")

