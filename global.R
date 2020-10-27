# Workspace ---------------------------------------------------------------
library(data.table)
library(dplyr)
library(tidyr)
library(igraph)
library(visNetwork)
library(ggplot2)
library(ggrepel)
library(DT)
#https://github.com/oncogenetics/oncofunco
library(oncofunco)

set.seed(12)

# Data prep ---------------------------------------------------------------
load("Data/20171111.RData")
setDT(hitsLD)

hitsType <- merge(hitsType, MAP[ , list(rsid, BP)], by.x = "oncoID", by.y = "rsid")

# update hits' SNP names
metaRegionClean <- merge(metaRegionClean, 
                         unique(hitsType[ , .(SNPid, temp_SNP = SNP)]),
                         by = "SNPid", all.x = TRUE)
metaRegionClean[, SNP := ifelse(!is.na(metaRegionClean$temp_SNP),
                              metaRegionClean$temp_SNP,
                              metaRegionClean$SNP) ]
metaRegionClean[, temp_SNP := NULL ]

# for network plot, "hover mouse text"
metaRegionClean[, title := paste0(
  "<p><b>rsid: ", rsid, "</b><br></p>",
  "<p><b>SNP: ", SNP, "</b><br>",
  "<b>BP: ", BP, "</b><br>",
  "<b>A1: ", toupper(Allele1), "</b><br>",
  "<b>A2: ", toupper(Allele2), "</b><br>",
  "<b>Type: ", if_else(TYPED == 2, "typed", "imputed"), "</b><br>",
  "<b>MAF: ", Freq1,"</b><br>",
  "<b>Info: ", info, "</b><br>",
  "<b>Effect: ", Effect, "</b><br>",
  "<b>P: ", `P-value`, "</b><br>",
  "<b>MaxPostProb: ", MaxPostProb, "</b><br>",
  "<b>MaxBF: ", MaxBF, "</b><br>",
  "</p>") ]

# Clean MAP ---------------------------------------------------------------
MAP[, SNP := extractRsid(rsid) ]
MAP[, SNPcnt := .N, by = SNP ]
MAP[ SNPcnt > 1, SNP := rsid ]
MAP[, SNPid := paste0("SNP_", .I) ]

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
hitsLDclean <- xxxx[, .(SNP_A = ifelse(is.na(SNP_A), MAPrsid_A, SNP_A),
                        SNP_B = ifelse(is.na(SNP_B), MAPrsid_B, SNP_B),
                        SNPid_A = SNP_hit,
                        SNPid_B = SNP,
                        BP_A, BP_B,
                        R2) ]
rm(x, xx, xxx, xxxx)

