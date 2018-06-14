# 14/06/2018

# About -------------------------------------------------------------------
# Make LocusExplorer files:
# - LD
# CHR_A,BP_A,SNP_A,CHR_B,BP_B,SNP_B,R2
# 1,150283370,rs587636640,1,150162158,rs113833990,0.062
# 1,150283370,rs587636640,1,150163915,rs35987241,0.16
# ...
# - assoc
# CHR,SNP,BP,P,TYPED,EFFECT
# chr1,rs72694965,150158723,0.0001674,1,0.0416
# chr1,rs3118121,150160430,0.193,1,0.0122
# ...
# - stats
# SNP,rsid,BP,PP_best_tag,PP_tag_r2,PostProb,BF,JAM99
# rs587636640,1:150283370:C:AA,150283370,1:150283370:C:AA,1,0.0129,2.601395291,0
# rs570851450,1:150334248:AT:A,150334248,1:150283370:C:AA,0.916126571,0.0129,2.601395291,0
# ...


# Workspace ---------------------------------------------------------------
library(dplyr)
library(data.table)

setwd("C:/Users/tdadaev/Desktop/Work/GitHubProjects/8q24")
load("Data/20171111.RData")


# LD ----------------------------------------------------------------------
setDT(hitsLD)

# get hit rsids
metaRegionClean$SNPclean <-
  setNames(hitsType[ hitType == "final12", SNP],
           hitsType[ hitType == "final12", SNPid])[ metaRegionClean$SNPid ]
metaRegionClean$SNPclean <- ifelse(is.na(metaRegionClean$SNPclean),
                                   metaRegionClean$SNP, metaRegionClean$SNPclean)

outLD <- hitsLD[ SNP_hit %in% hitsType[ hitType == "final12", SNPid], ]
#add rsid and bp
outLD <- merge(outLD, metaRegionClean[, list(SNPid, SNP_B = SNPclean, BP_B = BP)],
               by.x = "SNP", by.y = "SNPid")
outLD <- merge(outLD, metaRegionClean[, list(SNPid , SNP_A = SNPclean, BP_A = BP)],
               by.x = "SNP_hit", by.y = "SNPid")

outLD <- outLD %>% 
  transmute(
    CHR_A = 8, BP_A, SNP_A,
    CHR_B = 8, BP_B, SNP_B,
    R2) %>% 
  arrange(BP_A, BP_B)

# Assoc -------------------------------------------------------------------
outAssoc <- 
  metaRegionClean %>% 
  transmute(CHR = "chr8",
            SNP = SNPclean,
            BP,
            P = `P-value`,
            TYPED,
            EFFECT = Effect) %>% 
  arrange(BP)

# Stats -------------------------------------------------------------------
outStats <- 
  merge(metaRegionClean, PP[ , list(rsid = name,
                                    PP_best_tag = best_tag,
                                    PP_tag_r2 = `r^2`)], by = "rsid", all.x = TRUE) %>% 
  filter(!(MaxPostProb == 0 | is.na(PP_tag_r2))) %>% 
  transmute(SNP = SNPclean,
            rsid, BP,
            PP_best_tag, PP_tag_r2,
            PostProb = MaxPostProb, BF = MaxBF,
            JAM99 = as.integer(SNPclean %in% hitsType[ hitType == "final12", SNP])) %>% 
  arrange(BP)



# Output ------------------------------------------------------------------
write.csv(outAssoc, "Data/LE/chr8_127333841_129040776_assoc.txt", row.names = FALSE, quote = FALSE)
write.csv(outLD, "Data/LE/chr8_127333841_129040776_LD.txt", row.names = FALSE, quote = FALSE)
write.csv(outStats, "Data/LE/chr8_127333841_129040776_stats.txt", row.names = FALSE, quote = FALSE)

