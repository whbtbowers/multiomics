setwd("/project/home17/whb17/Documents/project3/project_files/feature_selection/ex_12/")

library(ggplot2)

# Selected features for tb vs od
aucs.gp.tb_od <- read.csv("../../data/ex_12/feat_sel_2/gp_tb_od_BH_LFC_lasso_sig_factors_withaucs.csv", header=TRUE, row.names = 1)

# Selected features for tb vs ltbi
aucs.gp.tb_ltbi <- read.csv("../../data/ex_12/feat_sel_2/gp_tb_ltbi_BH_LFC_lasso_sig_factors_withaucs.csv", header=TRUE, row.names = 1)

# Selected features for tb vs non-tb
aucs.gp.tb_nontb <- read.csv("../../data/ex_12/feat_sel_2/gp_tb_nontb_BH_LFC_lasso_sig_factors_withaucs.csv", header=TRUE, row.names = 1)


ex_dir <- "ex_12/"
date <- "2018-09-01/"

############
# TB vs OD # 
############

aucs.tb_od <- as.data.frame(cbind(as.character(aucs.gp.tb_od$features), as.numeric(aucs.gp.tb_od$aucs), aucs.gp.tb_od$lower.ci, aucs.gp.tb_od$upper.ci))
colnames(aucs.tb_od) <- c("features", "aucs", "lower.ci", "upper.ci")
aucs.tb_od$aucs <-as.numeric(as.character(aucs.tb_od$aucs))
aucs.tb_od$lower.ci <-as.numeric(as.character(aucs.tb_od$lower.ci))
aucs.tb_od$upper.ci <-as.numeric(as.character(aucs.tb_od$upper.ci))

bars.tb_od <- ggplot(data=aucs.tb_od, aes(x=factor(features, levels= as.character(aucs.gp.tb_od$features)), y=aucs)) +
 geom_bar(stat='identity', fill=c("#048AFF", "#048AFF", "#048AFF", "#048AFF", "#048AFF", "#048AFF", "#048AFF", "#048AFF", "#FFDE04", "#048AFF", "#048AFF", "#048AFF", "#048AFF", "#048AFF", "#048AFF", "#048AFF", "#048AFF", "#048AFF", "#048AFF", "#048AFF")) +
  geom_errorbar(aes(ymin=lower.ci, ymax=upper.ci), width=0.3) +
  scale_x_discrete(name="Feature",
                   labels=c("CYB561", "Hs.131087", "MYL12A", "EBF1", "CD74", "DEFA1", "GBP6", "GBP5", "TTR", "CD93", "CALML4", "KLHL28", "HLA-DPB1", "HIST1H4C", "CCNJL", "RFX2", "MYL9", "MPZL1", "ALDH1A1", "AP2A1")) +
  scale_y_continuous(name="AUC",
                     limits=c(0, 1)) +
  theme_minimal() +
  theme(axis.text.x  = element_text(angle=90, size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20)
        )

ggsave(paste("../../img/", ex_dir, date, "tb_od_aucsbar.png", sep=""), 
       bars.tb_od,
       width=200,
       height = 150,
       units="mm",
       dpi=320
       )

##############
# TB vs LTBI # 
##############
  
aucs.tb_ltbi <- as.data.frame(cbind(as.character(aucs.gp.tb_ltbi$features), as.numeric(aucs.gp.tb_ltbi$aucs), aucs.gp.tb_ltbi$lower.ci, aucs.gp.tb_ltbi$upper.ci))
colnames(aucs.tb_ltbi) <- c("features", "aucs", "lower.ci", "upper.ci")

aucs.tb_ltbi$aucs <-as.numeric(as.character(aucs.tb_ltbi$aucs))
aucs.tb_ltbi$lower.ci <-as.numeric(as.character(aucs.tb_ltbi$lower.ci))
aucs.tb_ltbi$upper.ci <-as.numeric(as.character(aucs.tb_ltbi$upper.ci))

bars.tb_ltbi <- ggplot(data=aucs.tb_ltbi, aes(x=factor(features, levels=as.character(aucs.gp.tb_ltbi$features)), y=aucs)) +
  geom_bar(stat='identity', fill=c("#048AFF", "#048AFF", "#048AFF", "#048AFF", "#048AFF", "#048AFF", "#FFDE04", "#048AFF", "#048AFF", "#048AFF", "#048AFF", "#048AFF", "#FFDE04")) +
  geom_errorbar(aes(ymin=lower.ci, ymax=upper.ci), width=0.3) +
  scale_x_discrete(name="Feature",
                   labels=c("FCGR1A", "NDRG2", "SPIB", "ZNF296", "GBP5", "PITPNC1", "SAA", "ALML4", "CD79A", "CCR6", "OSBPL10", "NRG1", "TTR")) +
  scale_y_continuous(name="AUC",
                     limits=c(0, 1)) +
  theme_minimal() +
  theme(axis.text.x  = element_text(angle=90, size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20)
  )

ggsave(paste("../../img/", ex_dir, date, "tb_ltbi_aucsbar.png", sep=""), 
       bars.tb_ltbi,
       width=130,
       height = 150,
       units="mm",
       dpi=320
)

################
# TB vs non-TB # 
################

aucs.tb_nontb <- as.data.frame(cbind(as.character(aucs.gp.tb_nontb$features), as.numeric(aucs.gp.tb_nontb$aucs), aucs.gp.tb_nontb$lower.ci, aucs.gp.tb_nontb$upper.ci))
colnames(aucs.tb_nontb) <- c("features", "aucs", "lower.ci", "upper.ci")

aucs.tb_nontb$aucs <-as.numeric(as.character(aucs.tb_nontb$aucs))
aucs.tb_nontb$lower.ci <-as.numeric(as.character(aucs.tb_nontb$lower.ci))
aucs.tb_nontb$upper.ci <-as.numeric(as.character(aucs.tb_nontb$upper.ci))

bars.tb_nontb <- ggplot(data=aucs.tb_nontb, aes(x=factor(features, levels=as.character(aucs.gp.tb_nontb$features)
), y=aucs)) +
  geom_bar(stat='identity', fill=c("#048AFF", "#048AFF", "#048AFF", "#048AFF", "#048AFF", "#048AFF", "#FFDE04", "#048AFF", "#048AFF", "#048AFF", "#048AFF", "#048AFF", "#048AFF", "#048AFF", "#048AFF", "#048AFF", "#048AFF", "#048AFF", "#048AFF", "#048AFF", "#048AFF", "#048AFF", "#048AFF", "#048AFF", "#048AFF", "#048AFF", "#048AFF", "#FFDE04", "#FFDE04", "#048AFF")) +
  geom_errorbar(aes(ymin=lower.ci, ymax=upper.ci), width=0.3) +
  scale_x_discrete(name="Feature",
                   labels=c("CALML4", "GBP6", "C1QA", "LOC90925", "IL15", "ABHD12B", "TTR", "EBF1", "FZD2", "GBP5", "OSBPL10", "NRG1", "ALDH1A1", "COL9A2", "ELANE", "DEFA1", "DEFA1B", "CD151", "BG205162", "EOMES", "Hs.131087", "NT5C3", "GYPE", "IL1RN", "RPL10A", "LIMK2", "RN5S9", "IL-1RA", "CFH", "STK3")) +
  scale_y_continuous(name="AUC",
                     limits=c(0, 1)) +
  theme_minimal() +
  theme(axis.text.x  = element_text(angle=90, size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20)
  )

ggsave(paste("../../img/", ex_dir, date, "tb_nontb_aucsbar.png", sep=""), 
       bars.tb_nontb,
       width=300,
       height = 150,
       units="mm",
       dpi=320
)



