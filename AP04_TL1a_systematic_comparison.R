## AP04_TL1a AP

df_control <- read.table("./tl1a/output/tl1a_tnfa_followup/tl1a_tnfa_phosphoproteomics_dea.txt", sep = "\t", na.strings = c("NULL", "NaN", "NA"), header = T)
df_compare <- read.table("./tl1a/output/tl1a_tnfa_followup/tl1a_tnfa_phosphoproteomics_dea_tl1a_vs_tnfa.txt", sep = "\t", na.strings = c("NULL", "NaN", "NA"), header = T)

library(ggplot2)
library(plotly)

## Some notes about the experiment. 50 ng/mL and 100 ng/mL TNFa, which is 1.95 nM and 3.9 nM, respectively. We want to compare 20 and 200 nM TL1a to 50 ng TNFa, or tla1_high and tl1a_med vs. tnfa_low
## which condition comparing first?
table(df_compare$contrast)
## want these ones: tl1a_high_vs_tnfa_low and tl1a_med_vs_tnfa_low 
# df_compare_cond <- df_compare[which(df_compare$contrast == "tl1a_med_vs_tnfa_low"),]
# df_compare_cond <- df_compare_cond[order(df_compare_cond$logfc),]

## now let's take some of these "top" sites and examine the magnitude change compared to control
# hit <- df_compare_cond[1,c("id")]

# df_control_hit <- df_control[(which(df_control$id==hit & df_control$contrast %in% c("tl1a_020nM", "tnfa_050_ng_ml"))),]

## hmm
## let's reproduce the correlation plots Pierre showed
# df_control_corr <- df_control[which(df_control$contrast %in% c("tl1a_020nM", "tnfa_050_ng_ml")),]
## reshape to long format for ggplot
df_control_corr <- df_control[,c("id", "id_prot", "gene_symbol", "ptm_sites", "contrast","logfc", "logfc_se", "pval", "padj", "t")]
df_control_corr <- reshape(df_control_corr, idvar = c("id"), timevar = "contrast", direction = "wide")
# df_control_corr <- df_control_corr[,c("id", "id_prot.tl1a_000_2nM", "gene_symbol.tl1a_000_2nM", "ptm_sites.tl1a_000_2nM", "logfc.tl1a_000_2nM", "logfc_se.tl1a_000_2nM", "pval.tl1a_000_2nM", "padj.tl1a_000_2nM", "t.tl1a_000_2nM", "id_prot.tl1a_020nM", "gene_symbol.tl1a_020nM", "ptm_sites.tl1a_020nM", "logfc.tl1a_020nM", "logfc_se.tl1a_020nM", "pval.tl1a_020nM", "padj.tl1a_020nM", "t.tl1a_020nM", "id_prot.tl1a_200nM", "gene_symbol.tl1a_200nM", "ptm_sites.tl1a_200nM", "logfc.tl1a_200nM", "logfc_se.tl1a_200nM", "pval.tl1a_200nM", "padj.tl1a_200nM", "t.tl1a_200nM", "id_prot.tnfa_050_ng_ml", "gene_symbol.tnfa_050_ng_ml", "ptm_sites.tnfa_050_ng_ml", "logfc.tnfa_050_ng_ml", "logfc_se.tnfa_050_ng_ml", "pval.tnfa_050_ng_ml", "padj.tnfa_050_ng_ml", "t.tnfa_050_ng_ml", "id_prot.tnfa_100_ng_ml", "gene_symbol.tnfa_100_ng_ml", "ptm_sites.tnfa_100_ng_ml", "logfc.tnfa_100_ng_ml", "logfc_se.tnfa_100_ng_ml", "pval.tnfa_100_ng_ml", "padj.tnfa_100_ng_ml", "t.tnfa_100_ng_ml")]

df_control_corr <- df_control_corr[,c("id", "id_prot.tl1a_000_2nM", "gene_symbol.tl1a_000_2nM", "ptm_sites.tl1a_000_2nM", "logfc.tl1a_000_2nM", "logfc_se.tl1a_000_2nM", "pval.tl1a_000_2nM", "padj.tl1a_000_2nM", "t.tl1a_000_2nM",  "logfc.tl1a_020nM", "logfc_se.tl1a_020nM", "pval.tl1a_020nM", "padj.tl1a_020nM", "t.tl1a_020nM",  "logfc.tl1a_200nM", "logfc_se.tl1a_200nM", "pval.tl1a_200nM", "padj.tl1a_200nM", "t.tl1a_200nM",  "logfc.tnfa_050_ng_ml", "logfc_se.tnfa_050_ng_ml", "pval.tnfa_050_ng_ml", "padj.tnfa_050_ng_ml", "t.tnfa_050_ng_ml", "logfc.tnfa_100_ng_ml", "logfc_se.tnfa_100_ng_ml", "pval.tnfa_100_ng_ml", "padj.tnfa_100_ng_ml", "t.tnfa_100_ng_ml")]
df_control_corr$ratio_med <- df_control_corr$logfc.tl1a_020nM/df_control_corr$logfc.tnfa_050_ng_ml
df_control_corr$ratio_high <- df_control_corr$logfc.tl1a_200nM/df_control_corr$logfc.tnfa_050_ng_ml

## add a significance threshold
df_control_corr$significance_med <- ifelse(df_control_corr$padj.tl1a_020nM < alpha & df_control_corr$padj.tnfa_050_ng_ml < alpha, "yes", "no")
df_control_corr$significance_high <- ifelse(df_control_corr$padj.tl1a_200nM < alpha & df_control_corr$padj.tnfa_050_ng_ml < alpha, "yes", "no")

## what do these ratios look like?
alpha <- 0.05
ggplot(data = df_control_corr, aes(x = ratio_med)) +
  geom_histogram(binwidth = 0.5) +
  xlim(-10,10)
ggplot(data = df_control_corr, aes(x = ratio_high)) +
  geom_histogram(binwidth = 0.5) +
  xlim(-10,10)
ggplot(data = df_control_corr, aes(x = ratio_med, y = -log10(padj.tl1a_020nM), color = -log10(padj.tnfa_050_ng_ml))) +
  geom_point() +
  geom_hline(yintercept = -log10(alpha)) +
  xlim(-10,10)

ggplot(data = df_control_corr, aes(x = ratio_med, y = -log10(padj.tl1a_020nM), color = logfc.tl1a_020nM)) +
  geom_point() +
  geom_hline(yintercept = -log10(alpha)) +
  xlim(-10,10) +
  scale_color_viridis_c(option="magma")

ggplot(data = df_control_corr, aes(x = logfc.tl1a_020nM, y = ratio_med, shape = significance_med, color = -log10(padj.tl1a_020nM))) +
  geom_point() +
  # geom_hline(yintercept = -log10(alpha)) +
  scale_color_viridis_c(option="magma") +
  ylim(-10,10)

## This is the plot
corr_plot_med <- ggplot(data = df_control_corr[which(df_control_corr$significance_med=="yes"),], aes(x = logfc.tl1a_020nM, y = ratio_med, shape = significance_med, color = -log10(padj.tl1a_020nM))) +
  geom_point(aes(text=sprintf("Protein: %s<br>Site: %s", gene_symbol.tl1a_000_2nM, ptm_sites.tl1a_000_2nM))) +
  # geom_hline(yintercept = -log10(alpha)) +
  scale_color_viridis_c(option="magma") +
  ylim(0,2.5) +
  geom_hline(yintercept = 1) +
  geom_hline(yintercept = 1.5, linetype = "dashed") +
  ggtitle("Ratio LogFC(TL1a med/TNFa low) vs. LogFC Tl1a med")

ggplotly(corr_plot_med, tooltip = "text")
  
  
corr_plot_high<- ggplot(data = df_control_corr[which(df_control_corr$significance_high=="yes"),], aes(x = logfc.tl1a_200nM, y = ratio_high, shape = significance_high, color = -log10(padj.tl1a_200nM))) +
  geom_point(aes(text=sprintf("Protein: %s<br>Site: %s", gene_symbol.tl1a_000_2nM, ptm_sites.tl1a_000_2nM))) +
  # geom_hline(yintercept = -log10(alpha)) +
  scale_color_viridis_c(option="magma") +
  ylim(0,2.5) +
  geom_hline(yintercept = 1) +
  geom_hline(yintercept = 1.5, linetype = "dashed") +
  ggtitle("Ratio LogFC(TL1a high/TNFa low) vs. LogFC Tl1a high")

ggplotly(corr_plot_high, tooltip = "text")


# control_corr_plot_high <- ggplot(data = df_control_corr, aes(x = logfc.tl1a_200nM, y = logfc.tnfa_050_ng_ml, text = gene_symbol.tl1a_000_2nM)) +
#   geom_point(alpha = 0.2, aes(text=sprintf("Protein: %s<br>Site: %s", gene_symbol.tl1a_000_2nM, ptm_sites.tl1a_000_2nM))) +
#   geom_abline() +
#   coord_equal(xlim=c(-3,3),ylim=c(-3,3))+
#   ggtitle("TNFa 50 ng/mL (low) vs. TL1a 200 nM (high)")
# 
# 
# ggplotly(control_corr_plot_high, tooltip = c("x", "y","text"))
# 
# 
# control_corr_plot_med <- ggplot(data = df_control_corr, aes(x = logfc.tl1a_020nM, y = logfc.tnfa_050_ng_ml, text = gene_symbol.tl1a_000_2nM)) +
#   geom_point(alpha = 0.2, aes(text=sprintf("Protein: %s<br>Site: %s", gene_symbol.tl1a_000_2nM, ptm_sites.tl1a_000_2nM))) +
#   geom_abline() +
#   coord_equal(xlim=c(-3,3),ylim=c(-3,3))+
#   ggtitle("TNFa 50 ng/mL (low) vs. TL1a 20 nM (med)")
# 
# 
# 
# 
# df_control_corr$ratio <- df_control_corr$logfc.tl1a_200nM/df_control_corr$logfc.tnfa_050_ng_ml
# 
# df_temp <- df_control_corr[,c(1,2,5,7,27,30,42,45,57,60,ncol(df_control_corr))]
# 
# library(plotly)
# ggplotly(control_corr_plot_med, tooltip = c("x", "y","text"))




