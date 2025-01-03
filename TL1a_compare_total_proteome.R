### TL1a total proteome follow up. Can we capture any of the sites?

## Let's merge on previous dataset and see. 
library(MSstats)
library(MSstatsPTM)
library(ggplot2)
library(plotly)
library(ggrepel)
library(stringr)

## First attempt from scratch, ignore, skip to MSstats below
# df1_raw <- read.csv("tl1a_tnfa_phosphoproteomics_dea.txt", sep = "\t", na.strings = c(""))
# 
# df2_raw <- read.csv("24_12_09_total_proteome_from_Steve.csv", na.strings = c(""))
# ## pull out only phosphosites
# df2 <- df2_raw[grepl("phospho", df2_raw$Modifications.in.Master.Proteins, ignore.case = TRUE),]
# nrow(df2)/nrow(df2_raw)
# ## 2.4% peptides are phospho
# df2_raw <- NULL
# 
# ## how many proteins are shared?
# # sum(!is.na(match(unique(df1_raw$id_prot), unique(df2$Master.Protein.Accessions))))
# # ## 887 proteins are shared
# # vec_overlap_coord <- match(unique(df1_raw$id_prot), unique(df2$Master.Protein.Accessions))
# # ## matches in df2 below
# # vec_overlap_coord <- vec_overlap_coord[which(!is.na(vec_overlap_coord))]
# # vec_overlap <- unique(df2$Master.Protein.Accessions)[vec_overlap_coord]
# ## sanity check 
# # sum(unique(df1_raw$id_prot) %in% vec_overlap)
# 
# # df1_comp <- df1[which(df1$id_prot %in% vec_overlap),]
# df1_comp <- df1_raw
# df1_comp$ptm_sites <- gsub("\\|.*$", "", df1_comp$ptm_sites)
# df1_comp$ptm <- df1_comp$ptm_sites
# df1_comp$ptm <- gsub("^.*_p", "", df1_comp$ptm)
# df1_comp$ptm <- as.numeric(df1_comp$ptm)
# df1_comp <- unique(df1_comp[,c("id_prot", "ptm")])
# nrow(df1_comp)
# ## 17,921
# nrow(df1_comp[which(!is.na(df1_comp$ptm)),])
# ## 15,330 have quant site
# nrow(df1_comp[which(!is.na(df1_comp$ptm)),])/nrow(df1_comp)
# ## 85.5% of df1
# 
# # df2_comp <- df2[which(df2$Master.Protein.Accessions %in% vec_overlap),]
# df2_comp <- df2
# df2_comp$ptm <- gsub("^.*\\[S|^.*\\[T|^.*\\[Y", "", df2_comp$Modifications.in.Master.Proteins)
# df2_comp$ptm <- gsub("\\]$", "", df2_comp$ptm)
# df2_comp$ptm <- as.numeric(df2_comp$ptm)
# df2_comp <- unique(df2_comp[,c("Master.Protein.Accessions", "ptm")])
# nrow(df2_comp)
# ## 1735
# nrow(df2_comp[which(!is.na(df2_comp$ptm)),])
# ## 1462 have quant site
# nrow(df2_comp[which(!is.na(df2_comp$ptm)),])/nrow(df2_comp)
# ## 84.3% of df2
# 
# ## Now that we have unique protein ID + ptm site, let's use that as unique matching. 
# # df1_comp <- unique(df1_comp[,c("id_prot", "ptm")])
# names(df1_comp)[which(names(df1_comp)=="id_prot")] <- "Master.Protein.Accessions"
# df1_comp$Experiment_phospho <- "Yes"
# # df2_comp <- unique(df2_comp[,c("Master.Protein.Accessions", "ptm")])
# df2_comp$Experiment_totalproteome <- "Yes"
# 
# df_join <- merge(df1_comp[which(!is.na(df1_comp$ptm)),], df2_comp[which(!is.na(df2_comp$ptm)),], all.x = T, all.y = T)
# df_join <- df_join[which(!is.na(df_join$ptm)),]
# sum(!is.na(df_join$Experiment_phospho)& !is.na(df_join$Experiment_totalproteome))
# 
# df_both <- df_join[which(!is.na(df_join$Experiment_phospho) & !is.na(df_join$Experiment_totalproteome)),]
# df_both$Entry <- gsub("-.*$", "", df_both$Master.Protein.Accessions)
# 
# df_info <- read.csv("uniprotkb_Human_AND_model_organism_9606_2024_10_21.tsv", sep = "\t")
# 
# df_both <- merge(df_both, df_info, by = "Entry", all.x = T)
# ## matches in df2 below
# vec_overlap_coord <- vec_overlap_coord[which(!is.na(vec_overlap_coord))]
# vec_overlap <- unique(df2$Master.Protein.Accessions)[vec_overlap_coord]
# ## sanity check 
# # sum(unique(df1_raw$id_prot) %in% vec_overlap)


#### Trying with MSstatsPTM ####
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("MSstatsPTM")
library(MSstatsPTM)
df_raw_from_PDP <- read.csv("NG_09062024_TL1A_TNFa_PSM_.csv", header = T, na.strings = c("", "NA", "NaN", "NULL", "Null", "null"))
## Need one more column:
# Error in MSstatsPTMSiteLocator(input, protein_name_col = which_proteinid,  : 
# FASTA file not provided and `Start` column missing from data. 
# MSstatsPTMSiteLocator requires one of these to be provided to 
# identify modification site number.

# df_annot <- read.csv("PD_annotation.csv", header = T)
df_annot <- read.csv("PD_annot_Steve.csv", header = T)
df_annot <- do.call("rbind", replicate(length(unique(df_raw_from_PDP$Spectrum.File)), df_annot, simplify = FALSE))
df_annot$TechRepMixture <- 1
df_annot$Fraction <- rep(c(1:length(unique(df_raw_from_PDP$Spectrum.File))), each = length(unique(df_annot$Channel)))
vec_spectrum_files <- unique(df_raw_from_PDP$Spectrum.File)
vec_spectrum_files <- vec_spectrum_files[order(nchar(vec_spectrum_files), vec_spectrum_files)]
df_annot$Run <- rep(vec_spectrum_files, each = length(unique(df_annot$Channel)))
df_annot$Condition <- gsub("-.*$", "", df_annot$Condition)


  ## Run the MSstats converter. "Master.Protein.Accessions" is used for "ProteinName". "Sequence" and "Modifications" are used for "PeptideSequence". "Charge"
  list_quant <- MSstatsPTM::PDtoMSstatsPTMFormat(input = df_raw_from_PDP, 
                                                 annotation = df_annot,
                                                 fasta_path = "/Users/PANOVA02/Library/CloudStorage/OneDrive-Pfizer/projects/AP04_TL1a_phospho/NG_09062024_TL1A_TNFa_Proteome_24FX_[24].fasta",
                                                 labeling_type = "TMT",
                                                 which_proteinid = "Master.Protein.Accessions",
                                                 ## phospho is default, can omit if you want
                                                 mod_id = "\\(Phospho\\)")
  ## save memory if needed
  # df_raw_from_PDP <- NULL
  
  ## dataProcess step with the output from the converter
  list_PD_proposed <- MSstatsPTM::dataSummarizationPTM_TMT(
    data = list_quant,
    method = "msstats",
    ## only 2k sites, global norm of PTMs might derail, but not likely ALL 2k or even majority of 2k sites depend on TL1a or TNFa
    global_norm = TRUE,
    global_norm.PTM = TRUE,
    reference_norm = FALSE,
    reference_norm.PTM = FALSE,
    MBimpute = FALSE,
    MBimpute.PTM = FALSE,
    use_log_file = FALSE
  )
  ## save memory, compost
  # list_quant <- NULL
  ## need to hold onto the entire list_spectronaut_proposed for the groupComparison step later
  
  save.image("./24_12_10_AP04_steve_total_proteome.RData")
 
  dataProcessPlotsPTM(list_PD_proposed,
                      type = 'QCPLOT',
                      which.PTM = "allonly",
                      address = FALSE)

  #### Assigning Control and Experimental Groups ####
  ## This code is the same for any experiment
  ## make the controls, aka denominator
  ## PDP SOURCE
  levels(list_PD_proposed$PTM$ProteinLevelData$Condition)
  vec_control_conditions <- c("control")
  ## create a list to collapse at end of building, one matrix for each control condition
  ## each element of this list will become a matrix, which we will rbind together at end
  list_mat_comparison_msstats <- vector("list", length= length(vec_control_conditions))
  ## make the experimental conditions. These have to match what PDP assigns as vec_control_conditions, have to match what is in the dataset. For Spectronaut, this is R.Condition column. 
  vec_experimental_conditions <- levels(list_PD_proposed$PTM$ProteinLevelData$Condition)[!(levels(list_PD_proposed$PTM$ProteinLevelData$Condition) %in% vec_control_conditions)]
  ## how to automate this matrix, hmm
  
  ## First create matrices that are length(vec_experimental_conditions) rows long with rownames == vec_experimental_conditions and colnames == all.
  for (i in 1:length(list_mat_comparison_msstats)) {
    ## for now, the rownames do not contain the name of the control condition. We will replace later after grepl step. 
    list_mat_comparison_msstats_names01 <- list(mat_comparison_rows = paste0(vec_experimental_conditions),
                                                mat_comparison_cols = levels(list_PD_proposed$PTM$ProteinLevelData$Condition))
    list_mat_comparison_msstats[[i]] <- matrix(ncol = length(levels(list_PD_proposed$PTM$ProteinLevelData$Condition)),
                                               nrow = length(vec_experimental_conditions),
                                               dimnames = list_mat_comparison_msstats_names01)
  }
  
  ## This is the major matrix creation step
  for (j in 1:length(list_mat_comparison_msstats)) {
    for (i in 1:length(vec_experimental_conditions)) {
      a <- vec_experimental_conditions[i]
      b <- list_mat_comparison_msstats[[j]][i,]
      list_mat_comparison_msstats[[j]][i,] <- as.integer(grepl(a, names(b)))
    }
  }
  
  ## Now that we have that, we can assign the controls a value of -1
  for (i in 1:length(vec_control_conditions)) {
    list_mat_comparison_msstats[[i]][,grep(paste0("^", vec_control_conditions[i], "$"), colnames(list_mat_comparison_msstats[[i]]))] <- -1
  }
  
  ## last step is to add the control to the rownames
  for (i in 1:length(list_mat_comparison_msstats)) {
    list_mat_comparison_msstats_names02 <- paste0(vec_experimental_conditions, "-", vec_control_conditions[i])
    rownames(list_mat_comparison_msstats[[i]]) <-  list_mat_comparison_msstats_names02
    
  }
  
  ## Collapse into comparison matrix
  mat_comparison_msstats <- do.call(rbind, list_mat_comparison_msstats)
  ## BOOM

  MSstatsPTM.model <- MSstatsPTM::groupComparisonPTM(list_PD_proposed,
                                      data.type = "TMT",
                                      contrast.matrix = mat_comparison_msstats,
                                      use_log_file = FALSE)
  
  ## df from MSstats final
  df_p <- MSstatsPTM.model$PTM.Model
# save.image("./24_12_10_msstats_alldfs.RData")
## add col to have uniprot ID aka Entry
  df_p$Entry <- gsub("_.*$|-[[:digit:]]_.*$|-[[:digit:]][[:digit:]]_.*$", "", df_p$Protein)
  df_up <- read.csv("uniprotkb_Human_AND_model_organism_9606_2024_10_21.tsv" , sep = "\t", header = T, na.strings = c("NA", "NULL", "Null", "null", "NaN", "Na"))
  
  df_all <- merge(df_p, df_up, by = "Entry", all.x = T)
  # df_all <- df_all[which(is.na(df_all$issue)),]
length(unique(df_all$Protein)) 
## 2338 sites, not bad

list_of_genes_bryce <- c("TNFRSF25", "TRADD", "FADD", "TRAF1","TRAF2", "TRAF3", "TRAF6", "RIPK1", "RIPK3", "CIAP", "ASK1", "CASP8", "TANK", "NEMO", "CYLD", "LUBAC")
list_of_genes_bryce_regex <- paste0(list_of_genes_bryce, collapse = "|")
df_all$custom <- ifelse(grepl(list_of_genes_bryce_regex, df_all$Gene.Names), "Custom", NA)

## Need to compare this to the top hits of other experiment. Correlation? LogFC to control on x and y axis? 
df_og <- read.table("./tl1a/output/tl1a_tnfa_followup/tl1a_tnfa_phosphoproteomics_dea.txt", sep = "\t", na.strings = c("NULL", "NaN", "NA"), header = T)
df_og$Protein <- gsub("\\|.*$", "", df_og$ptm_sites)
df_og$PTM <- gsub("_p", "_", df_og$Protein)
df_og[df_og=="tl1a_000_2nM"] <- vec_experimental_conditions[1]
df_all$PTM <- gsub("_[[:alpha:]]", "_", df_all$Protein)
df_all$PTM_helper <- str_split_fixed(df_all$PTM, "_", n = 2)[,2]
df_all$PTM_helper <- gsub("_.*$", "", df_all$PTM_helper)
df_all$PTM <- paste0(
  str_split_fixed(df_all$PTM, "_", n = 2)[,1], 
  "_",
  df_all$PTM_helper
)

i <- 5
df_plot <- df_all[grep(vec_experimental_conditions[i], df_all$Label),]
df_plot_og <- df_og[which(df_og$contrast==vec_experimental_conditions[i]),]



## Plot individual contrasts/labels
df_m <- merge(df_plot, df_plot_og, by = "PTM")
df_m$site <- gsub("^.*_", "", df_m$Protein.x)
num_log_threshold <- 0.5
df_m$plot_label <- ifelse(abs(df_m$log2FC)> num_log_threshold & abs(df_m$logfc) > num_log_threshold, paste0(df_m$Gene.Names..primary., "\n", df_m$site), NA)
df_m$LogFC_threshold <- ifelse(!is.na(df_m$plot_label), "yes", "no")

num_cor <- cor(df_m$log2FC[which(is.finite(df_m$log2FC))], df_m$logfc[which(is.finite(df_m$log2FC))], use = "everything", method = "pearson")

ggplot(data = df_m[which(is.finite(df_m$log2FC)),], aes(x = log2FC, y = logfc, label = plot_label, color = LogFC_threshold)) +
  geom_abline(alpha = 0.5, linetype = "dotted") +
  geom_hline(yintercept = 0, alpha = 0.5) +
  geom_vline(xintercept = 0, alpha = 0.5) +
  geom_point() +
  scale_x_continuous(breaks = seq(-2,3,by =1), limits = c(-2,3)) +
  scale_y_continuous(breaks = seq(-2,3,by =1), limits = c(-2,3)) +
  ggrepel::geom_label_repel(show.legend = FALSE, max.overlaps = 30, force = 5, alpha = 0.8) +
  scale_color_viridis_d(option = "magma", begin = 0.1, end = 0.7) +
  xlab("Total proteome experiment Log2FC") +
  ylab("Phosphoproteomic experiment Log2FC") +
  ggtitle(paste0(vec_experimental_conditions[i])) +
  annotate(geom = "text", x = -1.5, y = 3, label = paste0("Pearson's cor. = ", round(num_cor, 2)))


## Some stats.
length(unique(df_og$id))
## 27,739 total peptides, whether site is present or not
length(unique(df_og$Protein))
## 15,328 ("Protein" is with isoforms removed, best we can do to merge)
length(unique(df_p$Protein))
## 2,338 total peptides?

df_unique_all <- unique(df_all[,c("Protein","PTM")])
## 2,338
df_unique_og <- unique(df_og[,c("id", "ptm_sites", "PTM")])
## 27,739
length(unique(df_unique_og$ptm_sites))
## 16,076 actual unique sites with isoforms removed
sum(is.na(df_unique_og$ptm_sites))
## 7,020 that cannot be assigned
df_unique_og_edit <- unique(df_og[,c("ptm_sites", "PTM")])
length(unique(df_unique_og_edit$ptm_sites))
length(unique(df_unique_og_edit$PTM))
## of those, 15,328 are the same site
vec_unique_all <- unique(df_unique_all$PTM)
vec_unique_og <- unique(df_unique_og_edit$PTM)
sum(!is.na(match(vec_unique_all, vec_unique_og)))
## 1047 total


## Let's also examine the venn diagram of PTM sites captured for each experiment. How many of the top sites that we see in the phos experiment do we see in the total proteome experiment?
## Merge in the site data left join phos, TP, all = T for global rankings

## First reshape the dfs for eventual merge. df_p and df_og
df_p_w <- reshape(df_p[,c("Protein", "Label", "log2FC", "adj.pvalue")], idvar = "Protein", timevar = "Label", direction = "wide")
df_p_w$PTM <- gsub("_[[:alpha:]]", "_", df_p_w$Protein)
df_p_w$PTM_helper <- str_split_fixed(df_p_w$PTM, "_", n = 2)[,2]
df_p_w$PTM_helper <- gsub("_.*$", "", df_p_w$PTM_helper)
df_p_w$PTM <- paste0(
  str_split_fixed(df_p_w$PTM, "_", n = 2)[,1], 
  "_",
  df_p_w$PTM_helper
)
df_p_w$PTM_helper <- NULL


# df_og_na <- df_og[which(!is.na(df_og$PTM)),]
# df_og_1 <- group_by(df_og_na, PTM, contrast)
# df_og_2 <- summarize(df_og_1, avg = mean(logfc))
df_og_w <- reshape(df_og[,c("id", "contrast", "logfc", "padj")], idvar = "id", timevar = "contrast", direction = "wide")
df_og_w <- merge(df_unique_og, df_og_w, by = "id")
## drop the NA. since can't do anything with those
df_og_w <- df_og_w[which(!is.na(df_og_w$ptm_sites)),]
## 15,328 unique peptides. Unique peptides, nonunique sites. Ask if duplicated, if yes, then do mean(padj of all), take peptide with lowest pval rbind back to unique. do this before the reshape.

## Went back to read Bryce's email. He asked for the top 50 p-sites from expt 1, ranked by TL1a 20 nM treatment. Let's make that df, then merge with TP. 
df_og_w <- df_og_w[order(-df_og_w$logfc.tl1a_020nM),]
df_og_w$id <- NULL
df_og_w$ptm_sites <- NULL
df_og_w$dup <- NULL
df_og_top <- df_og_w[1:50,]

names(df_og_top) <- gsub("logfc", "log2FC", names(df_og_top))
names(df_og_top) <- gsub("padj", "adj.pvalue", names(df_og_top))
names(df_p_w) <- gsub("-control", "", names(df_p_w))
names(df_og_top)[which(!(names(df_og_top) %in% c("Protein", "PTM")))] <- paste0("expt1_", names(df_og_top)[which(!(names(df_og_top) %in% c("Protein", "PTM")))])
names(df_og_top)
names(df_p_w)[which(!(names(df_p_w) %in% c("Protein", "PTM")))] <- paste0("expt2_", names(df_p_w)[which(!(names(df_p_w) %in% c("Protein", "PTM")))])
names(df_p_w)

df_top50_merged <- merge(df_og_top, df_p_w, by = "PTM")
df_top50_merged$Entry <- gsub("-.*$|_.*$", "", df_top50_merged$PTM)
df_top50_merged <- merge(df_top50_merged, df_up, by = "Entry", all.x = T)
df_top50_merged <- df_top50_merged[,c(1:2, 24:28, 3:23)]
if(!file.exists("./df_top50_merged.csv")) {
  write.csv(df_top50_merged, "./df_top50_merged.csv")
}

sum(is.na(match(unique(df_og_top$PTM) ,unique(df_top50_merged$PTM))))
# save.image("./24_12_13_df_top50_merged.RData")


## Interactive plot
# i <- 5
# df_plot <- df_all[grep(vec_experimental_conditions[i], df_all$Label),]
# g <- ggplot() +
#   geom_point(data = df_plot, aes(x = log2FC, y = -log10(adj.pvalue),  text=sprintf(
#     "Protein: %s<br>Site: %s<br>LogFC: %s<br>-log10(pval): %s", 
#     Gene.Names, 
#     Protein, 
#     round(log2FC, 2),
#     round(-log10(adj.pvalue),2)
#   )), 
#   alpha = 0.5, shape = 1, color = "gray") +
#   geom_point(data = df_plot[which(!is.na(df_plot$custom)),], aes(x = log2FC, y = -log10(adj.pvalue),  text=sprintf(
#     "Protein: %s<br>Site: %s<br>LogFC: %s<br>-log10(pval): %s", 
#     Gene.Names, 
#     Protein, 
#     round(log2FC, 2),
#     round(-log10(adj.pvalue),2)
#   ))) +
#   scale_color_viridis_c(begin = 0.1, end = 0.9, option = "magma") +
#   # geom_hline(yintercept= -log10(alpha), linetype = "dashed", alpha = 0.2) +
#   # labs(color = "Log2FC(TL1a 20nM/control)") +
#   xlab("Log2FC(experimental/control)") +
#   ylab("-log10(pval)") +
#   ggtitle(paste0(unique(df_plot$Label))) 
# 
# ggplotly(g, tooltip = c("text"))