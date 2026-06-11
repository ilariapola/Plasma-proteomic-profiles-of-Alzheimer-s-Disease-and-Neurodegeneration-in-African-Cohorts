#### LOAD PACKAGES ####

library(dplyr)
library(tidyr)
library(purrr)
library(broom)
library(ggplot2)
library(ggrepel)
library(table1)

#### DEMOGRAPHICS TABLE ####

final$DX_AB <- factor(final$DX_AB, levels = c("CU-", "CU+", "MCI+", "Dem+", "MCI-", "Dem-"))

render_mean_sd <- function(x) {
  with(stats.apply.rounding(stats.default(x), digits=2), 
       sprintf("%s (&plusmn; %s)", MEAN, SD))
}

table1(~ sex + age + years_edu + AB_status| DX_AB, 
       data = final, 
       render.continuous = render_mean_sd)


              
#### VOLCANO PLOT ####

##### VOLCANO PLOT amyloid positive versus amyloid negative #####

# Columns 85:207 correspond to protein measurements in the `final` dataset.
ex <- final[, c(85:207)]

# Run one linear model per protein.
# Main exposure: AB_pos
# AB_pos compares amyloid-positive participants versus amyloid-negative participants.
# Note: if AB_pos is coded as 0/1 or factor with 0 as reference,
# the term AB_pos1 represents amyloid-positive vs amyloid-negative.


ex <- final[, c(85:207)] #85:207 represent the columns of the proteomic data 
results_list <- list()

for (protein in colnames(ex)) {
  model <- lm(as.formula(paste(protein, "~ AB_status + age + sex + average_protein")), data = final)
  
  model_results <- tidy(model)
  
  model_results_filtered <- model_results %>%
    dplyr::filter(term == "AB_status") %>%
    dplyr::mutate(Protein = protein,        
                  logFC = estimate,        
                  log10P = -log10(p.value)) 
  
  results_list[[protein]] <- model_results_filtered
}


final_results <- bind_rows(results_list)

final_results <- final_results %>%
  dplyr::mutate(adj_p_value = p.adjust(p.value, method = "BH")) 

final_results_P <- final_results

tT2_plot <- final_results_P %>%
  dplyr::mutate(diffexpressed = case_when(
    adj_p_value < 0.05 ~ "FDR_P", 
    p.value < 0.05 & adj_p_value >= 0.05 ~ "Pval_P",
    p.value >= 0.05 ~ "NO"    
  ))
colnames(tT2_plot)


tT2_plot$delabel <- NA
tT2_plot$delabel[tT2_plot$diffexpressed != "NO"] <- tT2_plot$Protein[tT2_plot$diffexpressed != "NO"]

tT2_plot <- tT2_plot %>%
  mutate(
    delabel = ifelse(is.na(delabel), "", delabel),
    sig = case_when(
      diffexpressed == "FDR_P"  ~ "FDR",
      diffexpressed == "Pval_P" ~ "nominal",
      TRUE                      ~ "none"
    ),
    dir = case_when(
      logFC > 0 ~ "pos",
      logFC < 0 ~ "neg",
      TRUE      ~ "zero"
    ),
    color_group = case_when(
      sig == "FDR"     & dir == "pos" ~ "FDR_pos",      
      sig == "nominal" & dir == "pos" ~ "nominal_pos", 
      sig == "FDR"     & dir == "neg" ~ "FDR_neg",     
      sig == "nominal" & dir == "neg" ~ "nominal_neg",  
      TRUE                            ~ "none"          
    ),
    color_group = factor(color_group,
                         levels = c("FDR_pos","nominal_pos","FDR_neg","nominal_neg","none"))
  )

vc <- ggplot(tT2_plot, aes(x = logFC, y = log10P, color = color_group, label = delabel)) +
  geom_point() +
  scale_color_manual(values = c(
    "FDR_pos"     = "#a63603", 
    "nominal_pos" = "#e6550d", 
    "FDR_neg"     = "#542788",  
    "nominal_neg" = "#756bb1",  
    "none"        = "#bdbdbd"  
  ))+
  geom_text_repel(family = "sans", size = 3) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#949395") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#949395") +
  labs(x = "Standardized Beta", y = "-log10 (P-value)") +
  theme(
    text = element_text(family = "sans", size = 7),
    axis.title.x = element_text(size = 11, margin = margin(t = 10), family = "sans"),
    axis.title.y = element_text(size = 11, margin = margin(r = 10), family = "sans"),
    panel.background = element_rect(fill = "white"),
    axis.line = element_line(colour = "#565555", linewidth = 0.5, linetype = "solid"),
    panel.grid.major = element_line(linewidth = 0.15, linetype = "solid", colour = "#E1DFDF"),
    panel.grid.minor = element_line(linewidth = 0.15, linetype = "solid", colour = "#E1DFDF"),
    legend.position = "none"
  )

print(vc)



##### VOLCANO PLOT cognitive unimpaired versus cognitive impaired (amyloid positive or amyloid negative) #####

# separately within amyloid-negative or amyloid-positive participants.

# Amyloid-negative subgroup
df_sub <- final %>%
  dplyr::filter(AB_status == 0)

# Amyloid-positive subgroup
# df_sub <- final %>%
#   dplyr::filter(AB_status == 1)

# Select proteomic variables.
# Columns 85:207 correspond to protein measurements.
ex <- df_sub[, c(85:207)]

# Run one linear model per protein.
# Main exposure: CU_CI (as 0 and 1 for cognitively unimparied and 1 as cognitively impaired)
# CU_CI compares cognitively impaired versus cognitively unimpaired participants
# within the selected amyloid subgroup.
results_list <- list()

for (protein in colnames(ex)) {
  model <- lm(
    as.formula(paste(protein, "~ CU_CI + age + sex + average_protein")),
    data = df_sub
  )
  
  model_results <- tidy(model)
  
  # Extract only the CUCI group comparison
  model_results_filtered <- model_results %>%
    dplyr::filter(term == "CU_CI") %>%
    dplyr::mutate(
      Protein = protein,
      logFC = estimate,
      log10P = -log10(p.value)
    )
  
  results_list[[protein]] <- model_results_filtered
}
 
final_results <- bind_rows(results_list)

final_results <- final_results %>%
  dplyr::mutate(adj_p_value = p.adjust(p.value, method = "BH")) 

final_results_P <- final_results

tT2_plot <- final_results_P %>%
  dplyr::mutate(diffexpressed = case_when(
    adj_p_value < 0.05 ~ "FDR_P", 
    p.value < 0.05 & adj_p_value >= 0.05 ~ "Pval_P",
    p.value >= 0.05 ~ "NO"    
  ))
colnames(tT2_plot)


tT2_plot$delabel <- NA
tT2_plot$delabel[tT2_plot$diffexpressed != "NO"] <- tT2_plot$Protein[tT2_plot$diffexpressed != "NO"]

tT2_plot <- tT2_plot %>%
  mutate(
    delabel = ifelse(is.na(delabel), "", delabel),
    sig = case_when(
      diffexpressed == "FDR_P"  ~ "FDR",
      diffexpressed == "Pval_P" ~ "nominal",
      TRUE                      ~ "none"
    ),
    dir = case_when(
      logFC > 0 ~ "pos",
      logFC < 0 ~ "neg",
      TRUE      ~ "zero"
    ),
    color_group = case_when(
      sig == "FDR"     & dir == "pos" ~ "FDR_pos",      
      sig == "nominal" & dir == "pos" ~ "nominal_pos", 
      sig == "FDR"     & dir == "neg" ~ "FDR_neg",     
      sig == "nominal" & dir == "neg" ~ "nominal_neg",  
      TRUE                            ~ "none"          
    ),
    color_group = factor(color_group,
                         levels = c("FDR_pos","nominal_pos","FDR_neg","nominal_neg","none"))
  )

vc <- ggplot(tT2_plot, aes(x = logFC, y = log10P, color = color_group, label = delabel)) +
  geom_point() +
  scale_color_manual(values = c(
    "FDR_pos"     = "#a63603", 
    "nominal_pos" = "#e6550d", 
    "FDR_neg"     = "#542788",  
    "nominal_neg" = "#756bb1",  
    "none"        = "#bdbdbd"  
  ))+
  geom_text_repel(family = "sans", size = 3) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#949395") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#949395") +
  labs(x = "Standardized Beta", y = "-log10 (P-value)") +
  theme(
    text = element_text(family = "sans", size = 7),
    axis.title.x = element_text(size = 11, margin = margin(t = 10), family = "sans"),
    axis.title.y = element_text(size = 11, margin = margin(r = 10), family = "sans"),
    panel.background = element_rect(fill = "white"),
    axis.line = element_line(colour = "#565555", linewidth = 0.5, linetype = "solid"),
    panel.grid.major = element_line(linewidth = 0.15, linetype = "solid", colour = "#E1DFDF"),
    panel.grid.minor = element_line(linewidth = 0.15, linetype = "solid", colour = "#E1DFDF"),
    legend.position = "none"
  )

print(vc)

#### STACKED HISTOGRAM ####

#based on the categotization
neg_slices <- c(5, 1, 4, 4)
names(neg_slices) <- c("Synaptic","Vascular","Neurodegeneration","Inflammation")

#based on the categotization
pos_slices <- c(8, 2, 1, 4, 7)
names(pos_slices) <- c("Amyloid & Tau","Synaptic","Vascular","Neurodegeneration","Inflammation")


cols <- c(
  "Amyloid & Tau"     = "#a796bd",  
  "Synaptic"          = "#abb695",  
  "Vascular"          = "#a0d6cb",  
  "Neurodegeneration" = "#9b896b",  
  "Inflammation"      = "#edb63f" 
)

df <- bind_rows(
  data.frame(group = "Negative", category = names(neg_slices), value = as.numeric(neg_slices)),
  data.frame(group = "Positive", category = names(pos_slices), value = as.numeric(pos_slices))
) %>%
  group_by(group) %>%
  mutate(pct = value / sum(value) * 100) %>%
  ungroup()

cat_levels <- c("Amyloid & Tau","Synaptic","Vascular","Neurodegeneration","Inflammation")
df$category <- factor(df$category, levels = cat_levels)


p <- ggplot(df, aes(x = group, y = pct, fill = category)) +
  geom_col(width = 0.4, color = "white", linewidth = 0.5) +
  scale_fill_manual(values = cols, drop = FALSE) +
  scale_y_continuous(expand = c(0,0), limits = c(0,100), labels = function(x) paste0(x, "%")) +
  labs(x = NULL, y = "Percent", fill = NULL,
       title = "Category Composition: Positive vs Negative") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    plot.title = element_text(face = "bold", hjust = 0.5)
  )


#### COMORBIDITIES ####

# Proteins included in the analysis 

# proteins are the significant proteins from the differential expression analysis for amyloid positive versus negative
proteins_keep <- c(
  "pTau231","pTau181","IL5","IL13","pTau217","MAPT","FABP3","IGF1R",
  "FOLR1","IL33","VEGFD","CALB2","VCAM1","KLK6","REST","MSLN","CX3CL1",
  "GOT1","CD63","IGFBP7","ICAM1","CST3","NGF","ARSA","NRGN","HTT","ENO2",
  "BDNF","PGK1","CD40LG","IL6R","VGF","ANXA5","CNTN2","TIMP3","FCN2"
)
# proteins are the significant proteins from the differential expression analysis for CU vs CI amyloid positive
# proteins_keep <- c( "TAFA5","CCL11","NPTX1","GDI1","pTau217","KLK6", "UCHL1","TNF", "TIMP3", "MME", "FCN2")

# proteins are the significant proteins from the differential expression analysis for CU vs CI amyloid negative
#proteins_keep <- c("CHIT1", "NEFL", "GFAP", "TREM1", "VEGFD", "CALB2", "MSLN", "SFTPD", "PDGFRB", 
                   #"NPTX2", "SNCB", "GDI1", "SLIT2", "IGFBP7", "SNAP25", "UCHL1", "APOE", "VCAM1", 
                   #"NGF", "HTT", "S100B", "CSF2", "CCL26", "IFNG", "BASP1")


# Comorbidity exposures tested
exposures <- c(
  "heart_disease",
  "lipid_metabolism",
  "comorbidities", #indicating other pathologies as described in the methods 
  "inflammatory",
  "infectious",
  "cholesterol"
)

# Covariates included in each model
covariates <- c("AB_status", "age", "sex", "average_protein")

# Full model predictors
predictors <- c(exposures, covariates)

# For each protein, the following linear model was fitted:
# protein ~ comorbidity exposure + AB_status + age + sex + average_protein
#
# The model includes all comorbidity exposures simultaneously.
# AB_pos compares amyloid-positive versus amyloid-negative participants.
# average_protein was included as a technical covariate.

# Fit linear models for each protein

fit_one <- function(protein) {
  
  model <- lm(
    reformulate(predictors, response = protein),
    data = final
  )
  
  tidy(model, conf.int = TRUE) %>%
    filter(term %in% exposures) %>%
    mutate(
      Protein = protein,
      Exposure = term
    )
}

all_results <- map_dfr(proteins_keep, fit_one) %>%
  group_by(Exposure) %>%
  mutate(
    p_fdr = p.adjust(p.value, method = "BH")
  ) %>%
  ungroup()

# Plot comorbidity associations #

exposure_labels <- c(
  cholesterol         = "Cholesterol",
  heart_disease    = "Heart diseases",
  lipid_metabolism = "Metabolic diseases",
  comorbidities    = "Other conditions",
  inflammatory     = "Respiratory diseases",
  infectious       = "Infectious diseases"
)

pal <- c(
  cholesterol        = "#3b8b4f",
  heart_disease    = "#db4260",
  lipid_metabolism = "#555876",
  comorbidities    = "#cbb8a9",
  inflammatory     = "#e7a35b",
  infectious       = "#8c9e9a"
)

plot_beta <- all_results %>%
  mutate(
    sig_label = case_when(
      p_fdr < 0.05 ~ "**",
      p.value < 0.05 ~ "*",
      TRUE ~ ""
    ),
    Protein = factor(Protein, levels = rev(proteins_keep)),
    Exposure = factor(
      Exposure,
      levels = c(
        "cholesterol",
        "heart_disease",
        "infectious",
        "inflammatory",
        "lipid_metabolism",
        "comorbidities"
      )
    )
  ) %>%
  group_by(Exposure) %>%
  mutate(
    nudge = 0.03 * (max(conf.high, na.rm = TRUE) - min(conf.low, na.rm = TRUE))
  ) %>%
  ungroup()

p_beta <- ggplot(plot_beta, aes(y = Protein)) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.25,
    colour = "grey55"
  ) +
  geom_errorbarh(
    aes(xmin = conf.low, xmax = conf.high, colour = Exposure),
    height = 0,
    linewidth = 0.6
  ) +
  geom_point(
    aes(x = estimate, colour = Exposure),
    size = 2
  ) +
  geom_text(
    aes(x = estimate + nudge, label = sig_label),
    size = 5,
    colour = "black"
  ) +
  facet_wrap(
    ~ Exposure,
    nrow = 1,
    scales = "free_x",
    labeller = as_labeller(exposure_labels)
  ) +
  scale_colour_manual(values = pal, guide = "none") +
  scale_x_continuous(
    expand = expansion(mult = c(0.05, 0.10))
  ) +
  labs(
    x = "β (95% CI)",
    y = NULL,
    subtitle = "* p < 0.05; ** FDR q < 0.05 (BH)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major.x = element_line(colour = "grey88", linewidth = 0.3),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "none"
  )

print(p_beta)


#### LOESS PLOTS ####

##### LOESS PLOTS AMYLOID POSITIVE #####
# Z-score proteins using CU- as reference group #

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

proteins_to_scale <- c(
  "ANXA5", "BDNF", "CALB2", "CD63", "CX3CL1", "FCN2", "GOT1",
  "HTT", "IL13", "IL5", "KLK6", "NGF", "NRGN", "VEGFD"
)

reference_group <- final %>%
  filter(DX_AB == "CU-")

final_zscored <- final

for (protein in proteins_to_scale) {
  
  reference_mean <- mean(reference_group[[protein]], na.rm = TRUE)
  reference_sd   <- sd(reference_group[[protein]], na.rm = TRUE)
  
  final_zscored[[protein]] <- 
    (final[[protein]] - reference_mean) / reference_sd
}

final_zscored <- final_zscored %>%
  filter(
    DX_AB %in% c("CU-", "CU+", "MCI+", "Dem+"),
    !is.na(DX_AB),
    !is.na(sex)
  ) %>%
  mutate(
    DX_AB = factor(DX_AB, levels = c("CU-", "CU+", "MCI+", "Dem+"))
  )

df_long <- final_zscored %>%
  select(DX_AB, all_of(proteins_to_scale)) %>%
  tidyr::pivot_longer(
    cols = all_of(proteins_to_scale),
    names_to = "Protein",
    values_to = "Z_score"
  )

df_summary <- df_long %>%
  group_by(DX_AB, Protein) %>%
  summarise(
    Z_score_mean = mean(Z_score, na.rm = TRUE),
    .groups = "drop"
  )

# Plot protein trajectories across disease groups #

p_trajectory <- ggplot() +
  geom_smooth(
    data = df_long,
    aes(x = DX_AB, y = Z_score, color = Protein, group = Protein),
    method = "loess",
    span = 0.8,
    se = FALSE,
    size = 0.8
  ) +
  scale_color_hue() +  
  theme_minimal() +
  labs(
    x = "Disease Groups",
    y = "Z-scored Proteins",
    title = "LOESS Plot across Disease",
    color = "Protein"
  ) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.position = "right"
  ) +
  coord_cartesian(ylim = c(-2, 1.5))


##### LOESS PLOTS AMYLOID NEGATIVE  #####
# Z-score proteins using CU- as reference group #

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

proteins_to_scale <- c(
  "ANXA5", "BDNF", "CALB2", "CD63", "CX3CL1", "FCN2", "GOT1",
  "HTT", "IL13", "IL5", "KLK6", "NGF", "NRGN", "VEGFD"
)

reference_group <- final %>%
  filter(DX_AB == "CU-")

final_zscored <- final

for (protein in proteins_to_scale) {
  
  reference_mean <- mean(reference_group[[protein]], na.rm = TRUE)
  reference_sd   <- sd(reference_group[[protein]], na.rm = TRUE)
  
  final_zscored[[protein]] <- 
    (final[[protein]] - reference_mean) / reference_sd
}

final_zscored <- final_zscored %>%
  filter(
    DX_AB %in% c("CU-", "MCI-", "Dem-"),
    !is.na(DX_AB),
    !is.na(sex)
  ) %>%
  mutate(
    DX_AB = factor(DX_AB, levels = c("CU-", "MCI-", "Dem-"))
  )

df_long <- final_zscored %>%
  select(DX_AB, all_of(proteins_to_scale)) %>%
  tidyr::pivot_longer(
    cols = all_of(proteins_to_scale),
    names_to = "Protein",
    values_to = "Z_score"
  )

df_summary <- df_long %>%
  group_by(DX_AB, Protein) %>%
  summarise(
    Z_score_mean = mean(Z_score, na.rm = TRUE),
    .groups = "drop"
  )

# Plot protein trajectories across disease groups #

p_trajectory <- ggplot() +
  geom_smooth(
    data = df_long,
    aes(x = DX_AB, y = Z_score, color = Protein, group = Protein),
    method = "loess",
    span = 0.8,
    se = FALSE,
    size = 0.8
  ) +
  scale_color_hue() +  
  theme_minimal() +
  labs(
    x = "Disease Groups",
    y = "Z-scored Proteins",
    title = "LOESS Plot across Disease",
    color = "Protein"
  ) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.position = "right"
  ) +
  coord_cartesian(ylim = c(-2, 1.5))
