# Phase I analysis adapted from lca_work_Helen_v6.R.
# Usage: Rscript analysis/lca_phase1.R /absolute/path/to/cleaned.csv /absolute/output/dir
# The preconstructed indicators and weights must be supplied; this script
# does not reconstruct or audit the upstream NHANES cleaning pipeline.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("Provide cleaned CSV and output directory.")
input_file <- normalizePath(args[[1]], mustWork = TRUE)
dir.create(args[[2]], recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(args[[2]], mustWork = TRUE)
setwd(output_dir)

library(tidyverse)
library(data.table)
library(poLCA)
library(corrplot)
library(ggcorrplot)
library(FactoMineR)
library(factoextra)
library(gtsummary)
library(xtable)
library(gridExtra)
library(survey)

theme_gtsummary_journal(journal = "jama")

# read data and recode
df_raw = fread(file = input_file)

df = df_raw %>%
  mutate(
    real_age = age,
    # --- V4 CHANGE START:age groups ---
    age_group = factor(
      case_when(
        age <= 26 ~ "≤26 years",
        age >= 27 & age <= 45 ~ "27-45 years",
        age >= 46 ~ "≥46 years"
      ),
      levels = c("≤26 years", "27-45 years", "≥46 years")
    ),
    # --- V4 CHANGE END ---

    gender_label    = as.factor(gender),
    education_label = as.factor(recode(education,
                                       "Less than 9th grade" = "High school or less",
                                       "High school graduate/GED or equivalent" = "High school or less",
                                       "9-11th grade (Includes 12th grade with no diploma)" = "High school or less",
                                       "Some college or AA degree" = "Some college or more",
                                       "College graduate or above" = "Some college or more")),
    marital_label   = as.factor(recode(marital,
                                       "Married/Partnered" = "Married/Partnered",
                                       "Never Married" = "Never/Previously Married",
                                       "Previously Married" = "Never/Previously Married")),

    # --- recode 0/1 to 1/2 ---
    high_depression = ifelse(high_depression == 1, 2, 1),
    hpv_vaccinated  = ifelse(hpv_vaccinated == 1, 2, 1),
    mets_waist      = ifelse(mets_waist == 1, 2, 1),
    mets_bp         = ifelse(mets_bp == 1, 2, 1),
    mets_glu        = ifelse(mets_glu == 1, 2, 1),
    mets_hdl        = ifelse(mets_hdl == 1, 2, 1),
    mets_trig       = ifelse(mets_trig == 1, 2, 1),

    high_depression_f = factor(ifelse(high_depression == 2, "Yes", "No"), levels = c("Yes", "No")),
    hpv_vaccinated_f  = factor(ifelse(hpv_vaccinated == 2, "Yes", "No"),  levels = c("Yes", "No")),
    mets_waist_f      = factor(ifelse(mets_waist == 2, "Yes", "No"),     levels = c("Yes", "No")),
    mets_bp_f         = factor(ifelse(mets_bp == 2, "Yes", "No"),        levels = c("Yes", "No")),
    mets_glu_f        = factor(ifelse(mets_glu == 2, "Yes", "No"),       levels = c("Yes", "No")),
    mets_hdl_f        = factor(ifelse(mets_hdl == 2, "Yes", "No"),       levels = c("Yes", "No")),
    mets_trig_f       = factor(ifelse(mets_trig == 2, "Yes", "No"),      levels = c("Yes", "No"))

  )

# --- V6 CHANGE START: preserve full sample and define primary sample ---

df = df %>%
  filter(wtsaf2yrad > 0)
# --- V6 CHANGE END ---

lca_vars = c("high_depression", "hpv_vaccinated",
             "mets_waist", "mets_bp", "mets_glu", "mets_hdl", "mets_trig")

df_lca = df %>% dplyr::select(all_of(lca_vars)) %>% na.omit()

# ==============================================================================
# Overall Participant Characteristics
# ==============================================================================
# --- V6 CHANGE START: Table 1 unweighted and weighted comparison ---
nhanes_design = svydesign(
  ids = ~sdmvpsu,
  strata = ~sdmvstra,
  weights = ~wtsaf2yrad,
  nest = TRUE,
  data = df
)

make_table1_rows = function(variable, label, design) {
  values = design$variables[[variable]]
  response_levels = levels(values)

  unweighted_n = table(factor(values, levels = response_levels))
  unweighted_prop = prop.table(unweighted_n)
  weighted_prop = coef(
    svymean(
      as.formula(paste0("~", variable)),
      design = design,
      na.rm = TRUE
    )
  )

  bind_rows(
    tibble(
      Characteristic = label,
      `Unweighted n` = NA_integer_,
      `Unweighted %` = NA_real_,
      `Weighted %` = NA_real_
    ),
    tibble(
      Characteristic = paste0("  ", response_levels),
      `Unweighted n` = as.integer(unweighted_n),
      `Unweighted %` = round(as.numeric(unweighted_prop) * 100, 1),
      `Weighted %` = round(as.numeric(weighted_prop) * 100, 1)
    )
  )
}

table1_variables = c(
  age_group = "Age group",
  gender_label = "Gender",
  education_label = "Education level",
  marital_label = "Marital status",
  high_depression_f = "High depression symptoms",
  hpv_vaccinated_f = "HPV vaccinated",
  mets_waist_f = "Elevated waist circumference",
  mets_bp_f = "Elevated blood pressure",
  mets_glu_f = "Elevated fasting glucose",
  mets_hdl_f = "Low HDL cholesterol",
  mets_trig_f = "Elevated triglycerides"
)

tbl0 = map2_dfr(
  names(table1_variables),
  unname(table1_variables),
  ~ make_table1_rows(.x, .y, nhanes_design)
)

print(tbl0)
# --- V6 CHANGE END ---

tbl0_export = tbl0
names(tbl0_export)[names(tbl0_export) == "Unweighted n"] = paste0(
  "Unweighted n (N = ",
  format(nrow(df), big.mark = ","),
  ")"
)

write.csv(tbl0_export, "Table_Participant_Characteristics.csv", row.names = FALSE, na = "")
# print(xtable(as.data.frame(tbl0), caption = "Overall Participant Characteristics"),
#       file = "Table_Participant_Characteristics.tex", include.rownames = FALSE)

# ==============================================================================
#  Model Fit (1 to 5 Classes)
# ==============================================================================
df_lca = df %>%
  dplyr::select(high_depression, hpv_vaccinated, mets_waist, mets_trig, mets_hdl, mets_bp, mets_glu)

# Formula
f = as.formula(paste(paste("cbind(", paste(colnames(df_lca), collapse = ","), ")"), "~1"))

models = list()
nk = 5
set.seed(123)

for (k in 1:nk) {
  models[[k]] = poLCA(f, nclass=k, nrep = 10, data=df_lca, verbose = F, na.rm = F)
  ps = models[[k]]$probs.start
  models[[k]] = poLCA(f, nclass=k, nrep = 10, data=df_lca, verbose = F, na.rm = F, probs.start = ps)
}

# AIC/BIC
measures = map_df(models, magrittr::extract, c("aic", "bic")) %>%
  rename("AIC" = "aic", "BIC" = "bic") %>%
  rowid_to_column("Number of Classes")

# Calculate Entropy
calc_entropy = function(mod) {
  p = mod$posterior
  if (ncol(p) == 1) return(1.0)
  1 - (-sum(p * log(p + 1e-12)) / (nrow(p) * log(ncol(p))))
}

# Add Entropy and Smallest Class Size：
measures$Entropy = sapply(models, calc_entropy)
measures$`Smallest Class Size (%)` = sapply(models, function(m) round(min(prop.table(table(m$predclass)) * 100), 2))

print(measures)

write.csv(measures, "Table_Model_Fit.csv", row.names = FALSE)
# print(xtable(measures, caption = "LCA Model Fit Statistics (1-5 Classes)"),
#       file = "Table_Model_Fit.tex", include.rownames = FALSE)

# ==============================================================================
# Latent Class Analysis (PHASE I)
# Model Selection Graphs
# ==============================================================================
# AIC & BIC
fig1_data = measures %>%
  dplyr::select(`Number of Classes`, AIC, BIC) %>%
  pivot_longer(cols = c(AIC, BIC), names_to = "Metric", values_to = "Value")

p1 = ggplot(fig1_data, aes(x = as.factor(`Number of Classes`), y = Value, group = Metric, color = Metric)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  labs(title = "Figure 1: LCA Model Selection (AIC / BIC)", x = "Number of Classes", y = "Criterion Value") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("Figure1_Model_Selection_AIC_BIC.png", p1, width = 16, height = 10, units = "cm", dpi = 300)

# AIC/BIC + Class size pie
class_sizes_all = data.frame()
for (k in 1:5) {
  sizes = as.data.frame(prop.table(table(models[[k]]$predclass)) * 100)
  colnames(sizes) = c("Class", "Pct")
  sizes$Model = paste0(k, "-Class Model")
  class_sizes_all = rbind(class_sizes_all, sizes)
}

p2 = ggplot(class_sizes_all, aes(x = "", y = Pct, fill = Class)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar(theta = "y", start = 0) +
  facet_wrap(~ Model, nrow = 1) +
  theme_void() +
  labs(title = "Class Size Distribution across Models", fill = "Class") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    strip.text = element_text(face = "bold", size = 10),
    legend.position = "bottom"
  )

fig1_1 = grid.arrange(p1, p2, ncol = 1, heights = c(2, 1))

ggsave("Model_Selection_with_Pie.png", fig1_1, width = 20, height = 18, units = "cm", dpi = 300)

# ==============================================================================
# Model Diagnostics & Select Best Model (k = 3)
# ==============================================================================
selected_k = 3

final_model = models[[selected_k]]
df$Classes = final_model$predclass
df_lca$Classes = final_model$predclass

# 1. Posterior Probabilities Matrix
post_probs = final_model$posterior
ave_post_matrix = matrix(NA, nrow = selected_k, ncol = selected_k)
rownames(ave_post_matrix) = paste0("Assigned Class ", 1:selected_k)
colnames(ave_post_matrix) = paste0("Prob Class ", 1:selected_k)

for (i in 1:selected_k) {
  members = which(final_model$predclass == i)
  if (length(members) > 0) {
    ave_post_matrix[i, ] = colMeans(post_probs[members, , drop = FALSE])
  }
}

formatted_matrix = matrix("", nrow = selected_k, ncol = selected_k)
rownames(formatted_matrix) = rownames(ave_post_matrix)
colnames(formatted_matrix) = colnames(ave_post_matrix)

for (i in 1:selected_k) {
  for (j in 1:selected_k) {
    val = ave_post_matrix[i, j]

    if (!is.na(val)) {
      if (val < 0.001) {
        formatted_matrix[i, j] = "<0.001"
      } else {
        formatted_matrix[i, j] = sprintf("%.3f", val)
      }
    }
  }
}

write.csv(formatted_matrix, "Table_Model_Diagnostics_Posterior.csv")
# print(xtable(formatted_matrix, caption = "Average Posterior Probabilities by Assigned Class"),
#       file = "Table_Model_Diagnostics_Posterior.tex")

# 2. Local Independence (Bivariate Residuals)
# --- V4 BVR CHANGE START: calculate expected counts from the fitted LCA model ---
bvr_calc = function(model, data) {
  data = as.data.frame(data)
  vars = colnames(data)
  nvar = length(vars)

  # Store the pairwise Pearson BVR statistic in the lower triangle.
  bvr_mat = matrix(NA_real_, nrow = nvar, ncol = nvar,
                   dimnames = list(vars, vars))

  for (i in 1:(nvar - 1)) {
    for (j in (i + 1):nvar) {
      var_i = vars[i]
      var_j = vars[j]

      # poLCA stores P(response category | class) for every indicator.
      prob_i = model$probs[[var_i]]
      prob_j = model$probs[[var_j]]

      # Use participants observed on both indicators for this pair.
      pair_complete = !is.na(data[[var_i]]) & !is.na(data[[var_j]])
      n_pair = sum(pair_complete)

      obs = table(
        factor(data[[var_i]][pair_complete], levels = seq_len(ncol(prob_i))),
        factor(data[[var_j]][pair_complete], levels = seq_len(ncol(prob_j)))
      )

      # Under local independence, within each class:
      # P(Y_i, Y_j | class) = P(Y_i | class) * P(Y_j | class).
      expected_prob = matrix(0, nrow = ncol(prob_i), ncol = ncol(prob_j))

      for (k in seq_along(model$P)) {
        expected_prob = expected_prob +
          model$P[k] * outer(prob_i[k, ], prob_j[k, ])
      }

      # Convert the model-implied joint probabilities to expected counts.
      expected_count = n_pair * expected_prob

      # Pearson bivariate residual statistic: larger values indicate that
      # the fitted LCA model does not reproduce this pair well.
      valid_cells = expected_count > 0
      bvr_mat[j, i] = sum(
        (obs[valid_cells] - expected_count[valid_cells])^2 /
          expected_count[valid_cells]
      )
    }
  }

  return(round(bvr_mat, 3))
}
# --- V4 BVR CHANGE END ---

bvr_res = bvr_calc(final_model, as.data.frame(df_lca)[, lca_vars])
write.csv(bvr_res, "Table_Local_Independence_BVR.csv")
# print(xtable(bvr_res, caption = "Bivariate Residuals for Local Independence Check"),
#       file = "Table_Local_Independence_BVR.tex")

# ==============================================================================
# Class-Specific Conditional Probabilities
# ==============================================================================
cond_probs = matrix(NA, nrow = length(lca_vars), ncol = selected_k)
rownames(cond_probs) = lca_vars
colnames(cond_probs) = paste0("Class ", 1:selected_k)

for (i in seq_along(lca_vars)) {
  var_name = lca_vars[i]
  cond_probs[i, ] = final_model$probs[[var_name]][, 2]
}

cond_probs_rounded = round(cond_probs, 3)

write.csv(cond_probs_rounded, "Table_Conditional_Probabilities.csv")
# print(xtable(cond_probs_rounded, caption = "Class-Specific Conditional Probabilities"),
#       file = "Table_Conditional_Probabilities.tex")

# ==============================================================================
# Latent class profiles (Heatmap & MCA Biplot)
# ==============================================================================
class_prop = round(prop.table(table(final_model$predclass)) * 100, 1)

base_names = c(
  "Relatively lower metabolic burden",
  "Lipid-related and multicomponent burden",
  "Blood pressure/glucose burden; low HPV vaccination"
)

heatmap_labels = paste0("Class ", 1:selected_k, " (", class_prop, "%)")
class_labels = paste0("Class ", 1:selected_k, " (", class_prop, "%)", " - ", base_names)

# --- Heatmap ---
heatmap_data = cond_probs_rounded
colnames(heatmap_data) = heatmap_labels

# Plot Heatmap
gg2 = ggcorrplot(heatmap_data,
                 method = "circle",
                 ggtheme = theme_bw(),
                 lab = TRUE,
                 lab_size = 3.5) +
  scale_fill_gradient2(limit = c(0, 1), low = "white", mid = "grey80", high = "firebrick3", midpoint = 0.5) +
  labs(title = "Figure 2: Latent Class Profiles", fill = "Pr(Yes)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
        axis.text.y = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("Figure2_Class_Profiles_Heatmap.png", gg2, width = 25, height = 12, units = "cm", dpi = 300)


# --- MCA Biplot ---
dfs = df_lca %>%
  dplyr::select(all_of(lca_vars)) %>%
  mutate_all(as.character)

dfs[dfs == "2"] = "yes"
dfs[dfs == "1"] = "no"
dfs = dfs %>% mutate_all(as.factor)

mcaQ = MCA(dfs, graph = FALSE, ncp = 2)
grp = as.factor(final_model$predclass)

# plot MCA Biplot
gg3 = fviz_mca_biplot(mcaQ,
                      repel = TRUE,
                      label = "var",
                      arrows = c(FALSE, FALSE),
                      col.var = "black",
                      habillage = grp,
                      addEllipses = TRUE,
                      title = "Figure 3: MCA Biplot of Latent Classes") +
  theme_bw() +
  scale_color_brewer(palette = "Set1", name = "LCA Classes", labels = class_labels) +
  scale_fill_brewer(palette = "Set1", name = "LCA Classes", labels = class_labels) +
  theme(legend.position = "bottom",
        legend.direction = "vertical",
        axis.text = element_text(face = "bold", size = 10),
        axis.title = element_text(size = 10, face = "bold"),
        legend.text = element_text(face = "bold", size = 10),
        legend.title = element_text(face = "bold", size = 10))
ggsave("Figure3_MCA_Biplot.png", gg3, width = 25, height = 20, units = "cm", dpi = 300)

print(gg2)
print(gg3)

# ==============================================================================
# Class Demographic Comparison
# ==============================================================================
# Calculate the sample size (N) and participation rate (%) for each group
class_counts = table(final_model$predclass)
class_prop = round(prop.table(class_counts) * 100, 1)

class_headers = paste0("Class ", 1:selected_k, " (", class_prop, "%)\n")

df$Classes_factor = factor(
  final_model$predclass,
  levels = 1:selected_k,
  labels = class_headers
)

tbl1_class = df %>%
  dplyr::select(
    age_group, gender_label, education_label, marital_label,
    high_depression_f, hpv_vaccinated_f,
    mets_waist_f, mets_bp_f, mets_glu_f, mets_hdl_f, mets_trig_f,
    Classes_factor
  ) %>%
  tbl_summary(
    by = Classes_factor,
    type = list(all_categorical() ~ "categorical"),
    label = list(
      age_group ~ "Age group",
      gender_label ~ "Gender",
      education_label ~ "Education level",
      marital_label ~ "Marital status",
      high_depression_f ~ "High depression symptoms",
      hpv_vaccinated_f ~ "HPV vaccinated",
      mets_waist_f ~ "Elevated waist circumference",
      mets_bp_f ~ "Elevated blood pressure",
      mets_glu_f ~ "Elevated fasting glucose",
      mets_hdl_f ~ "Low HDL cholesterol",
      mets_trig_f ~ "Elevated triglycerides"
    ),
    statistic = list(all_categorical() ~ "{p}%"),
    missing = "no"
  ) %>%
  bold_labels()

print(tbl1_class)
write.csv(as.data.frame(tbl1_class), "Table_Class_Demographics.csv", row.names = FALSE)
# print(xtable(as.data.frame(tbl1_class), caption = "Participant Demographic Characteristics Across Latent Classes"),
#       file = "Table_Class_Demographics.tex", include.rownames = FALSE)


# Keep fit objects for likelihood-replication and reproducibility checks.
saveRDS(list(models = models, final_model = final_model), "lca_phase1_fits.rds")
