# ==============================================================================
# Electrochemical Impedance Spectroscopy (EIS) - Equivalent Circuit Analysis
# 3-Way ANOVA, Outlier Diagnostics (Cook's Distance), and Orthogonal Contrasts
# ==============================================================================

library(dplyr)
library(emmeans)
library(stats)

# ------------------------------------------------------------------------------
# 1. Data Import and Preparation
# ------------------------------------------------------------------------------
data_path <- "./data/EIS_Fitted_Parameters.csv"
eis_data <- read.csv(data_path, stringsAsFactors = TRUE)

# Translate column names to English for international consistency
colnames(eis_data)[colnames(eis_data) == "Lavagem"] <- "Washing"
colnames(eis_data)[colnames(eis_data) == "Solucao"] <- "Solution"

# Ensure all parameters are numeric
eis_data$R1 <- as.numeric(as.character(eis_data$R1))

# Define the Equivalent Circuit Parameters to analyze
# C1: Ca, C2: Cdl, R1: Rs, R2: Ra, R3: Rct
parameters <- c("C1", "C2", "R2", "R3", "R1") 
param_labels <- c("Ca (Adsorption Capacitance)", 
                  "Cdl (Double-Layer Capacitance)", 
                  "Ra (Adsorption Resistance)", 
                  "Rct (Charge Transfer Resistance)", 
                  "Rs (Solution Resistance)")

# ------------------------------------------------------------------------------
# 2. Define Custom Orthogonal Contrasts
# ------------------------------------------------------------------------------
# NOTE: The length of these vectors (12) corresponds to the specific combinations 
# of the 3-way interaction (Washing * Plasma * Solution).
custom_contrasts <- list(
  "Au-P"   = c(0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0),
  "Au-SP"  = c(0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0),
  "PEG-P"  = c(0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0),
  "PEG-SP" = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1),
  "Ab-P"   = c(1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
  "Ab-SP"  = c(0, 0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0)
)

# ------------------------------------------------------------------------------
# 3. Core Analytical Function
# ------------------------------------------------------------------------------
analyze_parameter <- function(data, param_name, label) {
  
  cat(sprintf("\n======================================================\n"))
  cat(sprintf(" ANALYZING PARAMETER: %s (%s)\n", param_name, label))
  cat(sprintf("======================================================\n"))
  
  # Construct and fit the 3-Way ANOVA model dynamically
  formula_str <- paste(param_name, "~ Plasma * Washing * Solution")
  model <- aov(as.formula(formula_str), data = data)
  
  # --- 3.1 Outlier Diagnostics (Cook's Distance) ---
  cooksD <- cooks.distance(model)
  cat("\n--- Top 5 Cook's Distances (Potential Outliers) ---\n")
  print(head(sort(cooksD, decreasing = TRUE), 5))
  
  # Note for reproducibility: In standard workflows, outliers identified here 
  # (e.g., Cook's D > 4/n) should be investigated and potentially filtered 
  # before finalizing the emmeans contrasts.
  
  # --- 3.2 Estimated Marginal Means (emmeans) ---
  em <- emmeans(model, ~ Washing * Plasma * Solution)
  
  # --- 3.3 Apply Orthogonal Contrasts with Holm Correction ---
  contrast_result <- contrast(em, custom_contrasts)
  contrast_summary <- summary(contrast_result)
  
  # Apply Holm-Bonferroni correction for multiple comparisons
  p_values <- contrast_summary$p.value
  contrast_summary$adjusted_p <- p.adjust(p_values, method = "holm")
  
  cat("\n--- Orthogonal Contrasts Results ---\n")
  print(contrast_summary)
  
  return(contrast_summary)
}

# ------------------------------------------------------------------------------
# 4. Execute Pipeline
# ------------------------------------------------------------------------------
# Run the analysis for all parameters automatically
results_list <- list()

for (i in seq_along(parameters)) {
  results_list[[parameters[i]]] <- analyze_parameter(eis_data, parameters[i], param_labels[i])
}

cat("\n\nPipeline execution complete. Review Cook's distances above for data integrity.\n")