## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  echo = TRUE,
  message = FALSE,
  warning = FALSE,
  fig.align = "center",
  fig.retina = 2,
  dpi = 96,
  out.width = "100%",
  fig.width = 8,
  fig.height = 5
)

## ----setup--------------------------------------------------------------------
library(CytoProfile)
library(dplyr)   # used only for filter() in the examples below

## ----load-data----------------------------------------------------------------
data("ExampleData1")
data_df <- ExampleData1

# Inspect structure
dim(data_df)
names(data_df)[1:6]
head(data_df[, 1:5])
table(data_df$Group)
table(data_df$Treatment)

## ----boxplots, fig.height=5---------------------------------------------------
# Ungrouped boxplots with log2 transformation
# All numeric columns are plotted; up to bin_size = 25 per page
cyt_bp(data_df[, -c(1:3)], output_file = NULL, scale = "log2")

## ----boxplots-grouped, fig.height=4-------------------------------------------
# Grouped boxplots: one plot per cytokine, colored by Group
# Only passing Group and two cytokines for a concise display
cyt_bp(
  data_df[, c("Group", "IL-10", "CCL-20/MIP-3A")],
  group_by = "Group",
  scale = "zscore"
)

## ----violins, fig.height=5----------------------------------------------------
# Ungrouped violin plots with z-score scaling
cyt_violin(data_df[, -c(1:3)], output_file = NULL, scale = "zscore")

## ----violins-grouped, fig.height=4--------------------------------------------
# Grouped violin plots with boxplot overlay and log2 scaling
cyt_violin(
  data_df[, c("Group", "IL-10", "CCL-20/MIP-3A")],
  group_by = "Group",
  boxplot_overlay = TRUE,
  scale = "log2"
)

## ----skku, fig.height=6-------------------------------------------------------
# Overall distributional assessment (no grouping)
cyt_skku(data_df[, -c(1:3)], output_file = NULL, group_cols = NULL)

## ----skku-grouped, fig.height=6-----------------------------------------------
# Grouped assessment by "Group"
cyt_skku(data_df[, -c(2:3)], output_file = NULL, group_cols = c("Group"))

## ----errbp, fig.height=6, fig.width=9-----------------------------------------
# Basic error bar plot: default mean +/- SE
df_err <- ExampleData1[, c("Group", "CCL-20/MIP-3A", "IL-10")]
cyt_errbp(
  df_err,
  group_col = "Group",
  x_lab = "Group",
  y_lab = "Concentration"
)

## ----errbp-symbols, fig.height=6, fig.width=9---------------------------------
# Mean +/- SD with log2 transformation, symbols, and t-test
cyt_errbp(
  df_err,
  group_col = "Group",
  stat = "mean",
  error = "sd",
  scale = "log2",
  class_symbol = TRUE,
  method = "ttest",
  x_lab = "Group",
  y_lab = "Concentration (log2)"
)

## ----univariate---------------------------------------------------------------
data_uni <- ExampleData1[, -c(3)]
data_uni <- dplyr::filter(data_uni, Group != "ND", Treatment != "Unstimulated")

# Tidy output with log2 transformation and automatic test selection
cyt_univariate(
  data_uni[, c("Group", "Treatment", "IL-10", "CCL-20/MIP-3A")],
  scale    = "log2",
  method   = "auto",
  format_output = TRUE,
  p_adjust_method = "BH"
)

## ----univariate-multi---------------------------------------------------------
# Kruskal-Wallis with pairwise Wilcoxon for multi-group comparison
cyt_univariate_multi(
  ExampleData1[, c("Group", "IL-10", "CCL-20/MIP-3A")],
  method = "kruskal",
  format_output = TRUE
)

## ----volcano, fig.height=5, fig.width=7---------------------------------------
data_volc <- ExampleData1[, -c(2:3)]

volc_plots <- cyt_volc(
  data_volc,
  group_col = "Group",
  cond1 = "T2D",
  cond2 = "ND",
  fold_change_thresh = 2.0,
  p_value_thresh = 0.05,
  top_labels = 15,
  method = "ttest"
)

# The function returns a named list of ggplot objects; print the comparison
volc_plots[["T2D vs ND"]]

## ----heatmap, fig.height=6----------------------------------------------------
cyt_heatmap(
  data = data_df[, -c(2:3)],
  scale = "log2",
  annotation_col = "Group",
  annotation_side = "auto",
  show_row_names = FALSE,
  title = NULL
)

## ----dualflash, fig.height=6, fig.width=7-------------------------------------
data_dfp <- ExampleData1[, -c(2:3)]

dfp <- cyt_dualflashplot(
  data_dfp,
  group_var    = "Group",
  group1       = "T2D",
  group2       = "ND",
  ssmd_thresh  = 0.2,
  log2fc_thresh = 1,
  top_labels   = 10,
  verbose      = FALSE
)
dfp

## ----dualflash-table----------------------------------------------------------
# Inspect the underlying statistics
head(dfp$data[order(abs(dfp$data$ssmd), decreasing = TRUE), ], 10)

## ----pca, fig.show='hide'-----------------------------------------------------
data_pca <- ExampleData1[, -c(3, 23)]
data_pca <- dplyr::filter(data_pca, Group != "ND" & Treatment != "Unstimulated")

pca_results <- cyt_pca(
  data_pca,
  output_file  = NULL,
  colors       = c("black", "red2"),
  scale        = "log2",
  comp_num     = 2,
  pch_values   = c(16, 4),
  group_col    = "Group",
  ellipse      = TRUE
)

## ----pca-plots, fig.show="hold", fig.width=9, fig.height=5--------------------
pca_results$overall_indiv_plot
pca_results$scree_plot

## ----pca-loadings, fig.show="hold", fig.width=9, fig.height=5-----------------
pca_results$loadings$Comp1()
pca_results$correlation_circle()

## ----splsda, fig.show='hide'--------------------------------------------------
data_spls <- ExampleData1[, -c(3)]
data_spls <- dplyr::filter(data_spls, Group != "ND" & Treatment == "CD3/CD28")

spls_results <- cyt_splsda(
  data_spls,
  output_file  = NULL,
  colors       = c("black", "purple"),
  scale        = "log2",
  ellipse      = TRUE,
  var_num      = 25,
  cv_opt       = "loocv",
  comp_num     = 2,
  pch_values   = c(16, 4),
  group_col    = "Group",
  group_col2   = "Treatment",
  roc          = FALSE,
  verbose      = FALSE
)

## ----splsda-plots, fig.show="hold", fig.width=9, fig.height=5-----------------
spls_results$overall_indiv_plot
spls_results$vip_indiv_plot

## ----splsda-loadings, fig.show="hold", fig.width=9, fig.height=5--------------
spls_results$loadings$Comp1()
spls_results$vip_scores$Comp1()

## ----mint, fig.show='hide'----------------------------------------------------
data_mint <- ExampleData5[, -c(2, 4)]
data_mint <- dplyr::filter(data_mint, Group != "ND")

mint_results <- cyt_mint_splsda(
  data_mint,
  group_col  = "Group",
  batch_col  = "Batch",
  colors     = c("black", "purple"),
  ellipse    = TRUE,
  var_num    = 25,
  comp_num   = 2,
  scale      = "log2",
  verbose    = FALSE
)

## ----mint-plots, fig.show="hold", fig.width=9, fig.height=5-------------------
mint_results$global_indiv_plot
mint_results$partial_indiv_plot

## ----mint-loadings, fig.show="hold", fig.width=9, fig.height=5----------------
mint_results$global_loadings_plots$Comp1()

## ----xgb, fig.height=5--------------------------------------------------------
data_ml <- data.frame(ExampleData1[, 1:3], log2(ExampleData1[, -c(1:3)]))
data_ml <- data_ml[, -c(2:3)]
data_ml <- dplyr::filter(data_ml, Group != "ND")

cyt_xgb(
  data          = data_ml,
  group_col     = "Group",
  nrounds       = 250,
  max_depth     = 4,
  learning_rate = 0.05,
  nfold         = 5,
  cv            = FALSE,
  objective     = "multi:softprob",
  eval_metric   = "mlogloss",
  top_n_features = 10,
  verbose       = 0,
  plot_roc      = FALSE,
  print_results = FALSE,
  seed          = 123
)

## ----rf, fig.height=5---------------------------------------------------------
cyt_rf(
  data       = data_ml,
  group_col  = "Group",
  ntree      = 500,
  mtry       = 4,
  k_folds    = 5,
  run_rfcv   = TRUE,
  plot_roc   = FALSE,
  verbose    = FALSE,
  seed       = 123
)

## ----session-info-------------------------------------------------------------
sessionInfo()

