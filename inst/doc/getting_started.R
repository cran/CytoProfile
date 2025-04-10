## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  out.width = "100%",
  fig.show = "asis",
  fig.keep = "all"
)

## ----setup, echo=TRUE, warning=FALSE, message=FALSE, cache=TRUE---------------
# Loading all packages required
# Data manipulation and reshaping
library(dplyr)       # For data filtering, grouping, and summarising.
library(tidyr)       # For reshaping data (e.g., pivot_longer, pivot_wider).

# Plotting and visualization
library(ggplot2)     # For creating all the ggplot-based visualizations.
library(gridExtra)   # For arranging multiple plots on a single page.
library(ggrepel)     # For improved label placement in plots (e.g., volcano plots).
library(gplots)      # For heatmap.2, which is used to generate heatmaps.
library(plot3D)      # For creating 3D scatter plots in PCA and sPLS-DA analyses.
library(reshape2)    # For data transformation (e.g., melt) in cross-validation plots.

# Statistical analysis
library(mixOmics)    # For multivariate analyses (PCA, sPLS-DA, etc.).
library(e1071)     # For computing skewness and kurtosis.
library(pROC)        # For ROC curve generation in machine learning model evaluation.

# Machine learning
library(xgboost)     # For building XGBoost classification models.
library(randomForest) # For building Random Forest classification models.
library(caret)       # For cross-validation and other machine learning utilities.

# Package development and document rendering
library(knitr)       # For knitting RMarkdown files and setting chunk options.
library(devtools)    # For installing the development version of the package from GitHub.

# Load in the CytoProfile package
library(CytoProfile)

# Loading in data
data("ExampleData1")
data_df <- ExampleData1

## ----EDA1,  echo=TRUE, warning=FALSE, message=FALSE, cache=TRUE---------------
# Generating boxplots to check for outliers for raw values
cyt_bp(data_df[, -c(1:3)], 
        pdf_title = NULL)  
# Removing the first 3 columns to retain only continuous variables.

# Generating boxplots to check for outliers for log2 values
cyt_bp(data_df[, -c(1:3)], 
       pdf_title = NULL,
       scale = "log2")  


## ----EDA2, echo=TRUE, warning=FALSE, message=FALSE, cache=TRUE----------------
# Raw values for group-specific boxplots
cyt_bp2(data_df[, -c(3, 5:28)], 
        pdf_title = NULL, 
        scale = NULL)

# Log2-transformed group-specific boxplots
cyt_bp2(data_df[, -c(3, 5:28)], 
        pdf_title = NULL, 
        scale = "log2")



## ----EDA3, echo=TRUE, warning=FALSE, message=FALSE, cache=TRUE----------------
data_df <- ExampleData1
# Histogram for overall raw data
cyt_skku(data_df[, -c(1:3)], 
         pdf_title = NULL, 
         group_cols = NULL)

# Histogram with grouping (e.g., "Group")
cyt_skku(ExampleData1[, -c(2:3)], 
         pdf_title = NULL, 
         group_cols = c("Group"))


## ----EDA4,  echo=TRUE, warning=FALSE, message=FALSE, cache=TRUE---------------
# Generating basic error bar plots
data_df <- ExampleData1
cyt_errbp(data_df[,c("Group", "CCL.20.MIP.3A", "IL.10")], group_col = "Group", p_lab = FALSE, 
es_lab = FALSE, class_symbol = FALSE, x_lab = "Cytokines", y_lab = "Concentrations in log2 scale", log2 = TRUE)

## ----EDA5,  echo=TRUE, warning=FALSE, message=FALSE, cache=TRUE---------------
# Generating Error Bar Plot enriched with p-value and effect size 
data_df <- ExampleData1
cyt_errbp(data_df[,c("Group", "CCL.20.MIP.3A", "IL.10")], group_col = "Group", p_lab = TRUE, 
es_lab = TRUE, class_symbol = TRUE, x_lab = "Cytokines", y_lab = "Concentrations in log2 scale", log2 = TRUE)


## ----Univariate1,  echo=TRUE, warning=FALSE, message=FALSE, cache=TRUE--------
# Performing Two Sample T-test and Mann Whitney U Test
data_df <- ExampleData1[, -c(3)]
data_df <- dplyr::filter(data_df, Group != "ND", Treatment != "Unstimulated")
# Two sample T-test
cyt_ttest(data_df[, c(1:2, 5:6)], scale = "log2", verbose = TRUE, format_output = TRUE)
# Mann-Whitney U Test
cyt_ttest(data_df[, c(1:2, 5:6)], verbose = TRUE)

## ----Univariate2,  echo=TRUE, warning=FALSE, message=FALSE, cache=TRUE--------
# Perform ANOVA comparisons test (example with 2 cytokines)
data_df <- ExampleData1[, -c(3)]
cyt_anova(data_df[, c(1:2, 5:6)], format_output = TRUE)

## ----Multivariate1, echo=TRUE, warning=FALSE, message=FALSE, cache=TRUE, fig.show = "hold", out.width = "50%"----
# cyt_plsda function. 
data <- ExampleData1[, -c(3)]
data_df <- dplyr::filter(data, Group != "ND" & Treatment == "CD3/CD28")
cyt_splsda(data_df, 
           pdf_title = NULL, 
           colors = c("black", "purple"),
           bg = FALSE, scale = "log2", ellipse = TRUE,
           conf_mat = FALSE, var_num = 25, 
           cv_opt = "loocv", comp_num = 2, 
           pch_values = c(16, 4), group_col = "Group", group_col2 = "Treatment", 
           roc = TRUE)

## ----Multivariate2, echo=TRUE, warning=FALSE, message=FALSE, cache=TRUE, fig.show = "hold", out.width = "50%"----
data <- ExampleData1[, -c(3,23)]
data_df <- dplyr::filter(data, Group != "ND" & Treatment != "Unstimulated")
# cyt_pca(data_df, 
#         pdf_title = file.path(path, "example_pca_analysis.pdf"), 
#         colors = c("black", "red2"), scale = "log2", 
#         comp_num = 3, pch_values = c(16, 4), 
#         style = "3D", group_col = "Group", group_col2 = "Treatment")
cyt_pca(data_df, 
        pdf_title = NULL, 
        colors = c("black", "red2"), scale = "log2", 
        comp_num = 2, pch_values = c(16, 4), group_col = "Group")


## ----EDA6, echo=TRUE, warning=FALSE, message=FALSE, cache=TRUE----------------
# Generating Volcano Plot
data_df <- ExampleData1[, -c(2:3)]
cyt_volc(data_df, group_col = "Group", 
                      cond1 = "T2D", cond2 = "ND", 
                      fold_change_thresh = 2.0, 
                      top_labels = 15)


## ----EDA7, echo=TRUE, warning=FALSE, message=FALSE, cache=TRUE, fig.width=5, fig.height=4----
# Generating Heat map
cyt_heatmap(data = data_df,
            scale = "log2",        # Optional scaling
            annotation_col_name = "Group",
            title = NULL)

## ----EDA8, echo=TRUE, warning=FALSE, message=FALSE, cache=TRUE----------------
# Generating dual flashlights plot
data_df <- ExampleData1[, -c(2:3)]
dfp <- cyt_dualflashplot(data_df, group_var = "Group", 
                         group1 = "T2D", group2 = "ND", 
                         ssmd_thresh = -0.2, log2fc_thresh = 1, 
                         top_labels = 10)
# Print the plot
dfp
# Print the table data used for plotting
print(dfp$data, n=25)

## ----ML1, echo=TRUE, warning=FALSE, message=FALSE, cache=TRUE-----------------
# Using XGBoost for classification
data_df0 <- ExampleData1
data_df <- data.frame(data_df0[, 1:3], log2(data_df0[, -c(1:3)]))
data_df <- data_df[, -c(2:3)]
data_df <- dplyr::filter(data_df, Group != "ND")

cyt_xgb(data = data_df, group_col = "Group",
                       nrounds = 500, max_depth = 4, eta = 0.05,
                       nfold = 5, cv = TRUE, eval_metric = "mlogloss",
                       early_stopping_rounds = NULL, top_n_features = 10,
                       verbose = 0, plot_roc = TRUE, print_results = FALSE)

## ----ML2, echo=TRUE, warning=FALSE, message=FALSE, cache=TRUE-----------------
# Using Random Forest for classification
cyt_rf(data = data_df, group_col = "Group", k_folds = 5,
                     ntree = 1000, mtry = 4, run_rfcv = TRUE,
                     plot_roc = TRUE, verbose = FALSE)

