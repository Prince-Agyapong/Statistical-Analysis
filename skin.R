library(ucimlrepo) 
library(tidyverse)
library(janitor)
data <- fetch_ucirepo(id = 33)

data <- data$data$original |>
  mutate(class = factor(class)) |>
  clean_names()




# EDA ---------------------------------------------------------------------

ggplot(data, aes(x = class, fill = class)) +
  geom_bar(show.legend = FALSE) +
  geom_text(stat='count', aes(label=..count..), vjust=-0.5) +
  labs(title = "Class Distribution", x = "Disease Class", y = "Count") +
  theme_bw()




gg_miss_var(data) +
  labs(title = "Missing Values per Feature", x = "Feature", y = "Number Missing")

# Age could be a differentiating factor across diseases.
ggplot(data, aes(x = class, y = age, fill = class)) +
  geom_boxplot(show.legend = FALSE) +
  labs(title = "Age Distribution by Disease Class", x = "Disease Class", y = "Age") +
  theme_minimal()

# Why: Features with large variance across classes are more likely to help classification
# Calculate variance by class for each feature
# --- helpers ---------------------------------------------------------------
eta2_oneway <- function(x, g) {
  # x: numeric vector (feature), g: factor (class)
  ok <- complete.cases(x, g)
  x <- x[ok]; g <- g[ok]
  if (length(unique(g)) < 2L || length(x) < 3L) return(NA_real_)
  fit <- aov(x ~ g)
  ss <- anova(fit)[["Sum Sq"]]
  ss_between <- ss[1]; ss_total <- sum(ss)
  if (ss_total == 0) return(NA_real_)
  as.numeric(ss_between / ss_total)  # η² in [0,1]
}

kruskal_eps2 <- function(x, g) {
  ok <- complete.cases(x, g)
  x <- x[ok]; g <- g[ok]
  k <- length(unique(g)); n <- length(x)
  if (k < 2L || n < 3L) return(NA_real_)
  H <- suppressWarnings(kruskal.test(x ~ g)$statistic)
  as.numeric((H - k + 1) / (n - k))  # ε² in [0,1]
}

# --- compute per-feature metrics ------------------------------------------
feat_names <- setdiff(names(data), "class")

rank_tbl <- tibble(
  Feature = feat_names,
  eta2 = map_dbl(feat_names, ~ eta2_oneway(data[[.x]], data$class)),
  eps2 = map_dbl(feat_names, ~ kruskal_eps2(data[[.x]], data$class))
) %>%
  arrange(desc(eta2))

# Top 10 by η²
top_eta <- rank_tbl %>% filter(eta2 >= 0.7)

ggplot(top_eta,
       aes(x = reorder(Feature, eta2), y = eta2)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(title = "Features by Between-Class Variance (η²) of At least 0.7",
       subtitle = "Proportion of variance explained by disease class (0–1)",
       x = "Feature", y = "η²") +
  theme_minimal()


set.seed(123)
library(missForest)
library(randomForest)
library(ranger)
library(ggplot2)
library(dplyr)
library(caret)

set.seed(123)

# -----------------------------------
# 1. Data prep & imputation
# -----------------------------------
data <- data %>%
  mutate(class = factor(class),
         age = as.numeric(age))

idx <- createDataPartition(data$class, p = 0.8, list = FALSE)
train_data <- data[idx, ]
test_data  <- data[-idx, ]

# Impute training
train_imp <- missForest(train_data)$ximp

# Impute test (keep same structure)
tmp <- bind_rows(train_data %>% select(-class),
                 test_data %>% select(-class))
tmp_imp <- missForest(tmp)$ximp
test_imp <- cbind(tmp_imp[(nrow(train_data)+1):nrow(tmp_imp), ], class = test_data$class)

# -----------------------------------
# 2. Random Forest (Gini impurity)
# -----------------------------------
rf_gini <- randomForest(
  class ~ ., data = train_imp,
  ntree = 500,
  mtry = floor(sqrt(ncol(train_imp) - 1)),
  importance = TRUE
)

imp_gini <- as.data.frame(randomForest::importance(rf_gini)) %>%
  rownames_to_column("Feature") %>%
  arrange(desc(MeanDecreaseGini)) %>%
  mutate(Method = "Gini Impurity")

# -----------------------------------
# 3. Random Forest (Permutation importance)
# -----------------------------------
rf_perm <- ranger(
  class ~ ., data = train_imp,
  num.trees = 500,
  mtry = floor(sqrt(ncol(train_imp) - 1)),
  importance = "permutation",
  probability = TRUE
)

imp_perm <- as.data.frame(rf_perm$variable.importance) %>%
  rownames_to_column("Feature") %>%
  rename(Importance = 2) %>%
  arrange(desc(Importance)) %>%
  mutate(Method = "Accuracy Reduction")

# -----------------------------------
# 4. Combine & plot
# -----------------------------------
# Standardize column name for plotting
imp_gini <- imp_gini %>%
  rename(Importance = MeanDecreaseGini)

# Top 15 per method
top10_gini <- imp_gini %>% slice(1:15)
top10_perm <- imp_perm %>% slice(1:10)

imp_all <- bind_rows(top10_perm, top10_gini)

# Split into two datasets
imp_acc_plot  <- top10_perm
imp_gini_plot <- top10_gini

# Accuracy reduction plot
p1 <- ggplot(imp_gini %>% slice(1:15), 
              aes(x = reorder(Feature, Importance), y = Importance, fill = Method)) +
  geom_col(show.legend = FALSE) +
  coord_flip() +
  labs(title = "(a) Accuracy Reduction",
       x = "Feature", y = "Importance") +
  theme_minimal()

# Gini impurity plot
p2 <- ggplot(imp_perm %>% 
               slice(1:15), aes(x = reorder(Feature, Importance), y = Importance, fill = Method)) +
  geom_col(show.legend = FALSE) +
  coord_flip() +
  labs(title = "(b) Gini Impurity",
       x = "Feature", y = "Importance") +
  theme_minimal()

p1
p2



# Borda -------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(ggplot2)
set.seed(123)

# ---------- helper: turn importance scores into descending ranks (1 = best)
rank_features <- function(scores_named) {
  s <- scores_named[is.finite(scores_named)]
  r <- rank(-s, ties.method = "average")
  out <- rep(NA_real_, length(scores_named))
  names(out) <- names(scores_named)
  out[names(r)] <- as.numeric(r)
  out
}

# ---------- compute the two RF rankings on THIS training set ----------
rank_on_train <- function(df) {
  stopifnot("class" %in% names(df))
  df$class <- droplevels(df$class)
  
  rf <- randomForest::randomForest(
    class ~ ., data = df,
    ntree = 500, importance = TRUE
  )
  imp <- randomForest::importance(rf)
  
  # two vectors of importances
  mda  <- imp[, "MeanDecreaseAccuracy"]; names(mda)  <- rownames(imp)
  gini <- imp[, "MeanDecreaseGini"];    names(gini) <- rownames(imp)
  
  # convert to ranks (1 = most important)
  ranks <- list(
    RF_Accuracy = rank_features(mda),
    RF_Gini     = rank_features(gini)
  )
  
  # long tibble: Method / Feature / Rank
  dplyr::bind_rows(lapply(names(ranks), function(m) {
    tibble::tibble(Method = m, Feature = names(ranks[[m]]), Rank = as.numeric(ranks[[m]]))
  }))
}

# ---------- run once on your imputed training set ----------
ranks_long <- rank_on_train(train_imp)

# ---------- Borda aggregation across *methods only* (no resampling) ----------
borda <- ranks_long %>%
  group_by(Feature) %>%
  summarise(
    BordaScore = sum(Rank, na.rm = TRUE),   
    Coverage   = sum(!is.na(Rank))
  ) %>%
  ungroup() %>%
  arrange(BordaScore) %>%
  mutate(BordaRank = rank(BordaScore, ties.method = "first")) %>%
  select(BordaRank, Feature, BordaScore, Coverage)

# inspect
print(borda, n = 20)

# top-k vector to use downstream (e.g., Top10)
top10_feats <- borda %>% slice_min(BordaScore, n = 10) %>% pull(Feature)

# ---------- plot: smallest (best) scores on top ----------
borda %>%
  slice_min(BordaScore, n = 15) %>%
  ggplot(aes(x = reorder(Feature, BordaScore), y = BordaScore)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(
    title = "Borda aggregation (RF MeanDecreaseAccuracy + RF MeanDecreaseGini)",
    subtitle = "Lower BordaScore = better (sum of ranks across the two methods)",
    x = "Feature", y = "Borda score"
  ) +
  theme_minimal()




############################################################################
# Evaluate RF with top‑k features from Borda list
############################################################################
library(dplyr)
library(ranger)
library(pROC)
library(ggplot2)
library(purrr)
set.seed(123)

# --- 0) Preconditions ------------------------------------------------------
# Assumes:
# - train_imp: data.frame with predictors + factor column `class`
# - test_imp : data.frame with predictors + factor column `class`
# - borda    : tibble with columns Feature, BordaScore (lower = better)

stopifnot("class" %in% names(train_imp), is.factor(train_imp$class))
stopifnot("class" %in% names(test_imp),  is.factor(test_imp$class))
stopifnot(all(c("Feature","BordaScore") %in% names(borda)))

# Keep only features that actually exist in the data
borda <- borda %>%
  filter(Feature %in% setdiff(names(train_imp), "class")) %>%
  arrange(BordaScore)

# --- 1) Helpers ------------------------------------------------------------
macro_f1_from_cm <- function(cm_tbl) {
  # cm_tbl: confusion matrix as table (Reference x Prediction)
  classes <- rownames(cm_tbl)
  f1s <- sapply(classes, function(cl) {
    tp <- cm_tbl[cl, cl]
    fp <- sum(cm_tbl[, cl]) - tp
    fn <- sum(cm_tbl[cl, ]) - tp
    precision <- if ((tp + fp) == 0) 0 else tp / (tp + fp)
    recall    <- if ((tp + fn) == 0) 0 else tp / (tp + fn)
    if ((precision + recall) == 0) 0 else 2 * precision * recall / (precision + recall)
  })
  mean(f1s, na.rm = TRUE)
}

macro_auc_ovr <- function(y_true_factor, prob_df) {
  # one-vs-rest macro AUC
  classes <- levels(y_true_factor)
  aucs <- map_dbl(classes, function(cl) {
    y_bin  <- as.numeric(y_true_factor == cl)
    roc(y_bin, prob_df[[cl]], quiet = TRUE)$auc |> as.numeric()
  })
  mean(aucs, na.rm = TRUE)
}

fit_eval_rf <- function(train_df, test_df, feat_vec, num.trees = 800) {
  # Build formula with selected features
  fml <- as.formula(paste("class ~", paste(feat_vec, collapse = " + ")))
  # mtry ~ sqrt(p)
  mtry_val <- max(1, floor(sqrt(length(feat_vec))))
  
  fit <- ranger(
    formula = fml,
    data = train_df,
    num.trees = num.trees,
    mtry = mtry_val,
    probability = TRUE,
    importance = "permutation",
    seed = 123
  )
  
  probs <- predict(fit, data = test_df)$predictions
  pred_class <- colnames(probs)[max.col(probs, ties.method = "first")]
  pred_class <- factor(pred_class, levels = levels(test_df$class))
  
  acc <- mean(pred_class == test_df$class)
  
  # Confusion matrix (Reference x Prediction)
  cm <- table(Reference = test_df$class, Prediction = pred_class)
  macro_f1 <- macro_f1_from_cm(cm)
  macro_auc <- macro_auc_ovr(test_df$class, as.data.frame(probs))
  
  list(
    model = fit,
    metrics = tibble(
      k = length(feat_vec),
      Accuracy = acc,
      MacroF1 = macro_f1,
      MacroAUC = macro_auc
    ),
    conf_mat = cm
  )
}

# --- 2) Define subsets -----------------------------------------------------
ks <- 1:10
topk_list <- lapply(ks, function(k) head(borda$Feature, k))
names(topk_list) <- paste0("Top", ks)

# All features (baseline)
all_feats <- setdiff(names(train_imp), "class")
topk_list[["All"]] <- all_feats

# --- 3) Run evaluations ----------------------------------------------------
results <- map(topk_list, ~ fit_eval_rf(train_imp, test_imp, .x))
perf_tbl <- bind_rows(lapply(names(results), function(nm) {
  results[[nm]]$metrics %>% mutate(Subset = nm)
}))

# Order rows
perf_tbl <- perf_tbl %>% arrange(match(Subset, c(paste0("Top", ks), "All")))
print(perf_tbl)


#----- Plot: Borda list, top two circled in red
borda <- borda %>%
  mutate(Feature = factor(Feature, levels = Feature),
         Top2 = row_number() <= 7)
p_borda <- ggplot(borda, aes(x = Feature, y = BordaScore)) +
  geom_point(size = 2) +
  coord_flip() +
  labs(title = "Aggregated Feature Rankings (Borda list)",
       subtitle = "Lower score = higher average rank",
       x = "Feature", y = "Borda count score (average rank)") +
  theme_minimal()

# add red circles around top 2
p_borda <- p_borda +
  geom_point(data = subset(borda, Top2),
             aes(x = Feature, y = BordaScore),
             shape = 21, size = 4, color = "red", fill = NA, stroke = 1.2) +
  theme_minimal()

p_borda 


# --- Gather confusion matrices from results ---
cm_list <- lapply(names(results), function(nm) {
  as.data.frame.matrix(results[[nm]]$conf_mat) %>%
    mutate(Reference = rownames(results[[nm]]$conf_mat),
           Subset = nm)
})
cm_df <- bind_rows(cm_list)

# Convert to tidy long format for ggplot
cm_long <- cm_df %>%
  pivot_longer(cols = -c(Reference, Subset),
               names_to = "Prediction",
               values_to = "Count") |>
  filter(Subset == "Top7")

# Ensure factor ordering matches numeric class labels
cm_long$Reference <- factor(cm_long$Reference, levels = levels(train_imp$class))
cm_long$Prediction <- factor(cm_long$Prediction, levels = levels(train_imp$class))


# --- Plot heatmaps for each subset ---
ggplot(cm_long, aes(x = Prediction, y = Reference, fill = Count)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Count), color = "black") +
  scale_fill_gradient(low = "white", high = "steelblue") +
  labs(title = "Confusion Matrices by Feature Subset",
       subtitle = "Top-7 features from Borda ranking vs All features",
       x = "Predicted Class", y = "Actual Class") +
  theme_bw()




# 1) Take the Borda top-10 as the anchor
borda_top7 <- borda %>%
  arrange(BordaScore) %>%
  slice(1:7)

# 2) Merge in Gini and Accuracy importances (keep NAs if not in that top-10)
top7_table <- borda_top7 %>%
  mutate(
    BordaRank = row_number(),
    BordaScore = round(BordaScore, 2)
  ) %>%
  left_join(borda |>
              pivot_wider(names_from = Method,
                          values_from = Rank),
            by = "Feature") |>
  inner_join(imp_perm |>
               select(Feature, Accuracy_decrease = Importance),
             by = "Feature") |>
  inner_join(imp_gini |>
               select(Feature, Gini_impurity = Importance),
             by = "Feature") |>
  select(BordaRank, Feature, BordaScore, RF_Accuracy, RF_Gini, Accuracy_decrease, Gini_impurity)

       # 3) Pretty table
top7_table %>%
  kable("html", caption = "Top 7 Features by Borda Score (with Gini & Accuracy importances)") %>%
  kable_styling(full_width = FALSE, position = "center")


