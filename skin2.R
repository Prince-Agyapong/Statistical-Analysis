library(ucimlrepo) 
library(tidyverse)
library(janitor)
library(naniar)
library(missForest)
library(randomForest)
library(ranger)
library(ggplot2)
library(dplyr)
library(caret)
library(partykit)
library(ggparty)
library(patchwork)
library(pROC)
library(purrr)

set.seed(123)
data <- fetch_ucirepo(id = 33)

data <- data$data$original |>
  mutate(
    class = factor(class),
    class_names = recode_factor(
      class,
      `1` = "Psoriasis",
      `2` = "Seboreic Dermatitis",
      `3` = "Lichen Planus",
      `4` = "Pityriasis Rosea",
      `5` = "Cronic Dermatitis",
      `6` = "Pityriasis Rubra Pilaris"
    )
  )


gg_miss_var(data |>
              select(-class_names)) +
  labs(title = "Missing Values per Feature", x = "Feature", y = "Number Missing")

cls_cols <- c("#F8766D","#00BA38","#619CFF","#E76BF3","#FFD700","#00C2A0")
p <-  ggplot(data, aes(x = class_names, fill = class_names)) +
  geom_bar(show.legend = FALSE) +
  scale_fill_manual(values = cls_cols) +
  geom_text(stat='count', aes(label=..count..), 
            vjust=-0.5,
            size = 4.5,
            fontface = "bold") +
  labs(title = "Class Distribution", x = "Disease Class", y = "Count") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    axis.text = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 14, face = "bold")
  )





# Age could be a differentiating factor across diseases.
ggplot(data, aes(x = class_names, y = age, fill = class)) +
  geom_boxplot(show.legend = FALSE) +
  labs(title = "Age Distribution by Disease Class", x = "Disease Class", y = "Age") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    axis.text = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 14, face = "bold")
  )


# Why: Features with large variance across classes are more likely to help classification
# Calculate variance by class for each feature
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


# compute per-feature metrics
feat_names <- setdiff(names(data), c("class", "class_names"))

rank_tbl <- tibble(
  Feature = feat_names,
  eta2 = map_dbl(feat_names, ~ eta2_oneway(data[[.x]], data$class))
) |>
  arrange(desc(eta2))

# Top 10 by η²
top_eta <- rank_tbl |> filter(eta2 >= 0.7)

ggplot(top_eta,
       aes(x = reorder(Feature, eta2), y = eta2)) +
  geom_col(fill = "forestgreen") +
  coord_flip() +
  labs(
    title = "Features by Between-Class Variance (\u03B7\u00B2) of At least 0.7",
    subtitle = "Proportion of variance explained by disease class (0–1)",
    x = "Feature", 
    y = "\u03B7\u00B2"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    axis.text = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 14, face = "bold")
  )


# Data prep & imputation
mod_data <- data |>
  mutate(class = factor(class),
         age = as.numeric(age)) |>
  clean_names() |>
  select(-class_names)

idx <- createDataPartition(mod_data$class, p = 0.8, list = FALSE)
train_data <- mod_data[idx, ]
test_data  <- mod_data[-idx, ]

# Impute training
train_imp <- missForest(train_data)$ximp

# Impute test
tmp <- bind_rows(train_data |> select(-matches("class")),
                 test_data |> select(-matches("class")))
tmp_imp <- missForest(tmp)$ximp
test_imp <- cbind(tmp_imp[(nrow(train_data)+1):nrow(tmp_imp), ], class = test_data$class)

# Random Forest (Gini impurity)
rf_gini <- randomForest(
  class ~ ., data = train_imp,
  ntree = 500,
  mtry = floor(sqrt(ncol(train_imp) - 1)),
  importance = TRUE
)

imp_gini <- as.data.frame(randomForest::importance(rf_gini)) |>
  rownames_to_column("Feature") |>
  arrange(desc(MeanDecreaseGini)) |>
  mutate(Method = "Gini Impurity")

# Random Forest (Permutation importance)
rf_perm <- ranger(
  class ~ ., data = train_imp,
  num.trees = 500,
  mtry = floor(sqrt(ncol(train_imp) - 1)),
  importance = "permutation",
  probability = TRUE
)

imp_perm <- as.data.frame(rf_perm$variable.importance) |>
  rownames_to_column("Feature") |>
  rename(Importance = 2) |>
  arrange(desc(Importance)) |>
  mutate(Method = "Accuracy Reduction")

# Combine & plot
imp_gini <- imp_gini |>
  rename(Importance = MeanDecreaseGini)


# Accuracy reduction plot
p1 <- ggplot(imp_gini |> 
               slice(1:15) |>
               mutate(Feature = str_replace_all(Feature, "_", " ")), 
             aes(x = reorder(Feature, Importance), y = Importance, fill = Method)) +
  geom_col(show.legend = FALSE,
           color = "black", fill = "tomato") +
  coord_flip() +
  labs(title = "(a) Accuracy Reduction",
       x = "Feature", y = "Importance") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    axis.text = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 14, face = "bold")
  )

# Gini impurity plot
p2 <- ggplot(imp_perm |> 
               slice(1:15) |>
               mutate(Feature = str_replace_all(Feature, "_", " ")), 
             aes(x = reorder(Feature, Importance), y = Importance, fill = Method)) +
  geom_col(show.legend = FALSE,
           color = "black", fill = "forestgreen") +
  coord_flip() +
  labs(title = "(b) Gini Impurity",
       x = "Feature", y = "Importance") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    axis.text = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 14, face = "bold")
  )

p1
p2


# Borda -------------------------------------------------------------------
gini_ranked <- imp_gini |>
  mutate(rank_gini = rank(-Importance, ties.method = "average")) |>
  rename(Importance_Gini = Importance) |>
  select(Feature, Importance_Gini, rank_gini)

perm_ranked <- imp_perm |>
  mutate(rank_acc = rank(-Importance, ties.method = "average")) |>
  rename(Importance_Acc = Importance) |>
  select(Feature, Importance_Acc, rank_acc)

# Combine and compute Borda
borda_tbl <- gini_ranked |>
  inner_join(perm_ranked, by = "Feature") |>
  mutate(BordaScore = rank_gini + rank_acc) |>
  arrange(BordaScore, rank_gini, rank_acc) |>
  mutate(BordaRank = row_number())

# Plot: smallest (best) scores on top
borda_tbl |>
  arrange(BordaScore) |>
  slice(1:15) |>
  mutate(Feature = str_replace_all(Feature, "_", " ")) |>
  ggplot(aes(x = reorder(Feature, -BordaScore), y = BordaScore)) +
  geom_col(fill = "#3abeff", color = "black") +
  coord_flip() +
  labs(
    title = "Borda aggregation: permutation (accuracy) + Gini",
    subtitle = "Lower Borda score = better",
    x = "Feature",
    y = "Borda score"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    axis.text = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 14, face = "bold")
  )



# Evaluate RF with top-k features from Borda list
set.seed(123)

# Preconditions
stopifnot(
  "class" %in% names(train_imp), is.factor(train_imp$class),
  "class" %in% names(test_imp),  is.factor(test_imp$class),
  all(c("Feature","BordaScore") %in% names(borda_tbl))
)

# Filter Borda features to those in data
borda_tbl <- borda_tbl |>
  filter(Feature %in% setdiff(names(train_imp), "class")) |>
  arrange(BordaScore)

# --- Metric helpers ---
macro_f1 <- function(cm) {
  mean(sapply(rownames(cm), function(cl) {
    tp <- cm[cl, cl]
    fp <- sum(cm[, cl]) - tp
    fn <- sum(cm[cl, ]) - tp
    p  <- if ((tp+fp)==0) 0 else tp/(tp+fp)
    r  <- if ((tp+fn)==0) 0 else tp/(tp+fn)
    if ((p+r)==0) 0 else 2*p*r/(p+r)
  }), na.rm = TRUE)
}

macro_auc <- function(y, probs) {
  mean(sapply(levels(y), function(cl) {
    y_bin <- as.numeric(y == cl)
    as.numeric(roc(y_bin, probs[[cl]], quiet = TRUE)$auc)
  }), na.rm = TRUE)
}

# --- RF evaluation ---
fit_eval_rf <- function(train_df, test_df, feats, num.trees = 800) {
  fml <- reformulate(feats, "class")
  fit <- ranger(
    fml, data = train_df, num.trees = num.trees,
    mtry = max(1, floor(sqrt(length(feats)))), probability = TRUE,
    importance = "permutation", seed = 123
  )
  probs <- predict(fit, data = test_df)$predictions
  pred  <- factor(colnames(probs)[max.col(probs)], levels = levels(test_df$class))
  cm    <- table(test_df$class, pred)
  
  tibble(
    k        = length(feats),
    Accuracy = mean(pred == test_df$class),
    MacroF1  = macro_f1(cm),
    MacroAUC = macro_auc(test_df$class, as.data.frame(probs))
  )
}

# subsets 
ks <- 1:15
topk_list <- setNames(lapply(ks, \(k) head(borda_tbl$Feature, k)), paste0("Top", ks))

# evaluations 
perf_tbl <- imap_dfr(topk_list, ~ fit_eval_rf(train_imp, test_imp, .x) |> mutate(Subset = .y)) |>
  arrange(match(Subset, paste0("Top", ks)))

print(perf_tbl)


perf_tbl_long <- perf_tbl |>
  pivot_longer(
    cols = c(Accuracy, MacroF1, MacroAUC),
    names_to = "Metric",
    values_to = "Value"
  )

ggplot(perf_tbl_long, aes(x = k, y = Value, color = Metric)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  geom_vline(xintercept = c(6), linetype = "dashed", color = "grey50") +
  scale_x_continuous(breaks = 1:15) +
  labs(
    title = "RF Performance vs Number of Top Features",
    x = "Number of Top Features (k)",
    y = "Metric Value",
    color = "Metric"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    axis.text = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 14, face = "bold")
  )



borda_plot_df <- borda_tbl |>
  mutate(Feature = str_replace_all(Feature, "_", " ")) |>
  arrange(BordaScore) |>       
  mutate(
    Feature = factor(Feature, levels = Feature),  
    Top6 = row_number() <= 6
  )

ggplot(borda_plot_df, aes(x = BordaScore, y = Feature)) +
  geom_point(size = 2.6, color = "#3abeff") +
  geom_point(data = dplyr::filter(borda_plot_df, Top6),
             shape = 21, size = 5, stroke = 1.2, color = "red", fill = NA) +
  labs(
    title = "Aggregated Feature Rankings (Borda)",
    subtitle = "Lower Borda score = higher overall importance",
    x = "Borda score (sum of ranks)",
    y = "Feature"
  ) +
  expand_limits(x = 0) +
  theme_bw() +
  theme(
    axis.ticks.y = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    axis.text = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 14, face = "bold")
  )

# Select top 6 features from Borda scores
top_features <- borda_plot_df |>
  arrange(BordaScore) |>
  slice(1:6) |>
  pull(Feature)

# Keep only target + top 6 features
df_top6 <- data |>
  select(class, all_of(top_features))

# Tree
ctree_pk <- partykit::ctree(class ~ ., data = df_top6)

# Class mapping 
cls_levels <- levels(df_top6$class)
cls_names <- c(
  `1` = "Psoriasis",
  `2` = "Seboreic Dermatitis",
  `3` = "Lichen Planus",
  `4` = "Pityriasis Rosea",
  `5` = "Cronic Dermatitis",
  `6` = "Pityriasis Rubra Pilaris"
)[cls_levels]

names(cls_cols) <- cls_levels

p_text <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "p < 0.001", sprintf("p = %.3f", p)))

# Tree
p_tree  <-  ggparty(ctree_pk) +
  geom_edge() +
  geom_edge_label(size = 3) +
  geom_node_splitvar(
    ids = "inner",
    aes(label = paste0(splitvar, "\n", p_text(p.value))),
    fill = "lightblue", color = "black", size = 3
  ) +
  geom_node_plot(
    ids = "terminal",
    gglist = list(
      geom_bar(aes(x = factor(1), fill = class),
               stat = "count", position = "fill", width = 0.9),
      scale_fill_manual(values = cls_cols, labels = cls_names, name = "Class"),
      scale_x_discrete(labels = NULL, breaks = NULL),
      labs(x = NULL, y = NULL),
      theme_minimal(base_size = 9),
      theme(
        legend.position = "none",
        panel.grid = element_blank()
      )
    )
  ) +
  geom_node_label(
    aes(label = paste0("n = ", nodesize)),
    nudge_y = -0.04, size = 3,
    label.r = unit(0.15, "lines"),
    fill = "white"
  ) +
  theme_void(base_size = 12) +
  labs(title = "Pathways (Top 6 Features)",
       x = "",
       y = "") +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 14, face = "bold")
  )

leg_df <- data.frame(
  class = factor(cls_levels, levels = cls_levels),
  name  = unname(cls_names[cls_levels])  
)

legend_strip <-
  ggplot(leg_df, aes(x = class, y = 1, fill = class)) +
  geom_tile(height = 0.7, width = 0.9) +
  geom_text(aes(label = name), vjust = 0.5, hjust = 0.5, size = 3.3,
            color = "black", fontface = "bold") +
  scale_fill_manual(values = cls_cols, guide = "none") +
  scale_x_discrete(expand = expansion(mult = c(0.02, 0.02))) +
  scale_y_continuous(limits = c(0.5, 1.5)) +
  theme_void() 


final_plot <- p_tree / legend_strip + plot_layout(heights = c(1, 0.12))
final_plot
