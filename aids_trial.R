library(ucimlrepo) 
library(tidyverse)
aids_data <- fetch_ucirepo(id = 890)

# Extract features and targets if structured like a list
X <- aids_data$data$original |>
  select(-pidnum)

X <- aids_data$data$features 
y <- aids_data$data$targets

# Load survival package
library(survival)

# Create survival object
surv_obj <- Surv(time = X$time, event = y$cid)

# Univariate screening (p < 0.20)
vars <- setdiff(names(X), "time")

univ_results <- lapply(vars, function(var) {
  formula <- as.formula(paste0("surv_obj ~ ", var))
  model <- try(coxph(formula, data = X), silent = TRUE)
  
  if (inherits(model, "try-error")) return(NULL)
  
  p_val <- try(summary(model)$coefficients[,"Pr(>|z|)"], silent = TRUE)
  
  if (inherits(p_val, "try-error") || is.na(p_val)) return(NULL)
  
  data.frame(var = var, p = p_val)
})

univ_table <- do.call(rbind, univ_results)
selected_vars <- univ_table$var[univ_table$p < 0.20]

# Drop NA if any slipped through
selected_vars <- selected_vars[!is.na(selected_vars)]

# Multivariable model (keep p < 0.10)
multi_formula <- as.formula(paste("surv_obj ~", paste(selected_vars, collapse = " + ")))
multi_model <- coxph(multi_formula, data = X)

# Extract significant variables
final_vars <- summary(multi_model)$coefficients
final_vars <- final_vars[final_vars[,"Pr(>|z|)"] < 0.10, , drop = FALSE]

# Final Cox model
final_formula <- as.formula(paste("surv_obj ~", paste(rownames(final_vars), collapse = " + ")))
final_model <- coxph(final_formula, data = X)
summary(final_model)

# Check proportional hazards assumption
cox.zph(final_model) 








library(survival)

X <- aids_data$data$original |>
  select(-pidnum)

# Define survival object
surv_obj <- Surv(time = X$time, event = X$cid)

# Full model with all covariates (excluding 'time')
full_vars <- setdiff(names(X), c("time", "cid"))
full_formula <- as.formula(paste("surv_obj ~", paste(full_vars, collapse = " + ")))

# Fit full Cox model
full_model <- coxph(full_formula, data = X)

# Backward elimination using stepwise AIC
back_model <- step(full_model, direction = "backward", trace = 0)

# Keep only variables with p < 0.05 in final model
final_summary <- summary(back_model)
signif_vars <- rownames(final_summary$coefficients)[final_summary$coefficients[,"Pr(>|z|)"] < 0.05]

# Refit final model with significant variables
final_formula <- as.formula(paste("surv_obj ~", paste(signif_vars, collapse = " + ")))
final_model <- coxph(final_formula, data = X, method = "breslow", ties = "breslow")

# Output summary
summary(final_model)

library(survminer)
library(survival)
# Create survival object
surv_obj <- Surv(time = X$time, event = X$cid)

# Fit Kaplan-Meier curve by treatment group
km_fit <- survfit(surv_obj ~ trt, data = X)

# Plot
ggsurvplot(
  km_fit,
  data = X,
  palette = "Set1",
  legend.title = "Treatment",
  title = "Kaplan-Meier Survival Curves by Treatment",
  xlab = "Time (days)",
  ylab = "Survival Probability"
)

# For Cox proportional hazards model
ggforest(final_model, data = X)




ggplot(X, aes(x = time)) +
  geom_histogram(bins = 30, fill = "#0072B2", color = "white") +
  scale_x_log10() +
  labs(title = "Distribution of Time to Event", x = "Time (days)", y = "Count")

ggplot(X, aes(x = factor(trt), fill = factor(cid))) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = c("gray70", "#D55E00"), labels = c("Censored", "Event")) +
  labs(title = "Event Proportion by Treatment", x = "Treatment Group", y = "Proportion", fill = "Status")




## CD4 plots
library(survival)
library(survminer)

dat0 <- X |>
  mutate(
    treatment_group = case_when(
      trt == 3 ~ "With ddI only",
      TRUE ~ "Without ddI only"
    ),
    cd4_group = case_when(
      cd40 <= 50 ~ "≤50",
      cd40 > 50 & cd40 <= 200 ~ "51–200",
      TRUE ~ "99"
    )
  )

# --- 1. All Patients ---
surv_all <- Surv(dat0$time, dat0$cid)
plot_all <- ggsurvplot(
  survfit(surv_all ~ treatment_group, data = dat0),
  data = dat0,
  conf.int = TRUE,
  risk.table = TRUE,
  title = "All Patients",
  legend.title = "Treatment",
  legend.labs = c("Without ddI only", "With ddI only")
)



# --- 3. Patients with CD4 51–200 ---
sub_200 <- dat0[dat0$cd4_group == "51–200", ]
surv_200 <- Surv(sub_200$time, sub_200$cid)
plot_200 <- ggsurvplot(
  survfit(surv_200 ~ treatment_group, data = sub_200),
  data = sub_200,
  conf.int = TRUE,
  risk.table = TRUE,
  title = "Patients with CD4 51–200",
  legend.title = "Treatment",
  legend.labs = c("Without ddI only", "With ddI only")
)

# --- Display or save the plots ---
# Print them one by one
print(plot_all)
print(plot_50)
print(plot_200)


## from paper
library(dplyr)
library(survival)
library(broom)


dat <- X |>  
  mutate(
    treatment = case_when(
      trt == 0 ~ "Zidovudine",
      trt == 1 ~ "Zidovudine + Didanosine",
      trt == 2 ~ "Zidovudine + Zalcitabine",
      trt == 3 ~ "Didanosine"
    ) |> factor(levels = c("Zidovudine", "Zidovudine + Didanosine", "Zidovudine + Zalcitabine", "Didanosine")),
    cd4_decline50 = if_else(cd420 < 0.5 * cd40, 1, 0)
  )


summarize_model <- function(df, outcome_var, endpoint_label, population_label) {
  # N and %
  n_tbl <- df |>
    group_by(treatment) |>
    summarise(
      N = n(),
      Events = sum(.data[[outcome_var]], na.rm = TRUE),
      Percent = round(100 * Events / N, 1),
      .groups = "drop"
    )  |>
    mutate(n_pct = sprintf("%.0f (%.1f%%)", N, Percent))
  
  # Cox model
  model <- coxph(as.formula(paste0("Surv(time, ", outcome_var, ") ~ treatment")), data = df)
  hr_tbl <- tidy(model, exponentiate = TRUE, conf.int = TRUE) |>
    mutate(treatment = gsub("^treatment", "", term)) |>
    select(treatment, HR = estimate, conf.low, conf.high, p.value) |>
    mutate(
      HR = sprintf("%.2f", HR),
      CI = sprintf("%.2f–%.2f", conf.low, conf.high),
      p.value = if_else(p.value < 0.001, "<0.001", sprintf("%.3f", p.value))
    )
  
  # Merge and format
  full_tbl <- n_tbl |>
    left_join(hr_tbl, by = "treatment") |>
    mutate(
      HR = if_else(is.na(HR), "Ref", HR),
      CI = if_else(is.na(CI), "", CI),
      p.value = if_else(is.na(p.value), "", p.value),
      Endpoint = endpoint_label,
      Population = population_label
    ) |>
    relocate(Population, Endpoint, treatment, N, Events, Percent, HR, CI, p.value) |>
    rename(Treatment = treatment) |>
    select(-c(matches("conf")))
  
  return(full_tbl)
}

results <- bind_rows(
  summarize_model(dat, "cd4_decline50", "CD4 decline", "All Patients"),
  summarize_model(dat, "cid", "Death", "All Patients")
)

results_wide <- results |>
  select(-c(N, Events)) |>
  pivot_longer(
    cols = c(n_pct, HR, CI, p.value),
    names_to = "STAT",
    values_to = "Value"
  ) |>
  pivot_wider(
    id_cols = c(Population, Endpoint, STAT),
    names_from = Treatment,
    values_from = Value
  ) |>
  relocate(Population, Endpoint, STAT) 



dat |> 
  ggplot(aes(x = cd40, fill = treatment)) +
  geom_density(alpha = 0.5) +
  labs(title = "Baseline CD4 Count by Treatment Group", x = "CD4 at Baseline", y = "Density") +
  theme_minimal()


dat |> 
  ggplot(aes(x = cd40, y = cd420, color = treatment)) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(title = "CD4 Week 20 vs Baseline", x = "CD4 at Baseline", y = "CD4 at Week 20") +
  theme_minimal()


dat |> 
  group_by(treatment) |> 
  summarise(
    pct_decline = mean(cd4_decline50 == 1, na.rm = TRUE) * 100
  ) |> 
  ggplot(aes(x = treatment, y = pct_decline, fill = treatment)) +
  geom_col() +
  labs(title = "Proportion with ≥50% CD4 Decline", y = "% of Patients") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")


dat |>
  group_by(treatment) |> 
  summarise(
    death = mean(cid == 1, na.rm = TRUE) * 100,
    composite = mean(composite_event == 1, na.rm = TRUE) * 100
  ) |> 
  pivot_longer(cols = c(death, composite), names_to = "event", values_to = "percent") |>
  ggplot(aes(x = treatment, y = percent, fill = event)) +
  geom_col(position = position_dodge()) +
  labs(title = "Event Rates by Treatment", y = "% of Patients", x = "Treatment") +
  scale_fill_manual(values = c("death" = "red", "composite" = "purple"),
                    labels = c("Death", "Composite (CD4↓ or Death)")) +
  theme_bw()


ggplot(dat, aes(x = cd40)) +
  geom_histogram(binwidth = 50, fill = "skyblue", color = "black") +
  labs(title = "Distribution of Baseline CD4 Counts", x = "CD4 at Baseline", y = "Count") +
  theme_bw()

ggplot(dat, aes(x = treatment, y = cd40, fill = treatment)) +
  geom_boxplot() +
  labs(title = "Baseline CD4 by Treatment Group", x = "Treatment", y = "CD4 at Baseline") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(dat, aes(x = cd40, y = cd420, color = treatment)) +
  geom_point(alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(title = "Week 20 CD4 vs Baseline", x = "Baseline CD4", y = "CD4 at Week 20") +
  theme_minimal()



######################################
plot_km <- function(df, outcome_var, endpoint_label, population_label) {
  surv_obj <- Surv(df$time, df[[outcome_var]])
  fit <- survfit(surv_obj ~ treatment, data = df)
  
  ggsurvplot(
    fit,
    data = df,
    risk.table = TRUE,
    conf.int = FALSE,
    pval = TRUE,
    legend.title = "Treatment",
    legend.labs = levels(df$treatment),
    palette = "Dark2",
    title = paste0("Kaplan–Meier: ", endpoint_label, " (", population_label, ")"),
    xlab = "Time (days)",
    ylab = "Survival Probability",
    ggtheme = theme_minimal()
  )
}

km_all_composite <- plot_km(dat, "cd4_decline50", "CD4 Decline", "All Patients")
km_all_death <- plot_km(dat, "cid", "Death", "All Patients")

# 2. Prior ART
km_prior_composite <- plot_km(dat |> filter(str2 == 1), "cd4_decline50", "CD4 Decline", "Prior ART")
km_prior_death <- plot_km(dat |> filter(str2 == 1), "cid", "Death", "Prior ART")

# 3. No Prior ART
km_noart_composite <- plot_km(dat |> filter(str2 == 0), "cd4_decline50", "CD4 Decline", "No Prior ART")
km_noart_death <- plot_km(dat |> filter(str2 == 0), "cid", "Death", "No Prior ART")
