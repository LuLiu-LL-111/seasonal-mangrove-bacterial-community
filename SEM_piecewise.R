
rm(list = ls())
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(piecewiseSEM)
  library(car)
  library(lmtest)
})

setwd("your_path")

FILE_ICAMP <- "icamp.txt"
FILE_REP   <- "rep.txt"
FILE_ASV   <- "asv_table.txt"
FILE_ENV   <- "env.txt"
OUTPUT_DIR <- "SEM_output"

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

required_files <- c(FILE_ICAMP, FILE_REP, FILE_ASV, FILE_ENV)
if (!all(file.exists(required_files))) {
  stop(
    paste0(
      "The following required files do not exist: ",
      paste(required_files[!file.exists(required_files)], collapse = ", ")
    )
  )
}

cat("Current working directory: ", getwd(), "\n")

## 2. Utility functions
clean_sample_id <- function(x) {
  x <- toupper(trimws(as.character(x)))
  gsub("[^A-Z0-9]", "", x)
}

make_pair_key <- function(x, y) {
  x <- clean_sample_id(x)
  y <- clean_sample_id(y)
  paste(pmin(x, y), pmax(x, y), sep = "__")
}

as_num_checked <- function(x, var_name) {
  old <- x
  out <- suppressWarnings(as.numeric(as.character(x)))
  bad <- is.na(out) & !is.na(old) & trimws(as.character(old)) != ""
  if (any(bad)) {
    stop(
      paste0(
        var_name,
        " contains values that cannot be converted to numeric, for example: ",
        paste(head(unique(as.character(old[bad])), 10), collapse = ", ")
      )
    )
  }
  out
}

zscore <- function(x) {
  x <- as.numeric(x)
  sx <- sd(x, na.rm = TRUE)
  if (!is.finite(sx) || sx == 0) stop("A variable has a standard deviation of 0 and cannot be standardized.")
  as.numeric(scale(x))
}

safe_mean <- function(x) {
  if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
}

safe_median <- function(x) {
  if (all(is.na(x))) NA_real_ else median(x, na.rm = TRUE)
}

check_columns <- function(dat, needed, file_name) {
  miss <- setdiff(needed, colnames(dat))
  if (length(miss) > 0) {
    stop(
      paste0(file_name, " is missing required columns: ", paste(miss, collapse = ", "))
    )
  }
}

repair_names <- function(x) {
  x <- as.data.frame(x, check.names = FALSE, stringsAsFactors = FALSE)
  nm <- names(x)
  if (length(nm) == 0 && ncol(x) > 0) nm <- rep("", ncol(x))
  bad <- is.na(nm) | trimws(nm) == ""
  if (any(bad)) nm[bad] <- paste0("Unnamed_", seq_len(sum(bad)))
  names(x) <- make.unique(nm, sep = "_")
  x
}

## 3. Read iCAMP data
icamp_raw <- read.delim(
  FILE_ICAMP,
  sep = "\t",
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

check_columns(
  icamp_raw,
  c("Var1", "Var2", "Homogeneous.Selection"),
  FILE_ICAMP
)

icamp_pair <- icamp_raw %>%
  transmute(
    Sample1 = pmin(clean_sample_id(Var1), clean_sample_id(Var2)),
    Sample2 = pmax(clean_sample_id(Var1), clean_sample_id(Var2)),
    pair_key = make_pair_key(Var1, Var2),
    HoS = as_num_checked(Homogeneous.Selection, "Homogeneous.Selection")
  ) %>%
  filter(Sample1 != "", Sample2 != "", Sample1 != Sample2) %>%
  group_by(pair_key, Sample1, Sample2) %>%
  summarise(HoS = safe_mean(HoS), .groups = "drop")

cat("Number of iCAMP sample pairs: ", nrow(icamp_pair), "\n")

## 4. Read richness-difference data
rep_raw <- read.delim(
  FILE_REP,
  sep = "\t",
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

check_columns(
  rep_raw,
  c("Var1", "Var2", "replacement", "richness"),
  FILE_REP
)

rep_pair <- rep_raw %>%
  transmute(
    Sample1 = pmin(clean_sample_id(Var1), clean_sample_id(Var2)),
    Sample2 = pmax(clean_sample_id(Var1), clean_sample_id(Var2)),
    pair_key = make_pair_key(Var1, Var2),
    Replacement = as_num_checked(replacement, "replacement"),
    RichnessDifference = as_num_checked(richness, "richness")
  ) %>%
  filter(Sample1 != "", Sample2 != "", Sample1 != Sample2) %>%
  group_by(pair_key, Sample1, Sample2) %>%
  summarise(
    Replacement = safe_mean(Replacement),
    RichnessDifference = safe_mean(RichnessDifference),
    .groups = "drop"
  )

cat("Number of sample pairs in rep.txt: ", nrow(rep_pair), "\n")

## 5. Match the two pairwise files
pair_data <- icamp_pair %>%
  inner_join(
    rep_pair %>% select(pair_key, Replacement, RichnessDifference),
    by = "pair_key"
  ) %>%
  drop_na(HoS, RichnessDifference)

if (nrow(pair_data) == 0) stop("No matching sample pairs were found between icamp.txt and rep.txt.")

cat("Number of successfully matched sample pairs: ", nrow(pair_data), "\n")

write.csv(
  anti_join(icamp_pair, rep_pair, by = "pair_key"),
  file.path(OUTPUT_DIR, "unmatched_pairs_in_icamp.csv"),
  row.names = FALSE
)

write.csv(
  anti_join(rep_pair, icamp_pair, by = "pair_key"),
  file.path(OUTPUT_DIR, "unmatched_pairs_in_rep.csv"),
  row.names = FALSE
)

## 6. Read environmental data
env_raw <- read.delim(
  FILE_ENV,
  sep = "\t",
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

check_columns(
  env_raw,
  c("SampleID", "pH", "Salinity", "season", "site"),
  FILE_ENV
)

env <- env_raw %>%
  transmute(
    Sample = clean_sample_id(SampleID),
    pH = as_num_checked(pH, "pH"),
    Salinity = as_num_checked(Salinity, "Salinity"),
    season = tolower(trimws(as.character(season))),
    site = toupper(trimws(as.character(site)))
  ) %>%
  filter(Sample != "", season != "", site != "")

if (anyDuplicated(env$Sample) > 0) {
  duplicated_env <- env %>% count(Sample) %>% filter(n > 1)
  write.csv(
    duplicated_env,
    file.path(OUTPUT_DIR, "duplicated_environment_samples.csv"),
    row.names = FALSE
  )
  stop("Duplicate sample names were found in env.txt.")
}

cat("Number of samples in the environmental table: ", nrow(env), "\n")
cat("Sites: ", paste(sort(unique(env$site)), collapse = ", "), "\n")
cat("Seasons: ", paste(sort(unique(env$season)), collapse = ", "), "\n")

## 7. Calculate observed richness from the ASV table
ASV_raw <- read.delim(
  FILE_ASV,
  sep = "\t",
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

if (ncol(ASV_raw) < 2) stop("Invalid ASV table format: at least two columns are required.")

ASV_sample_names <- colnames(ASV_raw)[-1]

## Avoid using [[ ]] to prevent the previous bracket-related error
ASV_numeric_list <- Map(
  function(x, sample_name) {
    as_num_checked(x, paste0("ASV sample column: ", sample_name))
  },
  ASV_raw[, -1, drop = FALSE],
  ASV_sample_names
)

ASV_matrix <- do.call(cbind, ASV_numeric_list)
colnames(ASV_matrix) <- clean_sample_id(ASV_sample_names)

if (anyDuplicated(colnames(ASV_matrix)) > 0) {
  dup <- unique(colnames(ASV_matrix)[duplicated(colnames(ASV_matrix))])
  stop(
    paste0("Duplicate ASV sample names after cleaning: ", paste(dup, collapse = ", "))
  )
}

ASV_richness <- data.frame(
  Sample = colnames(ASV_matrix),
  Richness = colSums(ASV_matrix > 0, na.rm = TRUE),
  stringsAsFactors = FALSE
)

cat("Number of samples in the ASV table: ", nrow(ASV_richness), "\n")

## 8. Filter valid samples
valid_samples <- intersect(env$Sample, ASV_richness$Sample)
if (length(valid_samples) == 0) stop("No shared samples were found between env.txt and ASVr.txt.")

pair_data_valid <- pair_data %>%
  filter(Sample1 %in% valid_samples, Sample2 %in% valid_samples)

if (nrow(pair_data_valid) == 0) stop("No valid sample pairs were found.")

cat("Number of valid samples: ", length(valid_samples), "\n")
cat("Number of valid sample pairs: ", nrow(pair_data_valid), "\n")

write.csv(
  data.frame(
    MissingSample = setdiff(
      unique(c(pair_data$Sample1, pair_data$Sample2)),
      valid_samples
    )
  ),
  file.path(OUTPUT_DIR, "pairwise_samples_missing_env_or_ASV.csv"),
  row.names = FALSE
)

## 9. Add site and season metadata
metadata_1 <- env %>%
  select(Sample, site, season) %>%
  rename(Sample1 = Sample, site1 = site, season1 = season)

metadata_2 <- env %>%
  select(Sample, site, season) %>%
  rename(Sample2 = Sample, site2 = site, season2 = season)

pair_with_meta <- pair_data_valid %>%
  inner_join(metadata_1, by = "Sample1") %>%
  inner_join(metadata_2, by = "Sample2")

if (nrow(pair_with_meta) == 0) stop("No data remained after matching sample pairs with the environmental table.")

write.csv(
  pair_with_meta,
  file.path(OUTPUT_DIR, "matched_pairwise_data_with_metadata.csv"),
  row.names = FALSE
)

## 10. Aggregate pairwise metrics to the sample level
aggregate_pair_metrics <- function(pair_df, method = c("mean", "median")) {
  method <- match.arg(method)
  if (nrow(pair_df) == 0) stop("The pairwise dataset is empty.")

  pair_long <- bind_rows(
    pair_df %>%
      transmute(
        Sample = Sample1,
        Partner = Sample2,
        HoS = HoS,
        RichnessDifference = RichnessDifference
      ),
    pair_df %>%
      transmute(
        Sample = Sample2,
        Partner = Sample1,
        HoS = HoS,
        RichnessDifference = RichnessDifference
      )
  ) %>%
    drop_na(HoS, RichnessDifference)

  summary_fun <- if (method == "mean") safe_mean else safe_median

  pair_long %>%
    group_by(Sample) %>%
    summarise(
      HoS_sample = summary_fun(HoS),
      RichDiff_sample = summary_fun(RichnessDifference),
      HoS_SD = sd(HoS, na.rm = TRUE),
      RichDiff_SD = sd(RichnessDifference, na.rm = TRUE),
      n_pairs = n(),
      n_unique_partners = n_distinct(Partner),
      .groups = "drop"
    )
}

## 11. Build the SEM dataset
build_sem_data <- function(pair_summary) {
  dat <- env %>%
    select(Sample, pH, Salinity, season, site) %>%
    inner_join(ASV_richness, by = "Sample") %>%
    inner_join(pair_summary, by = "Sample") %>%
    drop_na(
      pH,
      Salinity,
      season,
      site,
      Richness,
      HoS_sample,
      RichDiff_sample
    ) %>%
    filter(Richness > 0) %>%
    mutate(
      season = factor(season),
      site = factor(site),
      log_Richness = log(Richness)
    )

  if (nrow(dat) < 10) stop(paste0("SEM sample size is too small: ", nrow(dat)))
  if (nlevels(dat$season) < 2) stop("The season factor has fewer than two levels.")
  if (nlevels(dat$site) < 2) stop("The site factor has fewer than two levels.")

  dat %>%
    mutate(
      pH_z = zscore(pH),
      Salinity_z = zscore(Salinity),
      Richness_z = zscore(log_Richness),
      HoS_z = zscore(HoS_sample),
      RichDiff_z = zscore(RichDiff_sample)
    )
}

## 12. Fit the piecewise SEM
fit_piecewise_sem <- function(dat) {
  dat <- droplevels(dat)

  model_pH <- lm(
    pH_z ~ season + Salinity_z + site,
    data = dat
  )

  model_richness <- lm(
    Richness_z ~ pH_z + season + Salinity_z + site,
    data = dat
  )

  model_HoS <- lm(
    HoS_z ~ pH_z + Richness_z + Salinity_z + season + site,
    data = dat
  )

  model_RichDiff <- lm(
    RichDiff_z ~ pH_z + Richness_z + HoS_z + season + site,
    data = dat
  )

  sem_model <- piecewiseSEM::psem(
    model_pH,
    model_richness,
    model_HoS,
    model_RichDiff,
    data = dat
  )

  list(
    sem = sem_model,
    model_pH = model_pH,
    model_richness = model_richness,
    model_HoS = model_HoS,
    model_RichDiff = model_RichDiff
  )
}

## 13. Save SEM results
save_sem_results <- function(fit, analysis_name) {
  coefficients <- repair_names(
    piecewiseSEM::coefs(fit$sem, standardize = "scale")
  )
  coefficients$Analysis <- analysis_name

  r_squared <- repair_names(piecewiseSEM::rsquared(fit$sem))
  r_squared$Analysis <- analysis_name

  fisher_c <- tryCatch(
    {
      out <- repair_names(piecewiseSEM::fisherC(fit$sem))
      out$Analysis <- analysis_name
      out
    },
    error = function(e) {
      data.frame(
        Analysis = analysis_name,
        Error = conditionMessage(e),
        stringsAsFactors = FALSE
      )
    }
  )

  d_sep <- tryCatch(
    {
      out <- repair_names(piecewiseSEM::dSep(fit$sem))
      out$Analysis <- analysis_name
      out
    },
    error = function(e) {
      data.frame(
        Analysis = analysis_name,
        Error = conditionMessage(e),
        stringsAsFactors = FALSE
      )
    }
  )

  write.csv(
    coefficients,
    file.path(OUTPUT_DIR, paste0(analysis_name, "_path_coefficients.csv")),
    row.names = FALSE
  )

  write.csv(
    r_squared,
    file.path(OUTPUT_DIR, paste0(analysis_name, "_R2.csv")),
    row.names = FALSE
  )

  write.csv(
    fisher_c,
    file.path(OUTPUT_DIR, paste0(analysis_name, "_FisherC.csv")),
    row.names = FALSE
  )

  write.csv(
    d_sep,
    file.path(OUTPUT_DIR, paste0(analysis_name, "_dSep.csv")),
    row.names = FALSE
  )

  capture.output(
    summary(fit$sem),
    file = file.path(OUTPUT_DIR, paste0(analysis_name, "_summary.txt"))
  )

  list(
    coefficients = coefficients,
    r_squared = r_squared,
    fisher_c = fisher_c,
    d_sep = d_sep
  )
}

## 14. General SEM analysis runner
run_sem_analysis <- function(pair_df, aggregation_method, analysis_name) {
  pair_summary <- aggregate_pair_metrics(pair_df, aggregation_method)
  sem_data <- build_sem_data(pair_summary)
  fit <- fit_piecewise_sem(sem_data)
  results <- save_sem_results(fit, analysis_name)

  write.csv(
    sem_data,
    file.path(OUTPUT_DIR, paste0(analysis_name, "_sample_level_data.csv")),
    row.names = FALSE
  )

  cat("\n========================================\n")
  cat("Analysis: ", analysis_name, "\n")
  cat("SEM sample size: ", nrow(sem_data), "\n")
  cat(
    "Range of pair counts per sample: ",
    min(sem_data$n_pairs),
    "to",
    max(sem_data$n_pairs),
    "\n"
  )
  cat("Median number of pairs per sample: ", median(sem_data$n_pairs), "\n")
  cat("========================================\n")
  print(summary(fit$sem))

  list(
    pair_summary = pair_summary,
    sem_data = sem_data,
    fit = fit,
    results = results
  )
}

## 15. Main analysis
main_analysis <- run_sem_analysis(
  pair_df = pair_with_meta,
  aggregation_method = "mean",
  analysis_name = "Main_all_pairs_mean"
)

sem_data_main <- main_analysis$sem_data
fit_main <- main_analysis$fit
result_main <- main_analysis$results

## 16. Variance inflation factor (VIF)
extract_vif <- function(model, model_name) {
  tryCatch(
    {
      v <- car::vif(model)
      if (is.matrix(v)) {
        out <- data.frame(
          Variable = rownames(v),
          v,
          row.names = NULL,
          check.names = FALSE
        )
      } else {
        out <- data.frame(
          Variable = names(v),
          VIF = as.numeric(v),
          row.names = NULL
        )
      }
      out$Model <- model_name
      out
    },
    error = function(e) {
      data.frame(
        Variable = NA_character_,
        Model = model_name,
        Error = conditionMessage(e),
        stringsAsFactors = FALSE
      )
    }
  )
}

vif_results <- bind_rows(
  extract_vif(fit_main$model_pH, "pH"),
  extract_vif(fit_main$model_richness, "Richness"),
  extract_vif(fit_main$model_HoS, "Homogeneous_selection"),
  extract_vif(fit_main$model_RichDiff, "Richness_difference")
)

write.csv(
  vif_results,
  file.path(OUTPUT_DIR, "Main_SEM_VIF.csv"),
  row.names = FALSE
)

## 17. Breusch-Pagan test
extract_bptest <- function(model, model_name) {
  tryCatch(
    {
      test <- lmtest::bptest(model)
      data.frame(
        Model = model_name,
        Statistic = unname(test$statistic),
        Df = unname(test$parameter),
        P_value = test$p.value,
        stringsAsFactors = FALSE
      )
    },
    error = function(e) {
      data.frame(
        Model = model_name,
        Statistic = NA_real_,
        Df = NA_real_,
        P_value = NA_real_,
        Error = conditionMessage(e),
        stringsAsFactors = FALSE
      )
    }
  )
}

bptest_results <- bind_rows(
  extract_bptest(fit_main$model_pH, "pH"),
  extract_bptest(fit_main$model_richness, "Richness"),
  extract_bptest(fit_main$model_HoS, "Homogeneous_selection"),
  extract_bptest(fit_main$model_RichDiff, "Richness_difference")
)

write.csv(
  bptest_results,
  file.path(OUTPUT_DIR, "Main_SEM_Breusch_Pagan_test.csv"),
  row.names = FALSE
)

## 18. Cook's distance and influential observations
extract_influence <- function(model, model_name, sample_names) {
  if (length(sample_names) != length(residuals(model))) {
    sample_names <- rownames(model$model)
  }

  data.frame(
    Sample = sample_names,
    Model = model_name,
    StandardizedResidual = rstandard(model),
    CookDistance = cooks.distance(model),
    Leverage = hatvalues(model),
    stringsAsFactors = FALSE
  )
}

influence_results <- bind_rows(
  extract_influence(fit_main$model_pH, "pH", sem_data_main$Sample),
  extract_influence(
    fit_main$model_richness,
    "Richness",
    sem_data_main$Sample
  ),
  extract_influence(
    fit_main$model_HoS,
    "Homogeneous_selection",
    sem_data_main$Sample
  ),
  extract_influence(
    fit_main$model_RichDiff,
    "Richness_difference",
    sem_data_main$Sample
  )
) %>%
  mutate(
    Cook_threshold = 4 / nrow(sem_data_main),
    Potentially_influential = CookDistance > Cook_threshold
  )

write.csv(
  influence_results,
  file.path(OUTPUT_DIR, "Main_SEM_influence_diagnostics.csv"),
  row.names = FALSE
)

## 19. Sensitivity analysis 1: median aggregation
median_analysis <- run_sem_analysis(
  pair_df = pair_with_meta,
  aggregation_method = "median",
  analysis_name = "Sensitivity_all_pairs_median"
)

result_median <- median_analysis$results

## 20. Sensitivity analysis 2: within-site sample pairs only
pair_same_site <- pair_with_meta %>% filter(site1 == site2)
if (nrow(pair_same_site) == 0) stop("No within-site sample pairs were found.")

cat("\nNumber of within-site sample pairs: ", nrow(pair_same_site), "\n")

write.csv(
  pair_same_site,
  file.path(OUTPUT_DIR, "same_site_pairwise_data.csv"),
  row.names = FALSE
)

same_site_analysis <- run_sem_analysis(
  pair_df = pair_same_site,
  aggregation_method = "mean",
  analysis_name = "Sensitivity_same_site_pairs"
)

result_same_site <- same_site_analysis$results

## 21. Combine results from the three analyses
all_coefficients <- bind_rows(
  result_main$coefficients,
  result_median$coefficients,
  result_same_site$coefficients
)

all_fisher_c <- bind_rows(
  result_main$fisher_c,
  result_median$fisher_c,
  result_same_site$fisher_c
)

all_r_squared <- bind_rows(
  result_main$r_squared,
  result_median$r_squared,
  result_same_site$r_squared
)

write.csv(
  all_coefficients,
  file.path(OUTPUT_DIR, "All_SEM_path_coefficients.csv"),
  row.names = FALSE
)

write.csv(
  all_fisher_c,
  file.path(OUTPUT_DIR, "All_SEM_FisherC_results.csv"),
  row.names = FALSE
)

write.csv(
  all_r_squared,
  file.path(OUTPUT_DIR, "All_SEM_R2_results.csv"),
  row.names = FALSE
)

## 22. Extract key pH paths
if (all(c("Predictor", "Response") %in% colnames(all_coefficients))) {
  important_pH_paths <- all_coefficients %>%
    filter(
      Predictor == "pH_z",
      Response %in% c("Richness_z", "HoS_z", "RichDiff_z")
    )
} else {
  important_pH_paths <- data.frame(
    Message = paste0(
      "The path coefficient table does not contain Predictor or Response columns; actual column names: ",
      paste(colnames(all_coefficients), collapse = ", ")
    ),
    stringsAsFactors = FALSE
  )
}

write.csv(
  important_pH_paths,
  file.path(OUTPUT_DIR, "Important_pH_paths_across_analyses.csv"),
  row.names = FALSE
)

## 23. Data summary
data_summary <- data.frame(
  Item = c(
    "Environmental samples",
    "ASV table samples",
    "Valid samples",
    "Matched pairwise comparisons",
    "Main SEM samples",
    "Main minimum pairs per sample",
    "Main maximum pairs per sample",
    "Main median pairs per sample",
    "Median sensitivity SEM samples",
    "Same-site pairwise comparisons",
    "Same-site sensitivity SEM samples"
  ),
  Value = c(
    nrow(env),
    nrow(ASV_richness),
    length(valid_samples),
    nrow(pair_with_meta),
    nrow(main_analysis$sem_data),
    min(main_analysis$sem_data$n_pairs),
    max(main_analysis$sem_data$n_pairs),
    median(main_analysis$sem_data$n_pairs),
    nrow(median_analysis$sem_data),
    nrow(pair_same_site),
    nrow(same_site_analysis$sem_data)
  ),
  stringsAsFactors = FALSE
)

write.csv(
  data_summary,
  file.path(OUTPUT_DIR, "SEM_data_summary.csv"),
  row.names = FALSE
)

capture.output(
  sessionInfo(),
  file = file.path(OUTPUT_DIR, "R_sessionInfo.txt")
)

## 24. Completion message
cat("\n========================================\n")
cat("All SEM analyses completed successfully.\n")
cat("Output directory: ", file.path(getwd(), OUTPUT_DIR), "\n")
cat("Main analysis sample size: ", nrow(main_analysis$sem_data), "\n")
cat("Median-aggregation sensitivity analysis sample size: ", nrow(median_analysis$sem_data), "\n")
cat("Within-site sensitivity analysis sample size: ", nrow(same_site_analysis$sem_data), "\n")
cat("\nKey output files:\n")
cat("1. All_SEM_path_coefficients.csv\n")
cat("2. Important_pH_paths_across_analyses.csv\n")
cat("3. All_SEM_FisherC_results.csv\n")
cat("4. All_SEM_R2_results.csv\n")
cat("5. Main_SEM_VIF.csv\n")
cat("6. Main_SEM_Breusch_Pagan_test.csv\n")
cat("7. Main_SEM_influence_diagnostics.csv\n")
cat("8. SEM_data_summary.csv\n")
cat("========================================\n")
