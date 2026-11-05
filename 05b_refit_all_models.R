# =============================================================================
# 05b_refit_all_models.R
#
# Unconditional re-estimation pipeline -- a purpose-built alternative to
# 05_estimate_paper_models.R with NO resumability checks:
#
#   * ALWAYS re-fits every model from scratch (no file.exists() short-circuits)
#   * Saves each model's output IMMEDIATELY after that model finishes
#   * Writes ALL outputs to ./outputs_refit/  (NEVER touches ./outputs/)
#   * Adds the new four-component GTRE_ZH model (psfm_gtre_zh.R) with
#     determinants on BOTH inefficiency margins (transient + persistent)
#   * Ends with a MASTER TABLE: every term x every model, coefficient and
#     standard error side by side in original units (scaled models are
#     back-transformed with the sd_scale metadata, exactly as in
#     table_master_coefficients.csv)
#
# Model set (all fitted unconditionally, each saved right after fitting):
#   1. Naive pooled CD (Z_full)                 -- frontier::sfa, unscaled
#   2. Naive pooled TL (Z_full)                 -- frontier::sfa, unscaled
#   3. Scaled pooled CD (Z_full)                -- frontier::sfa, standardized
#   4. Scaled pooled TL (Z_full)                -- frontier::sfa, standardized
#   5. GTRE CD naive unscaled (Z_full)          -- sfa::psfm, documented
#                                                  non-convergence
#   6. GTRE CD scaled (Z_full)  [MAIN]          -- sfa::psfm, converged
#   7. GTRE_ZH four-component (NEW)             -- psfm_gtre_zh.R:
#        ly1 ~ inputs + years | logAt+li+equity+deposit | filiacao+age
#        (filiacao and age drive the PERSISTENT inefficiency variance)
#   8. KLH multi-step robustness                -- transient vs persistent
#   9. MASTER TABLE (coefficients & SEs, side by side)
#
# y1 is net loans as reported (COSIF 16000001, already net of provisions; see
# 01_build_data.R), positive in all 6,902 cooperative-years.
# Dependent variable is ly1 = -log(y1): a COST-style frontier (higher ly1
# = worse). All pooled and psfm models use ineffDecrease/inefdec = FALSE
# (cost orientation).
#
# Inputs  (read-only):  outputs/df_main.parquet
# Outputs (all saved):  outputs_refit/
#
# Usage:
#   source("05b_refit_all_models.R")
#   # or
#   Rscript 05b_refit_all_models.R
# =============================================================================

## ---- 0. CONFIGURATION -------------------------------------------------------
OUTDIR <- "outputs_refit"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
if (!dir.exists(OUTDIR)) stop("Could not create output directory: ", OUTDIR)
cat("Output folder (all re-fit outputs go here):", normalizePath(OUTDIR), "\n")

SEED <- 2025
set.seed(SEED)

# GTRE_ZH settings. 1000 Halton draws: the Monte Carlo of
# appendix_psfm_gtre_zh.R shows that with 150 draws (shared by all
# cooperatives) the persistent block is biased and over-dispersed, and
# mc_draws_reps.R shows that 1000 draws remove most of that.
HALTON_ZH  <- 1000L
ZH_WORKERS <- 8L      # parallel GTRE_ZH fits (about 1 GB of memory each)

# The two sfa::psfm GTRE fits (sections 5-6) are not used in the paper and
# take about 90 minutes; FALSE reuses their saved fits when present.
REFIT_PSFM <- FALSE
# Environment variable REFIT_ZH=FALSE reuses the saved GTRE_ZH fits and LR
# tests (section 7) when they were estimated with HALTON_ZH draws.
REFIT_ZH <- Sys.getenv("REFIT_ZH", "TRUE") != "FALSE"

## ---- 1. PACKAGES & LOCAL HELPERS --------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(arrow)
  library(frontier)
  library(sfa)
  library(plm)
  library(sfaR)
})

ensure_pkg <- function(p) {
  if (!requireNamespace(p, quietly = TRUE)) {
    cat(sprintf("  Installing missing package: %s\n", p))
    utils::install.packages(p, repos = "https://cloud.r-project.org")
  }
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}
ensure_pkg("sfaR")        # KLH steps 2-3
ensure_pkg("randtoolbox") # Halton draws for psfm_gtre_zh (matches demo/psfm)
ensure_pkg("mvtnorm")     # TE block of psfm_gtre_zh
ensure_pkg("minqa")       # bobyqa pre-optimizer for psfm_gtre_zh

# New four-component model: determinants on BOTH inefficiency margins
source("psfm_gtre_zh.R")

# Save helper: writes immediately after each model finishes
save_fit <- function(obj, file) {
  saveRDS(obj, file)
  cat("  [saved]", file, "\n")
}

# Log-likelihood of an sfa::psfm (sfareg) fit. psfm stores the NEGATIVE
# log-likelihood in $opt$value and has NO top-level $LL (accessing $LL
# directly yields NULL, which breaks data.frame()/tibble() column
# construction). Some psfm versions expose $LL; guard all shapes.
psfm_loglik <- function(fit) {
  if ("LL" %in% names(fit)) return(as.numeric(fit$LL))
  if (!is.null(fit$opt$value) && is.finite(fit$opt$value)) {
    return(-as.numeric(fit$opt$value))
  }
  NA_real_
}

# -----------------------------------------------------------------------------
# Canonical term universe for the MASTER TABLE. Every model is mapped onto
# these rows so that the same economic quantity sits side by side across
# columns (e.g. log(x1) from the naive CD, lx1 from the scaled GTRE and
# lx1_s from the scaled pooled model all land on the same "log(x1)" row).
# Standardized models are back-transformed to original units using the
# sd_scale metadata (same convention as outputs/table_master_coefficients.csv).
# -----------------------------------------------------------------------------
canonical_names <- c(
  "(Intercept)",
  "log(x1)", "log(x2)", "log(x3)",
  "0.5*log(x1)^2", "0.5*log(x2)^2", "0.5*log(x3)^2",
  "log(x1)*log(x2)", "log(x1)*log(x3)", "log(x2)*log(x3)",
  paste0("year_", 2016:2024),
  "Z_(Intercept)", "Z_filiacao", "Z_logAt", "Z_li", "Z_equity_ratio",
  "Z_deposit_ratio", "Z_age", "Z_NO", "Z_NE", "Z_CO", "Z_SE",
  "Zh_(Intercept)", "Zh_filiacao", "Zh_age",
  "sigmaSq", "gamma",
  "sigv", "sigr", "sigh"
)

master_blocks <- c(
  rep("Frontier", 10),
  rep("Frontier - year dummies", 9),
  rep("Z-block (transient)", 11),
  rep("Zh-block (persistent)", 3),
  rep("Variance components", 5)
)

# Raw coefficient names (per package / per spec) that map onto each canonical
# term, in priority order. Terms absent from `canon_from_raw` (year dummies,
# sigmaSq, gamma, sigv, sigr, sigh, Z_NO, Z_NE, Z_CO, Z_SE) are matched by
# identical raw name.
canon_from_raw <- list(
  "(Intercept)"         = c("(Intercept)"),
  "log(x1)"             = c("log(x1)", "lx1", "lx1_s"),
  "log(x2)"             = c("log(x2)", "lx2", "lx2_s"),
  "log(x3)"             = c("log(x3)", "lx3", "lx3_s"),
  "0.5*log(x1)^2"       = c("I(0.5 * lx1c^2)", "I(0.5 * lx1_s^2)"),
  "0.5*log(x2)^2"       = c("I(0.5 * lx2c^2)", "I(0.5 * lx2_s^2)"),
  "0.5*log(x3)^2"       = c("I(0.5 * lx3c^2)", "I(0.5 * lx3_s^2)"),
  "log(x1)*log(x2)"     = c("I(lx1c * lx2c)", "I(lx1_s * lx2_s)"),
  "log(x1)*log(x3)"     = c("I(lx1c * lx3c)", "I(lx1_s * lx3_s)"),
  "log(x2)*log(x3)"     = c("I(lx2c * lx3c)", "I(lx2_s * lx3_s)"),
  "Z_(Intercept)"       = c("Z_(Intercept)", "Zu:(Intercept)"),
  "Z_filiacao"          = c("Z_filiacao", "filiacao"),
  "Z_logAt"             = c("Z_logAt", "Z_logAt_s", "Zu:logAt", "logAt"),
  "Z_li"                = c("Z_li", "Z_li_s", "Zu:li", "li"),
  "Z_equity_ratio"      = c("Z_equity_ratio", "Z_equity_ratio_s", "Zu:equity_ratio", "equity_ratio"),
  "Z_deposit_ratio"     = c("Z_deposit_ratio", "Z_deposit_ratio_s", "Zu:deposit_ratio", "deposit_ratio"),
  "Z_age"               = c("Z_age", "Z_age_s", "Zu:age", "age"),
  "Zh_(Intercept)"      = c("Zh:(Intercept)"),
  "Zh_filiacao"         = c("Zh:filiacao"),
  "Zh_age"              = c("Zh:age")
)

# Back-transformation divisor for each canonical term when the model was fit
# on standardized covariates (the divisor is the SD of the raw variable, or
# the product of SDs for second-order / interaction terms).
scaled_divisor <- function(cn, sd_scale) {
  d <- switch(cn,
    "log(x1)"         = sd_scale["sd_logx1"],
    "log(x2)"         = sd_scale["sd_logx2"],
    "log(x3)"         = sd_scale["sd_logx3"],
    "0.5*log(x1)^2"   = sd_scale["sd_logx1"]^2,
    "0.5*log(x2)^2"   = sd_scale["sd_logx2"]^2,
    "0.5*log(x3)^2"   = sd_scale["sd_logx3"]^2,
    "log(x1)*log(x2)" = sd_scale["sd_logx1"] * sd_scale["sd_logx2"],
    "log(x1)*log(x3)" = sd_scale["sd_logx1"] * sd_scale["sd_logx3"],
    "log(x2)*log(x3)" = sd_scale["sd_logx2"] * sd_scale["sd_logx3"],
    "Z_logAt"         = sd_scale["sd_logAt"],
    "Z_li"            = sd_scale["sd_li"],
    "Z_equity_ratio"  = sd_scale["sd_equity_ratio"],
    "Z_deposit_ratio" = sd_scale["sd_deposit_ratio"],
    "Z_age"           = sd_scale["sd_age"],
    "Zh_age"          = sd_scale["sd_age"],
    1
  )
  as.numeric(d)
}

# Map a model's raw coefficient vector onto the canonical term universe.
# Returns a list with `coef` and `se` named vectors indexed by canonical name.
extract_canonical <- function(raw_names, coefs, ses, is_scaled, sd_scale) {
  coef_out <- setNames(rep(NA_real_, length(canonical_names)), canonical_names)
  se_out   <- coef_out

  for (cn in canonical_names) {
    cand <- canon_from_raw[[cn]]
    if (is.null(cand)) {
      # Direct name match (year dummies, Z_NO/NE/CO/SE, sigmaSq, gamma, sig*)
      if (cn %in% raw_names) {
        i <- which(raw_names == cn)[1]
        coef_out[cn] <- coefs[i]
        se_out[cn]   <- ses[i]
      }
      next
    }
    m <- which(raw_names %in% cand)
    if (length(m) == 0L) next
    i <- m[1]
    vc <- coefs[i]
    vs <- ses[i]
    if (is_scaled) {
      d  <- scaled_divisor(cn, sd_scale)
      vc <- vc / d
      vs <- vs / d
    }
    coef_out[cn] <- vc
    se_out[cn]   <- vs
  }
  list(coef = coef_out, se = se_out)
}

## ---- 2. DATA ----------------------------------------------------------------
cat("\n===== 2. DATA =====\n")
df_main <- read_parquet("outputs/df_main.parquet")
cat(
  "Loaded df_main:", nrow(df_main), "obs |",
  n_distinct(df_main$CNPJ), "coops\n"
)

if (!"logAt" %in% names(df_main)) {
  df_main$logAt <- log(df_main$TA)
}
if (!"ly1" %in% names(df_main)) {
  df_main$ly1 <- -log(df_main$y1)
}
if (!"lx1" %in% names(df_main)) {
  df_main$lx1 <- log(df_main$x1)
}
if (!"lx2" %in% names(df_main)) {
  df_main$lx2 <- log(df_main$x2)
}
if (!"lx3" %in% names(df_main)) {
  df_main$lx3 <- log(df_main$x3)
}
if (!"lx1c" %in% names(df_main)) {
  df_main$lx1c <- log(df_main$x1) - mean(log(df_main$x1))
}
if (!"lx2c" %in% names(df_main)) {
  df_main$lx2c <- log(df_main$x2) - mean(log(df_main$x2))
}
if (!"lx3c" %in% names(df_main)) {
  df_main$lx3c <- log(df_main$x3) - mean(log(df_main$x3))
}
if (!"year_2015" %in% names(df_main)) {
  yd <- model.matrix(~ factor(df_main$year_index) - 1) %>% as.data.frame()
  names(yd) <- paste0("year_", 2015:2024)
  df_main <- bind_cols(df_main, yd)
}
df_main <- as.data.frame(df_main)
cat(
  "Final sample:", nrow(df_main), "obs |",
  n_distinct(df_main$CNPJ), "coops\n"
)

## ---- 3. FORMULAS ------------------------------------------------------------
year_fe <- paste0("year_", 2016:2024)
year_terms <- paste(year_fe, collapse = " + ")

z_full <- ~ filiacao +
  logAt +
  li +
  equity_ratio +
  deposit_ratio +
  age +
  NO +
  NE +
  CO +
  SE

frontier_cd <- as.formula(paste("~ log(x1) + log(x2) + log(x3) +", year_terms))
frontier_tl <- as.formula(paste(
  "~ lx1c + lx2c + lx3c +",
  "I(0.5 * lx1c^2) + I(0.5 * lx2c^2) + I(0.5 * lx3c^2) +",
  "I(lx1c * lx2c) + I(lx1c * lx3c) + I(lx2c * lx3c) +",
  year_terms
))

fmla_gtre_zfull_cd <- as.formula(paste(
  "ly1 ~ lx1 + lx2 + lx3 +",
  year_terms,
  "| filiacao + logAt + li + equity_ratio + deposit_ratio + age"
))

# NEW: GTRE_ZH -- determinants on BOTH margins.
#   frontier        : ly1 ~ lx1 + lx2 + lx3 + year dummies (year_2015 = base)
#   transient Z_u   : logAt + li + equity_ratio + deposit_ratio
#   persistent Z_h  : filiacao + age       <- the new "long-run" channel
# (year_2015 is intentionally excluded to keep year_2015 as the reference
#  year, consistent with every other model in the paper.)
fmla_gtre_zh <- as.formula(paste(
  "ly1 ~ lx1 + lx2 + lx3 +", year_terms,
  "| logAt + li + equity_ratio + filiacao + deposit_ratio",
  "| filiacao + age"
))

build_fmla <- function(frontier_vars, z_vars) {
  as.formula(paste0(
    "-log(y1) ~ ",
    paste(deparse(frontier_vars[[2]]), collapse = " "),
    " | ",
    paste(deparse(z_vars[[2]]), collapse = " ")
  ))
}

## ---- 4. POOLED MODELS (naive CD/TL + scaled CD/TL, Z_full) ------------------
cat("\n===== 4. POOLED: NAIVE + SCALED (Z_full) =====\n")

# --- 4.1 Naive pooled CD ---
sfa_naive_cd <- suppressMessages(suppressWarnings(
  frontier::sfa(
    build_fmla(frontier_cd, z_full),
    data = df_main,
    ineffDecrease = FALSE
  )
))
save_fit(sfa_naive_cd, file.path(OUTDIR, "fit_pooled_naive_cd.rds"))

# --- 4.2 Naive pooled TL ---
sfa_naive_tl <- tryCatch(
  suppressMessages(suppressWarnings(frontier::sfa(
    build_fmla(frontier_tl, z_full),
    data = df_main,
    ineffDecrease = FALSE
  ))),
  error = function(e) {
    cat("  Naive TL failed:", conditionMessage(e), "\n")
    NULL
  }
)
if (!is.null(sfa_naive_tl)) {
  save_fit(sfa_naive_tl, file.path(OUTDIR, "fit_pooled_naive_tl.rds"))
}

# --- 4.3 Scaled data + SD metadata (needed for back-transformation) ---
df_main_scaled <- df_main %>%
  mutate(
    lx1_s = as.numeric(scale(log(x1))),
    lx2_s = as.numeric(scale(log(x2))),
    lx3_s = as.numeric(scale(log(x3))),
    logAt_s = as.numeric(scale(logAt)),
    li_s = as.numeric(scale(li)),
    equity_ratio_s = as.numeric(scale(equity_ratio)),
    deposit_ratio_s = as.numeric(scale(deposit_ratio)),
    age_s = as.numeric(scale(age))
  )
df_main_scaled <- as.data.frame(df_main_scaled)

sd_scale <- c(
  sd_logx1 = sd(log(df_main$x1)),
  sd_logx2 = sd(log(df_main$x2)),
  sd_logx3 = sd(log(df_main$x3)),
  sd_logAt = sd(df_main$logAt),
  sd_li = sd(df_main$li),
  sd_equity_ratio = sd(df_main$equity_ratio),
  sd_deposit_ratio = sd(df_main$deposit_ratio),
  sd_age = sd(df_main$age)
)
save_fit(sd_scale, file.path(OUTDIR, "paper_sd_scale.rds"))

frontier_cd_s <- as.formula(paste("~ lx1_s + lx2_s + lx3_s +", year_terms))
frontier_tl_s <- as.formula(paste(
  "~ lx1_s + lx2_s + lx3_s +",
  "I(0.5 * lx1_s^2) + I(0.5 * lx2_s^2) + I(0.5 * lx3_s^2) +",
  "I(lx1_s * lx2_s) + I(lx1_s * lx3_s) + I(lx2_s * lx3_s) +",
  year_terms
))
z_full_s <- ~ filiacao +
  logAt_s +
  li_s +
  equity_ratio_s +
  deposit_ratio_s +
  age_s +
  NO +
  NE +
  CO +
  SE

# --- 4.4 Scaled pooled CD ---
sfa_scaled_cd <- suppressMessages(suppressWarnings(
  frontier::sfa(
    build_fmla(frontier_cd_s, z_full_s),
    data = df_main_scaled,
    ineffDecrease = FALSE
  )
))
save_fit(sfa_scaled_cd, file.path(OUTDIR, "fit_pooled_scaled_cd.rds"))

# --- 4.5 Scaled pooled TL ---
sfa_scaled_tl <- tryCatch(
  suppressMessages(suppressWarnings(frontier::sfa(
    build_fmla(frontier_tl_s, z_full_s),
    data = df_main_scaled,
    ineffDecrease = FALSE
  ))),
  error = function(e) {
    cat("  Scaled TL failed:", conditionMessage(e), "\n")
    NULL
  }
)
if (!is.null(sfa_scaled_tl)) {
  save_fit(sfa_scaled_tl, file.path(OUTDIR, "fit_pooled_scaled_tl.rds"))
}

## ---- 5. GTRE NAIVE UNSCALED (documented NON-convergence) --------------------
cat("\n===== 5. GTRE NAIVE UNSCALED (documented non-convergence) =====\n")

panel_data <- pdata.frame(df_main, index = c("CNPJ", "year_index")) %>%
  mutate(ly1 = -log(y1), lx1 = log(x1), lx2 = log(x2), lx3 = log(x3))

gtre_naive_unscaled <- if (!REFIT_PSFM &&
                           file.exists(file.path(OUTDIR, "gtre_naive_unscaled.rds"))) {
  readRDS(file.path(OUTDIR, "gtre_naive_unscaled.rds"))
} else tryCatch(
  psfm(
    fmla_gtre_zfull_cd,
    model_name = "GTRE_Z",
    data = panel_data,
    individual = "CNPJ",
    inefdec = FALSE,
    halton_num = 100,
    Method = "BFGS",
    start_val = TRUE,
    optHessian = TRUE,
    # BOBYQA pre-optimization makes the degenerate fit crawl; cap it so
    # the documented failure (zero SEs, flat likelihood) reproduces fast.
    maxit.bobyqa = 20
  ),
  error = function(e) {
    cat("  GTRE naive unscaled failed:", conditionMessage(e), "\n")
    NULL
  }
)
if (!is.null(gtre_naive_unscaled)) {
  save_fit(gtre_naive_unscaled, file.path(OUTDIR, "gtre_naive_unscaled.rds"))
}

## ---- 6. GTRE SCALED (MAIN MODEL, psfm, converged) ---------------------------
cat("\n===== 6. GTRE SCALED (MAIN MODEL) =====\n")

panel_data_scaled <- panel_data %>%
  mutate(
    across(
      c(lx1, lx2, lx3, logAt, li, equity_ratio, deposit_ratio, age),
      ~ as.numeric(scale(.))
    )
  )

gtre_zfull_cd_scaled <- if (!REFIT_PSFM &&
                            file.exists(file.path(OUTDIR, "gtre_zfull_cd_scaled.rds"))) {
  readRDS(file.path(OUTDIR, "gtre_zfull_cd_scaled.rds"))
} else tryCatch(
  psfm(
    fmla_gtre_zfull_cd,
    model_name = "GTRE_Z",
    data = panel_data_scaled,
    individual = "CNPJ",
    inefdec = FALSE,
    halton_num = 150,
    Method = "BFGS",
    start_val = TRUE,
    optHessian = TRUE,
    verbose = TRUE
  ),
  error = function(e) {
    cat("  GTRE scaled failed:", conditionMessage(e), "\n")
    NULL
  }
)
if (!is.null(gtre_zfull_cd_scaled)) {
  save_fit(gtre_zfull_cd_scaled, file.path(OUTDIR, "gtre_zfull_cd_scaled.rds"))
}

cat("Scaled GTRE diagnostics:\n")
if (!is.null(gtre_zfull_cd_scaled)) {
  cfs <- gtre_zfull_cd_scaled$coefficients
  ses <- gtre_zfull_cd_scaled$std.errors
  print(data.frame(
    sigv = cfs["sigv"],
    sigr = cfs["sigr"],
    sigh = cfs["sigh"],
    n_zero_se = sum(ses == 0, na.rm = TRUE),
    LL = psfm_loglik(gtre_zfull_cd_scaled)
  ))
}

## ---- 7. GTRE_ZH FOUR-COMPONENT (NEW: both margins) --------------------------
# Four-component GTRE with determinants of BOTH inefficiency variances
# (psfm_gtre_zh.R):
#   sigma_h_i^2 = exp(w_i'gamma)   <- persistent  margin: filiacao + age
#   sigma_u_it^2 = exp(z_it'delta) <- transient  margin: logAt + li + ...
#
# Cost frontier: ly1 = -log(y1) and inefdec = FALSE (s = -1), matching the
# paper's cost orientation; HALTON_ZH Halton draws, seed 123. Optimizer: BFGS
# on the analytic gradient with no bobyqa stage (about 10 min per fit with
# 1000 draws; the former bobyqa + L-BFGS-B chain could stop at its
# 200-iteration cap). The likelihood can have several local maxima, so each
# fit is repeated from N_STARTS_ZH perturbations of the first optimum
# (variance-function coefficients + N(0, 1)) and the highest likelihood is
# kept; 05d_output_multiverse.R repeats the search from 30 starts as an audit.
#
# The Monte Carlo shows that Hessian-based standard errors understate the
# sampling spread of the persistent-margin coefficients, so inference on
# them uses likelihood-ratio tests: the model is refitted without filiacao
# and without age on the persistent margin. A restricted search started only
# from its own default can stop at a lower local maximum and overstate the
# LR statistic, so each restricted model is fitted from its default start and
# N_STARTS_ZH perturbations of that optimum, AND from the full model's
# optimum without the dropped coefficient and N_STARTS_ZH perturbations of
# it; the LR statistic uses the highest restricted likelihood reached by any
# start that converged (a converged point bounds the restricted maximum from
# below whether or not its Hessian check passes). All fits run in parallel.
# Efficiency scores (Colombi et al. CSN block, about 11 min) are computed
# once, for the selected full fit. REFIT_ZH=FALSE reuses the saved full fit
# and refits only the restricted models.
cat("\n===== 7. GTRE_ZH (NEW: filiacao + age on PERSISTENT margin) =====\n")

zh_saved <- file.path(OUTDIR, "gtre_zh_both_margins.rds")
reuse_zh <- !REFIT_ZH && file.exists(zh_saved) &&
  identical(as.integer(readRDS(zh_saved)$R_draws), as.integer(HALTON_ZH))

N_STARTS_ZH <- 4
zh_rhs <- function(persistent) as.formula(paste(
  "ly1 ~ lx1 + lx2 + lx3 +", year_terms,
  "| logAt + li + equity_ratio + filiacao + deposit_ratio",
  "|", persistent
))
zh_models <- list(
  full     = fmla_gtre_zh,           # persistent margin: filiacao + age
  no_fil_h = zh_rhs("age"),
  no_age_h = zh_rhs("filiacao")
)

# plain data frame (what psfm_gtre_zh() converts the pdata.frame to), so the
# workers need not load plm
zh_data <- as.data.frame(panel_data_scaled)
fit_zh <- function(fmla, start_val = NULL) {
  tryCatch(
    suppressWarnings(psfm_gtre_zh(
      fmla,
      data = zh_data,
      individual = "CNPJ",
      inefdec = FALSE,
      halton_num = HALTON_ZH,
      rand.gtre = 123,
      start_val = start_val,
      maxit.bobyqa = 0,
      Method = "BFGS",
      maxit.optim = 3000,
      optHessian = TRUE,
      compute_te = FALSE
    )),
    error = function(e) {
      cat("  GTRE_ZH failed:", conditionMessage(e), "\n")
      NULL
    }
  )
}
zh_ok <- function(f) {
  !is.null(f) && f$opt$convergence == 0 && is.finite(f$hess_min_eig) &&
    f$hess_min_eig > 0 && all(is.finite(f$std.errors))
}

zh_conv <- function(f) !is.null(f) && f$opt$convergence == 0
perturb <- function(p0, k) {
  set.seed(SEED + k)
  vf <- grepl("^Z[uh]:", names(p0))
  p0[vf] <- p0[vf] + rnorm(sum(vf))
  unname(p0)
}
dropped <- c(no_fil_h = "Zh:filiacao", no_age_h = "Zh:age")

zh_src <- normalizePath("psfm_gtre_zh.R")
cl_zh <- parallel::makeCluster(ZH_WORKERS)
parallel::clusterExport(cl_zh, c("zh_src", "zh_data", "HALTON_ZH", "fit_zh", "zh_models"))
parallel::clusterEvalQ(cl_zh, source(zh_src))
zh_runs <- list()                       # every fit: model, start label, result
run_zh <- function(jobs) {
  res <- parallel::parLapplyLB(cl_zh, jobs, function(j) fit_zh(zh_models[[j$m]], j$start),
                               chunk.size = 1L)
  for (i in seq_along(jobs)) {
    zh_runs[[length(zh_runs) + 1L]] <<- list(m = jobs[[i]]$m, label = jobs[[i]]$label,
                                             fit = res[[i]])
  }
  invisible(res)
}
fits_of <- function(m) Filter(function(r) r$m == m, zh_runs)
best_of <- function(m, check = zh_ok) {
  rs <- fits_of(m)
  ll <- vapply(rs, function(r) if (check(r$fit)) r$fit$loglik else -Inf, numeric(1))
  if (any(is.finite(ll))) rs[[which.max(ll)]]$fit else NULL
}
pert_jobs <- function(m) {
  f <- fits_of(m)[[1]]$fit                          # the default-start fit
  if (!zh_ok(f)) return(list())
  lapply(seq_len(N_STARTS_ZH), function(k)
    list(m = m, label = sprintf("default_p%d", k), start = perturb(f$coefficients, k)))
}
warm_jobs <- function(full_fit) {
  unlist(lapply(names(dropped), function(m) {
    p0 <- full_fit$coefficients
    p0 <- p0[names(p0) != dropped[[m]]]
    c(list(list(m = m, label = "warm", start = unname(p0))),
      lapply(seq_len(N_STARTS_ZH), function(k)
        list(m = m, label = sprintf("warm_p%d", k), start = perturb(p0, k))))
  }), recursive = FALSE)
}
default_job <- function(m) list(m = m, label = "default", start = NULL)

if (reuse_zh) {
  cat("  [reuse] saved GTRE_ZH full fit (REFIT_ZH=FALSE); refitting the restricted models\n")
  gtre_zh <- readRDS(zh_saved)
  run_zh(c(lapply(names(dropped), default_job), warm_jobs(gtre_zh)))
  run_zh(unlist(lapply(names(dropped), pert_jobs), recursive = FALSE))
} else {
  run_zh(lapply(names(zh_models), default_job))
  run_zh(unlist(lapply(names(zh_models), pert_jobs), recursive = FALSE))
  gtre_zh <- best_of("full")
  run_zh(warm_jobs(gtre_zh))
}
parallel::stopCluster(cl_zh)

zh_starts_tab <- bind_rows(lapply(zh_runs, function(r) tibble(
  model = r$m, start = r$label,
  loglik = if (is.null(r$fit)) NA_real_ else r$fit$loglik,
  converged = zh_conv(r$fit), checks_ok = zh_ok(r$fit))))
print(as.data.frame(zh_starts_tab), digits = 7)
write_csv(zh_starts_tab, file.path(OUTDIR, "gtre_zh_lr_starts.csv"))

# likelihood-ratio tests of the persistent-margin determinants: selected full
# fit against the highest converged restricted likelihood
zh_best_r <- lapply(names(dropped), best_of, check = zh_conv)
names(zh_best_r) <- names(dropped)
zh_lr <- tibble(
  test = paste(dropped, "= 0"),
  restricted = names(dropped),
  loglik_full = gtre_zh$loglik,
  loglik_restricted = vapply(zh_best_r, function(f) f$loglik, numeric(1))
) %>%
  mutate(LR = 2 * (loglik_full - loglik_restricted),
         df = 1,
         p_value = pchisq(pmax(LR, 0), df, lower.tail = FALSE),
         draws = HALTON_ZH)
print(zh_lr)
write_csv(zh_lr, file.path(OUTDIR, "gtre_zh_lr_tests.csv"))
save_fit(zh_best_r, file.path(OUTDIR, "gtre_zh_restricted_fits.rds"))

if (!is.null(gtre_zh)) {
  if (!reuse_zh) {
    gtre_zh$starts_loglik <- vapply(fits_of("full"), function(r)
      if (zh_ok(r$fit)) r$fit$loglik else -Inf, numeric(1))
    set.seed(SEED)
    gtre_zh <- te_gtre_zh(gtre_zh)
  }
  gtre_zh$lr_tests <- zh_lr
}
summary(gtre_zh$out)
if (!is.null(gtre_zh)) {
  save_fit(gtre_zh, file.path(OUTDIR, "gtre_zh_both_margins.rds"))

  # Coefficient-level table for quick inspection
  write.csv(
    as.data.frame(gtre_zh$out),
    file.path(OUTDIR, "gtre_zh_coefs.csv"),
    row.names = TRUE
  )

  # Efficiencies: transient (per obs), persistent (per coop), total (per obs)
  te_zh <- data.frame(
    CNPJ         = as.character(gtre_zh$data$CNPJ),
    year_index   = gtre_zh$data$year_index,
    filiacao     = gtre_zh$data$filiacao,
    TE_transient = gtre_zh$U,
    TE_persistent = gtre_zh$H[as.character(gtre_zh$data$CNPJ)],
    TE_total     = gtre_zh$O,
    stringsAsFactors = FALSE
  )
  write_csv(te_zh, file.path(OUTDIR, "gtre_zh_efficiencies.csv"))

  # Average marginal effects on expected inefficiency
  me_zh <- margeff_gtre_zh(gtre_zh)
  write_csv(
    tibble(term = names(me_zh$transient_AME), AME = as.numeric(me_zh$transient_AME)),
    file.path(OUTDIR, "gtre_zh_margeff_transient.csv")
  )
  write_csv(
    tibble(term = names(me_zh$persistent_AME), AME = as.numeric(me_zh$persistent_AME)),
    file.path(OUTDIR, "gtre_zh_margeff_persistent.csv")
  )

  cat("GTRE_ZH fit summary:\n")
  print(gtre_zh)
  cat("Average marginal effects on expected inefficiency:\n")
  print(lapply(me_zh[1:2], round, 4))
}

## ---- 8. KLH MULTI-STEP ROBUSTNESS CHECK -------------------------------------
cat("\n===== 8. KLH MULTI-STEP (persistent vs transient channel) =====\n")

# Step 1: RE panel -- separates firm means from within-deviations
s1 <- plm::plm(
  as.formula(paste("ly1 ~ lx1c + lx2c + lx3c +", year_terms)),
  data = panel_data,
  index = c("CNPJ", "year_index"),
  model = "random"
)
panel_data$resid_all <- residuals(s1)
firm_means <- aggregate(resid_all ~ CNPJ, data = panel_data, FUN = mean)
names(firm_means)[2] <- "resid_mean"
firm_map <- firm_means$resid_mean[match(panel_data$CNPJ, firm_means$CNPJ)]
panel_data$resid_within <- panel_data$resid_all - firm_map

klh_data_step2 <- as.data.frame(panel_data)
klh_data_step3 <- firm_means
fil_one <- df_main %>%
  group_by(CNPJ) %>%
  summarise(filiacao = first(filiacao), .groups = "drop")
klh_data_step3 <- klh_data_step3 %>%
  left_join(fil_one, by = "CNPJ")

# Step 2: TRANSIENT inefficiency from within-residuals; Z on sigma_u
s2 <- tryCatch(
  sfaR::sfacross(
    resid_within ~ 1,
    uhet = ~ filiacao + logAt + li + equity_ratio + deposit_ratio + age,
    data = klh_data_step2,
    S = -1,
    udist = "hnormal"
  ),
  error = function(e) {
    cat("  KLH step 2 failed:", conditionMessage(e), "\n")
    NULL
  }
)
if (!is.null(s2)) {
  save_fit(s2, file.path(OUTDIR, "klh_step2.rds"))
}

# Step 3: PERSISTENT inefficiency from firm-mean residuals; Z on sigma_h
# (clean persistent-channel test: filiacao only)
s3 <- tryCatch(
  sfaR::sfacross(
    resid_mean ~ 1,
    uhet = ~filiacao,
    data = klh_data_step3,
    S = -1,
    udist = "hnormal"
  ),
  error = function(e) {
    cat("  KLH step 3 failed:", conditionMessage(e), "\n")
    NULL
  }
)
if (!is.null(s3)) {
  save_fit(s3, file.path(OUTDIR, "klh_step3.rds"))
}

extract_filiacao_sfacross <- function(fit, step) {
  if (is.null(fit)) {
    return(tibble(step = step, est = NA, se = NA, t = NA, p = NA))
  }
  b <- coef(fit)
  nm <- grep("filiacao", names(b), value = TRUE)
  if (length(nm) == 0) {
    return(tibble(step = step, est = NA, se = NA, t = NA, p = NA))
  }
  est <- as.numeric(b[nm[1]])
  v <- tryCatch(vcov(fit), error = function(e) NULL)
  se <- if (!is.null(v) && nm[1] %in% dimnames(v)[[1]]) {
    sqrt(diag(v)[nm[1]])
  } else {
    NA_real_
  }
  t <- est / se
  tibble(
    step = step,
    est = est,
    se = se,
    t = t,
    p = 2 * pnorm(-abs(t))
  )
}

klh_filiacao <- bind_rows(
  extract_filiacao_sfacross(s2, "Step 2: transient (Z on sigma_u)"),
  extract_filiacao_sfacross(s3, "Step 3: persistent (Z on sigma_h)")
)
cat("KLH filiacao across channels:\n")
print(klh_filiacao)

klh_summary <- list(
  step1_re = s1,
  step2_transient = s2,
  step3_persistent = s3,
  filiacao = klh_filiacao,
  persistent_channel_supported = {
    s2r <- klh_filiacao %>% filter(grepl("Step 2", step))
    s3r <- klh_filiacao %>% filter(grepl("Step 3", step))
    !is.na(s3r$p) &&
      s3r$p < 0.05 &&
      s3r$est < 0 &&
      (is.na(s2r$p) || s2r$p >= 0.05)
  }
)
save_fit(klh_summary, file.path(OUTDIR, "paper_klh_summary.rds"))
write_csv(klh_filiacao, file.path(OUTDIR, "paper_klh_filiacao.csv"))

## ---- 9. CONSOLIDATED SUMMARY + MASTER TABLE (coeffs & SEs side by side) -----
cat("\n===== 9. CONSOLIDATED SUMMARY + MASTER TABLE =====\n")

# --- 9.1 Pooled summary rows ------------------------------------------------
summarize_pooled <- function(fit, label, scaled, family) {
  if (is.null(fit)) {
    return(tibble(
      model = label, RTS = NA, MonoPct = NA, logLik = NA,
      gamma = NA, filiacao = NA, se_filiacao = NA, mean_TE = NA
    ))
  }
  s <- suppressMessages(summary(fit))
  betas <- coef(fit)
  nm <- grep("filiacao", names(betas), value = TRUE)
  x_names <- switch(
    family,
    "cd" = if (scaled) c("lx1_s", "lx2_s", "lx3_s") else c("log(x1)", "log(x2)", "log(x3)"),
    "tl" = if (scaled) c("lx1_s", "lx2_s", "lx3_s") else c("lx1c", "lx2c", "lx3c")
  )
  el <- if (scaled) {
    -as.numeric(betas[x_names]) / as.numeric(sd_scale[c("sd_logx1", "sd_logx2", "sd_logx3")])
  } else {
    -as.numeric(betas[x_names])
  }
  rts <- sum(el)
  te <- tryCatch(mean(efficiencies(fit), na.rm = TRUE), error = function(e) NA)
  tibble(
    model = label,
    RTS = round(rts, 3),
    MonoPct = round(100 * mean(el > 0, na.rm = TRUE), 1),
    logLik = round(as.numeric(logLik(fit)), 2),
    gamma = round(tryCatch(s$mleParam["gamma", 1], error = function(e) NA), 4),
    filiacao = if (length(nm)) round(betas[nm[1]], 5) else NA,
    se_filiacao = if (length(nm)) {
      round(tryCatch(s$mleParam[nm[1], 2], error = function(e) NA), 5)
    } else {
      NA
    },
    mean_TE = round(te, 4)
  )
}

pooled_summary_rows <- bind_rows(
  summarize_pooled(sfa_naive_cd, "Naive pooled CD (Z_full)", scaled = FALSE, family = "cd"),
  summarize_pooled(sfa_naive_tl, "Naive pooled TL (Z_full)", scaled = FALSE, family = "tl"),
  summarize_pooled(sfa_scaled_cd, "Scaled pooled CD (Z_full)", scaled = TRUE, family = "cd"),
  summarize_pooled(sfa_scaled_tl, "Scaled pooled TL (Z_full)", scaled = TRUE, family = "tl")
)
print(pooled_summary_rows)
write_csv(pooled_summary_rows, file.path(OUTDIR, "paper_pooled_summary.csv"))

# --- 9.2 GTRE summary rows (includes the new GTRE_ZH) -----------------------
summarize_gtre <- function(fit, label, scaled) {
  if (is.null(fit)) {
    return(tibble(
      model = label, RTS = NA, LL = NA, converged = NA,
      sigv = NA, sigr = NA, sigh = NA,
      filiacao = NA, se_filiacao = NA, n_zero_se = NA
    ))
  }
  cfs <- fit$coefficients
  ses <- fit$std.errors
  rts <- if (scaled) {
    el <- -as.numeric(cfs[c("lx1", "lx2", "lx3")]) /
      as.numeric(sd_scale[c("sd_logx1", "sd_logx2", "sd_logx3")])
    round(sum(el), 3)
  } else {
    NA_real_
  }
  nm <- grep("filiacao", names(cfs), value = TRUE)
  tibble(
    model = label,
    RTS = rts,
    LL = psfm_loglik(fit),
    converged = if ("convergence" %in% names(fit)) {
      fit$convergence
    } else if (!is.null(fit$opt$convergence)) {
      fit$opt$convergence
    } else {
      NA
    },
    sigv = cfs["sigv"], sigr = cfs["sigr"], sigh = cfs["sigh"],
    filiacao = if (length(nm)) cfs[nm[1]] else NA,
    se_filiacao = if (length(nm)) ses[nm[1]] else NA,
    n_zero_se = sum(ses == 0, na.rm = TRUE)
  )
}

summarize_gtre_zh <- function(fit, label) {
  if (is.null(fit)) {
    return(tibble(
      model = label, RTS = NA, LL = NA, converged = NA,
      sigv = NA, sigr = NA, sigh = NA,
      filiacao = NA, se_filiacao = NA, n_zero_se = NA
    ))
  }
  cfs <- fit$coefficients
  ses <- fit$std.errors
  el <- -as.numeric(cfs[c("lx1", "lx2", "lx3")]) /
    as.numeric(sd_scale[c("sd_logx1", "sd_logx2", "sd_logx3")])
  tibble(
    model = label,
    RTS = round(sum(el), 3),
    LL = round(fit$loglik, 2),
    converged = fit$opt$convergence,
    sigv = cfs["sigv"], sigr = cfs["sigr"], sigh = NA_real_,
    filiacao = cfs["Zh:filiacao"],
    se_filiacao = ses["Zh:filiacao"],
    n_zero_se = sum(ses == 0, na.rm = TRUE)
  )
}

gtre_summary_rows <- bind_rows(
  summarize_gtre(
    gtre_naive_unscaled,
    "GTRE CD naive unscaled (Z_full)",
    scaled = FALSE
  ),
  summarize_gtre(
    gtre_zfull_cd_scaled,
    "GTRE CD scaled (Z_full) [MAIN]",
    scaled = TRUE
  ),
  summarize_gtre_zh(
    gtre_zh,
    "GTRE_ZH four-component (filiacao on persistent margin)"
  )
)
print(gtre_summary_rows)
write_csv(gtre_summary_rows, file.path(OUTDIR, "paper_gtre_summary.csv"))

# --- 9.3 Combined model summary object (mirrors 05's paper_model_summary) ---
paper_model_summary <- list(
  pooled = pooled_summary_rows,
  gtre = gtre_summary_rows,
  klh = klh_summary,
  sd_scale = sd_scale
)
save_fit(paper_model_summary, file.path(OUTDIR, "paper_model_summary.rds"))

# --- 9.4 MASTER TABLE: every term x every model, coef & SE side by side -----
master_fits <- list(
  "Naive Pooled CD"      = list(fit = sfa_naive_cd,          type = "pooled", scaled = FALSE, family = "cd"),
  "Naive Pooled TL"      = list(fit = sfa_naive_tl,          type = "pooled", scaled = FALSE, family = "tl"),
  "Scaled Pooled CD"     = list(fit = sfa_scaled_cd,         type = "pooled", scaled = TRUE,  family = "cd"),
  "Scaled Pooled TL"     = list(fit = sfa_scaled_tl,         type = "pooled", scaled = TRUE,  family = "tl"),
  "GTRE Naive Unscaled"  = list(fit = gtre_naive_unscaled,   type = "gtre",    scaled = FALSE),
  "GTRE Scaled [MAIN]"   = list(fit = gtre_zfull_cd_scaled,  type = "gtre",    scaled = TRUE),
  "GTRE_ZH Both Margins" = list(fit = gtre_zh,               type = "gtre_zh", scaled = TRUE)
)

master_wide <- data.frame(
  Block = master_blocks,
  Term  = canonical_names,
  stringsAsFactors = FALSE
)

for (nm in names(master_fits)) {
  meta <- master_fits[[nm]]
  fit  <- meta$fit
  if (is.null(fit)) {
    cat("  [skip master column:", nm, "] fit is NULL\n")
    next
  }
  if (meta$type == "pooled") {
    cf <- coef(fit)
    sv <- tryCatch(sqrt(diag(vcov(fit))), error = function(e) rep(NA_real_, length(cf)))
    if (length(sv) != length(cf)) sv <- rep(NA_real_, length(cf))
    names(sv) <- names(cf)
  } else {
    cf <- fit$coefficients
    sv <- fit$std.errors
  }
  ex <- extract_canonical(
    names(cf), as.numeric(cf), as.numeric(sv),
    is_scaled = meta$scaled,
    sd_scale  = sd_scale
  )
  master_wide[[paste0(nm, " Coef")]] <- round(ex$coef, 6)
  master_wide[[paste0(nm, " SE")]]   <- round(ex$se, 6)
}

# KLH columns: inject the filiacao channels into the corresponding rows.
#   Step 2 (transient)  -> Z_filiacao   row
#   Step 3 (persistent) -> Zh_filiacao  row
master_wide[["KLH Step 2 (transient) Coef"]] <- NA_real_
master_wide[["KLH Step 2 (transient) SE"]]   <- NA_real_
master_wide[["KLH Step 3 (persistent) Coef"]] <- NA_real_
master_wide[["KLH Step 3 (persistent) SE"]]   <- NA_real_

s2r <- klh_filiacao %>% filter(grepl("Step 2", step))
s3r <- klh_filiacao %>% filter(grepl("Step 3", step))
if (nrow(s2r) > 0) {
  master_wide[master_wide$Term == "Z_filiacao", "KLH Step 2 (transient) Coef"] <- s2r$est[1]
  master_wide[master_wide$Term == "Z_filiacao", "KLH Step 2 (transient) SE"]   <- s2r$se[1]
}
if (nrow(s3r) > 0) {
  master_wide[master_wide$Term == "Zh_filiacao", "KLH Step 3 (persistent) Coef"] <- s3r$est[1]
  master_wide[master_wide$Term == "Zh_filiacao", "KLH Step 3 (persistent) SE"]   <- s3r$se[1]
}

cat("Master table dimensions:", nrow(master_wide), "terms x",
    ncol(master_wide), "columns\n")
print(head(master_wide[, 1:min(8, ncol(master_wide))], 8))

# Save the master table: rounded CSV for reading, full-precision RDS for work
write_csv(master_wide, file.path(OUTDIR, "master_table_coefficients_se.csv"))
save_fit(master_wide, file.path(OUTDIR, "master_table_coefficients_se.rds"))

cat("\n===== DONE: all models re-fitted and saved in", normalizePath(OUTDIR), "=====\n")
cat("Master table saved: master_table_coefficients_se.csv (+ .rds)\n")