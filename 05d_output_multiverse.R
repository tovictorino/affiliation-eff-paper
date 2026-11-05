# =============================================================================
# 05d_output_multiverse.R
#
# Output-definition multiverse for the four-component GTRE_ZH model.
#
# The paper's frontier counts only (net) loans as output. Affiliated singular
# cooperatives place about 29% of assets with their central (COSIF 14500008),
# independents hold about 30% in securities (13000004): output the frontier
# does not count is correlated with affiliation by institutional design. The
# IDENTICAL model -- psfm_gtre_zh(), same frontier, same margins, same
# standardization, same Halton draws and seed -- is re-estimated under each
# output definition (all in 2024 reais, built in 01_build_data.R):
#
#   V0   y1             baseline, net loans as reported: 16000001 (already net
#                       of PCLD; 16000001 = 161 + 162 + 163 + 164 + 168 + 169)
#   V1   y_loans_gross  gross loans: 16000001 + |16900008|
#   V2   y_v2           net loans + securities (13000004)
#   V3   y_v3           earning assets: V2 + interfinancial liquidity
#                       (12000005) + central placements (14500008) [preferred]
#   V4   y_v3           V3, central-placement share added to the transient z
#
# Sample: every variant uses the same observations (positive output under all
# of them: the full 6,902), so only the output definition changes. No variant
# has an observation below R$0.01m, so no flooring is applied (counts in
# output_variants_validation.csv).
#
# Variants x starting values are independent jobs fitted in parallel by BFGS
# on the analytic gradient (10 starts per variant, section 4); the converged
# fit with the highest likelihood is kept (checks: optim code 0, positive-
# definite Hessian, remaining Newton step below STEP_TOL standard errors).
# Each variant is also fitted without filiacao on the persistent margin (6
# starts) for a likelihood-ratio test of gamma_fil = 0.
# Every fit is saved as soon as it finishes; REFIT = FALSE reuses saved fits
# and fits only the missing ones. Runtime with 1000 Halton draws on 8 workers:
# about 10-15 min per fit (80 fits, about 2.5 h) plus about 11 min for the
# efficiency scores. The 150-draw run (30 starts per variant) is kept in
# outputs_refit/_superseded/draws150_2310/multiverse/.
#
# Inputs  (read-only): outputs/df_main.parquet
#                      outputs_refit/gtre_zh_both_margins.rds (start, audit)
# Outputs (outputs_refit/multiverse/):
#   fit_<V>_<start>.rds              every fit, without its data copy
#                                    (start = ols | paper | rNN_sdS)
#   fit_<V>.rds                      selected fit per variant, with TE scores
#   output_variants_validation.csv   output/TA shares, zeros, floors, N
#   multiverse_affiliation.csv       affiliation coefficients, all variants
#   multiverse_diagnostics.csv       per fit: convergence, gradient, timing
#   multiverse_starts_summary.csv    per variant: starts at the optimum, and
#                                    coefficient ranges over local maxima
#   multiverse_mechanism.csv         transient scores vs central/securities
#                                    shares (affiliates, two-way FE)
#   multiverse_coefficients.csv      every coefficient x output, original units
#   table_output_multiverse.tex/.md  the referee table: every coefficient by
#                                    output, with a robustness count
#
# Usage: Rscript 05d_output_multiverse.R
# =============================================================================

## ---- 0. CONFIGURATION -------------------------------------------------------
OUTDIR <- file.path("outputs_refit", "multiverse")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
if (!dir.exists(OUTDIR)) stop("Could not create output directory: ", OUTDIR)
cat("Output folder (multiverse):", normalizePath(OUTDIR), "\n")

# environment variable REFIT=FALSE reuses saved fit_<V>_<start>.rds and fits
# only the missing ones (e.g. to regenerate the tables)
REFIT          <- Sys.getenv("REFIT", "TRUE") != "FALSE"
SEED           <- 2025
HALTON         <- 1000L   # Halton draws, as in 05b (see appendix_psfm_gtre_zh.R)
STEP_TOL       <- 0.05    # remaining Newton step accepted at an optimum, in SEs
N_RANDOM_SEEDS <- 4       # perturbed starts per SD and variant (+ ols + paper)
PERTURB_SD     <- c(1, 2) # every seed is used with each perturbation SD
N_RANDOM_R     <- 4       # perturbed starts (SD 1) of each restricted model
N_WORKERS      <- 8L      # parallel fits (about 1 GB of memory each)
set.seed(SEED)

suppressPackageStartupMessages({
  library(tidyverse)
  library(arrow)
  library(plm)
  library(parallel)
})
source("psfm_gtre_zh.R")

save_fit <- function(obj, file) {
  saveRDS(obj, file)
  cat("  [saved]", file, "\n")
}

## ---- 1. DATA (mirrors 05b exactly) -----------------------------------------
cat("\n===== 1. DATA =====\n")
df_main <- as.data.frame(read_parquet("outputs/df_main.parquet"))
needed <- c("y_loans_gross", "y_v2", "y_v3", "central_share")
if (!all(needed %in% names(df_main))) {
  stop("df_main lacks the output variants; re-run 01_build_data.R first.")
}

df_main$logAt <- log(df_main$TA)
yd <- model.matrix(~ factor(df_main$year_index) - 1) %>% as.data.frame()
names(yd) <- paste0("year_", 2015:2024)
df_main <- bind_cols(df_main, yd) %>% as.data.frame()

# Common estimation sample: positive output under every variant. V0 is the
# smallest variant (y1 <= gross, y1 <= V2 <= V3), so y1 > 0 suffices.
est <- df_main$y1 > 0
ly_of <- function(y) {
  out <- rep(NA_real_, length(y))
  out[est] <- -log(y[est])
  out
}
df_main <- df_main %>%
  mutate(
    ly_V0  = ly_of(y1),
    ly_V1  = ly_of(y_loans_gross),
    ly_V2  = ly_of(y_v2),
    ly_V3  = ly_of(y_v3)
  )
cat("Full panel:", nrow(df_main), "obs | common estimation sample:", sum(est),
    "obs |", n_distinct(df_main$CNPJ[est]), "coops\n")

# Same construction and standardization as 05b (scale() over the full panel);
# the central-placement share is standardized the same way.
panel_data <- pdata.frame(df_main, index = c("CNPJ", "year_index")) %>%
  mutate(ly1 = -log(y1), lx1 = log(x1), lx2 = log(x2), lx3 = log(x3))
panel_data_scaled <- panel_data %>%
  mutate(
    across(
      c(lx1, lx2, lx3, logAt, li, equity_ratio, deposit_ratio, age),
      ~ as.numeric(scale(.))
    ),
    central_share_s = as.numeric(scale(central_share))
  )
panel_data_scaled <- as.data.frame(panel_data_scaled)

## ---- 2. OUTPUT VARIANTS: VALIDATION ----------------------------------------
cat("\n===== 2. OUTPUT VARIANTS: VALIDATION =====\n")
variant_y <- c(V0 = "y1", V1 = "y_loans_gross", V2 = "y_v2", V3 = "y_v3")
FLOOR <- 1e4 # R$0.01m in 2024 reais

# Output/TA shares are taken over the common estimation sample.
validation <- imap_dfr(variant_y, function(yv, vid) {
  y   <- df_main[[yv]]
  sh  <- (y / df_main$TA)[est]
  aff <- (df_main$filiacao == 1)[est]
  tibble(
    variant           = vid,
    column            = yv,
    share_TA_all      = mean(sh),
    share_TA_aff      = mean(sh[aff]),
    share_TA_ind      = mean(sh[!aff]),
    n_nonpos_panel    = sum(y <= 0),
    n_below_floor_est = sum(y[est] < FLOOR),
    n_floored         = 0L,
    n_above_TA        = sum(y > df_main$TA * (1 + 1e-9)),
    n_est             = sum(est)
  )
})
print(as.data.frame(validation), digits = 3)
stopifnot(all(df_main$y_v3 >= 0))                # no negative earning assets
stopifnot(all(validation$n_below_floor_est == 0)) # no flooring needed
write_csv(validation, file.path(OUTDIR, "output_variants_validation.csv"))
cat("  [saved]", file.path(OUTDIR, "output_variants_validation.csv"), "\n")
cat("Estimation sample, share of TA (affiliated | independent):",
    "central placements",
    round(mean(df_main$central_share[est & df_main$filiacao == 1]), 3), "|",
    round(mean(df_main$central_share[est & df_main$filiacao == 0]), 3),
    "; securities",
    round(mean((df_main$securities / df_main$TA)[est & df_main$filiacao == 1]), 3), "|",
    round(mean((df_main$securities / df_main$TA)[est & df_main$filiacao == 0]), 3), "\n")

## ---- 3. SPECIFICATIONS -------------------------------------------------------
year_terms <- paste(paste0("year_", 2016:2024), collapse = " + ")
zu_base    <- "logAt + li + equity_ratio + filiacao + deposit_ratio"
fmla_zh <- function(lhs, zu, persistent = "filiacao + age") {
  as.formula(paste(lhs, "~ lx1 + lx2 + lx3 +", year_terms,
                   "|", zu, "|", persistent))
}

variants <- list(
  list(id = "V0",  label = "Net loans as reported, baseline (16000001)",
       lhs = "ly_V0", zu = zu_base),
  list(id = "V1",  label = "Gross loans (16000001 + |16900008|)",
       lhs = "ly_V1", zu = zu_base),
  list(id = "V2",  label = "Net loans + securities",
       lhs = "ly_V2", zu = zu_base),
  list(id = "V3",  label = "Earning assets (V2 + interbank + central placements)",
       lhs = "ly_V3", zu = zu_base),
  list(id = "V4",  label = "Earning assets, central share in z_u",
       lhs = "ly_V3", zu = paste(zu_base, "+ central_share_s"))
)
variants <- map(variants, function(v) c(v, list(fmla = fmla_zh(v$lhs, v$zu))))
names(variants) <- map_chr(variants, "id")

# Restricted models for the likelihood-ratio test of gamma_fil = 0 (the Monte
# Carlo shows that Hessian-based SEs understate the sampling spread of the
# persistent-margin coefficients): same specification without filiacao on
# the persistent margin. Ids end in "r" (V0r, ...).
variants_r <- map(variants, function(v) {
  list(id = paste0(v$id, "r"), base = v$id,
       label = paste(v$label, "[no affiliation on the persistent margin]"),
       fmla = fmla_zh(v$lhs, v$zu, persistent = "age"))
})
names(variants_r) <- map_chr(variants_r, "id")

## ---- 4. ESTIMATION (parallel over variants x starting values) ---------------
# Same model call as 05b's GTRE_ZH for every variant; only the response (and,
# for V4, one extra transient determinant) changes. Optimizer: BFGS on the
# analytic gradient (psfm_gtre_zh.R), which reaches 05b's V0 optimum in about
# 2 minutes instead of 2 hours; bobyqa is skipped and the iteration cap is
# raised so that the quasi-Newton stage is not cut off by maxit.
#
# The likelihood is multimodal for some outputs, so each variant is fitted
# from 2 + N_RANDOM_SEEDS * length(PERTURB_SD) starting values: the default
# OLS-based start, the paper's V0 optimum (the extra V4 coefficient starts at
# 0), and the paper optimum with the variance-function coefficients perturbed
# by N(0, sd) for every seed x sd, and sigv, sigr scaled by exp(N(0, 0.25)).
# The converged fit with the highest likelihood is kept; how many starts reach
# it, and the range of the affiliation coefficients over all converged local
# maxima, are reported.
paper_fit <- readRDS(file.path("outputs_refit", "gtre_zh_both_margins.rds"))

warm_start <- function(v) {
  p <- paper_fit$coefficients
  if (v$id == "V4") {
    k <- max(grep("^Zu:", names(p)))
    p <- c(p[seq_len(k)], `Zu:central_share_s` = 0, p[-seq_len(k)])
  }
  p
}
perturbed_start <- function(v, k, sd) {
  p  <- warm_start(v)
  set.seed(SEED + k)
  vf <- grepl("^Z[uh]:", names(p))
  p[vf]  <- p[vf] + rnorm(sum(vf), 0, sd)
  p[1:2] <- p[1:2] * exp(rnorm(2, 0, 0.25))
  unname(p)
}
rand_starts <- expand.grid(k = seq_len(N_RANDOM_SEEDS), sd = PERTURB_SD)
n_starts    <- 2L + nrow(rand_starts)

jobs <- unlist(lapply(variants, function(v) c(
  list(list(v = v, start = "ols",   start_val = NULL),
       list(v = v, start = "paper", start_val = unname(warm_start(v)))),
  lapply(seq_len(nrow(rand_starts)), function(r) {
    k <- rand_starts$k[r]; sd <- rand_starts$sd[r]
    list(v = v, start = sprintf("r%02d_sd%g", k, sd),
         start_val = perturbed_start(v, k, sd))
  })
)), recursive = FALSE)

# restricted models: OLS start, the warm start without Zh:filiacao, and
# N_RANDOM_R perturbations of it (SD 1)
drop_gamma_fil <- function(p) p[names(p) != "Zh:filiacao"]
n_starts_r <- 2L + N_RANDOM_R
jobs_r <- unlist(lapply(variants_r, function(vr) {
  v  <- variants[[vr$base]]
  p0 <- drop_gamma_fil(warm_start(v))
  c(list(list(v = vr, start = "ols",   start_val = NULL),
         list(v = vr, start = "paper", start_val = unname(p0))),
    lapply(seq_len(N_RANDOM_R), function(k) {
      list(v = vr, start = sprintf("r%02d_sd1", k),
           start_val = {
             set.seed(SEED + k)
             p  <- p0
             vf <- grepl("^Z[uh]:", names(p))
             p[vf]  <- p[vf] + rnorm(sum(vf), 0, 1)
             p[1:2] <- p[1:2] * exp(rnorm(2, 0, 0.25))
             unname(p)
           })
    }))
}), recursive = FALSE)
jobs <- c(jobs, jobs_r)

# Per-start fits are saved without their copy of the data (about 1.5 MB each);
# model_data() rebuilds it exactly as psfm_gtre_zh() does (complete cases on
# the model's variables, sorted by cooperative) when a selected fit needs it.
model_data <- function(fmla) {
  pf   <- .parse_formula3(fmla)
  used <- unique(c(pf$y, "CNPJ", all.vars(pf$fx), all.vars(pf$fzu), all.vars(pf$fzh)))
  d    <- panel_data_scaled[stats::complete.cases(panel_data_scaled[, used]), , drop = FALSE]
  d[order(d$CNPJ), , drop = FALSE]
}

fit_variant <- function(job) {
  set.seed(SEED)
  t0  <- Sys.time()
  fit <- tryCatch(
    psfm_gtre_zh(
      job$v$fmla,
      data         = panel_data_scaled,
      individual   = "CNPJ",
      inefdec      = FALSE,
      halton_num   = HALTON,
      rand.gtre    = 123,
      start_val    = job$start_val,
      maxit.bobyqa = 0,
      Method       = "BFGS",
      maxit.optim  = 3000,
      optHessian   = TRUE,
      compute_te   = FALSE   # efficiencies computed once per variant, below
    ),
    error = function(e) e
  )
  if (!inherits(fit, "error")) {
    fit$variant   <- job$v$id
    fit$start     <- job$start
    fit$wall_secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    slim <- fit; slim$data <- NULL
    saveRDS(slim, file.path(OUTDIR, sprintf("fit_%s_%s.rds", job$v$id, job$start)))
  }
  fit
}

# Remaining Newton step H^{-1} g in units of each parameter's SE: a converged
# optimum is within a small fraction of a standard error of the maximum.
newton_step_se <- function(fit) {
  st <- tryCatch(solve(fit$opt$hessian, fit$grad), error = function(e) NA)
  max(abs(st) / fit$std.errors)
}
fit_ok <- function(fit) {
  !inherits(fit, "error") &&
    fit$opt$convergence == 0 &&
    is.finite(fit$hess_min_eig) && fit$hess_min_eig > 0 &&
    all(is.finite(fit$std.errors)) &&
    isTRUE(newton_step_se(fit) < STEP_TOL)
}

n_workers <- min(length(jobs), max(1L, detectCores() - 1L), N_WORKERS)
cat("\n===== 4. ESTIMATION:", length(variants), "variants x",
    n_starts, "starts on", n_workers, "workers =====\n")
cl <- makePSOCKcluster(n_workers)
clusterCall(cl, setwd, getwd())
clusterExport(cl, c("panel_data_scaled", "OUTDIR", "SEED", "REFIT", "HALTON", "model_data"))
invisible(clusterEvalQ(cl, {
  source("psfm_gtre_zh.R")
  NULL
}))

job_names <- vapply(jobs, function(j) paste(j$v$id, j$start, sep = "_"), "")
job_files <- file.path(OUTDIR, paste0("fit_", job_names, ".rds"))
todo <- if (REFIT) seq_along(jobs) else which(!file.exists(job_files))
cat("Fitting", length(todo), "of", length(jobs), "jobs",
    if (!REFIT) "(saved fits reused)" else "", "\n")
t_start <- Sys.time()
all_fits <- vector("list", length(jobs))
if (length(todo)) all_fits[todo] <- parLapplyLB(cl, jobs[todo], fit_variant, chunk.size = 1L)
for (i in setdiff(seq_along(jobs), todo)) {
  all_fits[[i]] <- readRDS(job_files[i])
  all_fits[[i]]$start <- jobs[[i]]$start
}
names(all_fits) <- job_names
cat("Estimation wall time:",
    round(as.numeric(difftime(Sys.time(), t_start, units = "mins")), 1), "min\n")

starts <- map_dfr(names(all_fits), function(nm) {
  f <- all_fits[[nm]]
  if (inherits(f, "error")) {
    return(tibble(job = nm, error = conditionMessage(f)))
  }
  tibble(
    job            = nm,
    variant        = f$variant,
    start          = f$start,
    loglik         = f$loglik,
    delta_fil      = f$coefficients[["Zu:filiacao"]],
    gamma_fil      = unname(f$coefficients["Zh:filiacao"]),   # NA if restricted
    convergence    = f$opt$convergence,
    iterations     = f$opt$counts[["function"]],
    max_abs_grad   = f$max_abs_grad,
    newton_step_se = newton_step_se(f),
    hess_min_eig   = f$hess_min_eig,
    n_na_se        = sum(!is.finite(f$std.errors)),
    minutes        = f$wall_secs / 60,
    passes_checks  = fit_ok(f)
  )
})
print(as.data.frame(starts %>% select(-job)), digits = 4)

# Keep, per variant, the converged fit with the highest likelihood (if no
# start passes the checks, the highest likelihood is kept and flagged).
select_best <- function(id) {
  cand <- all_fits[startsWith(names(all_fits), paste0(id, "_"))]
  cand <- cand[!vapply(cand, inherits, logical(1), "error")]
  if (!length(cand)) return(simpleError(paste(id, ": every start failed")))
  ok   <- vapply(cand, fit_ok, logical(1))
  pool <- if (any(ok)) cand[ok] else cand
  pool[[which.max(vapply(pool, `[[`, numeric(1), "loglik"))]]
}
fits   <- map(variants,   function(v) select_best(v$id))
# V0 is the main model: when the multi-start search confirms 05b's fit (same
# sample, likelihood within 0.01), the V0 column reports that fit, so the
# baseline column matches the paper's main tables
if (!inherits(fits$V0, "error") && paper_fit$n == fits$V0$n &&
    abs(fits$V0$loglik - paper_fit$loglik) < 0.01) {
  cat(sprintf("V0: using the main fit (logLik %.4f; best of the search %.4f)\n",
              paper_fit$loglik, fits$V0$loglik))
  fits$V0 <- paper_fit
  fits$V0$variant <- "V0"; fits$V0$start <- "main_fit"
}
# restricted models: the highest likelihood reached by any converged start (a
# converged point bounds the restricted maximum from below even when its
# Newton-step or Hessian check fails, so the LR statistic is not overstated)
select_best_conv <- function(id) {
  cand <- all_fits[startsWith(names(all_fits), paste0(id, "_"))]
  cand <- cand[!vapply(cand, inherits, logical(1), "error")]
  cand <- cand[vapply(cand, function(f) f$opt$convergence == 0, logical(1))]
  if (!length(cand)) return(simpleError(paste(id, ": no converged start")))
  cand[[which.max(vapply(cand, `[[`, numeric(1), "loglik"))]]
}
fits_r <- map(variants_r, function(v) select_best_conv(v$id))

# likelihood-ratio test of gamma_fil = 0 for every output
lr_tests <- map_dfr(names(variants), function(id) {
  f <- fits[[id]]; fr <- fits_r[[paste0(id, "r")]]
  if (inherits(f, "error") || inherits(fr, "error")) {
    return(tibble(variant = id, LR = NA_real_, p_lr = NA_real_))
  }
  tibble(variant = id, loglik_full = f$loglik, loglik_restricted = fr$loglik,
         LR = 2 * (f$loglik - fr$loglik),
         p_lr = pchisq(pmax(2 * (f$loglik - fr$loglik), 0), 1, lower.tail = FALSE),
         restricted_ok = fit_ok(fr), restricted_start = fr$start)
})
cat("\nLikelihood-ratio tests of gamma_fil = 0:\n")
print(as.data.frame(lr_tests), digits = 4)
write_csv(lr_tests, file.path(OUTDIR, "multiverse_lr_gamma_fil.csv"))

# Efficiency scores (Colombi et al. CSN block) once per selected fit, in
# parallel; each completed fit is saved as fit_<V>.rds.
cat("\nEfficiency scores for the selected fits...\n")
t_te <- Sys.time()
fits <- parLapplyLB(cl, fits, function(f) {
  if (inherits(f, "error")) return(f)
  if (!is.null(f$U)) return(f)                # scores already computed (main fit)
  done <- file.path(OUTDIR, sprintf("fit_%s.rds", f$variant))
  if (!REFIT && file.exists(done)) {          # reuse scores of the same fit
    g <- readRDS(done)
    if (identical(g$start, f$start) && identical(g$coefficients, f$coefficients)) return(g)
  }
  if (is.null(f$data)) f$data <- model_data(f$formula)
  set.seed(SEED)
  f <- te_gtre_zh(f)
  saveRDS(f, done)
  f
})
names(fits) <- names(variants)
cat("Efficiency wall time:",
    round(as.numeric(difftime(Sys.time(), t_te, units = "mins")), 1), "min\n")
stopCluster(cl)

## ---- 5. V0 AUDIT: the main GTRE_ZH fit is the best optimum found ------------
# The V0 row is 05b's main model searched from 30 starts; 05b's fit should
# match it (same likelihood up to the optimizer tolerance).
cat("\n===== 5. V0 AUDIT (vs outputs_refit/gtre_zh_both_margins.rds) =====\n")
paper_fit <- readRDS(file.path("outputs_refit", "gtre_zh_both_margins.rds"))
if (!inherits(fits$V0, "error")) {
  d_coef <- fits$V0$coefficients - paper_fit$coefficients[names(fits$V0$coefficients)]
  cat(sprintf("logLik: main fit %.4f (n = %d) | V0 best of %d starts %.4f (n = %d)\n",
              paper_fit$loglik, paper_fit$n, n_starts, fits$V0$loglik, fits$V0$n))
  cat(sprintf("max |coef difference|: %.4f\n", max(abs(d_coef))))
  if (paper_fit$n == fits$V0$n && fits$V0$loglik > paper_fit$loglik + 0.01) {
    cat("WARNING: the multi-start search found a higher likelihood than the",
        "main fit; re-fit 05b from fit_V0.rds.\n")
  }
  for (term in c("Zu:filiacao", "Zh:filiacao")) {
    cat(sprintf("  %-12s main %+.4f (%.4f) | V0 %+.4f (%.4f)\n", term,
                paper_fit$coefficients[term], paper_fit$std.errors[term],
                fits$V0$coefficients[term], fits$V0$std.errors[term]))
  }
}

## ---- 6. AFFILIATION COEFFICIENTS ACROSS OUTPUT DEFINITIONS -------------------
cat("\n===== 6. MULTIVERSE: AFFILIATION COEFFICIENTS =====\n")
te_gap <- function(fit) {
  fil <- fit$data$filiacao
  H_o <- fit$H[as.character(fit$data$CNPJ)]
  c(total      = mean(fit$O[fil == 1]) - mean(fit$O[fil == 0]),
    transient  = mean(fit$U[fil == 1]) - mean(fit$U[fil == 0]),
    persistent = mean(H_o[fil == 1]) - mean(H_o[fil == 0]))
}

# Fitted inefficiency scales by group: sigma_eta per cooperative (affiliated =
# affiliated in most of its years), sigma_u per observation. A group mean
# sigma_eta below BOUNDARY_SIG means the persistent component of that group is
# at the zero boundary, where gamma and its SE are not individually meaningful.
BOUNDARY_SIG   <- 0.05 # both group-mean sigma_eta below this: the persistent
                       # component vanishes
BOUNDARY_RATIO <- 0.2  # one group-mean sigma_eta below this fraction of the
                       # other's: that group's persistent scale is at the boundary
NEAR_LL      <- 5    # log-points: local maxima this close compete with the best
group_scales <- function(fit) {
  fil_coop <- tapply(fit$data$filiacao, as.character(fit$data$CNPJ), mean)
  aff_coop <- fil_coop[names(fit$sigma_h)] > 0.5
  c(sig_eta_aff = mean(fit$sigma_h[aff_coop]),
    sig_eta_ind = mean(fit$sigma_h[!aff_coop]),
    sig_u_aff   = mean(fit$sigma_u[fit$data$filiacao == 1]),
    sig_u_ind   = mean(fit$sigma_u[fit$data$filiacao == 0]))
}

multiverse <- map_dfr(variants, function(v) {
  f <- fits[[v$id]]
  val <- validation %>% filter(variant == if (v$id == "V4") "V3" else v$id)
  if (inherits(f, "error")) {
    return(tibble(variant = v$id, label = v$label, ok = FALSE))
  }
  cf <- f$coefficients; se <- f$std.errors
  gap <- if (!is.null(f$O)) te_gap(f) else c(total = NA, transient = NA, persistent = NA)
  gs  <- group_scales(f)
  tibble(
    variant        = v$id,
    label          = v$label,
    share_TA_aff   = val$share_TA_aff,
    share_TA_ind   = val$share_TA_ind,
    delta_fil      = cf["Zu:filiacao"],
    se_delta       = se["Zu:filiacao"],
    p_delta        = 2 * pnorm(-abs(cf["Zu:filiacao"] / se["Zu:filiacao"])),
    ratio_u        = exp(cf["Zu:filiacao"] / 2),
    gamma_fil      = cf["Zh:filiacao"],
    se_gamma       = se["Zh:filiacao"],
    p_gamma_wald   = 2 * pnorm(-abs(cf["Zh:filiacao"] / se["Zh:filiacao"])),
    LR_gamma       = lr_tests$LR[lr_tests$variant == v$id],
    p_gamma        = lr_tests$p_lr[lr_tests$variant == v$id],   # LR test
    ratio_eta      = exp(cf["Zh:filiacao"] / 2),
    delta_central  = if ("Zu:central_share_s" %in% names(cf)) cf["Zu:central_share_s"] else NA_real_,
    se_central     = if ("Zu:central_share_s" %in% names(se)) se["Zu:central_share_s"] else NA_real_,
    dual_sign      = cf["Zu:filiacao"] > 0 & cf["Zh:filiacao"] < 0,
    sig_eta_aff    = gs[["sig_eta_aff"]],
    sig_eta_ind    = gs[["sig_eta_ind"]],
    sig_u_aff      = gs[["sig_u_aff"]],
    sig_u_ind      = gs[["sig_u_ind"]],
    # both groups at the boundary: the persistent component vanishes, so there
    # is no persistent contrast to identify
    eta_vanishes   = gs[["sig_eta_aff"]] < BOUNDARY_SIG & gs[["sig_eta_ind"]] < BOUNDARY_SIG,
    # one group's persistent scale near zero relative to the other's (group
    # means differ by more than BOUNDARY_RATIO): the contrast's direction is
    # identified (LR test) but gamma_fil's magnitude is not
    eta_boundary   = !(gs[["sig_eta_aff"]] < BOUNDARY_SIG & gs[["sig_eta_ind"]] < BOUNDARY_SIG) &
                     min(gs[["sig_eta_aff"]], gs[["sig_eta_ind"]]) <
                       BOUNDARY_RATIO * max(gs[["sig_eta_aff"]], gs[["sig_eta_ind"]]),
    te_gap_total   = gap[["total"]],
    te_gap_trans   = gap[["transient"]],
    te_gap_pers    = gap[["persistent"]],
    RTS            = sum(-cf[c("lx1", "lx2", "lx3")] /
                           c(sd(log(df_main$x1)), sd(log(df_main$x2)), sd(log(df_main$x3)))),
    sigv           = cf["sigv"],
    sigr           = cf["sigr"],
    loglik         = f$loglik,
    n              = f$n,
    N              = f$N,
    ok             = fit_ok(f)
  )
})
print(as.data.frame(multiverse %>% select(variant, delta_fil, se_delta, gamma_fil,
                                          se_gamma, dual_sign, sig_eta_aff,
                                          sig_eta_ind, eta_boundary, ok)), digits = 4)

selected <- vapply(fits, function(f) {
  if (inherits(f, "error")) NA_character_ else paste(f$variant, f$start, sep = "_")
}, "")
diagnostics <- starts %>%
  mutate(selected = job %in% selected) %>%
  group_by(variant) %>%
  mutate(
    gap_to_best = max(loglik[passes_checks]) - loglik,
    at_best     = passes_checks & gap_to_best < 0.01
  ) %>%
  ungroup()
write_csv(diagnostics, file.path(OUTDIR, "multiverse_diagnostics.csv"))

# Per variant: how many starts reach the selected optimum, and the best
# converged local optimum that is NOT the global one (if any).
start_summary <- diagnostics %>%
  filter(variant %in% names(variants)) %>%     # restricted fits: LR test only
  group_by(variant) %>%
  summarise(
    n_starts          = n(),
    n_converged       = sum(passes_checks),
    n_at_best         = sum(at_best),
    delta_spread_best = diff(range(delta_fil[at_best])),
    gamma_spread_best = diff(range(gamma_fil[at_best])),
    local_gap         = if (any(passes_checks & gap_to_best > 0.5))
                          min(gap_to_best[passes_checks & gap_to_best > 0.5]) else NA_real_,
    local_delta       = if (is.na(local_gap)) NA_real_ else
                          delta_fil[passes_checks & gap_to_best == local_gap][1],
    local_gamma       = if (is.na(local_gap)) NA_real_ else
                          gamma_fil[passes_checks & gap_to_best == local_gap][1],
    # every converged local maximum, whatever its likelihood
    delta_min_conv    = min(delta_fil[passes_checks]),
    delta_max_conv    = max(delta_fil[passes_checks]),
    gamma_min_conv    = min(gamma_fil[passes_checks]),
    gamma_max_conv    = max(gamma_fil[passes_checks]),
    delta_sign_stable = all(sign(delta_fil[passes_checks]) == sign(delta_fil[at_best][1])),
    gamma_sign_stable = all(sign(gamma_fil[passes_checks]) == sign(gamma_fil[at_best][1])),
    # local maxima within NEAR_LL log-points of the best: if they disagree on
    # the sign of gamma_fil, the persistent contrast is not identified
    n_near            = sum(passes_checks & gap_to_best < NEAR_LL),
    gamma_min_near    = min(gamma_fil[passes_checks & gap_to_best < NEAR_LL]),
    gamma_max_near    = max(gamma_fil[passes_checks & gap_to_best < NEAR_LL]),
    gamma_sign_stable_near = all(sign(gamma_fil[passes_checks & gap_to_best < NEAR_LL]) ==
                                   sign(gamma_fil[at_best][1])),
    .groups = "drop"
  )
cat("\nStarting values: starts reaching the selected optimum, the best other\n",
    "local optimum, and the range over all converged local maxima:\n")
print(as.data.frame(start_summary), digits = 3)
write_csv(start_summary, file.path(OUTDIR, "multiverse_starts_summary.csv"))

multiverse <- multiverse %>%
  left_join(start_summary %>%
              select(variant, n_converged, n_at_best, delta_sign_stable,
                     gamma_sign_stable_near, gamma_min_near, gamma_max_near),
            by = "variant") %>%
  mutate(gamma_fragile = !gamma_sign_stable_near)
write_csv(multiverse, file.path(OUTDIR, "multiverse_affiliation.csv"))
cat("  [saved]", file.path(OUTDIR, "multiverse_affiliation.csv"), "\n")
cat("  [saved]", file.path(OUTDIR, "multiverse_diagnostics.csv"), "\n")

## ---- 6b. MECHANISM: transient scores and asset composition ---------------
# Within an affiliated cooperative, are the years with a larger share of
# assets placed with the central (or held in securities) the years with
# larger measured transient inefficiency? Two-way (cooperative + year) FE
# regression of -log(TE_transient) on the two shares, affiliated
# cooperative-years, SEs clustered by cooperative. Under a loans-only output
# a reallocation from loans to central placements lowers measured output for
# given inputs, so a positive coefficient that vanishes once placements count
# as output identifies a portfolio-composition component of transient scores.
cat("\n===== 6b. MECHANISM: transient scores vs asset composition =====\n")
shares <- df_main %>%
  transmute(CNPJ = as.character(CNPJ), year_index = as.integer(year_index),
            central_share, securities_share = securities / TA)
mechanism <- map_dfr(names(fits), function(vid) {
  f <- fits[[vid]]
  if (inherits(f, "error") || is.null(f$U)) return(NULL)
  d <- tibble(CNPJ       = as.character(f$data$CNPJ),
              year_index = as.integer(as.character(f$data$year_index)),
              filiacao   = f$data$filiacao,
              ineff_u    = -log(f$U)) %>%
    left_join(shares, by = c("CNPJ", "year_index")) %>%
    filter(filiacao == 1)
  m  <- plm(ineff_u ~ central_share + securities_share,
            data = pdata.frame(as.data.frame(d), index = c("CNPJ", "year_index")),
            model = "within", effect = "twoways")
  ct <- lmtest::coeftest(m, vcov = vcovHC(m, cluster = "group", type = "HC1"))
  tibble(variant = vid, term = rownames(ct), estimate = ct[, 1],
         se = ct[, 2], p = ct[, 4], n = nobs(m))
})
print(as.data.frame(mechanism), digits = 3)
write_csv(mechanism, file.path(OUTDIR, "multiverse_mechanism.csv"))
cat("  [saved]", file.path(OUTDIR, "multiverse_mechanism.csv"), "\n")

## ---- 7. REFEREE TABLE: every coefficient across output definitions ---------
# One column per output definition, so a reader sees at a glance which
# coefficients survive every output. Coefficients are in original units (the
# standardized estimate divided by the SD of the covariate, as in 05b/09);
# frontier rows are the cost elasticities -beta_j / sd(log x_j). The last
# column counts the outputs in which a coefficient has its baseline (V0) sign
# and p < 0.05. Persistent-margin coefficients of an output whose persistent
# contrast is not identified (boundary solution or local maxima of both signs)
# are flagged and not counted.
stars <- function(p) {
  ifelse(is.na(p), "", ifelse(p < 0.01, "^{***}", ifelse(p < 0.05, "^{**}",
         ifelse(p < 0.1, "^{*}", ""))))
}
vids <- names(variants)
sd_of <- c(
  lx1 = sd(log(df_main$x1)), lx2 = sd(log(df_main$x2)), lx3 = sd(log(df_main$x3)),
  logAt = sd(df_main$logAt), li = sd(df_main$li),
  equity_ratio = sd(df_main$equity_ratio), deposit_ratio = sd(df_main$deposit_ratio),
  age = sd(df_main$age), central_share_s = sd(df_main$central_share)
)
coef_rows <- tribble(
  ~panel, ~label,                                   ~term,                 ~div,                       ~flip,
  "A",    "Transient scale, $\\hat\\delta_{fil}$",  "Zu:filiacao",         1,                          1,
  "A",    "Persistent scale, $\\hat\\gamma_{fil}$", "Zh:filiacao",         1,                          1,
  "B",    "log total assets",                       "Zu:logAt",            sd_of[["logAt"]],           1,
  "B",    "Labour intensity",                       "Zu:li",               sd_of[["li"]],              1,
  "B",    "Equity ratio",                           "Zu:equity_ratio",     sd_of[["equity_ratio"]],    1,
  "B",    "Deposit ratio",                          "Zu:deposit_ratio",    sd_of[["deposit_ratio"]],   1,
  "B",    "Central-placement share",                "Zu:central_share_s",  sd_of[["central_share_s"]], 1,
  "C",    "Age (per year)",                         "Zh:age",              sd_of[["age"]],             1,
  "D",    "Elasticity $x_1$ (costs)",                "lx1",                 sd_of[["lx1"]],            -1,
  "D",    "Elasticity $x_2$ (deposits)",             "lx2",                 sd_of[["lx2"]],            -1,
  "D",    "Elasticity $x_3$ (capital)",              "lx3",                 sd_of[["lx3"]],            -1
)

persistent_unidentified <- setNames(
  vapply(vids, function(v) {
    m <- multiverse[multiverse$variant == v, ]
    nrow(m) == 1 && isTRUE(m$eta_vanishes | m$gamma_fragile)
  }, logical(1)), vids)

coef_long <- map_dfr(vids, function(v) {
  f <- fits[[v]]
  if (inherits(f, "error")) return(NULL)
  cf <- f$coefficients; sev <- f$std.errors
  coef_rows %>%
    filter(term %in% names(cf)) %>%
    mutate(
      variant = v,
      est     = unname(flip * cf[term] / div),
      se      = unname(sev[term] / div),
      # Wald p-values, except on the persistent margin: gamma_fil is tested by
      # likelihood ratio, age is not tested (its Hessian SE is unreliable)
      p       = ifelse(term == "Zh:filiacao",
                       lr_tests$p_lr[lr_tests$variant == v],
                       ifelse(term == "Zh:age", NA_real_,
                              unname(2 * pnorm(-abs(cf[term] / sev[term]))))),
      flagged = startsWith(term, "Zh:") & persistent_unidentified[[v]]
    )
})
coef_long <- coef_long %>%
  group_by(term) %>%
  mutate(
    base_sign = sign(est[variant == "V0"][1]),
    robust    = sign(est) == base_sign & !is.na(p) & p < 0.05 & !flagged
  ) %>%
  ungroup()
write_csv(coef_long %>% select(-div, -flip),
          file.path(OUTDIR, "multiverse_coefficients.csv"))
cat("  [saved]", file.path(OUTDIR, "multiverse_coefficients.csv"), "\n")

robust_count <- coef_long %>%
  group_by(term) %>%
  summarise(k = sum(robust), n = n(), .groups = "drop")
cat("\nCoefficients with the baseline sign and p < 0.05, by output:\n")
print(as.data.frame(robust_count))

# --- cells -----------------------------------------------------------------
cell_coef <- function(trm, v) {
  r <- coef_long %>% filter(term == trm, variant == v)
  if (!nrow(r)) return("---")
  bnd  <- trm == "Zh:filiacao" && multiverse$eta_boundary[multiverse$variant == v]
  mark <- paste0(if (bnd) "\\dagger" else "", if (r$flagged) "\\ddagger" else "")
  sup  <- paste0(sub("^\\^\\{(.*)\\}$", "\\1", stars(r$p)), mark)
  sprintf("$%+.3f%s$", r$est, if (nzchar(sup)) paste0("^{", sup, "}") else "")
}
cell_se <- function(trm, v) {
  r <- coef_long %>% filter(term == trm, variant == v)
  if (!nrow(r)) "" else sprintf("(%.3f)", r$se)
}
cell_robust <- function(trm) {
  r <- robust_count %>% filter(term == trm)
  if (!nrow(r) || r$n < 2 || trm == "Zh:age") return("---")
  txt <- sprintf("%d/%d", r$k, r$n)
  if (r$k == r$n) paste0("\\textbf{", txt, "}") else txt
}
by_variant <- function(fun) vapply(vids, fun, character(1))
row_tex <- function(label, cells, last = "") {
  paste0(label, " & ", paste(cells, collapse = " & "), " & ", last, " \\\\")
}
mv <- function(v, col) multiverse[[col]][multiverse$variant == v]
val <- function(v, col) {
  validation[[col]][validation$variant == if (v == "V4") "V3" else v]
}

panel_rows <- function(pn) {
  unlist(lapply(which(coef_rows$panel == pn), function(i) {
    trm <- coef_rows$term[i]
    if (!any(coef_long$term == trm)) return(NULL)
    c(row_tex(coef_rows$label[i], by_variant(function(v) cell_coef(trm, v)), cell_robust(trm)),
      row_tex("", by_variant(function(v) cell_se(trm, v))))
  }))
}
ratio_row <- function(label, col) {
  row_tex(label, by_variant(function(v) {
    if (col == "ratio_eta") {
      if (isTRUE(mv(v, "eta_vanishes"))) return("n.i.")
      # range over the converged local maxima within NEAR_LL of the best
      lo <- exp(mv(v, "gamma_min_near") / 2); hi <- exp(mv(v, "gamma_max_near") / 2)
      if (is.finite(lo) && is.finite(hi) && hi - lo >= 0.005) {
        return(sprintf("%s--%.2f", if (lo < 0.005) "0.00" else sprintf("%.2f", lo), hi))
      }
    }
    sprintf("%.2f", mv(v, col))
  }))
}
sig_eta_row <- row_tex("\\quad $\\sigma_\\eta$, indep.\\ / affil.", by_variant(function(v) {
  sprintf("%.2f / %.2f", mv(v, "sig_eta_ind"), mv(v, "sig_eta_aff"))
}))
lr_row <- row_tex("\\quad LR test of $\\gamma_{fil}=0$", by_variant(function(v) {
  sprintf("%.1f", max(mv(v, "LR_gamma"), 0))
}))
rts_row <- row_tex("Returns to scale", by_variant(function(v) sprintf("%.3f", mv(v, "RTS"))))
sig_row <- function(label, col) {
  row_tex(label, by_variant(function(v) sprintf("%.3f", mv(v, col))))
}
ncol_all <- length(vids) + 2
panel_head <- function(txt) {
  sprintf("\\addlinespace\\multicolumn{%d}{l}{\\textit{%s}} \\\\", ncol_all, txt)
}

frag <- multiverse %>% filter(gamma_fragile)
van  <- multiverse %>% filter(eta_vanishes)
bnd  <- multiverse %>% filter(eta_boundary)
loc  <- start_summary %>% filter(variant %in% vids)
delta_stable <- all(loc$delta_sign_stable)
notes <- paste0(
  "\\item Each column re-estimates the identical model (",
  format(HALTON, big.mark = "{,}"), " Halton draws) on the same ",
  format(multiverse$n[1], big.mark = "{,}"), " cooperative-years, changing only ",
  "the output, and is the highest likelihood over ", n_starts, " starting values",
  if (delta_stable) "; $\\hat\\delta_{fil}$ keeps its sign at every converged local maximum" else "",
  ". Original units (per year for age); elasticities are $-\\beta_j$. ",
  "Hessian-based standard errors in parentheses. $^{*}$ $p<0.10$, $^{**}$ $p<0.05$, ",
  "$^{***}$ $p<0.01$: Wald tests, except $\\hat\\gamma_{fil}$, tested by likelihood ",
  "ratio against the highest likelihood without affiliation on the persistent margin (",
  n_starts_r, " starts, including the unrestricted optimum), and age, not tested ",
  "(persistent-margin standard errors are unreliable; Appendix~\\ref{app:mc}). ",
  "Robust count: outputs in which the coefficient has its V0 sign and $p<0.05$ (bold: ",
  "all five). $\\sigma_\\eta$ ratio: range of $e^{\\gamma/2}$ over the converged local ",
  "maxima within ", NEAR_LL, " log-points of the highest; next row: group means of ",
  "$\\hat\\sigma_\\eta$. Dual sign: $\\hat\\delta_{fil}>0$ and $\\hat\\gamma_{fil}<0$.",
  if (nrow(bnd)) paste0(
    " $^{\\dagger}$One group's persistent scale is near zero (",
    paste(bnd$variant, collapse = ", "),
    "): the sign of $\\hat\\gamma_{fil}$ is tested, its magnitude is not identified."
  ) else "",
  if (nrow(van) || nrow(frag)) paste0(
    " $^{\\ddagger}$Persistent margin not identified",
    if (nrow(van)) paste0(" (", paste(van$variant, collapse = ", "),
                          ": persistent scale near zero in both groups)") else "",
    if (nrow(frag)) paste0(" (", paste(frag$variant, collapse = ", "),
                           ": local maxima of both signs)") else "",
    "; not counted as robust."
  ) else ""
)

tex <- c(
  "\\begin{table}[!htbp]\\centering",
  "\\footnotesize",
  "\\caption{Output-definition multiverse: which coefficients survive every output.}",
  "\\label{tab:multiverse}",
  "\\setlength{\\tabcolsep}{3pt}",
  "\\renewcommand{\\arraystretch}{0.95}",
  "\\begin{threeparttable}",
  paste0("\\begin{tabular}{l", strrep("c", length(vids)), "c}"),
  "\\toprule",
  row_tex("", vids, "Robust"),
  row_tex("", c("Net loans", "Gross loans", "+ securities", "Earning assets",
                "V3 + central")[seq_along(vids)], "count"),
  row_tex("", c("(baseline)", "", "", "", "share in $z_u$")[seq_along(vids)], ""),
  "\\midrule",
  sub("\\addlinespace", "", panel_head("A. Affiliation"), fixed = TRUE),
  panel_rows("A")[1:2],
  ratio_row("\\quad $\\sigma_u$ ratio, $e^{\\delta/2}$", "ratio_u"),
  panel_rows("A")[3:4],
  lr_row,
  ratio_row("\\quad $\\sigma_\\eta$ ratio, $e^{\\gamma/2}$", "ratio_eta"),
  sig_eta_row,
  row_tex("Dual sign",
          by_variant(function(v) if (isTRUE(mv(v, "dual_sign"))) "yes" else "no")),
  panel_head("B. Other determinants of the transient scale"),
  panel_rows("B"),
  panel_head("C. Other determinant of the persistent scale"),
  panel_rows("C"),
  panel_head("D. Frontier and noise"),
  panel_rows("D"),
  rts_row,
  sig_row("$\\sigma_v$ (noise)", "sigv"),
  sig_row("$\\sigma_\\tau$ (heterogeneity)", "sigr"),
  panel_head("E. Output"),
  row_tex("Output/assets, affil.", by_variant(function(v) sprintf("%.2f", val(v, "share_TA_aff")))),
  row_tex("Output/assets, indep.", by_variant(function(v) sprintf("%.2f", val(v, "share_TA_ind")))),
  row_tex("Starts at the optimum", by_variant(function(v) {
    sprintf("%d/%d", mv(v, "n_at_best"), mv(v, "n_converged"))
  })),
  "\\bottomrule",
  "\\end{tabular}",
  "\\begin{tablenotes}\\footnotesize",
  notes,
  "\\end{tablenotes}",
  "\\end{threeparttable}",
  "\\end{table}"
)
writeLines(tex, file.path(OUTDIR, "table_output_multiverse.tex"))
cat("  [saved]", file.path(OUTDIR, "table_output_multiverse.tex"), "\n")

# Markdown: same matrix, coefficient and SE in one cell
md_cell <- function(trm, v) {
  r <- coef_long %>% filter(term == trm, variant == v)
  if (!nrow(r)) return("-")
  sprintf("%+.3f%s (%.3f)%s", r$est, gsub("[\\^{}]", "", stars(r$p)), r$se,
          if (r$flagged) " n.i." else "")
}
md <- c(
  "**Table. Output-definition multiverse: which coefficients survive every output (GTRE_ZH)**",
  "",
  paste0("| Coefficient | ", paste(vids, collapse = " | "), " | Same sign, p<0.05 |"),
  paste0("|---|", strrep("---:|", length(vids)), ":---:|"),
  unlist(lapply(seq_len(nrow(coef_rows)), function(i) {
    trm <- coef_rows$term[i]
    if (!any(coef_long$term == trm)) return(NULL)
    paste0("| ", gsub("\\$|\\\\hat|\\\\|\\{|\\}|\\^", "", coef_rows$label[i]), " | ",
           paste(by_variant(function(v) md_cell(trm, v)), collapse = " | "), " | ",
           gsub("\\\\textbf\\{|\\}", "", cell_robust(trm)), " |")
  })),
  paste0("| Dual sign | ", paste(by_variant(function(v) if (isTRUE(mv(v, "dual_sign"))) "yes" else "no"),
                                 collapse = " | "), " | |"),
  "",
  "V0 net loans (baseline), V1 gross loans, V2 net loans + securities, V3 earning",
  "assets, V4 V3 with the central-placement share in the transient scale. Original",
  "units; Hessian SEs; * p<0.10, ** p<0.05, *** p<0.01 (Wald; LR test for",
  "gamma_fil). n.i. = persistent margin not identified for that output (not",
  "counted as robust)."
)
writeLines(md, file.path(OUTDIR, "table_output_multiverse.md"))
cat("  [saved]", file.path(OUTDIR, "table_output_multiverse.md"), "\n")

## ---- 8. DECISION RULE --------------------------------------------------------
cat("\n===== 8. DECISION RULE =====\n")
alt <- multiverse %>% filter(variant %in% c("V1", "V2", "V3", "V4"))
n_alt   <- nrow(alt)
n_dual  <- sum(alt$dual_sign, na.rm = TRUE)
n_dsig  <- sum(alt$delta_fil > 0 & alt$p_delta < 0.05, na.rm = TRUE)
n_gsig  <- sum(alt$gamma_fil < 0 & alt$p_gamma < 0.05 & !alt$gamma_fragile & !alt$eta_vanishes,
               na.rm = TRUE)
n_gfrag <- sum(alt$gamma_fragile | alt$eta_vanishes, na.rm = TRUE)
cat(sprintf("Dual sign in %d of %d alternative outputs (V1-V4).\n", n_dual, n_alt))
cat(sprintf("Transient delta_fil > 0 and p < 0.05 in %d of %d.\n", n_dsig, n_alt))
cat(sprintf("Persistent gamma_fil < 0, p < 0.05 and identified in %d of %d (not identified in %d).\n",
            n_gsig, n_alt, n_gfrag))
v3 <- multiverse %>% filter(variant == "V3")
if (nrow(v3) && !is.na(v3$delta_fil)) {
  cat(sprintf(
    "Preferred output (V3): delta_fil %+.3f (p = %.3g) | gamma_fil %+.3f (LR p = %.3g)%s%s\n",
    v3$delta_fil, v3$p_delta, v3$gamma_fil, v3$p_gamma,
    if (v3$eta_boundary) " [persistent boundary]" else "",
    if (v3$gamma_fragile) " [persistent contrast not identified]" else ""
  ))
}
cat(if (n_dual == n_alt && n_dsig == n_alt) {
  "VERDICT: dual sign survives every output definition -> construct-validity threat answered.\n"
} else {
  "VERDICT: transient channel does not survive every output definition -> see H2b reframing.\n"
})

cat("\n===== DONE: multiverse saved in", normalizePath(OUTDIR), "=====\n")
