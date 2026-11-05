## psfm120_main_fit.R -- the paper's main GTRE_ZH model fitted with sfa::psfm
## (CRAN 1.2.0, which accepts y ~ x | z_u | z_h for GTRE_Z) instead of
## psfm_gtre_zh(). sfa 1.2.0 lives in r_lib_sfa120/ so the rest of the project
## keeps sfa 1.0.4 (mc_nested_check.R reads 1.0.4's internals).
##
## Usage: Rscript psfm120_main_fit.R default|warm
##   default : psfm's own start values and optimizer chain
##   warm    : started at psfm_gtre_zh()'s best maximum (logL -654.12)
## Same data, formula, orientation and 1000 Halton draws as the paper's fit.
## Output: outputs_refit/psfm120/fit_<mode>.rds

.libPaths(c(normalizePath("r_lib_sfa120"), .libPaths()))
suppressPackageStartupMessages(library(sfa))
stopifnot(packageVersion("sfa") >= "1.2.0")

mode <- commandArgs(trailingOnly = TRUE)[1]
stopifnot(mode %in% c("default", "warm"))
OUT <- file.path("outputs_refit", "psfm120"); dir.create(OUT, showWarnings = FALSE)

zh  <- readRDS("outputs_refit/gtre_zh_both_margins.rds")
dat <- plm::pdata.frame(zh$data, index = c("CNPJ", "year_index"))

t0  <- Sys.time()
fit <- psfm(zh$formula, model_name = "GTRE_Z", data = dat, individual = "CNPJ",
            inefdec = FALSE, halton_num = 1000, rand.gtre = 123,
            start_val = if (mode == "warm") unname(zh$coefficients) else FALSE,
            optHessian = TRUE)
fit$elapsed <- Sys.time() - t0
saveRDS(fit, file.path(OUT, sprintf("fit_%s.rds", mode)))

cat(sprintf("\nmode %s, elapsed %s, logL %.3f\n", mode, format(fit$elapsed),
            as.numeric(logLik(fit))))
print(round(rbind(psfm120 = fit$coefficients, se = fit$std.errors,
                  gtre_zh = zh$coefficients), 4))
