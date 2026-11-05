## refit_perfirm.R -- the paper's main GTRE_ZH model re-fitted with per-firm
## Halton draws (psfm_gtre_zh_pf.R, draws = "perfirm", as in sfa 1.2.0), same
## data, formula, 1,000 draws, BFGS settings and analytic gradient as fit_zh()
## in 05b_refit_all_models.R. Starts:
##   paper_best      the paper's best maximum (shared draws, logL -654.12)
##   paper_interior  the paper's best interior local maximum (-656.50)
##   psfm120         sfa 1.2.0 psfm()'s optimum from its own start (-648.05)
##   default         the routine's own OLS-based start
##   *_p1, *_p2      variance-function coefficients + N(0,1), seed 2025 + k
##                   (the perturbation rule of 05b_refit_all_models.R)
## Output: outputs_refit/psfm_pf/fits.rds, starts.csv

OUT <- file.path("outputs_refit", "psfm_pf"); dir.create(OUT, showWarnings = FALSE)
zhf <- readRDS("outputs_refit/gtre_zh_both_margins.rds")
P   <- readRDS("outputs_refit/psfm120/cross_ll.rds")$P
nm  <- names(zhf$coefficients)

perturb <- function(p0, k) {
  set.seed(2025 + k)
  vf <- grepl("^Z[uh]:", nm)
  p0[vf] <- p0[vf] + rnorm(sum(vf))
  unname(p0)
}
starts <- list(paper_best = P$gtre_zh_best, psfm120 = P$psfm120_default,
               paper_interior = P$gtre_zh_interior, default = NULL,
               paper_best_p1 = perturb(P$gtre_zh_best, 1), psfm120_p1 = perturb(P$psfm120_default, 1),
               paper_best_p2 = perturb(P$gtre_zh_best, 2), psfm120_p2 = perturb(P$psfm120_default, 2))

fit_pf <- function(st) {
  t0 <- Sys.time()
  f <- tryCatch(suppressWarnings(psfm_gtre_zh(zhf$formula, data = as.data.frame(zhf$data),
         individual = "CNPJ", inefdec = FALSE, halton_num = 1000, rand.gtre = 123,
         draws = "perfirm", start_val = st, maxit.bobyqa = 0, Method = "BFGS",
         maxit.optim = 3000, optHessian = TRUE, compute_te = FALSE)),
         error = function(e) { message(conditionMessage(e)); NULL })
  if (!is.null(f)) { f$elapsed <- Sys.time() - t0; f$data <- NULL }
  f
}

t0 <- Sys.time()
cl <- parallel::makeCluster(6)
parallel::clusterExport(cl, c("zhf", "fit_pf"))
parallel::clusterEvalQ(cl, source("psfm_gtre_zh_pf.R"))
fits <- parallel::parLapplyLB(cl, starts, fit_pf, chunk.size = 1L)
parallel::stopCluster(cl)
names(fits) <- names(starts)
saveRDS(fits, file.path(OUT, "fits.rds"))

tab <- do.call(rbind, lapply(names(fits), function(k) {
  f <- fits[[k]]
  if (is.null(f)) return(data.frame(start = k, loglik = NA))
  cf <- f$coefficients
  data.frame(start = k, loglik = f$loglik, conv = f$opt$convergence,
             max_abs_grad = f$max_abs_grad, hess_min_eig = f$hess_min_eig,
             sigr = cf[["sigr"]], zu_fil = cf[["Zu:filiacao"]],
             zh_int = cf[["Zh:(Intercept)"]], zh_fil = cf[["Zh:filiacao"]], zh_age = cf[["Zh:age"]],
             minutes = as.numeric(f$elapsed, units = "mins"))
}))
write.csv(tab, file.path(OUT, "starts.csv"), row.names = FALSE)
cat("total elapsed", format(Sys.time() - t0), "\n")
print(tab[order(-tab$loglik), ], digits = 5, row.names = FALSE)
