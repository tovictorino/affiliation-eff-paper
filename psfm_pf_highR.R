## psfm_pf_highR.R -- rank the maxima of refit_perfirm.R (and the three earlier
## points) by the paper's simulated logL at R = 10,000 per-firm draws, over
## four independent draw sets. Evaluation only; objective and draws from
## psfm_gtre_zh_pf.R. Chunks of 1,000 draws keep memory low.
## Output: outputs_refit/psfm_pf/highR.rds, highR.csv

source("psfm_gtre_zh_pf.R")
zhf  <- readRDS("outputs_refit/gtre_zh_both_margins.rds")
fits <- readRDS("outputs_refit/psfm_pf/fits.rds")
P0   <- readRDS("outputs_refit/psfm120/cross_ll.rds")$P
ok   <- vapply(fits, function(f) !is.null(f) && f$opt$convergence == 0, logical(1))
P    <- c(setNames(P0, c("old_paper_best", "old_paper_interior", "old_psfm120")),
          setNames(lapply(fits[ok], function(f) unname(f$coefficients)), paste0("pf_", names(fits)[ok])))
L1k  <- c(rep(NA, 3), vapply(fits[ok], `[[`, numeric(1), "loglik"))

pf <- .parse_formula3(zhf$formula)
d  <- as.data.frame(zhf$data); d <- d[order(d$CNPJ), ]
Yv <- d[[pf$y]]; X <- model.matrix(pf$fx, d); Zu <- model.matrix(pf$fzu, d); Zh <- model.matrix(pf$fzh, d)
indiv <- unique(d$CNPJ); gid <- match(d$CNPJ, indiv); N <- length(indiv)
W  <- rowsum(Zh, gid) / tabulate(gid, N)
CH <- 1000L
seeds <- c(123L, 7L, 2025L, 42L)

worker <- function(seed, RBIG = 10000L) {
  D <- .make_draws_perfirm(N, RBIG, seed = seed)
  kx <- ncol(X); pu <- ncol(Zu); qh <- ncol(W)
  vapply(P, function(p) {
    sigv <- p[1]; sigr <- p[2]
    b <- p[3:(2 + kx)]; dl <- p[(3 + kx):(2 + kx + pu)]; g <- p[(3 + kx + pu):(2 + kx + pu + qh)]
    sh  <- sqrt(exp(as.numeric(W %*% g)))
    su2 <- exp(as.numeric(Zu %*% dl)); sg <- sqrt(sigv^2 + su2); lam <- sqrt(su2) / sigv
    base <- -(Yv - as.numeric(X %*% b))                   # s = -1 (inefdec = FALSE)
    Lc <- sapply(seq_len(RBIG / CH), function(c) {         # per-firm log-sum per chunk
      cols <- ((c - 1) * CH + 1):(c * CH)
      E <- sh[gid] * D$h[gid, cols] + base + sigr * D$r[gid, cols]
      A <- E / sg
      lp <- log(2) - 0.5 * log(2 * pi) - log(sg) - 0.5 * A^2 + pnorm(-A * lam, log.p = TRUE)
      S <- rowsum(lp, gid); mx <- apply(S, 1, max)
      mx + log(rowSums(exp(S - mx)))
    })
    mx <- apply(Lc, 1, max)
    sum(mx + log(rowSums(exp(Lc - mx))) - log(RBIG))
  }, numeric(1))
}

## self-check: at R = 1,000 and seed 123 the evaluator must return the
## refit's own logL (and -649.4246 at the old psfm 1.2.0 point)
chk <- worker(123L, 1000L)
cat("self-check, R = 1,000 seed 123: max |diff| to refit logL =",
    signif(max(abs(chk - L1k), na.rm = TRUE), 3), "; old_psfm120 =", round(chk[["old_psfm120"]], 4), "\n")

t0 <- Sys.time()
cl <- parallel::makeCluster(length(seeds))
parallel::clusterExport(cl, c("P", "Yv", "X", "Zu", "W", "gid", "N", "CH"))
parallel::clusterEvalQ(cl, source("psfm_gtre_zh_pf.R"))
res <- parallel::parSapply(cl, seeds, worker)
parallel::stopCluster(cl)
colnames(res) <- paste0("perfirm_", seeds)
tab <- data.frame(point = names(P), logL_R1000_pf = L1k, res, mean_R10k = rowMeans(res),
                  sd_R10k = apply(res, 1, sd), check.names = FALSE)
tab <- tab[order(-tab$mean_R10k), ]
saveRDS(tab, "outputs_refit/psfm_pf/highR.rds")
write.csv(tab, "outputs_refit/psfm_pf/highR.csv", row.names = FALSE)
cat("elapsed", format(Sys.time() - t0), "\n")
print(tab, digits = 6, row.names = FALSE)
