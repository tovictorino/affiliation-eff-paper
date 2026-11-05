## Cross-evaluate the simulated logL of psfm_gtre_zh() and sfa 1.2.0 psfm()
## at each other's optima, separating the two differences:
##   agg   : "meanz"  sigma_h_i = exp(mean_t(z_it)'g/2)      (psfm_gtre_zh)
##           "meansig" sigma_h_i = mean_t exp(z_it'g/2)       (psfm 1.2.0)
##   draws : "shared" one Halton set for all firms            (psfm_gtre_zh)
##           "perfirm" own Halton block per firm              (psfm 1.2.0)
.libPaths(c(normalizePath("r_lib_sfa120"), .libPaths()))
source("psfm_gtre_zh.R")
zhf <- readRDS("outputs_refit/gtre_zh_both_margins.rds")
f12 <- readRDS("outputs_refit/psfm120/fit_default.rds")
v0  <- lapply(list.files("outputs_refit/multiverse", "^fit_V0_.*[.]rds$", full.names = TRUE), readRDS)
ok  <- vapply(v0, function(f) f$opt$convergence == 0 && is.finite(f$hess_min_eig) && f$hess_min_eig > 0, logical(1))
saff <- vapply(v0, function(f) exp(0.5 * (f$coefficients[["Zh:(Intercept)"]] + f$coefficients[["Zh:filiacao"]])), numeric(1))
fi  <- v0[ok & saff > 0.10][[which.max(vapply(v0[ok & saff > 0.10], `[[`, numeric(1), "loglik"))]]

pf <- .parse_formula3(zhf$formula)
d  <- as.data.frame(zhf$data); d <- d[order(d$CNPJ), ]
Yv <- d[[pf$y]]; X <- model.matrix(pf$fx, d); Zu <- model.matrix(pf$fzu, d); Zh <- model.matrix(pf$fzh, d)
id <- d$CNPJ; indiv <- unique(id); gid <- match(id, indiv); N <- length(indiv); n <- nrow(d)
ti <- tabulate(gid, N)
W  <- rowsum(Zh, gid) / ti
s  <- -1; R <- 1000L
Ds <- .make_draws(R, discard = 100L, seed = 123)
Dp <- sfa:::.gtre_halton_draws(N = N, R = R, rand.gtre = 123)
## psfm 1.2.0 orders firms as unique() of the pdata.frame
pd <- plm::pdata.frame(zhf$data, index = c("CNPJ", "year_index"))
ind12 <- as.character(unique(as.vector(pd[, "CNPJ"])))
Dp <- Dp[match(as.character(indiv), ind12)]
ZRp <- do.call(rbind, lapply(Dp, function(m) m[, 1]))   # N x R
ZHp <- abs(do.call(rbind, lapply(Dp, function(m) m[, 2])))
kx <- ncol(X); pu <- ncol(Zu); qh <- ncol(Zh)

ll <- function(p, agg, draws) {
  sigv <- p[1]; sigr <- p[2]
  b <- p[3:(2 + kx)]; dl <- p[(3 + kx):(2 + kx + pu)]; g <- p[(3 + kx + pu):(2 + kx + pu + qh)]
  su2 <- exp(as.numeric(Zu %*% dl)); sg <- sqrt(sigv^2 + su2); lam <- sqrt(su2) / sigv
  sh <- if (agg == "meanz") sqrt(exp(as.numeric(W %*% g))) else
        as.numeric(rowsum(sqrt(exp(as.numeric(Zh %*% g))), gid)) / ti
  base <- s * (Yv - as.numeric(X %*% b))
  E <- if (draws == "shared") outer(sh[gid], Ds[, 2]) + base - s * sigr * rep(Ds[, 1], each = n)
       else sh[gid] * ZHp[gid, ] + base - s * sigr * ZRp[gid, ]
  A <- E / sg
  lp <- log(2) - 0.5 * log(2 * pi) - log(sg) - 0.5 * A^2 + pnorm(-A * lam, log.p = TRUE)
  S <- rowsum(lp, gid); mx <- apply(S, 1, max)
  sum(mx + log(rowSums(exp(S - mx))) - log(R))
}
P <- list(gtre_zh_best = unname(zhf$coefficients), gtre_zh_interior = unname(fi$coefficients),
          psfm120_default = unname(f12$coefficients))
cat("reported logL: gtre_zh best", zhf$loglik, " interior", fi$loglik,
    " psfm120", as.numeric(logLik(f12)), "\n\n")
cells <- expand.grid(agg = c("meanz", "meansig"), draws = c("shared", "perfirm"), stringsAsFactors = FALSE)
tab <- sapply(P, function(p) mapply(function(a, dr) ll(p, a, dr), cells$agg, cells$draws))
rownames(tab) <- paste(cells$agg, cells$draws, sep = "/")
print(round(tab, 3))
cat("\nZh coefficients:\n")
print(round(sapply(P, function(p) p[(3 + kx + pu):(2 + kx + pu + qh)]), 3))
cat("sigr:", round(sapply(P, `[`, 2), 4), "\n")
saveRDS(list(tab = tab, P = P), file.path(dirname(normalizePath("outputs_refit/psfm120/fit_default.rds")), "cross_ll.rds"))
