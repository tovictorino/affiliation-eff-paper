## Re-rank the three optima (GTRE_ZH best, GTRE_ZH interior, psfm 1.2.0
## default) with simulation error pushed down: R = 10,000 draws, several
## independent draw sets of each kind, both persistent-margin aggregations.
## Evaluation only (no re-optimization). Chunks of 1,000 draws keep memory low.
.libPaths(c(normalizePath("r_lib_sfa120"), .libPaths()))
source("psfm_gtre_zh.R")
zhf <- readRDS("outputs_refit/gtre_zh_both_margins.rds")
P   <- readRDS("outputs_refit/psfm120/cross_ll.rds")$P

pf <- .parse_formula3(zhf$formula)
d  <- as.data.frame(zhf$data); d <- d[order(d$CNPJ), ]
Yv <- d[[pf$y]]; X <- model.matrix(pf$fx, d); Zu <- model.matrix(pf$fzu, d); Zh <- model.matrix(pf$fzh, d)
id <- d$CNPJ; indiv <- unique(id); gid <- match(id, indiv); N <- length(indiv); n <- nrow(d)
ti <- tabulate(gid, N); W <- rowsum(Zh, gid) / ti
RBIG <- 10000L; CH <- 1000L

sets <- list(perfirm_123 = c("perfirm", 123), perfirm_7 = c("perfirm", 7), perfirm_2025 = c("perfirm", 2025),
             shared_123 = c("shared", 123), shared_7 = c("shared", 7))

worker <- function(set) {
  .libPaths(c(normalizePath("r_lib_sfa120"), .libPaths()))
  kind <- set[1]; seed <- as.integer(set[2])
  if (kind == "perfirm") {
    Dp <- sfa:::.gtre_halton_draws(N = N, R = RBIG, rand.gtre = seed)
    pd <- plm::pdata.frame(zhf$data, index = c("CNPJ", "year_index"))
    ind12 <- as.character(unique(as.vector(pd[, "CNPJ"])))
    Dp <- Dp[match(as.character(indiv), ind12)]
    ZR <- do.call(rbind, lapply(Dp, function(m) m[, 1]))
    ZH <- abs(do.call(rbind, lapply(Dp, function(m) m[, 2])))
    rm(Dp)
  } else {
    Ds <- .make_draws(RBIG, discard = 100L, seed = seed)
    ZR <- matrix(Ds[, 1], 1, RBIG); ZH <- matrix(Ds[, 2], 1, RBIG)
  }
  kx <- ncol(X); pu <- ncol(Zu); qh <- ncol(Zh); s <- -1
  ll <- function(p, agg) {
    sigv <- p[1]; sigr <- p[2]
    b <- p[3:(2 + kx)]; dl <- p[(3 + kx):(2 + kx + pu)]; g <- p[(3 + kx + pu):(2 + kx + pu + qh)]
    su2 <- exp(as.numeric(Zu %*% dl)); sg <- sqrt(sigv^2 + su2); lam <- sqrt(su2) / sigv
    sh <- if (agg == "meanz") sqrt(exp(as.numeric(W %*% g))) else
          as.numeric(rowsum(sqrt(exp(as.numeric(Zh %*% g))), gid)) / ti
    base <- s * (Yv - as.numeric(X %*% b))
    Lc <- sapply(seq_len(RBIG / CH), function(c) {
      cols <- ((c - 1) * CH + 1):(c * CH)
      E <- if (kind == "perfirm") sh[gid] * ZH[gid, cols] + base - s * sigr * ZR[gid, cols] else
           outer(sh[gid], ZH[1, cols]) + base - s * sigr * rep(ZR[1, cols], each = n)
      A <- E / sg
      lp <- log(2) - 0.5 * log(2 * pi) - log(sg) - 0.5 * A^2 + pnorm(-A * lam, log.p = TRUE)
      S <- rowsum(lp, gid); mx <- apply(S, 1, max)
      mx + log(rowSums(exp(S - mx)))                  # per-firm log-sum over the chunk
    })
    mx <- apply(Lc, 1, max)
    sum(mx + log(rowSums(exp(Lc - mx))) - log(RBIG))
  }
  out <- expand.grid(point = names(P), agg = c("meanz", "meansig"), stringsAsFactors = FALSE)
  out$set <- paste(set, collapse = "_")
  out$logL <- mapply(function(pt, a) ll(P[[pt]], a), out$point, out$agg)
  out
}

t0 <- Sys.time()
cl <- parallel::makeCluster(length(sets))
parallel::clusterExport(cl, c("zhf", "P", "pf", "Yv", "X", "Zu", "Zh", "indiv", "gid", "N", "n", "ti", "W",
                              "RBIG", "CH", ".make_draws"))
res <- do.call(rbind, parallel::parLapply(cl, sets, worker))
parallel::stopCluster(cl)
cat("elapsed", format(Sys.time() - t0), "\n\n")
saveRDS(res, "outputs_refit/psfm120/highR_ll.rds")
for (a in c("meanz", "meansig")) {
  cat("aggregation:", a, "(R = 10,000)\n")
  tab <- xtabs(logL ~ set + point, data = res[res$agg == a, ])
  print(round(unclass(tab)[, names(P)], 3))
  cat("\n")
}
