## =============================================================================
## psfm_gtre_zh_pf.R  --  copy of psfm_gtre_zh.R with a choice of draw scheme
##
##   draws = "shared"  : one Halton set for all individuals (psfm_gtre_zh.R,
##                       sfa 1.0.4); reproduces psfm_gtre_zh() exactly.
##   draws = "perfirm" : each individual gets its own block of R consecutive
##                       2-D Halton points with a seeded random shift, as in
##                       sfa 1.2.0 (.gtre_halton_draws). The analytic gradient
##                       is kept. Same function names, so sourcing this file
##                       replaces psfm_gtre_zh().
##
## Four-component ("generalized true random effects", GTRE) panel stochastic
## frontier with determinants of BOTH the persistent and the transient
## inefficiency variances.
##
## Model (production case, inefdec = TRUE; s = +1):
##
##   y_it = x_it' beta + r_i - s*eta_i + v_it - s*u_it
##
##   r_i    ~ N(0, sigma_r^2)                    time-invariant heterogeneity
##   eta_i  ~ N+(0, sigma_h_i^2)                 PERSISTENT inefficiency
##            sigma_h_i^2  = exp(w_i' gamma)     <- NEW: determinants w_i
##   v_it   ~ N(0, sigma_v^2)                    noise
##   u_it   ~ N+(0, sigma_u_it^2)                TRANSIENT inefficiency
##            sigma_u_it^2 = exp(z_it' delta)    (as in sfa::psfm "GTRE_Z")
##
## For a cost frontier (inefdec = FALSE; s = -1) the signs on eta and u flip,
## handled exactly as in sfa::psfm via the inefdec_n multiplier.
##
## Estimation: maximum simulated likelihood over Halton draws for (r_i, eta_i),
## following the GTRE_Z branch of sfa::psfm (Filippini & Greene 2016, JPA).
## Conditional on a draw of (r_i, eta_i) the remaining composed error
## v_it - s*u_it has the Aigner-Lovell-Schmidt closed-form density with
## observation-specific (sigma_u_it, lambda_it). The likelihood is evaluated
## for all observations at once (no loop over individuals) and the
## quasi-Newton stage and the Hessian use its analytic gradient
## (.zh_make_objective); analytic_grad = FALSE restores finite differences.
##
## Efficiency: closed-skew-normal conditional expectations of Colombi,
## Kumbhakar, Martini & Vittadini (2014, JPA), as in psfm's TE block, with the
## prior covariance of (eta_i, u_i1, ..., u_iT) now
##      VEE_i = diag(sigma_h_i^2, sigma_u_i1^2, ..., sigma_u_iT^2).
##
## Formula interface (three parts):
##      y ~ x1 + x2 | zu1 + zu2 | zh1 + zh2
##                    ^ transient   ^ persistent determinants
## Use "| 1" for an intercept-only (homoskedastic) margin. Persistent
## determinants must be time-invariant within individual; individual means are
## taken otherwise (Mundlak-style), with a warning.
##
## Dependencies: base R + stats. Optional: minqa (bobyqa pre-optimizer),
## pso (particle swarm, if PSopt = TRUE), randtoolbox (Halton draws; falls
## back to pseudo-random draws with a warning), mvtnorm (required only for
## the efficiency block, compute_te = TRUE).
##
## References:
##   Colombi, Kumbhakar, Martini & Vittadini (2014), J. of Productivity Analysis.
##   Kumbhakar, Lien & Hardaker (2014), J. of Productivity Analysis.
##   Filippini & Greene (2016), J. of Productivity Analysis.
##   Badunenko & Kumbhakar (2017), European J. of Operational Research.
##   Wang (2002), J. of Productivity Analysis (exp(z'd) variance function).
##   Bernstein, Parmeter & Stead: R package 'sfa' (psfm), CRAN.
## =============================================================================

## ---------------------------------------------------------------- helpers ---

.zh_constants <- list(
  MIN_POSITIVE = 1e-300,
  SIG_LOWER    = 1e-6,    # lower bound for sigma_v, sigma_r
  LIN_CLIP     = 30,      # clip |z'delta|, |w'gamma| to avoid exp overflow
  BIG          = 1e10     # per-individual penalty for non-finite likelihood
)

.parse_formula3 <- function(formula) {
  fs <- paste(deparse(formula, width.cutoff = 500L), collapse = " ")
  lr <- strsplit(fs, "~", fixed = TRUE)[[1]]
  if (length(lr) != 2L)
    stop("formula must have the form  y ~ x | z_transient | z_persistent")
  parts <- trimws(strsplit(lr[2], "|", fixed = TRUE)[[1]])
  if (length(parts) != 3L)
    stop("formula needs THREE right-hand parts:  y ~ x | z_transient | z_persistent\n",
         "  (use '| 1' for an intercept-only margin)")
  list(y   = trimws(lr[1]),
       fx  = stats::as.formula(paste("~", parts[1])),
       fzu = stats::as.formula(paste("~", parts[2])),
       fzh = stats::as.formula(paste("~", parts[3])))
}

.make_draws <- function(R, discard = 100L, seed = NULL, n_perm = 500L) {
  ## Two quasi-random (Halton) columns as in psfm: col 1 -> N(0,1) draws for
  ## r_i, col 2 -> half-normal draws for eta_i. qnorm((u+1)/2) is identical to
  ## sqrt(2)*erfinv(u) used in psfm (drops the pracma dependency).
  halton_ok <- requireNamespace("randtoolbox", quietly = TRUE)
  if (halton_ok) {
    Hm <- randtoolbox::halton(R + discard, dim = 2)
    Hm <- Hm[-seq_len(discard), , drop = FALSE]
  } else {
    warning("randtoolbox not installed: using pseudo-random draws instead of ",
            "Halton sequences. Results are still valid MSL but noisier; ",
            "install 'randtoolbox' to reproduce psfm's quasi-random draws.")
    if (!is.null(seed)) set.seed(seed)
    Hm <- cbind(stats::runif(R), stats::runif(R))
  }
  z_r <- stats::qnorm(Hm[, 1])
  z_h <- stats::qnorm((Hm[, 2] + 1) / 2)          # >= 0, half-normal quantiles
  ## Decorrelate the two columns by permuting z_r (psfm searches 9999 random
  ## permutations for minimum |cor|; we keep the idea with n_perm candidates).
  if (!is.null(seed)) set.seed(seed)
  best <- z_r; bc <- abs(stats::cor(z_r, z_h))
  for (k in seq_len(n_perm)) {
    cand <- sample(z_r)
    cc   <- abs(stats::cor(cand, z_h))
    if (cc < bc) { bc <- cc; best <- cand }
  }
  cbind(r = best, h = z_h)
}

.make_draws_perfirm <- function(N, R, discard = 1000L, seed = NULL, clamp = 1e-6) {
  ## Per-individual draws as in sfa 1.2.0's .gtre_halton_draws(): N*R
  ## consecutive 2-D Halton points (after a burn-in), shifted modulo 1 by a
  ## seeded uniform, individual i taking rows (i-1)*R+1 .. i*R. Returns two
  ## N x R matrices: N(0,1) draws for r_i and half-normal draws for eta_i.
  if (!requireNamespace("randtoolbox", quietly = TRUE))
    stop("draws = \"perfirm\" needs the 'randtoolbox' package")
  U <- randtoolbox::halton(N * R + discard, dim = 2, start = 1, normal = FALSE)
  U <- U[-seq_len(discard), , drop = FALSE]
  if (!is.null(seed)) {
    set.seed(seed)
    sh <- stats::runif(2L)
    U[, 1] <- (U[, 1] + sh[1]) %% 1
    U[, 2] <- (U[, 2] + sh[2]) %% 1
  }
  U[] <- pmin(pmax(U, clamp), 1 - clamp)
  list(r = matrix(stats::qnorm(U[, 1]), N, R, byrow = TRUE),
       h = matrix(stats::qnorm((U[, 2] + 1) / 2), N, R, byrow = TRUE))
}

.orthant_p <- function(mean, sigma) {
  ## P( N(mean, sigma) >= 0 componentwise ) -- the CSN normalizing probability.
  ## Equivalent to tmvtnorm::ptmvnorm(lowerx = 0, upperx = Inf, ...) in psfm.
  p <- try(mvtnorm::pmvnorm(lower = rep(0, length(mean)),
                            upper = rep(Inf, length(mean)),
                            mean  = as.numeric(mean),
                            sigma = sigma), silent = TRUE)
  if (inherits(p, "try-error")) return(NA_real_)
  max(as.numeric(p), .zh_constants$MIN_POSITIVE)
}

.zh_make_objective <- function(Yv, X, Zu, W, gid, zr, zh, s) {
  ## Negative simulated log-likelihood and its analytic gradient, vectorized
  ## over all n observations: E is the n x R matrix of composed errors (one
  ## column per draw), sums over t within individual are taken with rowsum(),
  ## and each individual's simulated likelihood is a log-sum-exp over draws.
  ## Same objective as psfm's per-individual GTRE_Z loop, without the R-level
  ## loop over individuals.
  ##
  ## Gradient: with A = E/sigma_it, B = -lambda_it*A, M = phi(B)/Phi(B) and
  ## posterior draw weights w_ir = L_ir / sum_r L_ir,
  ##   dlp/dE = -(A + M*lambda)/sigma,  dlp/dsigma = (A^2 - 1 - M*B)/sigma,
  ##   dlp/dlambda = -M*A,
  ## chained through E(beta, sigr, gamma) and (sigma, lambda)(sigv, delta).
  ## Clipped linear indices (|.| >= LIN_CLIP) get zero derivative.
  ## zr, zh: length-R vectors (shared draws) or N x R matrices (per-firm draws)
  pfd <- is.matrix(zr)
  n  <- length(Yv); N <- nrow(W); R <- if (pfd) ncol(zr) else length(zr)
  kx <- ncol(X); pu <- ncol(Zu); qh <- ncol(W)
  npar <- 2L + kx + pu + qh
  i_b <- 3L:(2L + kx)
  i_d <- (3L + kx):(2L + kx + pu)
  i_g <- (3L + kx + pu):(2L + kx + pu + qh)
  CLIP <- .zh_constants$LIN_CLIP
  BIGN <- .zh_constants$BIG / N
  logR <- log(R)
  LOG_SQRT_2PI <- 0.5 * log(2 * pi)
  if (pfd) {                              # per-firm draws expanded to n x R
    zr_n <- zr[gid, , drop = FALSE]
    zh_n <- zh[gid, , drop = FALSE]
  } else {
    zr_n <- rep(zr, each = n)             # draw-specific shift, column-major
  }

  function(p, deriv = FALSE) {
    sigv <- p[1]; sigr <- p[2]
    if (!all(is.finite(p)) || sigv <= 0 || sigr <= 0)
      return(list(value = .zh_constants$BIG, grad = rep(0, npar)))
    beta <- p[i_b]; delta <- p[i_d]; gam <- p[i_g]

    lh_raw <- as.numeric(W %*% gam)
    sh     <- sqrt(exp(pmin(pmax(lh_raw, -CLIP), CLIP)))     # sigma_h_i, length N
    lu_raw <- as.numeric(Zu %*% delta)
    su2    <- exp(pmin(pmax(lu_raw, -CLIP), CLIP))
    sg     <- sqrt(sigv^2 + su2)                             # sigma_it, length n
    lam    <- sqrt(su2) / sigv                               # lambda_it
    base   <- s * (Yv - as.numeric(X %*% beta))
    ## eps (n x R): s*(y - x'b - sigr*z_r + s*sigma_h_i*z_h)
    ##           =  base - s*sigr*z_r + sigma_h_i*z_h         (since s^2 = 1)
    E  <- (if (pfd) sh[gid] * zh_n else outer(sh[gid], zh)) + base - s * sigr * zr_n
    A  <- E / sg
    B  <- -A * lam
    lPB <- stats::pnorm(B, log.p = TRUE)                     # log Phi(B), reused
    lp <- (log(2) - LOG_SQRT_2PI) - log(sg) - 0.5 * A * A + lPB  # dnorm(A, log) inline
    S  <- rowsum(lp, gid)                                    # N x R
    mx <- S[cbind(seq_len(N), max.col(S, ties.method = "first"))]
    ES <- exp(S - mx)
    sE <- rowSums(ES)
    v  <- -(mx + log(sE) - logR)                             # -log L_i
    bad <- !is.finite(v)
    v[bad] <- BIGN
    value <- sum(v)
    if (!deriv) return(list(value = value))

    Wn  <- (ES / sE)[gid, , drop = FALSE]                    # posterior draw weights
    M   <- exp(-0.5 * B * B - LOG_SQRT_2PI - lPB)          # phi(B) / Phi(B)
    GE  <- Wn * (-(A + M * lam) / sg)
    gSG <- rowSums(Wn * (A^2 - 1 - M * B)) / sg
    gLM <- -rowSums(Wn * M * A)
    if (any(bad)) {                                          # penalized individuals
      ob <- bad[gid]
      GE[ob, ] <- 0; gSG[ob] <- 0; gLM[ob] <- 0
    }
    GEh  <- if (pfd) rowSums(GE * zh_n) else as.numeric(GE %*% zh)
    d_lh <- -(sh / 2) * as.numeric(rowsum(GEh, gid)) *
            (abs(lh_raw) < CLIP)
    d_lu <- -(gSG * su2 / (2 * sg) + gLM * lam / 2) * (abs(lu_raw) < CLIP)
    grad <- c(-sum(gSG * sigv / sg - gLM * lam / sigv),     # sigv
              s * (if (pfd) sum(GE * zr_n) else sum(GE %*% zr)),  # sigr
              s * as.numeric(crossprod(X, rowSums(GE))),     # beta
              as.numeric(crossprod(Zu, d_lu)),               # delta
              as.numeric(crossprod(W, d_lh)))                # gamma
    list(value = value, grad = grad)
  }
}

.zh_efficiency <- function(e_all, idx, sigv, sigr, sh_vec, su_all,
                           verbose = FALSE) {
  ## Colombi et al. (2014) closed-skew-normal conditional expectations, as in
  ## psfm's TE block, with prior covariance of (eta_i, u_i1, ..., u_iT)
  ## diag(sigma_h_i^2, sigma_u_i1^2, ..., sigma_u_iT^2).
  ## e_all: s*(y - x'beta) for every observation, rows grouped by idx.
  N <- length(idx)
  U <- O <- rep(NA_real_, length(e_all))
  H <- rep(NA_real_, N); names(H) <- names(sh_vec)
  for (ii in seq_len(N)) {
    ti  <- length(idx[[ii]])
    e_i <- e_all[idx[[ii]]]
    A   <- -cbind(rep(1, ti), diag(ti))
    ## SIGinv by Woodbury: (sv^2 I + sr^2 J)^{-1}
    SIGinv <- (diag(ti) - (sigr^2 / (sigv^2 + ti * sigr^2)) *
                 matrix(1, ti, ti)) / sigv^2
    vee    <- c(sh_vec[ii]^2, su_all[idx[[ii]]]^2)   # <- both margins here
    LAM    <- try(solve(diag(1 / vee, ti + 1L) + t(A) %*% SIGinv %*% A),
                  silent = TRUE)
    if (inherits(LAM, "try-error")) next
    LAM  <- (LAM + t(LAM)) / 2
    ARR  <- LAM %*% t(A) %*% SIGinv
    m0   <- as.numeric(ARR %*% e_i)
    den  <- .orthant_p(m0, LAM)
    if (!is.finite(den)) next
    tilt <- function(cv) {
      num <- .orthant_p(m0 + as.numeric(LAM %*% cv), LAM)
      val <- (num / den) * exp(sum(cv * m0) +
                                 0.5 * as.numeric(t(cv) %*% LAM %*% cv))
      min(max(val, 0), 1)
    }
    ch <- c(-1, rep(0, ti))
    H[ii] <- tilt(ch)                                  # E[exp(-eta_i)|.]
    for (j in seq_len(ti)) {
      cu <- rep(0, ti + 1L); cu[1L + j] <- -1
      co <- cu; co[1L] <- -1
      U[idx[[ii]][j]] <- tilt(cu)                      # E[exp(-u_it)|.]
      O[idx[[ii]][j]] <- tilt(co)                      # E[exp(-eta-u_it)|.]
    }
    if (verbose && ii %% 50L == 0L) cat("TE:", ii, "/", N, "individuals\n")
  }
  list(U = U, H = H, O = O)
}

## ---------------------------------------------------------- main function ---

psfm_gtre_zh <- function(formula, data, individual,
                         inefdec       = TRUE,
                         halton_num    = NULL,
                         halton_discard= 100L,
                         draws         = c("shared", "perfirm"),
                         intercept_zu  = TRUE,
                         intercept_zh  = TRUE,
                         start_val     = NULL,
                         maxit.bobyqa  = 300,
                         maxit.psoptim = 10,
                         maxit.optim   = 200,
                         PSopt         = FALSE,
                         optHessian    = TRUE,
                         Method        = "L-BFGS-B",
                         analytic_grad = TRUE,
                         compute_te    = TRUE,
                         rand.gtre     = NULL,
                         rand.psoptim  = NULL,
                         verbose       = FALSE) {

  call  <- match.call()
  t0    <- Sys.time()
  s     <- if (isTRUE(inefdec)) 1 else -1        # psfm's inefdec_n
  draws <- match.arg(draws)

  ## ---- data & design matrices -----------------------------------------
  pf <- .parse_formula3(formula)
  if (!pf$y %in% names(data))       stop("response '", pf$y, "' not in data")
  if (!individual %in% names(data)) stop("individual '", individual, "' not in data")

  data <- as.data.frame(data)
  used <- unique(c(pf$y, individual, all.vars(pf$fx),
                   all.vars(pf$fzu), all.vars(pf$fzh)))
  cc   <- stats::complete.cases(data[, used, drop = FALSE])
  if (any(!cc)) {
    message(sum(!cc), " incomplete rows dropped.")
    data <- data[cc, , drop = FALSE]
  }
  ## sort by individual so block bookkeeping (and hence U/O row alignment)
  ## is safe even for unsorted input
  data <- data[order(data[[individual]]), , drop = FALSE]

  Yv <- as.numeric(data[[pf$y]])
  X  <- stats::model.matrix(pf$fx,  data)   # handles factors (e.g. filiacao)
  Zu <- stats::model.matrix(pf$fzu, data)
  Zh <- stats::model.matrix(pf$fzh, data)
  if (!intercept_zu) Zu <- Zu[, colnames(Zu) != "(Intercept)", drop = FALSE]
  if (!intercept_zh) Zh <- Zh[, colnames(Zh) != "(Intercept)", drop = FALSE]
  if (ncol(Zu) == 0L) stop("transient margin has no columns; use '| 1 |' or intercept_zu = TRUE")
  if (ncol(Zh) == 0L) stop("persistent margin has no columns; use '| 1' or intercept_zh = TRUE")

  id    <- data[[individual]]
  indiv <- unique(id)
  N     <- length(indiv)
  idx   <- split(seq_len(nrow(data)), factor(id, levels = indiv))
  t_i   <- vapply(idx, length, integer(1))
  n     <- nrow(data)
  kx <- ncol(X); pu <- ncol(Zu); qh <- ncol(Zh)

  ## persistent determinants must be time-invariant: take individual means
  W <- do.call(rbind, lapply(idx, function(r) colMeans(Zh[r, , drop = FALSE])))
  colnames(W) <- colnames(Zh)
  nonint <- setdiff(colnames(Zh), "(Intercept)")
  if (length(nonint)) {
    wv <- vapply(nonint, function(j)
      max(vapply(idx, function(r)
        if (length(r) > 1) stats::sd(Zh[r, j]) else 0, numeric(1))), numeric(1))
    if (any(wv > 1e-8))
      warning("persistent-margin variable(s) vary within individuals: ",
              paste(names(wv)[wv > 1e-8], collapse = ", "),
              ". Individual means were used (Mundlak-style).")
  }

  ## ---- starting values -------------------------------------------------
  npar   <- 2L + kx + pu + qh
  pnames <- c("sigv", "sigr", colnames(X),
              paste0("Zu:", colnames(Zu)), paste0("Zh:", colnames(W)))

  if (is.numeric(start_val)) {
    if (length(start_val) != npar)
      stop("start_val must have length ", npar, " (order: sigv, sigr, beta, delta, gamma)")
    start_v <- start_val
  } else {
    ols  <- stats::lm.fit(X, Yv)
    b0   <- ols$coefficients; b0[is.na(b0)] <- 0
    res  <- Yv - as.numeric(X %*% b0)
    am   <- vapply(idx, function(r) mean(res[r]), numeric(1))   # between
    ew   <- res - rep(am, times = t_i)                          # within
    s2w  <- max(stats::var(ew), 1e-4)
    s2b  <- max(stats::var(am), 1e-4)
    d0   <- if (intercept_zu) c(log(s2w / 2), rep(0, pu - 1L)) else rep(0.1, pu)
    g0   <- if (intercept_zh) c(log(s2b / 2), rep(0, qh - 1L)) else rep(0.1, qh)
    start_v <- unname(c(sqrt(s2w / 2), sqrt(s2b / 2), b0, d0, g0))
  }

  ## ---- Halton draws: shared (psfm 1.0.4) or per individual (psfm 1.2.0) --
  R <- if (is.numeric(halton_num)) as.integer(halton_num) else ceiling(sqrt(n)) + 100L
  if (draws == "shared") {
    D  <- .make_draws(R, discard = halton_discard, seed = rand.gtre)
    zr <- D[, 1]; zh <- D[, 2]
  } else {
    D  <- .make_draws_perfirm(N, R, seed = rand.gtre)
    zr <- D$r; zh <- D$h
  }

  ## ---- simulated log-likelihood ---------------------------------------
  i_b <- 3L:(2L + kx)
  i_d <- (3L + kx):(2L + kx + pu)
  i_g <- (3L + kx + pu):(2L + kx + pu + qh)
  CLIP <- .zh_constants$LIN_CLIP

  obj <- .zh_make_objective(Yv, X, Zu, W, gid = match(id, indiv),
                            zr = zr, zh = zh, s = s)
  fn  <- function(p) obj(p)$value

  ## Quasi-Newton stage: one pass returns value and analytic gradient; the
  ## gradient is cached for the gr() call optim makes at the same point.
  ## analytic_grad = FALSE falls back to optim's finite differences.
  last <- new.env()
  fn_g <- function(p) {
    r <- obj(p, deriv = TRUE)
    last$p <- p; last$g <- r$grad
    r$value
  }
  gr <- function(p) {
    if (identical(p, last$p)) return(last$g)
    obj(p, deriv = TRUE)$grad
  }

  ## ---- optimization chain: bobyqa -> (psoptim) -> optim ----------------
  lower <- c(rep(.zh_constants$SIG_LOWER, 2L), rep(-Inf, npar - 2L))
  p <- start_v
  if (verbose) cat("initial nll:", fn(p), "\n")

  if (maxit.bobyqa > 0 && requireNamespace("minqa", quietly = TRUE)) {
    bo <- try(minqa::bobyqa(p, fn, lower = lower,
                            control = list(maxfun = max(maxit.bobyqa,
                                                        10L * npar^2 + 50L))),
              silent = TRUE)
    if (!inherits(bo, "try-error") && is.finite(bo$fval) && bo$fval < fn(p)) p <- bo$par
    if (verbose) cat("after bobyqa, nll:", fn(p), "\n")
  }

  if (isTRUE(PSopt) && requireNamespace("pso", quietly = TRUE)) {
    if (!is.null(rand.psoptim)) set.seed(rand.psoptim)
    ps <- try(pso::psoptim(p, fn,
                           lower  = pmax(lower, p - 10),
                           upper  = p + 10,
                           control = list(maxit = maxit.psoptim)), silent = TRUE)
    if (!inherits(ps, "try-error") && is.finite(ps$value) && ps$value < fn(p)) p <- ps$par
    if (verbose) cat("after psoptim, nll:", fn(p), "\n")
  }

  f_opt <- if (isTRUE(analytic_grad)) fn_g else fn
  g_opt <- if (isTRUE(analytic_grad)) gr else NULL
  op <- if (Method == "L-BFGS-B") {
    stats::optim(p, f_opt, g_opt, method = Method, lower = lower,
                 hessian = optHessian,
                 control = list(maxit = maxit.optim, trace = if (verbose) 1 else 0))
  } else {
    stats::optim(p, f_opt, g_opt, method = Method, hessian = optHessian,
                 control = list(maxit = maxit.optim, trace = if (verbose) 1 else 0))
  }
  if (verbose) cat("final nll:", op$value, " convergence:", op$convergence, "\n")

  ## convergence diagnostics: analytic gradient and Hessian curvature at the
  ## optimum (a converged interior maximum has |grad| ~ 0 and a PD Hessian)
  grad_opt <- obj(op$par, deriv = TRUE)$grad
  names(grad_opt) <- pnames
  hess_min_eig <- if (isTRUE(optHessian) && !is.null(op$hessian) &&
                      all(is.finite(op$hessian))) {
    min(eigen((op$hessian + t(op$hessian)) / 2, symmetric = TRUE,
              only.values = TRUE)$values)
  } else NA_real_

  ## ---- coefficient table ----------------------------------------------
  se <- rep(NA_real_, npar)
  if (isTRUE(optHessian) && !is.null(op$hessian) && all(is.finite(op$hessian))) {
    Vc <- try(solve(op$hessian), silent = TRUE)
    if (!inherits(Vc, "try-error")) {
      dv <- diag(Vc)
      se <- ifelse(dv > 0, sqrt(dv), NA_real_)
    }
  }
  out <- rbind(par = op$par, st_err = se, `t-val` = op$par / se)
  colnames(out) <- pnames

  est   <- op$par
  sigv  <- est[1]; sigr <- est[2]
  beta  <- est[i_b]; delta <- est[i_d]; gam <- est[i_g]
  su_all <- sqrt(exp(pmin(pmax(as.numeric(Zu %*% delta), -CLIP), CLIP)))  # n
  sh_vec <- sqrt(exp(pmin(pmax(as.numeric(W  %*% gam),   -CLIP), CLIP)))  # N
  names(sh_vec) <- as.character(indiv)

  ## ---- efficiency: Colombi et al. (2014) CSN conditional expectations --
  U <- H <- O <- NULL
  if (isTRUE(compute_te)) {
    if (!requireNamespace("mvtnorm", quietly = TRUE)) {
      warning("mvtnorm not installed: skipping efficiency computation ",
              "(re-run with compute_te = TRUE after installing mvtnorm).")
    } else {
      e_all <- s * (Yv - as.numeric(X %*% beta))                # NOT clipped
      te <- .zh_efficiency(e_all, idx, sigv, sigr, sh_vec, su_all, verbose)
      U <- te$U; H <- te$H; O <- te$O
    }
  }

  total_time <- difftime(Sys.time(), t0, units = "secs")

  results <- list(out         = t(out),
                  opt         = op,
                  data        = data,
                  total_time  = total_time,
                  start_v     = start_v,
                  model_name  = "GTRE_ZH",
                  formula     = formula,
                  inefdec_n   = s,          # +1 production, -1 cost orientation
                  U           = U,          # transient TE, per observation
                  H           = H,          # persistent TE, per individual
                  O           = O,          # overall TE E[exp(-eta-u)], per obs
                  sigma_u     = su_all,     # fitted sigma_u_it
                  sigma_h     = sh_vec,     # fitted sigma_h_i
                  loglik      = -op$value,
                  grad        = grad_opt,   # analytic gradient of -logL at optimum
                  max_abs_grad= max(abs(grad_opt)),
                  hess_min_eig= hess_min_eig,
                  n           = n, N = N, t_i = t_i, R_draws = R,
                  draws       = draws,
                  individual  = individual,
                  index       = idx,
                  zu_names    = colnames(Zu),
                  zh_names    = colnames(W),
                  coefficients= out["par", ],
                  std.errors  = out["st_err", ],
                  t.values    = out["t-val", ],
                  call        = call)
  class(results) <- c("psfm_gtre_zh", "sfareg")
  results
}

## ------------------------------------------------------------- S3 methods --

print.psfm_gtre_zh <- function(x, digits = 4, ...) {
  cat("Four-component panel SFA with determinants on both margins (GTRE_ZH)\n")
  cat("  sigma_h_i^2 = exp(w_i'gamma)  [persistent]   ",
      "sigma_u_it^2 = exp(z_it'delta)  [transient]\n\n")
  cat("N =", x$N, " individuals,  n =", x$n, " observations,  R =",
      x$R_draws, "draws,  logLik =", round(x$loglik, 2), "\n\n")
  tab <- cbind(Estimate = x$coefficients,
               `Std. Error` = x$std.errors,
               `t value` = x$t.values,
               `Pr(>|t|)` = 2 * stats::pnorm(-abs(x$t.values)))
  stats::printCoefmat(tab, digits = digits, has.Pvalue = TRUE, na.print = "NA")
  invisible(x)
}

## Efficiency scores for a fit estimated with compute_te = FALSE (or to
## recompute them): fills U, H, O exactly as compute_te = TRUE would.
te_gtre_zh <- function(fit, verbose = FALSE) {
  stopifnot(inherits(fit, "psfm_gtre_zh"))
  if (!requireNamespace("mvtnorm", quietly = TRUE))
    stop("te_gtre_zh() needs the 'mvtnorm' package.")
  pf <- .parse_formula3(fit$formula)
  X  <- stats::model.matrix(pf$fx, fit$data)
  cf <- fit$coefficients
  s  <- if (!is.null(fit$inefdec_n)) fit$inefdec_n else
          if (isFALSE(eval(fit$call$inefdec))) -1 else 1
  e_all <- s * (as.numeric(fit$data[[pf$y]]) - as.numeric(X %*% cf[colnames(X)]))
  te <- .zh_efficiency(e_all, fit$index, cf[["sigv"]], cf[["sigr"]],
                       fit$sigma_h, fit$sigma_u, verbose)
  fit$U <- te$U; fit$H <- te$H; fit$O <- te$O
  fit
}

## Marginal effects of the determinants on EXPECTED inefficiency:
##   E[u_it]  = sqrt(2/pi) * sigma_u_it  =>  d/dz_k = sqrt(2/pi)*sigma_u_it*delta_k/2
##   E[eta_i] = sqrt(2/pi) * sigma_h_i   =>  d/dw_k = sqrt(2/pi)*sigma_h_i *gamma_k/2
margeff_gtre_zh <- function(fit) {
  stopifnot(inherits(fit, "psfm_gtre_zh"))
  cf <- fit$coefficients
  delta <- cf[paste0("Zu:", fit$zu_names)]
  gam   <- cf[paste0("Zh:", fit$zh_names)]
  k <- sqrt(2 / pi)
  me_u <- outer(k * fit$sigma_u / 2, delta)      # n  x pu
  me_h <- outer(k * fit$sigma_h / 2, gam)        # N  x qh
  colnames(me_u) <- fit$zu_names
  colnames(me_h) <- fit$zh_names
  list(transient_AME  = colMeans(me_u),
       persistent_AME = colMeans(me_h),
       transient_full = me_u,
       persistent_full= me_h)
}
