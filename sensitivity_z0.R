# ============================================================
#  Sensitivity / uncertainty budget for
#  flx_compute_HS_HL_fixed_ustar_prior()
#
#  Run:  Rscript sensitivity_z0.R          (demo on synthetic u*)
#  Or:   source("sensitivity_z0.R"); flx_sensitivity_report(my_ustar, ...)
#
#  Produces two publication-ready tables:
#
#   (A) ENSEMBLE SPREAD across roughness models, stratified by wind.
#       This is the "which z0 model?" systematic uncertainty.  It is NOT a
#       single number -- it is ~0 at low wind and grows with u*.
#
#   (B) ELASTICITIES  dln(HS)/dln(x)  for x in {dt, u*, z0, zot}.
#       Fractional error in HS per fractional error in each input.  This is
#       what tells you where accuracy is worth buying.
#
#  ------------------------------------------------------------
#  WHY z0 MATTERS SO LITTLE HERE
#  ------------------------------------------------------------
#  With u* PRESCRIBED, z0 never enters the flux directly:
#
#      HS   = -rho*Cp*u* * t*
#      t*   = -kappa*dt / (log(zt/zot) - psi_h(zeta))
#      zeta = C1*(t* + C2*q*)/u*^2          <- contains no z0
#
#  The only path is  z0 -> rr = z0*u*/nu -> zot -> log(zt/zot).  So the
#  dependence is logarithmic:
#
#      dln(HS)/dln(zot) = 1 / ln(zt/zot) ~ 0.07-0.11
#
#  A factor-of-10 error in z0 costs only ~10-16% in HS.  Moreover zot is
#  capped at 1.6e-4 for rr < 0.244 (u* < ~0.121, U10N ~ 4), and in that
#  regime HS is EXACTLY independent of z0 -- every physical roughness model
#  returns the same flux.
#
#  ------------------------------------------------------------
#  HOW THE PERTURBATIONS ARE DONE
#  ------------------------------------------------------------
#  z0 is computed internally, so it cannot be perturbed by an argument.  But
#  because rr = z0*u*/nu is linear in z0, scaling rr inside scalar_rough_fn
#  scales z0 EXACTLY, and it preserves the cap (so the elasticity correctly
#  goes to zero where zot is pinned).  zot is perturbed by scaling the output
#  of the same hook.  dt is perturbed through T_s; u* through the input
#  directly, so that z0 is recomputed consistently -- which is the realistic
#  case, and is why the u* elasticity is sub-linear (~0.7) rather than 1.
#
#  NOTE on d_dt: it is obtained by perturbing T_s so that dt changes by eps.
#  That also moves Qs, so the number is the TOTAL response of HS to an SST
#  error (including the weak feedback through q* and zeta), normalised by the
#  fractional change in dt.  That is the operationally useful quantity -- it
#  is what a bias in your SST field actually costs you -- rather than a
#  partial derivative holding humidity fixed.
# ============================================================

if (!exists("flx_compute_HS_HL_fixed_ustar_prior"))
  source("bulk_algo_no_circularity_with_prior_v2.R")

.COARE_ZOT <- function(rr) pmin(1.6e-4, 5.8e-5 / pmax(rr, 1e-12)^0.72)

# scale z0 by `fac` (exact: rr is linear in z0, and the cap is preserved)
.rough_scale_z0 <- function(fac)
  function(rr, visa, u_star) { z <- .COARE_ZOT(fac * rr); list(zot = z, zoq = z) }

# scale zot/zoq by `fac`
.rough_scale_zot <- function(fac)
  function(rr, visa, u_star) { z <- fac * .COARE_ZOT(rr); list(zot = z, zoq = z) }

.default_bins <- c(0, 4, 6, 8, 10, 12, 15, 20, Inf)

.bin_label <- function(br) {
  n <- length(br) - 1L
  lab <- character(n)
  for (k in seq_len(n))
    lab[k] <- if (is.infinite(br[k + 1L])) sprintf(">%g", br[k])
              else sprintf("%g-%g", br[k], br[k + 1L])
  lab
}


# ------------------------------------------------------------------
# (A) Ensemble spread across roughness models, stratified by wind
# ------------------------------------------------------------------
flx_ensemble_spread <- function(u_star, T_s, T_a, RH, P,
                                models = c("charnock_winddep", "charnock_const",
                                           "eq3", "yelland"),
                                Hs = NA_real_, Tp = NA_real_,
                                bins = .default_bins,
                                verbose = TRUE, ...) {

  f <- flx_compute_HS_HL_fixed_ustar_prior

  # wave models need wave data; drop them rather than erroring
  have_waves <- any(is.finite(Hs)) && (any(is.finite(Tp)) || any(is.finite(list(...)$cp)))
  if (!have_waves && any(models %in% c("eq3", "yelland"))) {
    if (verbose)
      message("No finite Hs/Tp supplied -- dropping the wave-state models ",
              "(eq3, yelland) from the ensemble. Supply per-observation Hs and ",
              "Tp to include them; a constant sea state makes them misleading.")
    models <- setdiff(models, c("eq3", "yelland"))
  }
  if (length(models) < 2L)
    stop("need at least two roughness models to form a spread.", call. = FALSE)

  base <- suppressWarnings(f(u_star = u_star, T_s = T_s, T_a = T_a, RH = RH,
                             P = P, Hs = Hs, Tp = Tp,
                             z0_model = "charnock_winddep", ...))

  HS <- HL <- Z0 <- matrix(NA_real_, nrow = length(u_star), ncol = length(models),
                           dimnames = list(NULL, models))
  for (m in models) {
    r <- suppressWarnings(f(u_star = u_star, T_s = T_s, T_a = T_a, RH = RH,
                            P = P, Hs = Hs, Tp = Tp, z0_model = m, ...))
    HS[, m] <- ifelse(r$flux_valid, r$HS, NA_real_)
    HL[, m] <- ifelse(r$flux_valid, r$HL, NA_real_)
    Z0[, m] <- r$z0
  }

  # A percentage spread is meaningless where the mean flux passes through
  # zero (near-neutral air-sea dT), so exclude those rows rather than let
  # them produce enormous spurious percentages.
  HS_floor <- 1.0                                   # W/m2
  keep <- stats::complete.cases(HS) & is.finite(base$U10N_proxy) &
    abs(rowMeans(HS)) > HS_floor
  n_null <- sum(stats::complete.cases(HS) & abs(rowMeans(HS)) <= HS_floor)
  if (verbose && n_null > 0)
    message(n_null, " row(s) excluded from the spread: |HS| < ", HS_floor,
            " W/m2, where a percentage spread is not meaningful.")
  if (!any(keep)) stop("no rows valid across all models.", call. = FALSE)

  U   <- base$U10N_proxy[keep]
  HSk <- HS[keep, , drop = FALSE]; HLk <- HL[keep, , drop = FALSE]
  Z0k <- Z0[keep, , drop = FALSE]

  rng <- function(M) apply(M, 1, function(v) max(v) - min(v))
  mu  <- function(M) rowMeans(M)

  spread_HS <- 100 * rng(HSk) / abs(mu(HSk))
  spread_HL <- 100 * rng(HLk) / abs(mu(HLk))
  z0_ratio  <- apply(Z0k, 1, function(v) max(v) / min(v))

  grp <- cut(U, breaks = bins, right = FALSE, include.lowest = TRUE)
  lab <- .bin_label(bins)

  agg <- function(x, fun) vapply(split(x, grp), function(v)
    if (length(v)) fun(v) else NA_real_, 0)

  out <- data.frame(
    U10N_bin      = lab,
    n             = as.integer(table(grp)),
    HS_mean       = round(agg(mu(HSk), mean), 2),
    HS_spread_pct = round(agg(spread_HS, mean), 2),
    HS_spread_p95 = round(agg(spread_HS, function(v) stats::quantile(v, 0.95)), 2),
    HL_spread_pct = round(agg(spread_HL, mean), 2),
    z0_ratio      = round(agg(z0_ratio, stats::median), 1),
    row.names = NULL
  )
  attr(out, "models") <- models
  attr(out, "n_used") <- sum(keep)
  attr(out, "per_model_mean_HS") <- round(colMeans(HSk), 2)
  out
}


# ------------------------------------------------------------------
# (B) Elasticities  dln(HS)/dln(x)
# ------------------------------------------------------------------
flx_elasticities <- function(u_star, T_s, T_a, RH, P,
                             z0_model = "charnock_winddep",
                             Hs = NA_real_, Tp = NA_real_,
                             eps = 1e-3, bins = .default_bins, ...) {

  f  <- flx_compute_HS_HL_fixed_ustar_prior
  gg <- function(...) suppressWarnings(f(...))

  base <- gg(u_star = u_star, T_s = T_s, T_a = T_a, RH = RH, P = P,
             Hs = Hs, Tp = Tp, z0_model = z0_model, ...)

  # central difference in log space:  (ln y+ - ln y-) / (2 eps)
  el <- function(rp, rm) {
    a <- ifelse(rp$flux_valid, rp$HS, NA_real_)
    b <- ifelse(rm$flux_valid, rm$HS, NA_real_)
    (log(abs(a)) - log(abs(b))) / (2 * eps)
  }

  # dt: dt = T_s - T_a - 0.0098*zt, so scaling dt by (1+eps) means adding
  # eps*dt to T_s
  dt0 <- base$dt
  e_dt <- el(gg(u_star = u_star, T_s = T_s + eps * dt0, T_a = T_a, RH = RH,
                P = P, Hs = Hs, Tp = Tp, z0_model = z0_model, ...),
             gg(u_star = u_star, T_s = T_s - eps * dt0, T_a = T_a, RH = RH,
                P = P, Hs = Hs, Tp = Tp, z0_model = z0_model, ...))

  # u*: perturb the input so z0 is recomputed consistently (realistic case)
  e_us <- el(gg(u_star = u_star * (1 + eps), T_s = T_s, T_a = T_a, RH = RH,
                P = P, Hs = Hs, Tp = Tp, z0_model = z0_model, ...),
             gg(u_star = u_star * (1 - eps), T_s = T_s, T_a = T_a, RH = RH,
                P = P, Hs = Hs, Tp = Tp, z0_model = z0_model, ...))

  # z0: scale rr (exactly equivalent to scaling z0; cap preserved)
  e_z0 <- el(gg(u_star = u_star, T_s = T_s, T_a = T_a, RH = RH, P = P,
                Hs = Hs, Tp = Tp, z0_model = z0_model,
                scalar_rough_fn = .rough_scale_z0(1 + eps), ...),
             gg(u_star = u_star, T_s = T_s, T_a = T_a, RH = RH, P = P,
                Hs = Hs, Tp = Tp, z0_model = z0_model,
                scalar_rough_fn = .rough_scale_z0(1 - eps), ...))

  # zot: scale the scalar roughness itself
  e_zot <- el(gg(u_star = u_star, T_s = T_s, T_a = T_a, RH = RH, P = P,
                 Hs = Hs, Tp = Tp, z0_model = z0_model,
                 scalar_rough_fn = .rough_scale_zot(1 + eps), ...),
              gg(u_star = u_star, T_s = T_s, T_a = T_a, RH = RH, P = P,
                 Hs = Hs, Tp = Tp, z0_model = z0_model,
                 scalar_rough_fn = .rough_scale_zot(1 - eps), ...))

  # a log-elasticity is meaningless where HS ~ 0, so drop near-null-flux rows
  keep <- base$flux_valid & is.finite(base$U10N_proxy) & abs(base$HS) > 1e-3 &
    is.finite(e_dt) & is.finite(e_us) & is.finite(e_z0) & is.finite(e_zot)
  if (!any(keep)) stop("no rows usable for elasticities.", call. = FALSE)

  grp <- cut(base$U10N_proxy[keep], breaks = bins, right = FALSE,
             include.lowest = TRUE)
  agg <- function(x) vapply(split(x[keep], grp),
                            function(v) if (length(v)) mean(v) else NA_real_, 0)

  # fraction of rows where zot is pinned at its cap -> z0 cannot matter
  capped <- base$zot >= 1.6e-4 - 1e-12

  out <- data.frame(
    U10N_bin   = .bin_label(bins),
    n          = as.integer(table(grp)),
    d_dt       = round(agg(e_dt),  3),
    d_ustar    = round(agg(e_us),  3),
    d_zot      = round(agg(e_zot), 4),
    d_z0       = round(agg(e_z0),  4),
    zot_capped_pct = round(agg(as.numeric(capped)) * 100, 1),
    row.names = NULL
  )
  attr(out, "eps") <- eps
  attr(out, "z0_model") <- z0_model
  out
}


# ------------------------------------------------------------------
# Report
# ------------------------------------------------------------------
flx_sensitivity_report <- function(u_star, T_s, T_a, RH, P,
                                   Hs = NA_real_, Tp = NA_real_,
                                   bins = .default_bins,
                                   csv_prefix = NULL, ...) {

  rule <- function(ch = "=") cat(strrep(ch, 78), "\n")

  cat("\n"); rule()
  cat("  SENSITIVITY / UNCERTAINTY BUDGET\n")
  cat("  n =", length(u_star), "input rows\n")
  rule()

  A <- flx_ensemble_spread(u_star, T_s, T_a, RH, P, Hs = Hs, Tp = Tp,
                           bins = bins, ...)
  cat("\n(A) ROUGHNESS-MODEL ENSEMBLE SPREAD, by wind bin\n")
  cat("    models:", paste(attr(A, "models"), collapse = ", "), "\n")
  cat("    rows valid across all models:", attr(A, "n_used"), "\n")
  cat("    ('fixed' is excluded by default: it is not a physical model and is\n")
  cat("     the low-wind outlier for the wrong reason.)\n\n")
  print(A, row.names = FALSE)
  cat("\n    mean HS by model:\n")
  print(attr(A, "per_model_mean_HS"))

  B <- flx_elasticities(u_star, T_s, T_a, RH, P, Hs = Hs, Tp = Tp,
                        bins = bins, ...)
  cat("\n(B) ELASTICITIES  dln(HS)/dln(x)  -- HS error per unit input error\n\n")
  print(B, row.names = FALSE)

  cat("\n"); rule("-")
  cat("  READING THIS\n")
  rule("-")
  ok <- is.finite(B$d_dt) & B$n > 0
  if (any(ok)) {
    rd <- mean(B$d_dt[ok]); ru <- mean(B$d_ustar[ok]); rz <- mean(abs(B$d_z0[ok]))
    cat(sprintf("  air-sea dT : %.2f   -> a 10%% dT error is ~%.0f%% in HS\n",
                rd, 10 * rd))
    cat(sprintf("  u*         : %.2f   -> a 10%% u* error is ~%.0f%% in HS\n",
                ru, 10 * ru))
    cat(sprintf("  z0         : %.3f  -> a 10%% z0 error is ~%.1f%% in HS  (%.0fx less)\n",
                -rz, 10 * rz, rd / max(rz, 1e-9)))
  }
  cat("  z0 enters only through zot, inside a logarithm, and not at all where\n")
  cat("  zot is capped.  Buy accuracy in dT and u* before roughness.\n")
  rule("-")

  if (!is.null(csv_prefix)) {
    utils::write.csv(A, paste0(csv_prefix, "_ensemble_spread.csv"), row.names = FALSE)
    utils::write.csv(B, paste0(csv_prefix, "_elasticities.csv"),   row.names = FALSE)
    cat("  wrote ", csv_prefix, "_ensemble_spread.csv and _elasticities.csv\n", sep = "")
  }
  invisible(list(ensemble = A, elasticities = B))
}


# ------------------------------------------------------------------
# Demo (synthetic u*; replace with your CYGNSS record)
# ------------------------------------------------------------------
flx_sensitivity_demo <- function(n = 20000, seed = 1) {
  set.seed(seed)
  # roughly Weibull-ish U10N, mapped through the Andreas prior to u*
  U  <- pmin(rweibull(n, shape = 2.2, scale = 8.5), 30)
  us <- 0.239 + 0.0433 * ((U - 8.271) + sqrt(0.120 * (U - 8.271)^2 + 0.181))
  us <- pmax(us, 0.02)
  flx_sensitivity_report(
    u_star = us,
    T_s = 20, T_a = 17, RH = 80, P = 1013,
    Hs  = pmax(0.3, 0.025 * U^2),      # crude fetch-limited Hs(U)
    Tp  = pmax(3.0, 0.6  * U + 2),     # crude Tp(U)
    u_star_min = 0.06                  # drop the ill-posed free-convection tail
  )
}

if (identical(environment(), globalenv()) && !interactive()) flx_sensitivity_demo()
