# ============================================================
#  HS / HL given a PRESCRIBED u*  (diagnostic bulk algorithm)
#  ------------------------------------------------------------
#  v2 : cleaned / hardened / vectorised rewrite
#
#  Physics (unchanged from v1):
#    - COARE-like stability iteration for t*, q*, and L
#    - u* is treated as a FIXED INPUT; momentum outputs (tau, U10N_diag,
#      Umean_zu) are DIAGNOSTICS ONLY and must not be regressed against
#      the inputs used to build them.
#    - Wind-dependent Charnock alpha is obtained WITHOUT circularity:
#        observed u*  --(independent prior u*(U10N), inverted)-->  U10N_proxy
#        U10N_proxy   --(COARE-like linear law + caps)-->           alpha
#      No externally observed wind is required, and alpha is never tuned
#      from a wind that was itself derived from u* and z0.
#    - Below ~4 m/s the smooth/viscous regime dominates and a wind-dependent
#      alpha is not well constrained, so alpha is held at alpha_min (or
#      smoothly ramped).
#
#  Wave-state z0 options (deep-water assumption):
#    Deep water requires h >~ Lp/2, i.e. h_min(m) ~ 0.7804 * Tp^2
#      Tp = 12 s -> h_min ~ 112 m ;  Tp = 15 s -> h_min ~ 176 m
#
#  ------------------------------------------------------------
#  CHANGES vs v1 (all behaviour-preserving unless noted):
#   B1  bucksat(): Pmb was not recycled to length(TC), so a scalar P with a
#       vector T silently produced NA saturation vapour pressure for every
#       sub-freezing element after the first.  FIXED (rep_len).
#   B2  Gustiness (Ug, gf) can NOT influence HS/HL when u* is prescribed --
#       it has no pathway into t*/q*.  It is now computed ONCE after
#       convergence and clearly labelled diagnostic; it no longer pollutes
#       the convergence test.
#   B3  eq3 branch hard-stopped on missing Tp even when cp was supplied
#       directly, contradicting its own error message.  FIXED.
#   B4  u* below the prior's floor (Andreas: u*(0) ~ 0.0063) had no root and
#       silently produced all-NA rows.  Now clamped to [U10N_lo, U10N_hi]
#       with a reported flag (named `U10N_out_of_range` as of D12), or NA if
#       clamp_inversion = FALSE.
#   B5  `sign(den) * 1e-6` guard failed at exactly zero (sign(0) == 0).  FIXED.
#   B6  Convergence used max(..., na.rm=TRUE): NA rows were ignored so the
#       loop could declare success while rows were unconverged, and an
#       all-NA input gave -Inf plus a warning.  Now per-row, with `converged`
#       / `n_iter` columns and a warning when maxit is exhausted.
#   B7  Convergence metric now mixes an absolute floor with the relative
#       term (defensive; the pure ratio is scale-invariant in dt so it was
#       not itself a failure mode).
#   B8  Plain Picard iteration does NOT converge for strongly stable, very
#       low-u* rows (e.g. u* = 0.03, SST 10 C, Ta 22 C: zeta ~ 14, oscillates
#       for all 60 iterations).  v1 hid this behind na.rm; v2 detects it and
#       automatically re-runs only the stalled rows with under-relaxation
#       (relax_retry = 0.5), which converges them in ~8-20 iterations and
#       leaves already-converged rows untouched.
#   R1  Five mutually exclusive z0 flags were resolved by if-chain order,
#       silently overriding the user.  Replaced by z0_model (match.arg);
#       legacy flags still accepted, validated, and warned about.
#   R2  Per-observation inputs are length-checked and recycled explicitly.
#   R3  Non-finite rows no longer abort the call; they return NA rows.
#   R4  Guards on zu/zt/zq vs z0/zot/zoq.
#   R5  `c <- 5` no longer shadows base::c inside psiu_26.
#   P1  Prior inversion: per-element uniroot() loop -> one monotone lookup
#       table + approx().  ~100-1000x faster on long records.
#   P2  log(zt/zot), log(zq/zoq), log(zu/z0) and the zeta constants are
#       loop-invariant and are now precomputed.
#   P3  Umean_zu is no longer recomputed inside the loop.
#   N1  (WITHDRAWN)  An earlier v2 draft used Cp = 1004.67*(1 + 0.00084*Q).
#       That was wrong on two counts: COARE 3.5 builds HS from DRY-air
#       cpa = 1004.67 (cpv is a different quantity), and the 0.00084
#       coefficient presumes Q in g/kg whereas Q is kg/kg here -- so the
#       correction was ~1000x too small and did nothing anyway.  Cp is back
#       to the dry constant, matching v1 exactly.  If you do want the
#       moisture-adjusted value, set cp_moist = TRUE, which applies the
#       correct cpv = 1004.67*(1 + 0.84*Q).
#   N2  The viscous add-on is now available on the eq3 branch too
#       (add_eq3_viscous, default FALSE = v1 behaviour).
#   N3  Constant scalars (prior_model, ...) are attributes, not replicated
#       character columns.
#
#   D1  u* must now be STRICTLY positive (u_star_min, default 0).  u* = 0 was
#       previously admitted and produced a z0 that existed only because of the
#       pmax(u, 1e-6) floor in the viscous term.
#   D2  New QC columns: used_relaxation, final_residual.
#   D3  New column Rr = z0*u*/nu, the roughness Reynolds number, so the flow
#       regime can be filtered on directly.  Note that on the Charnock
#       branches Rr = beta_viscous + alpha*u*^3/(g*nu), i.e. Rr is bounded
#       BELOW by 0.11 -- the smooth limit is Rr -> 0.11, not Rr < 0.11.
#       Rough flow is Rr > 2.  With the Andreas prior the flow is not
#       aerodynamically rough until U10N ~ 9 m/s (u* ~ 0.29), and the viscous
#       term dominates z0 below U10N ~ 4.3 m/s -- essentially where
#       U_alpha_start = 4.0 already sits.
#   D4  attr(out, "qc") roll-up + flx_qc_summary(out), including the clamped
#       out-of-range fraction (superseded by D12: see f_out_of_range_pct,
#       taken over n_prior_evaluated).
#   D5  Converged rows are now FROZEN.  Previously every row was advanced on
#       every sweep, so a row that met tol at iteration 8 kept being updated
#       until the slowest row in the batch finished, while still reporting
#       n_iter = 8 and the iteration-8 residual.  The returned (t*, q*) now
#       corresponds exactly to the state at which the reported convergence
#       criterion was met.
#   D6  A non-converged solve no longer becomes a flux.  HS/HL are set to NA
#       for those rows and a new logical column `flux_valid` marks them;
#       tstar/qstar/n_iter/final_residual are retained for debugging.
#       Downstream: subset(out, flux_valid).  tau is unaffected (it does not
#       depend on the solver).
#   D8  The blanking extends to EVERY solver-dependent diagnostic, not just
#       HS/HL: L, zeta_zu, Ug, gf and Umean_zu are all derived from the
#       iterated (t*, q*) and would otherwise look like usable numbers.
#       Retained on failed rows, explicitly as debug state: tstar, qstar,
#       n_iter, final_residual, used_relaxation.  Retained as genuinely
#       solver-independent: tau, z0, zot, zoq, Rr, alpha_used, U10N_proxy,
#       U10N_diag, rho_a, Q, Qs, dt, dq.
#   D9  qc$max_residual no longer returns -Inf when every residual is
#       non-finite; it returns NA.
#   D11 The convergence residual is now measured on the RAW fixed-point step
#       ||F(x) - x||, before under-relaxation.  Previously it was taken from
#       the relaxed update x_new - x = lambda*(F(x) - x), so with
#       relax_retry = 0.5 the reported residual was half the true one and the
#       convergence criterion depended on lambda -- final_residual was not
#       comparable between the undamped and damped passes.  Relaxation now
#       affects only the update, never the test.
#   D13 n_iter counts iterations PERFORMED and is now populated on failed rows
#       (previously NA on failure while final_residual was filled in, which is
#       inconsistent for the two fields whose whole job is diagnosing a
#       failure).  For a converged row the value is unchanged.  When
#       used_relaxation is TRUE, n_iter refers to the damped pass.
#   D12 U10N_clamped is renamed U10N_out_of_range and is NEVER erased.
#       Previously it was blanked wherever ok == FALSE, which with
#       clamp_inversion = FALSE deleted the flag on exactly the rows the flag
#       existed to identify (out-of-range -> U10N_proxy = NA -> row invalid ->
#       flag erased), so the QC count could report zero.  QC now separates
#       n_prior_evaluated / n_prior_out_of_range / n_prior_clamped, with
#       f_out_of_range_pct taken over rows the inversion was attempted on.
#   D10 zu, zt and zq accept length-1 OR length-n, like every other
#       per-observation input, for records with instrument-specific
#       measurement heights.  Scalars are bit-identical to before.  The
#       stability argument is now zeta_zu * (zt/zu) with per-row ratios, and
#       attr(out, "heights") reports min/max rather than a single triple.
#   D7  QC now separates n_solver_rows from n_flux_valid, with
#       n_flux_valid = n_solver_rows - n_nonconverged, so a failed solve
#       cannot be mistaken for a usable one in a summary table.  The Rr
#       counter is renamed n_near_smooth_limit (see D3: Rr < 0.11 is
#       unreachable, so "smooth flow" was the wrong label).
#
#  NOT included (deliberate omissions, as in v1):
#    no Webb correction on HL, no cool-skin / warm-layer, no rain heat flux.
#
#  Outputs: HS, HL upward positive (ocean -> atmosphere).
# ============================================================


# ------------------------------------------------------------------
# Internal helpers (top level: defined once, not re-created per call)
# ------------------------------------------------------------------

# Gravity at latitude (Moritz 1980 normal gravity)
.flx_grv <- function(lat) {
  gamma <- 9.7803267715
  c1 <- 0.0052790414; c2 <- 0.0000232718
  c3 <- 0.0000001262; c4 <- 0.0000000007
  x <- sin(lat * pi / 180)
  gamma * (1 + c1 * x^2 + c2 * x^4 + c3 * x^6 + c4 * x^8)
}

# Kinematic viscosity of air [m2/s], TC in degC
.flx_nu_air <- function(TC) {
  1.326e-5 * (1 + 6.542e-3 * TC + 8.301e-6 * TC^2 - 4.84e-9 * TC^3)
}

# Buck saturation vapour pressure [mb]; Tf = freezing point [degC]
# B1: every argument is recycled to a common length before indexing.
.flx_bucksat <- function(TC, Pmb, Tf) {
  n   <- max(length(TC), length(Pmb), length(Tf))
  TC  <- rep_len(TC,  n)
  Pmb <- rep_len(Pmb, n)
  Tf  <- rep_len(Tf,  n)

  exx <- 6.1121 * exp(17.502 * TC / (TC + 240.97)) * (1.0007 + 3.46e-6 * Pmb)

  ii <- which(TC < Tf)
  if (length(ii) > 0L) {
    exx[ii] <- (1.0003 + 4.18e-6 * Pmb[ii]) * 6.1115 *
      exp(22.452 * TC[ii] / (TC[ii] + 272.55))
  }
  exx
}

# Specific humidity of air [kg/kg] from RH [%]
.flx_qsat_air <- function(TC, Pmb, rh) {
  es <- .flx_bucksat(TC, Pmb, 0)
  em <- 0.01 * rh * es
  (622 * em / (Pmb - 0.378 * em)) / 1000
}

# Saturation specific humidity at the sea surface [kg/kg]
.flx_qsat_sea <- function(TC, Pmb, S) {
  Tf <- -0.0575 * S + 1.71052e-3 * S^1.5 - 2.154996e-4 * S^2
  es <- (1 - 0.02 * S / 35) * .flx_bucksat(TC, Pmb, Tf)
  (622 * es / (Pmb - 0.378 * es)) / 1000
}

# COARE 3.0 stability functions ------------------------------------
.flx_psit_26 <- function(zet) {
  psi <- rep(NA_real_, length(zet))
  ok  <- is.finite(zet)
  st  <- ok & zet >= 0
  un  <- ok & zet <  0

  if (any(st)) {
    z <- zet[st]
    dzet <- pmin(50, 0.35 * z)
    psi[st] <- -((1 + 0.6667 * z)^1.5 + 0.6667 * (z - 14.28) * exp(-dzet) + 8.525)
  }

  if (any(un)) {
    z    <- zet[un]
    x1   <- (1 - 15 * z)^0.5
    psik <- 2 * log((1 + x1) / 2)

    x2   <- (1 - 34.15 * z)^(1/3)
    psic <- 1.5 * log((1 + x2 + x2^2) / 3) -
      sqrt(3) * atan((1 + 2 * x2) / sqrt(3)) + pi / sqrt(3)

    f <- z^2 / (1 + z^2)
    psi[un] <- (1 - f) * psik + f * psic
  }
  psi
}

.flx_psiu_26 <- function(zet) {
  psi <- rep(NA_real_, length(zet))
  ok  <- is.finite(zet)
  st  <- ok & zet >= 0
  un  <- ok & zet <  0

  if (any(st)) {
    # R5: renamed from a,b,c,d -- `c` shadowed base::c
    pa <- 0.7; pb <- 3/4; pc <- 5; pd <- 0.35
    z <- zet[st]
    dzet <- pmin(50, 0.35 * z)
    psi[st] <- -(pa * z + pb * (z - pc / pd) * exp(-dzet) + pb * pc / pd)
  }

  if (any(un)) {
    z    <- zet[un]
    x1   <- (1 - 15 * z)^0.25
    psik <- 2 * log((1 + x1) / 2) + log((1 + x1^2) / 2) - 2 * atan(x1) + pi / 2

    x2   <- (1 - 10.15 * z)^(1/3)
    psic <- 1.5 * log((1 + x2 + x2^2) / 3) -
      sqrt(3) * atan((1 + 2 * x2) / sqrt(3)) + pi / sqrt(3)

    f <- z^2 / (1 + z^2)
    psi[un] <- (1 - f) * psik + f * psic
  }
  psi
}

# P1: invert a monotone 1-D function by lookup table instead of per-element
# uniroot().  `rule = 2` in approx() clamps outside the table, which is what
# fixes B4 (u* below the prior's floor).
.flx_invert_monotone <- function(y_target, fun, lo, hi, n_grid = 4001L,
                                 clamp = TRUE) {
  xg <- seq(lo, hi, length.out = n_grid)
  yg <- fun(xg)

  keep <- is.finite(xg) & is.finite(yg)
  xg <- xg[keep]; yg <- yg[keep]
  degenerate <- list(x = rep(NA_real_, length(y_target)),
                     out_of_range = rep(NA, length(y_target)),
                     range = c(NA_real_, NA_real_))
  if (length(xg) < 2L) return(degenerate)

  if (is.unsorted(yg, strictly = TRUE)) {
    ym    <- cummax(yg)
    keep2 <- !duplicated(ym)
    xg <- xg[keep2]; yg <- ym[keep2]
    warning("prior u*(U10N) is not strictly increasing on [", lo, ", ", hi,
            "]; inversion used its monotone envelope.", call. = FALSE)
    if (length(xg) < 2L) return(degenerate)
  }

  # D12: NA (= "not evaluated") is distinct from FALSE (= "evaluated, in
  # range").  Rows whose target was never finite must not be reported as
  # in-range.
  y_lo <- yg[1L]; y_hi <- yg[length(yg)]
  out_of_range <- ifelse(is.finite(y_target),
                         y_target < y_lo | y_target > y_hi,
                         NA)

  x <- approx(x = yg, y = xg, xout = y_target,
              rule = if (clamp) 2 else 1)$y
  list(x = x, out_of_range = out_of_range, range = c(y_lo, y_hi))
}


# ------------------------------------------------------------------
# Main
# ------------------------------------------------------------------
flx_compute_HS_HL_fixed_ustar_prior <- function(
  u_star, T_s, T_a, RH, P,
  # measurement heights [m]; scalar or length-n (D10).  Use length-n when the
  # record carries instrument-specific heights (buoy/ship); a scalar is the
  # right choice when observations are already standardized to fixed levels.
  zt = 10, zq = 10, zu = 10,
  S   = 35,
  lat = 45,
  zi  = 600,
  Beta = 1.2,

  # ---- roughness model (R1: one enum replaces five booleans) ----
  #   "charnock_winddep" : alpha(U10N_proxy) from inverted prior  [default]
  #   "charnock_const"   : alpha = alpha_charnock
  #   "fixed"            : z0 = z0_fixed
  #   "eq3"              : Drennan-type  z0 = 3.35 * Hs * (u*/cp)^3.4
  #   "yelland"          : Taylor-Yelland z0 = 1200 * Hs * (Hs/Lp)^4.5
  z0_model = c("charnock_winddep", "charnock_const", "fixed", "eq3", "yelland"),

  z0_fixed       = 1.0e-4,
  alpha_charnock = 0.011,
  beta_viscous   = 0.11,

  # ---- wave inputs (eq3 / yelland) ----
  Hs = NA_real_,   # significant wave height [m]
  Tp = NA_real_,   # peak period [s]
  Lp = NA_real_,   # peak wavelength [m]  (optional)
  cp = NA_real_,   # phase speed [m/s]    (optional; overrides Lp/Tp)
  eq3_coef = 3.35, eq3_exp = 3.4,
  YL_coef  = 1200, YL_exp  = 4.5,
  use_deepwater_cp_if_missing = TRUE,
  use_deepwater_Lp_if_missing = TRUE,
  add_YL_viscous  = TRUE,
  add_eq3_viscous = FALSE,   # N2: available now; FALSE reproduces v1

  # ---- prior u*(U10N) ----
  prior_model = c("andreas", "edson2013", "custom"),
  prior_custom_fun = NULL,
  andreas_ustar_fun = function(U10) {
    0.239 + 0.0433 * ((U10 - 8.271) + sqrt(0.120 * (U10 - 8.271)^2 + 0.181))
  },
  edson_U2 = 8.5, edson_blend_width = 0.6,
  edson_slope_low = 0.035, edson_slope_high = 0.062,
  edson_intercept_high = -0.28,

  # ---- lower bound on admissible u* ----
  # D1: u* must be STRICTLY positive.  At u* = 0 the viscous roughness term
  # 0.11*nu/u* is only finite because of the pmax(u, 1e-6) floor, so a zero
  # retrieval would yield a purely numerical z0 rather than a physical one.
  # Raise this to the valid lower bound of your retrieval if one exists.
  # For orientation, with the Andreas prior (nu ~ 1.5e-5):
  #     u* = 0.0063  -> U10N = 0     (the prior's own floor; below this the
  #                                   inversion can only clamp)
  #     u* = 0.063   -> U10N = 2     viscous z0 still ~6x the Charnock z0
  #     u* = 0.092   -> U10N = 3     viscous ~2x Charnock
  #     u* = 0.13    -> U10N ~ 4.3   viscous == Charnock  <-- smooth/rough
  #                                   crossover, and essentially the same
  #                                   threshold as U_alpha_start = 4.0
  #     u* = 0.35    -> U10N = 10    Rr ~ 3.5, first aerodynamically rough
  # Rows at or below u_star_min are returned as NA, not silently floored.
  #
  # SOLVABILITY, which is the sharper argument for a non-zero bound:
  # with u* prescribed, zeta = (kappa*g*zu/Ta)*(t* + 0.61*Ta*q*)/u*^2 scales
  # as 1/u*^2, so it diverges as u* -> 0.  Measured for SST 20 / Ta 17 C:
  #     u* = 0.004  |zeta| ~ 1.5e4   NEVER converges (tried relax 1.0 -> 0.02)
  #     u* = 0.010  |zeta| ~ 550     converges, 30 iterations
  #     u* = 0.030  |zeta| ~ 36      converges, 25 iterations
  #     u* = 0.060  |zeta| ~ 7.4     converges, 23 iterations
  #     u* = 0.130  |zeta| ~ 1.3     converges, 21 iterations
  # This is the free-convection limit, where MO similarity at fixed u* is
  # ill-posed -- not a solver deficiency, and not fixable by damping.  Such
  # rows are correctly rejected by D6 (converged = FALSE, flux_valid = FALSE).
  # Setting u_star_min ~ 0.06 (U10N ~ 2) removes them up front instead.
  u_star_min = 0,

  U10N_lo = 0, U10N_hi = 60,
  inversion_grid_n = 4001L,
  clamp_inversion  = TRUE,   # B4: FALSE -> out-of-range u* gives NA instead

  # ---- alpha(U10N) law ----
  alpha_slope = 0.0017, alpha_intercept = -0.0050,
  alpha_min = 0.011, alpha_max = 0.028,
  U10_cap_for_alpha = 19,
  alpha_lowwind_mode = c("hold", "ramp", "none"),
  U_alpha_start = 4.0, U_alpha_end = 6.0,

  # ---- scalar roughness ----
  use_coare_scalar = TRUE,
  scalar_rough_fn  = NULL,   # function(rr, visa, u_star) -> list(zot=, zoq=)

  # ---- specific heat ----
  # FALSE (default, = COARE 3.5 and = v1): dry cpa = 1004.67 J/kg/K.
  # TRUE: moisture-adjusted cpv = cpa * (1 + 0.84 * Q), Q in kg/kg (~ +0.8%).
  cp_moist = FALSE,

  # ---- gustiness (B2: DIAGNOSTIC ONLY when u* is prescribed) ----
  use_gustiness = TRUE,
  Ug_min = 0.20, Ug_max = 2.0,

  # ---- iteration ----
  maxit   = 60,
  tol     = 1e-6,
  tol_abs = 1e-8,            # absolute floor in the convergence metric
  relax_retry = 0.5,         # B8: damped second pass for rows that stall
                             #     (set NULL to disable)
  warn_nonconverged = TRUE,

  # ---- deprecated v1 flags (accepted, validated, warned) ----
  use_fixed_z0 = NULL, use_eq3_z0 = NULL, use_coare35_YL_z0 = NULL,
  use_charnock_constant_alpha = NULL, use_charnock_winddep_alpha = NULL
) {

  z0_model           <- match.arg(z0_model)
  prior_model        <- match.arg(prior_model)
  alpha_lowwind_mode <- match.arg(alpha_lowwind_mode)

  # -- R1: legacy flag shim -----------------------------------------
  legacy <- c(fixed            = isTRUE(use_fixed_z0),
              eq3              = isTRUE(use_eq3_z0),
              yelland          = isTRUE(use_coare35_YL_z0),
              charnock_const   = isTRUE(use_charnock_constant_alpha),
              charnock_winddep = isTRUE(use_charnock_winddep_alpha))
  if (any(legacy)) {
    if (sum(legacy) > 1L) {
      stop("More than one legacy use_* roughness flag is TRUE (",
           paste(names(legacy)[legacy], collapse = ", "),
           "). These are mutually exclusive -- use z0_model instead.",
           call. = FALSE)
    }
    z0_model <- names(legacy)[legacy]
    warning("The use_* roughness flags are deprecated; use ",
            "z0_model = \"", z0_model, "\".", call. = FALSE)
  }

  # -- R2: normalise every per-observation input --------------------
  n <- length(u_star)
  if (n == 0L) stop("u_star has length 0.", call. = FALSE)

  fit <- function(x, nm) {
    lx <- length(x)
    if (lx == n) return(as.numeric(x))
    if (lx == 1L) return(rep_len(as.numeric(x), n))
    stop(sprintf("`%s` has length %d; expected 1 or %d (= length(u_star)).",
                 nm, lx, n), call. = FALSE)
  }
  u_star <- as.numeric(u_star)
  T_s <- fit(T_s, "T_s"); T_a <- fit(T_a, "T_a")
  RH  <- fit(RH,  "RH");  P   <- fit(P,   "P")
  S   <- fit(S,   "S");   lat <- fit(lat, "lat")
  Hs  <- fit(Hs,  "Hs");  Tp  <- fit(Tp,  "Tp")
  Lp  <- fit(Lp,  "Lp");  cp  <- fit(cp,  "cp")

  # D10: measurement heights are per-observation inputs like everything else.
  # Passing scalars is unchanged and gives bit-identical results.
  zu  <- fit(zu,  "zu");  zt  <- fit(zt,  "zt");  zq <- fit(zq, "zq")

  stopifnot(maxit >= 1L, tol > 0)
  if (any(!is.finite(zu) | zu <= 0 |
          !is.finite(zt) | zt <= 0 |
          !is.finite(zq) | zq <= 0)) {
    stop("zu, zt and zq must all be finite and > 0.", call. = FALSE)
  }

  # -- R3/D1: valid-row mask; bad rows return NA rather than aborting ---
  stopifnot(length(u_star_min) == 1L, is.finite(u_star_min), u_star_min >= 0)

  ok <- is.finite(u_star) & u_star > u_star_min &
    is.finite(T_s) & is.finite(T_a) & is.finite(RH) & is.finite(P) &
    is.finite(S) & is.finite(lat)

  n_rejected_u <- sum(is.finite(u_star) & u_star <= u_star_min)
  if (n_rejected_u > 0L) {
    warning(n_rejected_u, " row(s) at or below the configured admissible u* ",
            "minimum (u_star_min = ", u_star_min, ") were returned as NA.",
            call. = FALSE)
  }

  na_out <- function() {
    data.frame(
      u_star = u_star,
      HS = NA_real_, HL = NA_real_, tau = NA_real_,
      tstar = NA_real_, qstar = NA_real_, L = NA_real_, zeta_zu = NA_real_,
      Ug = NA_real_, gf = NA_real_,
      z0 = NA_real_, zot = NA_real_, zoq = NA_real_,
      Rr = NA_real_,
      alpha_used = NA_real_, U10N_proxy = NA_real_, U10N_out_of_range = NA,
      Umean_zu = NA_real_, U10N_diag = NA_real_,
      rho_a = NA_real_, Q = NA_real_, Qs = NA_real_,
      dt = NA_real_, dq = NA_real_,
      converged = NA, n_iter = NA_integer_,
      used_relaxation = NA, final_residual = NA_real_,
      flux_valid = NA
    )
  }

  # ---------------- constants & thermodynamics (valid rows) --------
  kappa <- 0.4; tdk <- 273.15
  grav  <- .flx_grv(lat)

  Pmb <- ifelse(P > 2000, P / 100, P)      # mb
  Ppa <- ifelse(P > 2000, P,      P * 100) # Pa

  Q  <- .flx_qsat_air(T_a, Pmb, RH)
  Qs <- .flx_qsat_sea(T_s, Pmb, S)

  rho_a <- Ppa / (287.05 * (T_a + tdk) * (1 + 0.61 * Q))
  visa  <- .flx_nu_air(T_a)
  Le    <- (2.501 - 0.00237 * T_s) * 1e6

  # Specific heat for HS.  COARE 3.5 uses DRY-air cpa = 1004.67 J/kg/K in the
  # sensible heat flux; cpv = cpa*(1 + 0.84*Q) is a separate quantity and is
  # not what HS is built from.  Note the coefficient is 0.84 for Q in kg/kg
  # (0.00084 belongs to formulations carrying Q in g/kg) -- Q is kg/kg here.
  Cp <- if (isTRUE(cp_moist)) 1004.67 * (1 + 0.84 * Q) else rep(1004.67, n)

  ok <- ok & is.finite(Q) & is.finite(Qs) & is.finite(rho_a) & visa > 0

  # ---------------- z0 ---------------------------------------------
  z0                <- rep(NA_real_, n)
  alpha_used        <- rep(NA_real_, n)
  U10N_proxy        <- rep(NA_real_, n)
  U10N_out_of_range <- rep(NA, n)
  prior_ustar_range <- c(NA_real_, NA_real_)

  viscous <- function(u) beta_viscous * visa / pmax(u, 1e-6)

  # deep-water helpers ------------------------------------------------
  resolve_cp <- function() {
    v <- rep(NA_real_, n)
    i <- is.finite(cp) & cp > 0
    v[i] <- cp[i]
    i <- is.na(v) & is.finite(Lp) & Lp > 0 & is.finite(Tp) & Tp > 0
    v[i] <- Lp[i] / Tp[i]
    if (use_deepwater_cp_if_missing) {
      i <- is.na(v) & is.finite(Tp) & Tp > 0
      v[i] <- grav[i] * Tp[i] / (2 * pi)          # cp = g*Tp/(2*pi)
    }
    v
  }
  resolve_Lp <- function() {
    v <- rep(NA_real_, n)
    i <- is.finite(Lp) & Lp > 0
    v[i] <- Lp[i]
    if (use_deepwater_Lp_if_missing) {
      i <- is.na(v) & is.finite(Tp) & Tp > 0
      v[i] <- grav[i] * Tp[i]^2 / (2 * pi)        # Lp = g*Tp^2/(2*pi)
    }
    v
  }

  if (z0_model == "fixed") {
    stopifnot(length(z0_fixed) == 1L, is.finite(z0_fixed), z0_fixed > 0)
    z0[] <- z0_fixed

  } else if (z0_model == "eq3") {
    # B3: no premature Tp check -- cp alone is sufficient.
    if (all(!is.finite(Hs) | Hs <= 0))
      stop('z0_model = "eq3" requires Hs > 0 [m].', call. = FALSE)
    cp_use <- resolve_cp()
    if (all(!is.finite(cp_use) | cp_use <= 0))
      stop('z0_model = "eq3": need cp, or Lp and Tp, or Tp with ',
           'use_deepwater_cp_if_missing = TRUE.', call. = FALSE)

    good <- is.finite(Hs) & Hs > 0 & is.finite(cp_use) & cp_use > 0
    z0[good] <- eq3_coef * Hs[good] *
      (pmax(u_star[good], 1e-12) / cp_use[good])^eq3_exp
    if (isTRUE(add_eq3_viscous)) z0[good] <- z0[good] + viscous(u_star)[good]
    ok <- ok & good

  } else if (z0_model == "yelland") {
    if (all(!is.finite(Hs) | Hs <= 0))
      stop('z0_model = "yelland" requires Hs > 0 [m].', call. = FALSE)
    Lp_use <- resolve_Lp()
    if (all(!is.finite(Lp_use) | Lp_use <= 0))
      stop('z0_model = "yelland": need Lp [m], or Tp [s] with ',
           'use_deepwater_Lp_if_missing = TRUE.', call. = FALSE)

    good <- is.finite(Hs) & Hs > 0 & is.finite(Lp_use) & Lp_use > 0
    steep <- Hs[good] / Lp_use[good]
    z0[good] <- YL_coef * Hs[good] * pmax(steep, 1e-12)^YL_exp
    if (isTRUE(add_YL_viscous)) z0[good] <- z0[good] + viscous(u_star)[good]
    ok <- ok & good

  } else if (z0_model == "charnock_const") {
    alpha_used[] <- alpha_charnock
    z0 <- alpha_used * u_star^2 / grav + viscous(u_star)

  } else {  # "charnock_winddep"
    edson_ustar_fun <- function(U10N) {
      U <- pmax(U10N, 0)
      w <- 1 / (1 + exp(-(U - edson_U2) / max(edson_blend_width, 1e-6)))
      (1 - w) * (edson_slope_low * U) +
        w * (edson_slope_high * U + edson_intercept_high)
    }
    ustar_prior_fun <- switch(
      prior_model,
      andreas   = andreas_ustar_fun,
      edson2013 = edson_ustar_fun,
      custom    = {
        if (!is.function(prior_custom_fun))
          stop('prior_model = "custom" requires prior_custom_fun.', call. = FALSE)
        prior_custom_fun
      }
    )

    inv <- .flx_invert_monotone(
      y_target = ifelse(ok, u_star, NA_real_),
      fun = ustar_prior_fun, lo = U10N_lo, hi = U10N_hi,
      n_grid = inversion_grid_n, clamp = clamp_inversion
    )
    U10N_proxy        <- inv$x
    U10N_out_of_range <- inv$out_of_range
    prior_ustar_range <- inv$range

    if (any(U10N_out_of_range, na.rm = TRUE)) {
      warning(sum(U10N_out_of_range, na.rm = TRUE), " of ", sum(ok),
              " u* values fall outside the range the prior can represent on ",
              "U10N in [", U10N_lo, ", ", U10N_hi, "] (u* in [",
              signif(prior_ustar_range[1], 4), ", ",
              signif(prior_ustar_range[2], 4), "]); ",
              if (clamp_inversion)
                "U10N_proxy was clamped to the endpoints and the rows kept."
              else "U10N_proxy set to NA and the rows rejected.",
              call. = FALSE)
    }

    # alpha(U10N) with caps + low-wind handling
    Uc    <- pmin(U10N_proxy, U10_cap_for_alpha)
    a_law <- pmin(alpha_max, pmax(alpha_min, alpha_slope * Uc + alpha_intercept))

    alpha_used <- switch(
      alpha_lowwind_mode,
      none = a_law,
      hold = ifelse(U10N_proxy < U_alpha_start, alpha_min, a_law),
      ramp = {
        U1 <- U_alpha_start; U2 <- max(U_alpha_end, U1 + 1e-6)
        x  <- pmin(1, pmax(0, (U10N_proxy - U1) / (U2 - U1)))
        alpha_min + (x^2 * (3 - 2 * x)) * (a_law - alpha_min)
      }
    )

    z0 <- alpha_used * u_star^2 / grav + viscous(u_star)
    ok <- ok & is.finite(z0)
  }

  z0 <- pmax(z0, 1e-10)
  ok <- ok & is.finite(z0)      # keeps `ok` strictly TRUE/FALSE (never NA)

  # R4: the log-law is meaningless if the measurement height is inside z0
  if (any(ok & zu <= z0, na.rm = TRUE)) {
    warning(sum(ok & zu <= z0, na.rm = TRUE),
            " row(s) have zu <= z0; returned as NA.", call. = FALSE)
    ok <- ok & zu > z0
  }

  # ---------------- scalar roughness -------------------------------
  rr <- z0 * u_star / pmax(visa, 1e-12)

  if (is.function(scalar_rough_fn)) {
    tmp <- scalar_rough_fn(rr = rr, visa = visa, u_star = u_star)
    zot <- rep_len(as.numeric(tmp$zot), n)
    zoq <- rep_len(as.numeric(tmp$zoq), n)
  } else if (use_coare_scalar) {
    zoq <- pmin(1.6e-4, 5.8e-5 / pmax(rr, 1e-12)^0.72)
    zot <- zoq
  } else {
    zoq <- rep(1e-5, n); zot <- rep(1e-5, n)
  }
  zot <- pmax(zot, 1e-10); zoq <- pmax(zoq, 1e-10)
  ok  <- ok & zt > zot & zq > zoq & is.finite(zot) & is.finite(zoq)

  # ---------------- gradients --------------------------------------
  dt  <- T_s - T_a - 0.0098 * zt        # COARE potential-temperature lapse
  dq  <- Qs - Q
  taK <- T_a + tdk
  ok  <- ok & is.finite(dt) & is.finite(dq)

  if (!any(ok)) {
    warning("no valid rows; returning all NA.", call. = FALSE)
    out <- na_out()
    attr(out, "z0_model") <- z0_model
    attr(out, "prior_model") <- prior_model
    attr(out, "alpha_lowwind_mode") <- alpha_lowwind_mode
    attr(out, "u_star_min") <- u_star_min
    attr(out, "qc") <- list(
      n_input = n, n_rejected_ustar = n_rejected_u,
      n_solver_rows = 0L, n_flux_valid = 0L, n_nonconverged = 0L,
      n_relaxed = 0L,
      n_prior_evaluated = sum(!is.na(U10N_out_of_range)),
      n_prior_out_of_range = sum(U10N_out_of_range, na.rm = TRUE),
      n_prior_clamped = 0L, f_out_of_range_pct = NA_real_,
      prior_ustar_range = prior_ustar_range,
      max_residual = NA_real_,
      n_near_smooth_limit = 0L, n_rough_flow = 0L
    )
    return(out)
  }

  # ---------------- iterate on the valid subset --------------------
  i  <- which(ok)
  us <- u_star[i]

  # P2: loop-invariant quantities  (D10: heights are per-row)
  zu_i  <- zu[i]; zt_i <- zt[i]; zq_i <- zq[i]
  log_t <- log(zt_i / zot[i])
  log_q <- log(zq_i / zoq[i])
  log_u <- log(zu_i / z0[i])
  C1    <- kappa * grav[i] * zu_i / taK[i]  # zeta = C1*(t* + C2*q*)/u*^2
  C2    <- 0.61 * taK[i]
  rt    <- zt_i / zu_i                      # zeta at zt = zeta_zu * (zt/zu)
  rq    <- zq_i / zu_i
  us2   <- pmax(us^2, 1e-12)
  dt_i  <- dt[i]; dq_i <- dq[i]

  # B5: magnitude floor that also works at exactly zero
  guard <- function(x, eps = 1e-6) {
    s <- ifelse(x < 0, -1, 1)          # sign(0) == 0 would defeat the guard
    s * pmax(abs(x), eps)
  }

  # One fixed-point sweep over the rows selected by `sel`, with optional
  # under-relaxation.  relax = 1 is the plain Picard iteration of v1.
  iterate <- function(sel, relax, t0, q0) {
    lt <- log_t[sel]; lq <- log_q[sel]
    c1 <- C1[sel];    c2 <- C2[sel];   u2 <- us2[sel]
    dts <- dt_i[sel]; dqs <- dq_i[sel]
    r_t <- rt[sel];   r_q <- rq[sel]           # D10: per-row height ratios

    ts <- t0; qs <- q0
    cv <- rep(FALSE, length(sel))
    ni <- rep(NA_integer_, length(sel))
    rs <- rep(NA_real_, length(sel))   # D2: final residual

    for (it in seq_len(maxit)) {
      zeta <- c1 * (ts + c2 * qs) / u2

      # raw fixed-point map F(x)
      tn_raw <- -kappa * dts / guard(lt - .flx_psit_26(zeta * r_t))
      qn_raw <- -kappa * dqs / guard(lq - .flx_psit_26(zeta * r_q))

      # D11: the residual is ||F(x) - x||, measured on the RAW step, BEFORE
      # under-relaxation.  Measuring the relaxed step x_new - x = lambda*(F(x)-x)
      # would scale the residual by lambda, making the convergence criterion
      # depend on the relaxation coefficient and rendering final_residual
      # incomparable between the undamped and damped passes.
      # B6/B7: per-row, mixed absolute + relative
      d <- pmax(abs(tn_raw - ts) / (abs(ts) + tol_abs),
                abs(qn_raw - qs) / (abs(qs) + tol_abs))

      # under-relaxation affects only the UPDATE, never the test
      if (relax != 1) {
        tn <- ts + relax * (tn_raw - ts)
        qn <- qs + relax * (qn_raw - qs)
      } else {
        tn <- tn_raw; qn <- qn_raw
      }

      # D5: FREEZE converged rows.  Updating every row on every sweep meant a
      # row that met tol at iteration 8 kept being advanced until the slowest
      # row in the batch finished, while still reporting n_iter = 8 and the
      # residual from iteration 8.  Now the returned (t*, q*) is exactly the
      # state at which the stored convergence criterion was satisfied.
      active <- !cv
      ts[active] <- tn[active]
      qs[active] <- qn[active]
      rs[active] <- d[active]
      # D13: n_iter counts iterations PERFORMED, so it is populated for failed
      # rows too (= the iteration the pass stopped at).  Previously it was set
      # only on success, leaving n_iter = NA while final_residual was filled in
      # -- inconsistent for the two fields that exist to diagnose a failure.
      # For a converged row this is still the iteration it converged at, since
      # it stops being active at that point.
      ni[active] <- it

      done <- active & is.finite(d) & d < tol
      cv[done] <- TRUE

      if (all(cv | !is.finite(d))) break
    }
    list(tstar = ts, qstar = qs, conv = cv, n_iter = ni, resid = rs)
  }

  t0 <- -dt_i * kappa / log_t          # neutral first guess
  q0 <- -dq_i * kappa / log_q

  res    <- iterate(seq_along(i), 1, t0, q0)
  tstar  <- res$tstar; qstar <- res$qstar
  conv   <- res$conv;  n_iter <- res$n_iter
  resid  <- res$resid
  used_relax <- rep(FALSE, length(i))   # D2

  # B8: strongly stable rows at very low u* oscillate under plain Picard
  # iteration and never converge.  Damping fixes them; the converged answer
  # is unchanged (differences ~1e-4 W/m2), so only stalled rows are redone.
  if (!is.null(relax_retry) && any(!conv)) {
    sel <- which(!conv)
    r2  <- iterate(sel, relax_retry, t0[sel], q0[sel])
    tstar[sel] <- r2$tstar; qstar[sel] <- r2$qstar
    conv[sel]  <- r2$conv;  n_iter[sel] <- r2$n_iter
    resid[sel] <- r2$resid
    used_relax[sel] <- TRUE
  }

  zeta <- C1 * (tstar + C2 * qstar) / us2
  bad  <- !is.finite(tstar) | !is.finite(qstar)
  if (any(bad)) { tstar[bad] <- NA_real_; qstar[bad] <- NA_real_ }

  if (warn_nonconverged && any(!conv)) {
    warning(sum(!conv), " of ", length(i),
            " solver row(s) did not reach tol = ", tol, " within maxit = ",
            maxit, "; their HS/HL are set to NA (flux_valid = FALSE). ",
            "tstar/qstar/final_residual are retained for inspection.",
            call. = FALSE)
  }

  # ---------------- diagnostics (P3 / B2: computed once) -----------
  Umean <- pmax((us / kappa) * (log_u - .flx_psiu_26(zeta)), 1e-6)

  if (use_gustiness) {
    Bf   <- -grav[i] / taK[i] * us * (tstar + C2 * qstar)   # m2/s3
    Ug_i <- pmax(Ug_min, ifelse(is.finite(Bf) & Bf > 0,
                                Beta * (pmax(Bf, 0) * zi)^(1/3), Ug_min))
    if (!is.null(Ug_max)) Ug_i <- pmin(Ug_i, Ug_max)
    gf_i <- sqrt(Umean^2 + Ug_i^2) / Umean
  } else {
    Ug_i <- rep(NA_real_, length(i))
    gf_i <- rep(1, length(i))
  }

  # ---------------- assemble ---------------------------------------
  out <- na_out()

  out$z0[i]  <- z0[i];  out$zot[i] <- zot[i]; out$zoq[i] <- zoq[i]
  # only report roughness/proxy diagnostics on rows that actually produced a flux
  alpha_used[!ok] <- NA_real_
  U10N_proxy[!ok] <- NA_real_
  out$alpha_used <- alpha_used
  out$U10N_proxy <- U10N_proxy
  # D12: NEVER erase the out-of-range flag.  With clamp_inversion = FALSE an
  # out-of-range row is rejected (ok = FALSE), and blanking the flag alongside
  # it would erase the very reason for the rejection -- the QC count would then
  # report zero out-of-range rows precisely when they were most consequential.
  out$U10N_out_of_range <- U10N_out_of_range
  out$rho_a[i] <- rho_a[i]; out$Q[i] <- Q[i]; out$Qs[i] <- Qs[i]
  out$dt[i] <- dt[i]; out$dq[i] <- dq[i]

  out$tstar[i]   <- tstar
  out$qstar[i]   <- qstar
  out$zeta_zu[i] <- zeta
  out$L[i]       <- ifelse(abs(zeta) < 1e-12, Inf, zu_i / zeta)

  out$HS[i]  <- -rho_a[i] * Cp[i] * us * tstar     # upward positive
  out$HL[i]  <- -rho_a[i] * Le[i] * us * qstar
  out$tau[i] <- rho_a[i] * us^2                    # solver-independent

  # D6: a non-converged solve must not become a flux.  See the blanking block
  # below, after every column has been populated.
  flux_ok <- conv & is.finite(tstar) & is.finite(qstar)
  out$flux_valid[i] <- flux_ok

  out$Ug[i] <- Ug_i; out$gf[i] <- gf_i
  out$Umean_zu[i]  <- Umean
  out$U10N_diag[i] <- (us / kappa) * log(10 / z0[i])  # diagnostic only

  out$converged[i]       <- conv
  out$n_iter[i]          <- n_iter
  out$used_relaxation[i] <- used_relax    # D2
  out$final_residual[i]  <- resid         # D2
  out$Rr[i]              <- rr[i]         # D3: roughness Reynolds number

  # ---------------------------------------------------------------------
  # D6/D8: blank EVERY solver-dependent quantity on rows that failed, so a
  # non-converged state cannot be mistaken for a result further downstream.
  #
  #   blanked : HS, HL, L, zeta_zu, Ug, gf, Umean_zu
  #   kept    : tstar, qstar, n_iter, final_residual, used_relaxation
  #             (debugging only -- they are the UNCONVERGED state)
  #   kept    : tau, z0, zot, zoq, Rr, alpha_used, U10N_proxy, U10N_diag,
  #             rho_a, Q, Qs, dt, dq  -- none of these touch the iteration
  # ---------------------------------------------------------------------
  if (any(!flux_ok)) {
    j <- i[!flux_ok]
    out$HS[j]       <- NA_real_
    out$HL[j]       <- NA_real_
    out$L[j]        <- NA_real_
    out$zeta_zu[j]  <- NA_real_
    out$Ug[j]       <- NA_real_
    out$gf[j]       <- NA_real_
    out$Umean_zu[j] <- NA_real_
  }

  # N3: constants as attributes, not replicated columns
  attr(out, "z0_model")           <- z0_model
  attr(out, "prior_model")        <- prior_model
  attr(out, "alpha_lowwind_mode") <- alpha_lowwind_mode
  # D10: heights may now vary by row, so summarise rather than store n copies
  attr(out, "heights") <- c(
    zu_min = min(zu), zu_max = max(zu),
    zt_min = min(zt), zt_max = max(zt),
    zq_min = min(zq), zq_max = max(zq)
  )
  attr(out, "u_star_min")         <- u_star_min
  attr(out, "gustiness_note")     <-
    "Ug/gf are diagnostic only: with u* prescribed they cannot affect HS/HL."

  # D4: QC roll-up -- the clamped fraction in particular decides whether the
  # choice of prior range is a footnote or a scientific caveat.
  # D7: N_solver and N_flux_valid are NOT the same thing.  n_solver_rows counts
  # rows that reached the iteration; n_flux_valid counts rows that produced a
  # usable HS/HL, and n_flux_valid = n_solver_rows - n_nonconverged.  Reporting
  # only one of them is how a non-converged solve ends up in a mean.
  ns <- length(i)
  attr(out, "qc") <- list(
    n_input          = n,
    n_rejected_ustar = n_rejected_u,
    n_solver_rows    = ns,
    n_flux_valid     = sum(flux_ok),
    n_nonconverged   = sum(!conv),
    n_relaxed        = sum(used_relax),
    # D12: three distinct numbers, previously conflated into one.
    #   n_prior_evaluated    : rows the inversion was actually attempted on
    #   n_prior_out_of_range : of those, how many fell outside the range the
    #                          prior can represent  <-- the science number
    #   n_prior_clamped      : how many were clamped and KEPT (only possible
    #                          when clamp_inversion = TRUE; otherwise those
    #                          rows were rejected instead, and the difference
    #                          n_out_of_range - n_clamped counts the rejections)
    # The denominator is n_prior_evaluated, not n_solver_rows: with
    # clamp_inversion = FALSE the out-of-range rows never become solver rows,
    # so n_solver_rows would understate the fraction exactly when it matters.
    n_prior_evaluated    = sum(!is.na(out$U10N_out_of_range)),
    n_prior_out_of_range = sum(out$U10N_out_of_range, na.rm = TRUE),
    n_prior_clamped      = if (isTRUE(clamp_inversion))
                             sum(out$U10N_out_of_range & ok, na.rm = TRUE) else 0L,
    f_out_of_range_pct   = { ne <- sum(!is.na(out$U10N_out_of_range))
                             if (ne > 0) 100 * sum(out$U10N_out_of_range,
                                                   na.rm = TRUE) / ne
                             else NA_real_ },
    prior_ustar_range    = prior_ustar_range,
    # D9: max(x, na.rm = TRUE) returns -Inf (with a warning) when every
    # element is non-finite, which would print as a bogus residual.
    max_residual     = { rfin <- resid[is.finite(resid)]
                         if (length(rfin)) max(rfin) else NA_real_ },
    # Rr thresholds: with the viscous term in z0, Rr = beta_viscous +
    # alpha*u*^3/(g*nu), so Rr has a HARD FLOOR of 0.11.  Rr < 0.135 is
    # therefore "near the smooth limit" (Charnock contributing < ~20% of Rr),
    # NOT smooth flow in the strict sense -- named accordingly.
    n_near_smooth_limit = sum(out$Rr < 0.135, na.rm = TRUE),
    n_rough_flow        = sum(out$Rr > 2.0,   na.rm = TRUE)
  )
  out
}


# ------------------------------------------------------------------
#  QC summary -- print the fractions that matter for a publication table
# ------------------------------------------------------------------
flx_qc_summary <- function(out) {
  q <- attr(out, "qc")
  if (is.null(q)) stop("no `qc` attribute: not a flx_* result.", call. = FALSE)
  pc <- function(k) if (q$n_solver_rows > 0) 100 * k / q$n_solver_rows else NA_real_

  cat("\nflux QC summary\n"); cat(strrep("-", 58), "\n")
  cat(sprintf("  roughness model         %s (prior: %s)\n",
              attr(out, "z0_model"), attr(out, "prior_model")))
  cat(sprintf("  input rows              %d\n", q$n_input))
  cat(sprintf("  rejected u* <= %-8g %d\n", attr(out, "u_star_min"),
              q$n_rejected_ustar))
  cat(sprintf("  solver rows             %d  (%.2f%% of input)\n",
              q$n_solver_rows,
              if (q$n_input > 0) 100 * q$n_solver_rows / q$n_input else NA_real_))
  cat(sprintf("    of which NOT converged  %d  (%.3f%%)   HS/HL set to NA\n",
              q$n_nonconverged, pc(q$n_nonconverged)))
  cat(sprintf("    needed relaxation       %d  (%.3f%%)\n",
              q$n_relaxed, pc(q$n_relaxed)))
  cat(sprintf("  FLUX-VALID rows         %d  (%.2f%% of solver rows)\n",
              q$n_flux_valid, pc(q$n_flux_valid)))
  cat(sprintf("  prior inversion attempted on %d rows; u* range [%.4g, %.4g]\n",
              q$n_prior_evaluated, q$prior_ustar_range[1], q$prior_ustar_range[2]))
  cat(sprintf("    OUT OF PRIOR RANGE      %d  (%.3f%%)   <-- report this\n",
              q$n_prior_out_of_range, q$f_out_of_range_pct))
  cat(sprintf("      of which clamped+kept   %d\n", q$n_prior_clamped))
  cat(sprintf("      of which rejected       %d\n",
              q$n_prior_out_of_range - q$n_prior_clamped))
  cat(sprintf("  max final residual      %.3e  (||F(x)-x||, relaxation-free)\n",
              q$max_residual))
  cat(sprintf("  near-smooth limit Rr<0.135  %d  (%.2f%%)\n",
              q$n_near_smooth_limit, pc(q$n_near_smooth_limit)))
  cat(sprintf("  rough flow        Rr>2      %d  (%.2f%%)\n",
              q$n_rough_flow, pc(q$n_rough_flow)))
  cat(strrep("-", 58), "\n")
  cat("  N_flux_valid = N_solver - N_nonconverged; use subset(out, flux_valid)\n")
  if (isTRUE(q$f_out_of_range_pct > 1))
    cat("  NOTE: >1% of rows fall outside the prior's representable range --\n",
        "        that range is a scientific choice here, not a numerical\n",
        "        detail.  Report it.\n", sep = "")
  invisible(q)
}
