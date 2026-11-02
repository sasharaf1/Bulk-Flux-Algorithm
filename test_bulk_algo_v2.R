# ============================================================
#  Self-tests for bulk_algo_no_circularity_with_prior_v2.R
#  Run:  Rscript test_bulk_algo_v2.R
#  Base R only, no dependencies.
# ============================================================

source("bulk_algo_no_circularity_with_prior_v2.R")

pass <- 0L; fail <- 0L
chk <- function(label, cond) {
  if (isTRUE(cond)) { pass <<- pass + 1L; cat("  PASS  ", label, "\n") }
  else              { fail <<- fail + 1L; cat("  FAIL  ", label, "\n") }
}
hdr <- function(x) cat("\n", x, "\n", strrep("-", nchar(x) + 1), "\n", sep = "")

# round() coerces logical columns to 0/1 and flattens small z0 values to
# 0.0000, so print through this instead: signif for the numerics, logicals
# left alone.
show <- function(d, digits = 4) {
  for (nm in names(d)) if (is.numeric(d[[nm]])) d[[nm]] <- signif(d[[nm]], digits)
  print(d, row.names = FALSE)
}

f <- flx_compute_HS_HL_fixed_ustar_prior

# ------------------------------------------------------------------
hdr("1. B1 -- scalar P with sub-freezing temperatures (v1 gave NA)")
r <- f(u_star = c(0.15, 0.25, 0.35, 0.45),
       T_s = c(-1.5, -1.0, 2.0, 4.0),
       T_a = c(-8.0, -5.0, -2.0, 1.0),
       RH = 85, P = 1013)                    # <- P is a SCALAR: the v1 trap
chk("all HS finite", all(is.finite(r$HS)))
chk("all HL finite", all(is.finite(r$HL)))
chk("cold-air ocean warms atmosphere (HS > 0)", all(r$HS > 0))
show(r[, c("u_star", "z0", "HS", "HL", "L")])

# ------------------------------------------------------------------
hdr("2. B4 -- u* below the Andreas prior floor (u*(0) = 0.00629)")
r <- suppressWarnings(f(u_star = c(0.001, 0.004, 0.05, 0.30, 0.90),
                        T_s = 20, T_a = 17, RH = 80, P = 1013))
chk("no silent all-NA rows", all(is.finite(r$U10N_proxy)))
chk("out-of-range flag marks the two sub-floor rows",
    sum(r$U10N_out_of_range) == 2L)
chk("clamped rows sit at U10N_lo", all(r$U10N_proxy[r$U10N_out_of_range] == 0))
chk("in-range rows produce fluxes", all(is.finite(r$HS[3:5])))
# The two sub-floor rows are clamped (B4 works: no silent NA proxy) but they
# still fail to solve, and that is CORRECT, not a solver deficiency: with u*
# prescribed, zeta scales as 1/u*^2, so u* = 0.001 gives |zeta| ~ 1e4 -- the
# free-convection limit, where MO similarity at fixed u* is ill-posed.  It
# converges at no relaxation coefficient (checked from 1.0 down to 0.02).
# D6 is what makes this safe: the rows are flagged, not silently returned.
chk("ill-posed u* -> 0 rows are rejected, not silently solved",
    all(is.na(r$HS[1:2])) && all(!r$flux_valid[1:2]))
chk("...and they are flagged as non-converged, with a large residual",
    all(!r$converged[1:2]) && all(r$final_residual[1:2] > 1e-6))
chk("...while the proxy itself was still recovered (the B4 fix)",
    all(is.finite(r$U10N_proxy[1:2])))
show(r[, c("u_star", "U10N_proxy", "alpha_used", "converged",
           "flux_valid", "HS")])
cat("   NOTE: convergence returns around u* ~ 0.01 (|zeta| ~ 550) and is\n",
    "        comfortable by u* ~ 0.06 (U10N ~ 2).  This is the practical\n",
    "        argument for setting u_star_min > 0 on real data.\n", sep = "")

# ---- D12: the flag must survive when clamping is OFF ----------------
r_na <- suppressWarnings(f(u_star = c(0.001, 0.004, 0.30, 0.90),
                           T_s = 20, T_a = 17, RH = 80, P = 1013,
                           clamp_inversion = FALSE))
chk("clamp_inversion = FALSE gives NA proxy", all(is.na(r_na$U10N_proxy[1:2])))
chk("out-of-range rows are rejected", all(is.na(r_na$HS[1:2])))
chk("D12: the out-of-range flag is NOT erased on rejected rows",
    identical(r_na$U10N_out_of_range, c(TRUE, TRUE, FALSE, FALSE)))
qn_ <- attr(r_na, "qc")
chk("n_prior_out_of_range still counts them", qn_$n_prior_out_of_range == 2L)
chk("n_prior_clamped is 0 when clamping is off", qn_$n_prior_clamped == 0L)
chk("rejected = out_of_range - clamped",
    qn_$n_prior_out_of_range - qn_$n_prior_clamped == 2L)
chk("n_solver_rows excludes the rejected out-of-range rows",
    qn_$n_solver_rows == 2L)
chk("f_out_of_range_pct is over rows the inversion was attempted on",
    abs(qn_$f_out_of_range_pct - 100 * 2 / qn_$n_prior_evaluated) < 1e-12)

qc_ <- attr(r, "qc")
chk("with clamping ON the same rows are counted AND kept",
    qc_$n_prior_out_of_range == 2L && qc_$n_prior_clamped == 2L)
chk("prior u* range is reported", all(is.finite(qc_$prior_ustar_range)))
cat(sprintf("   Andreas prior spans u* [%.5f, %.4f] over U10N [0, 60]\n",
            qc_$prior_ustar_range[1], qc_$prior_ustar_range[2]))

# ------------------------------------------------------------------
hdr("3. B8 -- strongly stable, very low u* (plain Picard never converges)")
u <- c(0.02, 0.03, 0.05, 0.08, 0.20, 0.50)
r_off <- suppressWarnings(f(u_star = u, T_s = 10, T_a = 22, RH = 60, P = 1013,
                            relax_retry = NULL))
r_on  <- f(u_star = u, T_s = 10, T_a = 22, RH = 60, P = 1013)
cat("   without damped retry: ", sum(r_off$converged), "/", length(u),
    " converged\n", sep = "")
cat("   with    damped retry: ", sum(r_on$converged),  "/", length(u),
    " converged\n", sep = "")
chk("damped retry rescues the stalled rows", all(r_on$converged))
chk("retry does not disturb already-converged rows",
    max(abs(r_on$HS[r_off$converged] - r_off$HS[r_off$converged])) < 1e-3)
chk("stable stratification gives HS < 0 and L > 0",
    all(r_on$HS < 0) && all(r_on$L > 0))
show(r_on[, c("u_star", "zeta_zu", "L", "HS", "n_iter")])

# ------------------------------------------------------------------
hdr("3b. D2 -- used_relaxation / final_residual")
chk("used_relaxation flags exactly the stalled rows",
    identical(r_on$used_relaxation, !r_off$converged))
chk("final_residual finite on all valid rows", all(is.finite(r_on$final_residual)))
chk("converged rows have residual < tol", all(r_on$final_residual[r_on$converged] < 1e-6))
r_stall <- suppressWarnings(f(u_star = u, T_s = 10, T_a = 22, RH = 60, P = 1013,
                              relax_retry = NULL))
chk("stalled rows report a residual >= tol",
    all(r_stall$final_residual[!r_stall$converged] >= 1e-6))
show(r_on[, c("u_star", "converged", "n_iter", "used_relaxation",
              "final_residual")])

# ------------------------------------------------------------------
hdr("3c. D5 -- converged rows are frozen, not advanced with the batch")
# A batch mixing fast and slow rows must give each row exactly the same answer
# as solving it alone; if converged rows kept iterating, they would not.
ub <- c(0.05, 0.08, 0.20, 0.50)                 # row 1 slow, row 4 fast
batch <- f(u_star = ub, T_s = 10, T_a = 22, RH = 60, P = 1013)
solo  <- vapply(ub, function(x)
  f(u_star = x, T_s = 10, T_a = 22, RH = 60, P = 1013)$tstar, 0)
chk("batched t* == solo t* bit-for-bit", identical(batch$tstar, solo))
soloHS <- vapply(ub, function(x)
  f(u_star = x, T_s = 10, T_a = 22, RH = 60, P = 1013)$HS, 0)
chk("batched HS == solo HS bit-for-bit", identical(batch$HS, soloHS))
chk("fast rows report small n_iter despite slow rows in the batch",
    batch$n_iter[4] < batch$n_iter[1])
cat("   n_iter by row:", paste(batch$n_iter, collapse = ", "), "\n")

# ------------------------------------------------------------------
hdr("3d. D6 -- a non-converged solve must not become a flux")
chk("non-converged rows have HS = NA", all(is.na(r_stall$HS[!r_stall$converged])))
chk("non-converged rows have HL = NA", all(is.na(r_stall$HL[!r_stall$converged])))
chk("flux_valid is FALSE exactly there",
    identical(r_stall$flux_valid, r_stall$converged))
chk("tstar retained for debugging", all(is.finite(r_stall$tstar)))

# D8: no solver-dependent diagnostic may survive on a failed row
bad <- !r_stall$converged
for (v in c("HS", "HL", "L", "zeta_zu", "Ug", "gf", "Umean_zu"))
  chk(sprintf("solver-dependent %-9s blanked on failed rows", v),
      all(is.na(r_stall[[v]][bad])))
for (v in c("tstar", "qstar", "final_residual", "n_iter"))
  chk(sprintf("debug field %-15s retained on failed rows", v),
      all(!is.na(r_stall[[v]][bad])))
for (v in c("tau", "z0", "zot", "zoq", "Rr", "alpha_used", "U10N_proxy",
            "U10N_diag", "rho_a", "Q", "Qs", "dt", "dq"))
  chk(sprintf("solver-independent %-11s retained on failed rows", v),
      all(is.finite(r_stall[[v]][bad])))
chk("U10N_out_of_range retained on failed rows",
    all(!is.na(r_stall$U10N_out_of_range[bad])))

# ------------------------------------------------------------------
hdr("3e. D11 -- residual is the raw fixed-point norm, not the relaxed step")
# With relax != 1 the reported residual must NOT scale with lambda.
res_a <- suppressWarnings(f(u_star = u, T_s = 10, T_a = 22, RH = 60, P = 1013,
                            relax_retry = 0.5))
res_b <- suppressWarnings(f(u_star = u, T_s = 10, T_a = 22, RH = 60, P = 1013,
                            relax_retry = 0.25))
rel <- res_a$used_relaxation & res_a$converged &
  res_b$used_relaxation & res_b$converged
chk("relaxed rows exist to compare", any(rel))
chk("residual threshold is lambda-independent (both < tol)",
    all(res_a$final_residual[rel] < 1e-6) && all(res_b$final_residual[rel] < 1e-6))
# the undamped path must be untouched by this change
chk("relax = 1 rows unaffected by relax_retry choice",
    identical(res_a$HS[!res_a$used_relaxation],
              res_b$HS[!res_b$used_relaxation]))
cat("   residual at relax=0.50:", format(max(res_a$final_residual[rel]), digits = 3),
    "| at relax=0.25:", format(max(res_b$final_residual[rel]), digits = 3), "\n")
chk("converged rows keep their flux", all(is.finite(r_stall$HS[r_stall$converged])))
chk("subset(out, flux_valid) yields only finite fluxes",
    all(is.finite(subset(r_stall, flux_valid)$HS)))
qs_ <- attr(r_stall, "qc")
chk("n_flux_valid == n_solver_rows - n_nonconverged",
    qs_$n_flux_valid == qs_$n_solver_rows - qs_$n_nonconverged)
show(r_stall[, c("u_star", "converged", "flux_valid", "tstar", "HS", "tau")])

# ------------------------------------------------------------------
hdr("4. D1 -- NA / zero / negative u* return NA rows instead of aborting")
r <- suppressWarnings(f(u_star = c(0.25, NA, 0.30, -0.1, 0.0, 0.4),
                        T_s = 20, T_a = 17, RH = 80, P = 1013))
chk("nrow preserved", nrow(r) == 6L)
chk("NA, negative and ZERO u* all rejected", all(is.na(r$HS[c(2, 4, 5)])))
chk("good rows are finite", all(is.finite(r$HS[c(1, 3, 6)])))
chk("u* = 0 no longer yields a floor-generated z0", is.na(r$z0[5]))

r_cut <- suppressWarnings(f(u_star = c(0.05, 0.10, 0.20, 0.40),
                            T_s = 20, T_a = 17, RH = 80, P = 1013,
                            u_star_min = 0.13))   # smooth/rough crossover
chk("u_star_min rejects rows at or below the cutoff",
    all(is.na(r_cut$HS[1:2])) && all(is.finite(r_cut$HS[3:4])))

# ------------------------------------------------------------------
hdr("4b. D3/D4 -- Rr regime column and QC roll-up")
rq <- suppressWarnings(f(u_star = c(0.004, 0.03, 0.09, 0.20, 0.50, 0.90, NA),
                         T_s = 20, T_a = 17, RH = 80, P = 1013))
chk("Rr returned and finite on valid rows", all(is.finite(rq$Rr[1:6])))
chk("Rr increases with u*", all(diff(rq$Rr[2:6]) > 0))
chk("Rr is bounded below by beta_viscous = 0.11 (identity, not a threshold)",
    all(rq$Rr[1:6] >= 0.11 - 1e-12))
chk("Rr == 0.11 + alpha*u*^3/(g*nu) exactly",
    max(abs(rq$Rr[1:6] - (0.11 + rq$alpha_used[1:6] * rq$u_star[1:6]^3 /
                            (9.8062 * 1.326e-5 * (1 + 6.542e-3 * 17 +
                               8.301e-6 * 17^2 - 4.84e-9 * 17^3))))) < 2e-3)
chk("low u* sits at the smooth limit (Rr < 0.135)", rq$Rr[2] < 0.135)
chk("high u* is aerodynamically rough (Rr > 2)", rq$Rr[6] > 2)
show(rq[, c("u_star", "U10N_proxy", "Rr", "z0", "HS")])
q <- attr(rq, "qc")
chk("qc attribute present", !is.null(q))
chk("f_out_of_range_pct uses n_prior_evaluated as denominator",
    abs(q$f_out_of_range_pct -
          100 * q$n_prior_out_of_range / q$n_prior_evaluated) < 1e-12)
chk("n_solver_rows excludes the NA input row", q$n_solver_rows == 6L)
# D9: max_residual must be NA, never -Inf, when nothing finite is available
chk("max_residual is finite here", is.finite(q$max_residual))
q0_ <- attr(suppressWarnings(f(u_star = c(NA, NA), T_s = 20, T_a = 17,
                               RH = 80, P = 1013)), "qc")
chk("all-invalid input gives max_residual = NA, not -Inf",
    is.na(q0_$max_residual) && !is.infinite(q0_$max_residual))
chk("all-invalid input gives n_solver_rows = 0", q0_$n_solver_rows == 0L)
chk("QC label matches the threshold actually used (0.135)",
    q$n_near_smooth_limit == sum(rq$Rr < 0.135, na.rm = TRUE))
flx_qc_summary(rq)

# ------------------------------------------------------------------
hdr("4c. D10 -- per-observation measurement heights")
u4 <- c(0.15, 0.25, 0.35, 0.50)

# a scalar height must be bit-identical to the same value repeated
chk("scalar zt/zq/zu == length-n repeat, bit-for-bit",
    identical(f(u4, 20, 17, 80, 1013, zt = 10, zq = 10, zu = 10)$HS,
              f(u4, 20, 17, 80, 1013, zt = rep(10, 4), zq = rep(10, 4),
                zu = rep(10, 4))$HS))

# each row must match solving that row alone at its own heights
zt_v <- c(3.0, 4.5, 10.0, 18.0)
zq_v <- c(3.0, 4.5, 10.0, 18.0)
zu_v <- c(5.0, 5.0, 10.0, 20.0)
mixed <- f(u4, 20, 17, 80, 1013, zt = zt_v, zq = zq_v, zu = zu_v)
solo4 <- vapply(seq_along(u4), function(k)
  f(u4[k], 20, 17, 80, 1013, zt = zt_v[k], zq = zq_v[k], zu = zu_v[k])$HS, 0)
chk("per-row heights == per-row solo solves", identical(mixed$HS, solo4))
chk("all finite with mixed heights", all(is.finite(mixed$HS)))

# zt enters the lapse correction, so HS must respond to zt
chk("dt reflects the per-row lapse -0.0098*zt",
    max(abs(mixed$dt - (20 - 17 - 0.0098 * zt_v))) < 1e-12)
# zu != zt must not be silently ignored: zeta_zu is defined at zu
chk("L is defined at zu (L = zu / zeta_zu)",
    max(abs(mixed$L - zu_v / mixed$zeta_zu)) < 1e-9)
chk("heights attribute reports the range",
    identical(unname(attr(mixed, "heights")[c("zu_min", "zu_max")]),
              c(min(zu_v), max(zu_v))))
chk("bad height length is an error",
    inherits(try(f(u4, 20, 17, 80, 1013, zt = c(10, 10)), silent = TRUE),
             "try-error"))
chk("non-positive height is an error",
    inherits(try(f(u4, 20, 17, 80, 1013, zt = c(10, 10, 0, 10)), silent = TRUE),
             "try-error"))
show(cbind(mixed[, c("u_star", "dt", "zeta_zu", "L", "HS", "HL")],
           zt = zt_v, zu = zu_v))

# ------------------------------------------------------------------
hdr("5. R1 -- mutually exclusive roughness options")
chk("two legacy flags TRUE is an error",
    inherits(try(f(0.3, 20, 17, 80, 1013, use_fixed_z0 = TRUE,
                   use_eq3_z0 = TRUE), silent = TRUE), "try-error"))
r <- suppressWarnings(f(0.3, 20, 17, 80, 1013, use_fixed_z0 = TRUE))
chk("legacy flag maps onto z0_model", attr(r, "z0_model") == "fixed")
chk("fixed z0 honoured", abs(r$z0 - 1e-4) < 1e-12)

# ------------------------------------------------------------------
hdr("6. R2 -- input length validation")
chk("length-3 covariate against length-5 u* errors",
    inherits(try(f(rep(0.3, 5), c(20, 20, 20), 17, 80, 1013), silent = TRUE),
             "try-error"))
chk("length-1 recycles fine", nrow(f(rep(0.3, 5), 20, 17, 80, 1013)) == 5L)

# ------------------------------------------------------------------
hdr("7. B3 -- eq3 accepts cp alone (v1 wrongly demanded Tp)")
r <- f(u_star = c(0.2, 0.4), T_s = 20, T_a = 17, RH = 80, P = 1013,
       z0_model = "eq3", Hs = 2.5, cp = 15.6)
chk("cp-only call succeeds", all(is.finite(r$z0)))
r2 <- f(u_star = c(0.2, 0.4), T_s = 20, T_a = 17, RH = 80, P = 1013,
        z0_model = "eq3", Hs = 2.5, Tp = 10)
chk("deep-water cp from Tp matches cp = g*Tp/(2*pi)",
    max(abs(r$z0 - r2$z0)) / max(r$z0) < 0.02)
show(r2[, c("u_star", "z0", "HS", "HL")])

# ------------------------------------------------------------------
hdr("8. All five z0 models run and are ordered sensibly")
args <- list(u_star = c(0.15, 0.35, 0.60), T_s = 20, T_a = 17,
             RH = 80, P = 1013, Hs = 2.5, Tp = 10)
for (m in c("charnock_winddep", "charnock_const", "fixed", "eq3", "yelland")) {
  rr <- do.call(f, c(args, list(z0_model = m)))
  cat(sprintf("   %-17s z0 = %s   HS = %s\n", m,
              paste(format(signif(rr$z0, 3), width = 9), collapse = " "),
              paste(format(round(rr$HS, 1), width = 7), collapse = " ")))
  chk(paste(m, "produces finite fluxes"), all(is.finite(rr$HS)))
}

# ------------------------------------------------------------------
hdr("9. Physics sanity")
u  <- c(0.10, 0.20, 0.35, 0.50, 0.80)
ru <- f(u_star = u, T_s = 20, T_a = 17, RH = 80, P = 1013)   # unstable
rs <- f(u_star = u, T_s = 14, T_a = 18, RH = 80, P = 1013)   # stable
chk("warm ocean -> upward HS and HL", all(ru$HS > 0) && all(ru$HL > 0))
chk("warm ocean -> L < 0",            all(ru$L < 0))
chk("cool ocean -> downward HS",      all(rs$HS < 0))
chk("cool ocean -> L > 0",            all(rs$L > 0))
chk("HS increases monotonically with u*", all(diff(ru$HS) > 0))
chk("tau = rho * u*^2",  max(abs(ru$tau - ru$rho_a * u^2)) < 1e-12)
chk("alpha within [alpha_min, alpha_max]",
    all(ru$alpha_used >= 0.011 - 1e-12 & ru$alpha_used <= 0.028 + 1e-12))
chk("alpha held at alpha_min below 4 m/s",
    all(ru$alpha_used[ru$U10N_proxy < 4] == 0.011))
show(ru[, c("u_star", "U10N_proxy", "alpha_used", "z0", "L",
            "HS", "HL", "n_iter")])

# ------------------------------------------------------------------
hdr("10. B2 -- gustiness cannot change HS/HL when u* is prescribed")
a <- f(u_star = u, T_s = 20, T_a = 17, RH = 80, P = 1013, use_gustiness = TRUE)
b <- f(u_star = u, T_s = 20, T_a = 17, RH = 80, P = 1013, use_gustiness = FALSE)
chk("HS identical with and without gustiness", identical(a$HS, b$HS))
chk("HL identical with and without gustiness", identical(a$HL, b$HL))
cat("   (confirms Ug/gf are diagnostic only -- see attr(out,'gustiness_note'))\n")

# ------------------------------------------------------------------
hdr("10b. Cp -- HS must be built from DRY cpa = 1004.67 (COARE 3.5)")
r <- f(u_star = u, T_s = 20, T_a = 17, RH = 80, P = 1013)
chk("HS == -rho * 1004.67 * u* * t* exactly",
    max(abs(r$HS - (-r$rho_a * 1004.67 * r$u_star * r$tstar))) < 1e-10)

rm_ <- f(u_star = u, T_s = 20, T_a = 17, RH = 80, P = 1013, cp_moist = TRUE)
chk("cp_moist = TRUE uses cpv = cpa*(1 + 0.84*Q), Q in kg/kg",
    max(abs(rm_$HS - (-rm_$rho_a * 1004.67 * (1 + 0.84 * rm_$Q) *
                        rm_$u_star * rm_$tstar))) < 1e-10)
cat(sprintf("   Q = %.5f kg/kg -> moist Cp is %+.2f%% on HS\n",
            r$Q[1], 100 * (max(rm_$HS / r$HS) - 1)))
chk("moist option shifts HS by ~0.8%, not ~0.0008%",
    abs(max(rm_$HS / r$HS) - 1) > 0.002)

# ------------------------------------------------------------------
hdr("11. P1 -- lookup inversion agrees with uniroot")
andreas <- function(U) 0.239 + 0.0433 * ((U - 8.271) +
                                           sqrt(0.120 * (U - 8.271)^2 + 0.181))
us  <- seq(andreas(0) + 1e-4, andreas(60) - 1e-4, length.out = 300)
ref <- vapply(us, function(t)
  uniroot(function(U) andreas(U) - t, c(0, 60), tol = 1e-12)$root, 0)
got <- f(u_star = us, T_s = 20, T_a = 17, RH = 80, P = 1013)$U10N_proxy
chk("max |lookup - uniroot| < 5e-4 m/s", max(abs(got - ref)) < 5e-4)
cat("   max deviation:", format(max(abs(got - ref)), digits = 3), "m/s\n")

# ------------------------------------------------------------------
hdr("12. P1 -- speed on a realistic record")
set.seed(1)
N  <- 2e5
uu <- runif(N, 0.02, 1.2)
tm <- system.time(f(u_star = uu, T_s = 20, T_a = 17, RH = 80, P = 1013))
cat("   n =", N, "rows in", format(tm[["elapsed"]], digits = 3), "s\n")
chk("200k rows in under 10 s", tm[["elapsed"]] < 10)

# ------------------------------------------------------------------
hdr("13. Both priors, all alpha low-wind modes")
for (pm in c("andreas", "edson2013"))
  for (lw in c("hold", "ramp", "none")) {
    rr <- suppressWarnings(f(u_star = c(0.05, 0.15, 0.30, 0.60),
                             T_s = 20, T_a = 17, RH = 80, P = 1013,
                             prior_model = pm, alpha_lowwind_mode = lw))
    chk(sprintf("%s / %-4s finite", pm, lw), all(is.finite(rr$HS)))
  }

# ------------------------------------------------------------------
cat("\n", strrep("=", 46), "\n", sep = "")
cat(sprintf("  %d passed, %d failed\n", pass, fail))
cat(strrep("=", 46), "\n")
if (fail > 0L) quit(status = 1L)
