/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef GEMSTORE_MATH_SOLIDHARM
#define GEMSTORE_MATH_SOLIDHARM

#include <gemstore/math/integral.h>

#define SH_MAX_TERMS 24

/* One term in 𝒴_ℓm(α r + β R) = Σ coef r^{lr} R^{lR} Y_{lr mr}(r̂) Y_{lR mR}(R̂).
 * For ℓ=1 this is exactly α 𝒴_1(r) + β 𝒴_1(R). */
typedef struct {
    int lr;
    int mr;
    int lR;
    int mR;
    double coef;
} sh_term_t;

/* Solid-harmonic addition: 𝒴_ℓm(x) = |x|^ℓ Y_ℓm(x̂). Returns the number of terms. */
int solidharm_shift(int ell, int m, double alpha, double beta, sh_term_t *out);

/* ∫ Y_{l0 m0} ⋯ Y_{l_{n-1} m_{n-1}} dΩ. n=0 returns 4π. */
double ylm_angular_integral(int n, const int *l, const int *m);

/* ∫_0^∞ r^n exp(−a r²) dr */
double gauss_radial_power(double a, int n);

/* Central ME of pot(|r|) between two coupled Gaussians after the R-shift.
 * Caller supplies gem_pref for each (l,ν). Does not include GI OCent. */
double solidharm_central_me(
    int lrho_a, int llam_a, int L_a, int M,
    int lrho_b, int llam_b, int L_b,
    double al_a, double be_a, double ga_a, double de_a,
    double al_b, double be_b, double ga_b, double de_b,
    double b11, double aRR,
    potential_t pot, void *ctx);

#endif
