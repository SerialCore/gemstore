/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 *
 * Solid-harmonic addition and the 6D central matrix element after the
 * Jacobi R-shift. Angular integrals are products of Ylm reduced by
 * successive Clebsch–Gordan couplings.
 */

#include <gemstore/math/solidharm.h>
#include <gemstore/math/soc.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define YLM_MAX 48

typedef struct {
    int l;
    int m;
    double c;
} ylm_comp_t;

static double fact_int(int n)
{
    if (n < 0) {
        return 0.0;
    }
    return tgamma((double)n + 1.0);
}

int solidharm_shift(int ell, int m, double alpha, double beta, sh_term_t *out)
{
    /* 𝒴_ℓm(αr+βR) = √(4π) Σ_k α^k β^{ℓ-k} √[(2ℓ+1)!/((2k+1)!(2ℓ-2k+1)!)]
     *   × Σ CG(k m1, ℓ-k m2; ℓ m) Y_{k m1}(r̂) Y_{ℓ-k m2}(R̂) */
    int n = 0;

    if (ell < 0 || abs(m) > ell) {
        return 0;
    }

    for (int k = 0; k <= ell; k++) {
        int k2 = ell - k;
        double fac = sqrt(4.0 * M_PI) * pow(alpha, k) * pow(beta, k2)
            * sqrt(fact_int(2 * ell + 1) / (fact_int(2 * k + 1) * fact_int(2 * k2 + 1)));

        for (int m1 = -k; m1 <= k; m1++) {
            int m2 = m - m1;
            double cg;
            if (m2 < -k2 || m2 > k2) {
                continue;
            }
            cg = clebsch_gordan((double)k, (double)m1, (double)k2, (double)m2,
                (double)ell, (double)m);
            if (cg == 0.0) {
                continue;
            }
            if (n >= SH_MAX_TERMS) {
                fprintf(stderr, "Error: solidharm_shift overflow\n");
                return n;
            }
            out[n].lr = k;
            out[n].mr = m1;
            out[n].lR = k2;
            out[n].mR = m2;
            out[n].coef = fac * cg;
            n++;
        }
    }

    return n;
}

/* Multiply a Ylm expansion by one extra Y_{l2 m2}. CG(l1 0, l2 0; l 0)
 * vanishes unless l1+l2+l is even. */
static void ylm_mul_one(const ylm_comp_t *in, int nin, int l2, int m2,
    ylm_comp_t *out, int *nout)
{
    int n = 0;

    for (int i = 0; i < nin; i++) {
        int l1 = in[i].l;
        int m1 = in[i].m;
        int lmin = abs(l1 - l2);
        int lmax = l1 + l2;

        for (int l = lmin; l <= lmax; l++) {
            double cg_m, cg_0, pref, coef;
            int m = m1 + m2;
            int j, found;

            if (abs(m) > l) {
                continue;
            }
            cg_0 = clebsch_gordan((double)l1, 0.0, (double)l2, 0.0, (double)l, 0.0);
            if (cg_0 == 0.0) {
                continue;
            }
            cg_m = clebsch_gordan((double)l1, (double)m1, (double)l2, (double)m2,
                (double)l, (double)m);
            if (cg_m == 0.0) {
                continue;
            }
            pref = sqrt((2.0 * l1 + 1.0) * (2.0 * l2 + 1.0) / (4.0 * M_PI * (2.0 * l + 1.0)));
            coef = in[i].c * pref * cg_0 * cg_m;

            found = 0;
            for (j = 0; j < n; j++) {
                if (out[j].l == l && out[j].m == m) {
                    out[j].c += coef;
                    found = 1;
                    break;
                }
            }
            if (!found) {
                if (n >= YLM_MAX) {
                    fprintf(stderr, "Error: ylm_mul overflow\n");
                    *nout = n;
                    return;
                }
                out[n].l = l;
                out[n].m = m;
                out[n].c = coef;
                n++;
            }
        }
    }
    *nout = n;
}

double ylm_angular_integral(int n, const int *l, const int *m)
{
    ylm_comp_t cur[YLM_MAX];
    ylm_comp_t nxt[YLM_MAX];
    int ncur = 1;
    int i;

    if (n <= 0) {
        return 4.0 * M_PI;
    }

    /* Start from 1 = √(4π) Y_00; after n multiplies, pick the Y_00 piece
     * and convert back with another √(4π). */
    cur[0].l = 0;
    cur[0].m = 0;
    cur[0].c = sqrt(4.0 * M_PI);

    for (i = 0; i < n; i++) {
        int nnxt = 0;
        ylm_mul_one(cur, ncur, l[i], m[i], nxt, &nnxt);
        ncur = nnxt;
        for (int k = 0; k < ncur; k++) {
            cur[k] = nxt[k];
        }
    }

    for (i = 0; i < ncur; i++) {
        if (cur[i].l == 0 && cur[i].m == 0) {
            return cur[i].c * sqrt(4.0 * M_PI);
        }
    }
    return 0.0;
}

double gauss_radial_power(double a, int n)
{
    /* ∫_0^∞ r^n e^{−a r²} dr = (1/2) a^{−(n+1)/2} Γ((n+1)/2).
     * n=2 is √π/(4 a^{3/2}), matching the meson spectator measure. */
    if (a <= 0.0 || n < 0) {
        return 0.0;
    }
    return 0.5 * pow(a, -0.5 * (n + 1.0)) * tgamma(0.5 * (n + 1.0));
}

double solidharm_central_me(
    int lrho_a, int llam_a, int L_a, int M,
    int lrho_b, int llam_b, int L_b,
    double al_a, double be_a, double ga_a, double de_a,
    double al_b, double be_b, double ga_b, double de_b,
    double b11, double aRR,
    potential_t pot, void *ctx)
{
    double sum = 0.0;
    sh_term_t trhoa[SH_MAX_TERMS], tlama[SH_MAX_TERMS];
    sh_term_t trhob[SH_MAX_TERMS], tlamb[SH_MAX_TERMS];

    if (L_a != L_b) {
        return 0.0;
    }

    for (int mrho_a = -lrho_a; mrho_a <= lrho_a; mrho_a++) {
        for (int mlam_a = -llam_a; mlam_a <= llam_a; mlam_a++) {
            double cg_a = clebsch_gordan((double)lrho_a, (double)mrho_a,
                (double)llam_a, (double)mlam_a, (double)L_a, (double)M);
            int nra, nla;
            if (cg_a == 0.0) {
                continue;
            }
            /* Bra: 𝒴_{ℓm}^* = (−1)^m 𝒴_{ℓ,−m}. Use (m&1), not m%2:
             * C remainder of a negative odd m is −1. */
            nra = solidharm_shift(lrho_a, -mrho_a, al_a, be_a, trhoa);
            nla = solidharm_shift(llam_a, -mlam_a, ga_a, de_a, tlama);
            cg_a *= ((mrho_a + mlam_a) & 1) ? -1.0 : 1.0;

            for (int mrho_b = -lrho_b; mrho_b <= lrho_b; mrho_b++) {
                for (int mlam_b = -llam_b; mlam_b <= llam_b; mlam_b++) {
                    double cg_b = clebsch_gordan((double)lrho_b, (double)mrho_b,
                        (double)llam_b, (double)mlam_b, (double)L_b, (double)M);
                    int nrb, nlb;
                    if (cg_b == 0.0) {
                        continue;
                    }
                    nrb = solidharm_shift(lrho_b, mrho_b, al_b, be_b, trhob);
                    nlb = solidharm_shift(llam_b, mlam_b, ga_b, de_b, tlamb);

                    for (int ia = 0; ia < nra; ia++) {
                        for (int ja = 0; ja < nla; ja++) {
                            for (int ib = 0; ib < nrb; ib++) {
                                for (int jb = 0; jb < nlb; jb++) {
                                    int lr[4], mr[4], lR[4], mR[4];
                                    int nr, nR;
                                    double coef, ang_r, ang_R, rad_r, rad_R;

                                    coef = cg_a * cg_b
                                        * trhoa[ia].coef * tlama[ja].coef
                                        * trhob[ib].coef * tlamb[jb].coef;

                                    lr[0] = trhoa[ia].lr; mr[0] = trhoa[ia].mr;
                                    lr[1] = tlama[ja].lr; mr[1] = tlama[ja].mr;
                                    lr[2] = trhob[ib].lr; mr[2] = trhob[ib].mr;
                                    lr[3] = tlamb[jb].lr; mr[3] = tlamb[jb].mr;
                                    lR[0] = trhoa[ia].lR; mR[0] = trhoa[ia].mR;
                                    lR[1] = tlama[ja].lR; mR[1] = tlama[ja].mR;
                                    lR[2] = trhob[ib].lR; mR[2] = trhob[ib].mR;
                                    lR[3] = tlamb[jb].lR; mR[3] = tlamb[jb].mR;

                                    nr = lr[0] + lr[1] + lr[2] + lr[3];
                                    nR = lR[0] + lR[1] + lR[2] + lR[3];

                                    ang_R = ylm_angular_integral(4, lR, mR);
                                    if (ang_R == 0.0) {
                                        continue;
                                    }
                                    ang_r = ylm_angular_integral(4, lr, mr);
                                    if (ang_r == 0.0) {
                                        continue;
                                    }

                                    /* d³x = x² dx dΩ, so radial power is n+2 */
                                    rad_R = gauss_radial_power(aRR, nR + 2);
                                    rad_r = integral_exp_rn(pot, b11, nr + 2, ctx);
                                    sum += coef * ang_r * ang_R * rad_r * rad_R;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    return sum;
}
