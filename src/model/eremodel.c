/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 *
 * S-wave effective-range expansion.
 * The potential is a callback. Model projection lives in cscatter.
 *
 * Yuan–Liu Eq. (8) at l = 0, integrated with Numerov to 80 fm in steps
 * of 0.01 fm, using ħc = 197.3 MeV·fm. The outer wave is matched to
 * sin(kr + δ0). Miyamoto's window k ∈ (0, 0.5] MeV fixes a0 and r0
 * through k cot δ0 = 1/a0 + r0 k^2 / 2.
 * The match assumes V is negligible for r ≥ 40 fm, so Coulomb stays out.
 */

#include <gemstore/model/eremodel.h>

#include <math.h>
#include <stddef.h>

#define ERE_HBARC 197.3
#define ERE_RMAX 80.0
#define ERE_DR 0.01
#define ERE_NSTEP 8000
#define ERE_MATCH_I0 4000
#define ERE_NK 10
#define ERE_K_STEP 0.05

/* g = 2μV/(ħc)^2 − k^2, with r and k in fm. */
static long double ere_g(double mu_mev, long double k_fm, double r_fm, ere_pot_fn v, void *ctx)
{
    long double pot = (long double)v(r_fm, ctx);
    long double hbarc = ERE_HBARC;
    long double k2 = k_fm * k_fm;

    return 2.0L * (long double)mu_mev * pot / (hbarc * hbarc) - k2;
}

static int ere_outer_tail(ere_pot_fn v, void *ctx)
{
    int i;

    for (i = ERE_MATCH_I0; i <= ERE_NSTEP; i += 1000) {
        double r = i * ERE_DR;
        double pot = v(r, ctx);

        if (!isfinite(pot) || fabs(pot) > 1.0e-6) {
            return 0;
        }
    }
    return 1;
}

/* March u(0) = 0 and read off α, β in u → α sin(kr) + β cos(kr) on r ≥ 40 fm. */
static int ere_match(double mu_mev, double k_mev, ere_pot_fn v, void *ctx,
    long double *alpha, long double *beta)
{
    const long double h = ERE_DR;
    const long double h2 = h * h;
    long double k_fm;
    long double g0;
    long double u_im1;
    long double u_i;
    long double ss = 0.0L;
    long double sc = 0.0L;
    long double cc = 0.0L;
    long double su = 0.0L;
    long double cu = 0.0L;
    long double det;
    int i;

    if (k_mev <= 0.0) {
        return -1;
    }
    k_fm = (long double)k_mev / ERE_HBARC;
    g0 = ere_g(mu_mev, k_fm, 0.0, v, ctx);
    if (!isfinite((double)g0)) {
        return -2;
    }
    u_im1 = 0.0L;
    u_i = h + g0 * h * h * h / 6.0L;

    for (i = 1; i < ERE_NSTEP; i++) {
        double r_im1 = (i - 1) * ERE_DR;
        double r_i = i * ERE_DR;
        double r_ip1 = (i + 1) * ERE_DR;
        long double g_im1 = ere_g(mu_mev, k_fm, r_im1, v, ctx);
        long double g_i = ere_g(mu_mev, k_fm, r_i, v, ctx);
        long double g_ip1 = ere_g(mu_mev, k_fm, r_ip1, v, ctx);
        long double left = 1.0L - h2 * g_ip1 / 12.0L;
        long double right;
        long double u_ip1;
        int ip1 = i + 1;

        if (!isfinite((double)g_ip1) || fabsl(left) < 1.0e-18L) {
            return -2;
        }
        right = 2.0L * (1.0L + 5.0L * h2 * g_i / 12.0L) * u_i
            - (1.0L - h2 * g_im1 / 12.0L) * u_im1;
        u_ip1 = right / left;
        if (ip1 >= ERE_MATCH_I0) {
            long double s = sinl(k_fm * (long double)r_ip1);
            long double c = cosl(k_fm * (long double)r_ip1);

            ss += s * s;
            sc += s * c;
            cc += c * c;
            su += s * u_ip1;
            cu += c * u_ip1;
        }
        u_im1 = u_i;
        u_i = u_ip1;
    }

    det = ss * cc - sc * sc;
    if (!isfinite((double)det) || fabsl(det) < 1.0e-30L) {
        return -3;
    }
    *alpha = (cc * su - sc * cu) / det;
    *beta = (ss * cu - sc * su) / det;
    if (!isfinite((double)*alpha) || !isfinite((double)*beta)) {
        return -3;
    }
    return 0;
}

int ere_delta(double mu_mev, double k_mev, ere_pot_fn v, void *ctx, double *delta)
{
    long double alpha;
    long double beta;
    int rc;

    if (delta == NULL || mu_mev <= 0.0 || v == NULL || ctx == NULL || k_mev <= 0.0) {
        return -1;
    }
    if (!ere_outer_tail(v, ctx)) {
        return -3;
    }
    rc = ere_match(mu_mev, k_mev, v, ctx, &alpha, &beta);
    if (rc != 0) {
        return rc;
    }
    *delta = (double)atan2l(beta, alpha);
    return 0;
}

int ere_swave(double mu_mev, ere_pot_fn v, void *ctx, ere_t *out)
{
    long double sx = 0.0L;
    long double sy = 0.0L;
    long double sxx = 0.0L;
    long double sxy = 0.0L;
    long double den;
    long double slope;
    long double icept;
    int n_free = 0;
    int j;

    if (out == NULL || mu_mev <= 0.0 || v == NULL || ctx == NULL) {
        return -1;
    }
    if (!ere_outer_tail(v, ctx)) {
        return -3;
    }
    out->a0_fm = 0.0;
    out->r0_fm = 0.0;

    for (j = 1; j <= ERE_NK; j++) {
        double k_mev = j * ERE_K_STEP;
        long double k_fm = (long double)k_mev / ERE_HBARC;
        long double alpha;
        long double beta;
        long double x;
        long double y;
        int rc = ere_match(mu_mev, k_mev, v, ctx, &alpha, &beta);

        if (rc != 0) {
            return rc;
        }
        /* |δ| at the roundoff floor: the potential does not move the wave. */
        if (fabsl(beta) <= 1.0e-10L * (fabsl(alpha) + fabsl(beta))) {
            n_free++;
            continue;
        }
        x = k_fm * k_fm;
        y = k_fm * alpha / beta;
        sx += x;
        sy += y;
        sxx += x * x;
        sxy += x * y;
    }
    if (n_free == ERE_NK) {
        return 0;
    }
    if (n_free != 0) {
        return -3;
    }

    den = (long double)ERE_NK * sxx - sx * sx;
    if (!isfinite((double)den) || fabsl(den) < 1.0e-30L) {
        return -3;
    }
    slope = ((long double)ERE_NK * sxy - sx * sy) / den;
    icept = (sy - slope * sx) / (long double)ERE_NK;
    if (!isfinite((double)slope) || !isfinite((double)icept) || fabsl(icept) < 1.0e-30L) {
        return -3;
    }
    out->a0_fm = (double)(1.0L / icept);
    out->r0_fm = (double)(2.0L * slope);
    return 0;
}
