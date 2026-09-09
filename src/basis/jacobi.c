/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 *
 * Unweighted Jacobi coordinates and maps between the three pair frames.
 * The production central ME uses jacobi_gaussian_shift; raynal_revai is
 * kept for a possible coefficient-based path.
 */

#include <gemstore/basis/jacobi.h>
#include <gemstore/math/soc.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* (mi,mj) = the pair in channel c; mk = spectator. */
void jacobi_pair_mass(double m1, double m2, double m3, int c,
    double *mi, double *mj, double *mk)
{
    switch (c) {
        case 1: *mi = m1; *mj = m2; *mk = m3; break;
        case 2: *mi = m3; *mj = m1; *mk = m2; break;
        case 3: *mi = m2; *mj = m3; *mk = m1; break;
        default:
            fprintf(stderr, "Error: invalid Jacobi channel %d\n", c);
            *mi = *mj = *mk = 0.0;
            break;
    }
}

double jacobi_mu_rho(double mi, double mj)
{
    return mi * mj / (mi + mj);
}

double jacobi_mu_lam(double mi, double mj, double mk)
{
    return (mi + mj) * mk / (mi + mj + mk);
}

void jacobi_r_map(double m1, double m2, double m3, int from, int to,
    double *alpha, double *beta, double *gamma, double *delta)
{
    /* Unweighted map ρ_from = α ρ_to + β λ_to, λ_from = γ ρ_to + δ λ_to.
     * Rows = `from` (1..3), columns = `to`. |αδ − βγ| = 1. */
    double a[3][3] = {
        { 1.0, -m3 / (m1 + m3), -m3 / (m2 + m3) },
        { -m2 / (m1 + m2), 1.0, -m2 / (m2 + m3) },
        { -m1 / (m1 + m2), -m1 / (m1 + m3), 1.0 }
    };
    double b[3][3] = {
        { 0.0,  1.0, -1.0 },
        { -1.0, 0.0,  1.0 },
        {  1.0, -1.0, 0.0 }
    };
    double g[3][3] = {
        { 0.0,
          -1.0 + (m2 * m3) / ((m1 + m2) * (m1 + m3)),
           1.0 - (m1 * m3) / ((m1 + m2) * (m2 + m3)) },
        {  1.0 - (m2 * m3) / ((m1 + m2) * (m1 + m3)),
           0.0,
          -1.0 + (m1 * m2) / ((m1 + m3) * (m2 + m3)) },
        { -1.0 + (m1 * m3) / ((m1 + m2) * (m2 + m3)),
           1.0 - (m1 * m2) / ((m1 + m3) * (m2 + m3)),
           0.0 }
    };
    double d[3][3] = {
        { 1.0, -m2 / (m1 + m2), -m1 / (m1 + m2) },
        { -m3 / (m1 + m3), 1.0, -m1 / (m1 + m3) },
        { -m3 / (m2 + m3), -m2 / (m2 + m3), 1.0 }
    };

    if (from < 1 || from > 3 || to < 1 || to > 3) {
        fprintf(stderr, "Error: invalid Jacobi map %d -> %d\n", from, to);
        *alpha = *beta = *gamma = *delta = 0.0;
        return;
    }

    *alpha = a[from - 1][to - 1];
    *beta  = b[from - 1][to - 1];
    *gamma = g[from - 1][to - 1];
    *delta = d[from - 1][to - 1];
}

/* Mass-weighted Raynal–Revai angle; unused by the current central ME path. */
double jacobi_angle(double m1, double m2, double m3, int from, int to)
{
    double alpha, beta, gamma, delta;
    double mi_f, mj_f, mk_f, mi_t, mj_t, mk_t;
    double mur_f, mul_t;

    if (from == to) {
        return 0.0;
    }

    jacobi_pair_mass(m1, m2, m3, from, &mi_f, &mj_f, &mk_f);
    jacobi_pair_mass(m1, m2, m3, to, &mi_t, &mj_t, &mk_t);
    jacobi_r_map(m1, m2, m3, from, to, &alpha, &beta, &gamma, &delta);

    mur_f = jacobi_mu_rho(mi_f, mj_f);
    mul_t = jacobi_mu_lam(mi_t, mj_t, mk_t);

    /* Mass-weighted ρ_from ≈ cosφ ρ_to + sinφ λ_to. */
    return atan2(sqrt(mur_f) * beta / sqrt(mul_t),
                 sqrt(mur_f) * alpha / sqrt(jacobi_mu_rho(mi_t, mj_t)));
}

static double legendre_p(int k, double x)
{
    if (k == 0) return 1.0;
    if (k == 1) return x;

    double p0 = 1.0;
    double p1 = x;
    double p2 = 0.0;
    for (int n = 2; n <= k; n++) {
        p2 = ((2.0 * n - 1.0) * x * p1 - (n - 1.0) * p0) / n;
        p0 = p1;
        p1 = p2;
    }
    return p1;
}

/* ⟨lρ' lλ' L | lρ lλ L⟩_φ. Unused by solidharm_central_me. */
double raynal_revai(int lrho, int llam, int L, int lrhop, int llamp, double phi)
{
    if (L < 0 || lrho < 0 || llam < 0 || lrhop < 0 || llamp < 0) {
        return 0.0;
    }
    if (abs(lrho - llam) > L || lrho + llam < L) {
        return 0.0;
    }
    if (abs(lrhop - llamp) > L || lrhop + llamp < L) {
        return 0.0;
    }

    int power = lrho + llam - lrhop - llamp;
    if (power % 2 != 0) {
        return 0.0;
    }

    if (fabs(phi) < 1e-14) {
        return (lrho == lrhop && llam == llamp) ? 1.0 : 0.0;
    }

    double phase = pow(-1.0, power / 2);
    double hat = sqrt((2.0 * lrho + 1.0) * (2.0 * llam + 1.0)
        * (2.0 * lrhop + 1.0) * (2.0 * llamp + 1.0));
    double cphi = cos(phi);
    double sum = 0.0;
    int kmax = lrho + lrhop;
    int kmax2 = llam + llamp;
    if (kmax2 < kmax) kmax = kmax2;

    for (int k = 0; k <= kmax; k++) {
        double cg1 = clebsch_gordan(lrho, 0.0, k, 0.0, lrhop, 0.0);
        double cg2 = clebsch_gordan(llam, 0.0, k, 0.0, llamp, 0.0);
        if (cg1 == 0.0 || cg2 == 0.0) {
            continue;
        }
        sum += (2.0 * k + 1.0) * cg1 * cg2
            * sixJ_symbol(lrho, llam, L, llamp, lrhop, k)
            * legendre_p(k, cphi);
    }

    return phase * hat * sum;
}

double jacobi_gaussian_overlap(double m1, double m2, double m3,
    int ca, double nurho_a, double nulam_a,
    int cb, double nurho_b, double nulam_b)
{
    double alpha, beta, gamma, delta;
    double arr, arR, aRR, det;

    /* Express channel-b coordinates in channel a, then complete the square. */
    jacobi_r_map(m1, m2, m3, cb, ca, &alpha, &beta, &gamma, &delta);

    arr = nurho_a + nurho_b * alpha * alpha + nulam_b * gamma * gamma;
    arR = 2.0 * (nurho_b * alpha * beta + nulam_b * gamma * delta);
    aRR = nulam_a + nurho_b * beta * beta + nulam_b * delta * delta;
    det = arr * aRR - 0.25 * arR * arR;
    if (det <= 0.0) {
        return 0.0;
    }

    return pow(M_PI, 3.0) / pow(det, 1.5);
}

int jacobi_gaussian_reduce(
    double m1, double m2, double m3,
    int from_a, double nurho_a, double nulam_a,
    int from_b, double nurho_b, double nulam_b,
    int pair,
    double *b11, double *aRR)
{
    double aa, ba, ga, da;
    double ab, bb, gb, db;
    double arr, arR, aRRv;

    jacobi_r_map(m1, m2, m3, from_a, pair, &aa, &ba, &ga, &da);
    jacobi_r_map(m1, m2, m3, from_b, pair, &ab, &bb, &gb, &db);

    arr = nurho_a * aa * aa + nulam_a * ga * ga
        + nurho_b * ab * ab + nulam_b * gb * gb;
    arR = 2.0 * (nurho_a * aa * ba + nulam_a * ga * da
        + nurho_b * ab * bb + nulam_b * gb * db);
    aRRv = nurho_a * ba * ba + nulam_a * da * da
         + nurho_b * bb * bb + nulam_b * db * db;
    if (aRRv <= 0.0) {
        return 0;
    }

    *aRR = aRRv;
    *b11 = arr - arR * arR / (4.0 * aRRv);
    return (*b11 > 0.0);
}

int jacobi_gaussian_shift(
    double m1, double m2, double m3,
    int from_a, double nurho_a, double nulam_a,
    int from_b, double nurho_b, double nulam_b,
    int pair,
    jacobi_shift_t *out)
{
    double aa, ba, ga, da;
    double ab, bb, gb, db;
    double arr, arR, aRRv, kappa;

    jacobi_r_map(m1, m2, m3, from_a, pair, &aa, &ba, &ga, &da);
    jacobi_r_map(m1, m2, m3, from_b, pair, &ab, &bb, &gb, &db);

    arr = nurho_a * aa * aa + nulam_a * ga * ga
        + nurho_b * ab * ab + nulam_b * gb * gb;
    arR = 2.0 * (nurho_a * aa * ba + nulam_a * ga * da
        + nurho_b * ab * bb + nulam_b * gb * db);
    aRRv = nurho_a * ba * ba + nulam_a * da * da
         + nurho_b * bb * bb + nulam_b * db * db;
    if (aRRv <= 0.0) {
        return 0;
    }

    out->aRR = aRRv;
    out->b11 = arr - arR * arR / (4.0 * aRRv);
    if (out->b11 <= 0.0) {
        return 0;
    }

    /* R = R' − κ r  ⇒  ρ = (α − βκ) r + β R' */
    kappa = arR / (2.0 * aRRv);
    out->al_a = aa - ba * kappa;
    out->be_a = ba;
    out->ga_a = ga - da * kappa;
    out->de_a = da;
    out->al_b = ab - bb * kappa;
    out->be_b = bb;
    out->ga_b = gb - db * kappa;
    out->de_b = db;
    return 1;
}
