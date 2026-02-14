/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef GEMSTORE_BASIS_ORBIT
#define GEMSTORE_BASIS_ORBIT

#include <math.h>
#include <complex.h>

typedef enum orbit_type {
    ORBIT_GEM,
    ORBIT_SHO
} orbit_type_t;

typedef struct argsOrbit
{
    /* radial number & gaussian parameter */
    int n;
    /* orbital momentum */
    int l;
    /* scale factor, nu for GEM and beta for SHO */
    double scale;
} argsOrbit_t;

/* Scale ν[n, nmax, rmax, rmin] = 1/rmin^2 * (rmax/rmin)^((2 - 2n)/(nmax-1)) */
static inline double getnu(int n, int nmax, double rmax, double rmin)
{
    double fm = 5.0676896;
    double r1 = rmin * fm;

    if (nmax <= 1) {
        return 1.0 / (r1 * r1);
    }

    double exponent = (2.0 - 2.0 * n) / (nmax - 1.0);
    return pow(rmax / rmin, exponent) / (r1 * r1);
}

/* Associated Laguerre L_k^alpha(x) (exact recurrence from Mathematica LaguerreL) */
static inline double laguerrel(int k, double alpha, double x)
{
    if (k == 0) return 1.0;
    if (k == 1) return -x + alpha + 1.0;

    double L_prev2 = 1.0;
    double L_prev1 = -x + alpha + 1.0;
    double L_curr = 0.0;
    for (int m = 2; m <= k; m++) {
        L_curr = ((2.0 * m - 1.0 + alpha - x) * L_prev1 - (m - 1.0 + alpha) * L_prev2) / m;
        L_prev2 = L_prev1;
        L_prev1 = L_curr;
    }

    return L_curr;
}

/* define orbit wave function in coordinate space */
typedef double (*orbit_wfn_t)(double x, int n, int l, double scale);

/* define orbit wave function in momentum space */
typedef complex (*orbit_wfn_complex_t)(double x, int n, int l, double scale);

/* Gaussian basis in coordinate space */
double GRnlr(double r, int n, int l, double nu);

/* Gaussian basis in coordinate space without exponential */
double GRnlr_nonexp(double r, int n, int l, double nu);

/* Gaussian basis in momentum space */
complex GRnlp(double p, int n, int l, double nu);

/* Gaussian basis in momentum space without exponential */
complex GRnlp_nonexp(double p, int n, int l, double nu);

/* Spherical harmonic oscillator basis in coordinate space */
double SRnlr(double r, int n, int l, double beta);

/* Spherical harmonic oscillator basis in coordinate space without exponential */
double SRnlr_nonexp(double r, int n, int l, double beta);

/* Spherical harmonic oscillator basis in momentum space */
complex SRnlp(double p, int n, int l, double beta);

/* Spherical harmonic oscillator basis in momentum space without exponential */
complex SRnlp_nonexp(double p, int n, int l, double beta);

#endif