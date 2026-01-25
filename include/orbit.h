/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef ORBIT_H
#define ORBIT_H

#include <math.h>
#include <complex.h>

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
    double L_curr;
    for (int m = 2; m <= k; m++) {
        L_curr = ((2.0 * m - 1.0 + alpha - x) * L_prev1 - (m - 1.0 + alpha) * L_prev2) / m;
        L_prev2 = L_prev1;
        L_prev1 = L_curr;
    }

    return L_curr;
}

/* Gaussian basis in coordinate space */
static inline double GRnlr(double r, int n, int l, double nu)
{
    double pre_factor = pow(2.0, l/2.0 + 1.25) * pow(nu, l/2.0 + 0.75);
    double gamma_term = sqrt(1 / tgamma(l + 1.5));
    double exp_term = exp(-r * r * nu);

    return pre_factor * pow(r, l) * gamma_term * exp_term;
}

/* Gaussian basis in coordinate space without exponential */
static inline double GRnlr_nonexp(double r, int n, int l, double nu)
{
    double pre_factor = pow(2.0, l/2.0 + 1.25) * pow(nu, l/2.0 + 0.75);
    double gamma_term = sqrt(1 / tgamma(l + 1.5));

    return pre_factor * pow(r, l) * gamma_term; /* * exp(-r * r * nu); */
}

/* Gaussian basis in momentum space */
static inline complex GRnlp(double p, int n, int l, double nu)
{
    complex phase = pow(-I, l);
    double pre_factor = pow(2.0, -l/2.0 - 0.25) * pow(nu, -l/2.0 - 0.75);
    double gamma_term = sqrt(1 / tgamma(l + 1.5));
    double exp_term = exp(-p * p / (4.0 * nu));

    return phase * pre_factor * pow(p, l) * gamma_term * exp_term;
}

/* Gaussian basis in momentum space without exponential */
static inline complex GRnlp_nonexp(double p, int n, int l, double nu)
{
    complex phase = pow(-I, l);
    double pre_factor = pow(2.0, -l/2.0 - 0.25) * pow(nu, -l/2.0 - 0.75);
    double gamma_term = sqrt(1 / tgamma(l + 1.5));

    return phase * pre_factor * pow(p, l) * gamma_term; /* * exp(-p * p / (4.0 * nu)); */
}

/* Spherical harmonic oscillator basis in coordinate space */
static inline double SRnlr(double r, int n, int l, double beta)
{
    double pre_factor = pow(beta, 1.5);
    double r_beta = r * beta;
    double x = r_beta * r_beta;
    double gamma_term = sqrt((2.0 * tgamma(n + 1.0)) / tgamma(n + l + 1.5));
    double exp_term = exp(-0.5 * x);
    double lag = laguerrel(n, l + 0.5, x);

    return pre_factor * pow(r_beta, l) * gamma_term *  exp_term * lag;
}

/* Spherical harmonic oscillator basis in coordinate space without exponential */
static inline double SRnlr_nonexp(double r, int n, int l, double beta)
{
    double pre_factor = pow(beta, 1.5);
    double r_beta = r * beta;
    double x = r_beta * r_beta;
    double gamma_term = sqrt((2.0 * tgamma(n + 1.0)) / tgamma(n + l + 1.5));
    double lag = laguerrel(n, l + 0.5, x);

    return pre_factor * pow(r_beta, l) * gamma_term * lag; /* * exp(-r * r * beta * beta / 2) */
}

/* Spherical harmonic oscillator basis in momentum space */
static inline complex SRnlp(double p, int n, int l, double beta)
{
    complex phase = pow(-1.0, n) * pow(-I, l);
    double pre_factor = 1 / pow(beta, 1.5);
    double p_beta = p / beta;
    double x = p_beta * p_beta;
    double gamma_term = sqrt((2.0 * tgamma(n + 1.0)) / tgamma(n + l + 1.5));
    double exp_term = exp(-0.5 * x);
    double lag = laguerrel(n, l + 0.5, x);
    
    return phase * pre_factor * pow(p_beta, l) * gamma_term * exp_term * lag;
}

/* Spherical harmonic oscillator basis in momentum space without exponential */
static inline complex SRnlp_nonexp(double p, int n, int l, double beta)
{
    complex phase = pow(-1.0, n) * pow(-I, l);
    double pre_factor = 1 / pow(beta, 1.5);
    double p_beta = p / beta;
    double x = p_beta * p_beta;
    double gamma_term = sqrt((2.0 * tgamma(n + 1.0)) / tgamma(n + l + 1.5));
    double lag = laguerrel(n, l + 0.5, x);
    
    return phase * pre_factor * pow(p_beta, l) * gamma_term * lag; /* * exp(-p * p / (2 * beta * beta)) */
}

/* Normalized integrals for Gaussian basis in coordinate space
 * Return 1 */
double NormalizedGr();

/* Orthogonal integrals for Gaussian basis in coordinate space
 * Return 0.973337 */
double OrthogonalGr();

/* Normalized integrals for Gaussian basis in momentum space
 * Return 1 */
double NormalizedGp();

/* Orthogonal integrals for Gaussian basis in momentum space
 * Return 0.973337 */
double OrthogonalGp();

/* Normalized integrals for SHO basis in coordinate space
 * Return 1 */
double NormalizedSr();

/* Orthogonal integrals for SHO basis in coordinate space
 * Return 0 */
double OrthogonalSr();

/* Normalized integrals for SHO basis in momentum space
 * Return 1 */
double NormalizedSp();

/* Orthogonal integrals for SHO basis in momentum space
 * Return 0 */
double OrthogonalSp();

#endif