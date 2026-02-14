/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gemstore/basis/orbit.h>

#include <complex.h>

double GRnlr(double r, int n, int l, double nu)
{
    double pre_factor = pow(2.0, l/2.0 + 1.25) * pow(nu, l/2.0 + 0.75);
    double gamma_term = sqrt(1 / tgamma(l + 1.5));
    double exp_term = exp(-r * r * nu);

    return pre_factor * pow(r, l) * gamma_term * exp_term;
}

double GRnlr_nonexp(double r, int n, int l, double nu)
{
    double pre_factor = pow(2.0, l/2.0 + 1.25) * pow(nu, l/2.0 + 0.75);
    double gamma_term = sqrt(1 / tgamma(l + 1.5));

    return pre_factor * pow(r, l) * gamma_term; /* * exp(-r * r * nu); */
}

complex GRnlp(double p, int n, int l, double nu)
{
    complex phase = cpow(-I, l);
    double pre_factor = pow(2.0, -l/2.0 - 0.25) * pow(nu, -l/2.0 - 0.75);
    double gamma_term = sqrt(1 / tgamma(l + 1.5));
    double exp_term = exp(-p * p / (4.0 * nu));

    return phase * pre_factor * pow(p, l) * gamma_term * exp_term;
}

complex GRnlp_nonexp(double p, int n, int l, double nu)
{
    complex phase = cpow(-I, l);
    double pre_factor = pow(2.0, -l/2.0 - 0.25) * pow(nu, -l/2.0 - 0.75);
    double gamma_term = sqrt(1 / tgamma(l + 1.5));

    return phase * pre_factor * pow(p, l) * gamma_term; /* * exp(-p * p / (4.0 * nu)); */
}

double SRnlr(double r, int n, int l, double beta)
{
    double pre_factor = pow(beta, 1.5);
    double r_beta = r * beta;
    double x = r_beta * r_beta;
    double gamma_term = sqrt((2.0 * tgamma(n + 1.0)) / tgamma(n + l + 1.5));
    double exp_term = exp(-0.5 * x);
    double lag = laguerrel(n, l + 0.5, x);

    return pre_factor * pow(r_beta, l) * gamma_term *  exp_term * lag;
}

double SRnlr_nonexp(double r, int n, int l, double beta)
{
    double pre_factor = pow(beta, 1.5);
    double r_beta = r * beta;
    double x = r_beta * r_beta;
    double gamma_term = sqrt((2.0 * tgamma(n + 1.0)) / tgamma(n + l + 1.5));
    double lag = laguerrel(n, l + 0.5, x);

    return pre_factor * pow(r_beta, l) * gamma_term * lag; /* * exp(-r * r * beta * beta / 2) */
}

complex SRnlp(double p, int n, int l, double beta)
{
    complex phase = pow(-1.0, n) * cpow(-I, l);
    double pre_factor = 1 / pow(beta, 1.5);
    double p_beta = p / beta;
    double x = p_beta * p_beta;
    double gamma_term = sqrt((2.0 * tgamma(n + 1.0)) / tgamma(n + l + 1.5));
    double exp_term = exp(-0.5 * x);
    double lag = laguerrel(n, l + 0.5, x);
    
    return phase * pre_factor * pow(p_beta, l) * gamma_term * exp_term * lag;
}

complex SRnlp_nonexp(double p, int n, int l, double beta)
{
    complex phase = pow(-1.0, n) * cpow(-I, l);
    double pre_factor = 1 / pow(beta, 1.5);
    double p_beta = p / beta;
    double x = p_beta * p_beta;
    double gamma_term = sqrt((2.0 * tgamma(n + 1.0)) / tgamma(n + l + 1.5));
    double lag = laguerrel(n, l + 0.5, x);
    
    return phase * pre_factor * pow(p_beta, l) * gamma_term * lag; /* * exp(-p * p / (2 * beta * beta)) */
}
