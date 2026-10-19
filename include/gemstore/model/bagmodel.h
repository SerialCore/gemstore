/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef GEMSTORE_MODEL_BAGMODEL
#define GEMSTORE_MODEL_BAGMODEL

#include <gemstore/input.h>
#include <gemstore/math/matrix.h>
#include <gemstore/param/argset.h>

#include <math.h>

/* α_s(R) = 0.296 / Log[1 + 1/(0.281 R)] as defined */
static const double BAG_ALPHA_NUM = 0.296;
static const double BAG_ALPHA_DEN = 0.281;

/* Flavor order of n, x, w, and mu is n, s, c, b. */
typedef struct bag_mit
{
    double mass;
    double R;
    int n[4];
    double x[4];
    double w[4];
    double mu[4];
    double vac;
    double st;
    int nvec;
    double coeff[CMI_COLOR_MAX * CMI_SPIN_MAX];
    int color_id[CMI_COLOR_MAX * CMI_SPIN_MAX];
    cmi_spin_t spin_ch[CMI_COLOR_MAX * CMI_SPIN_MAX];
} bag_mit_t;

/* α_s(R) = 0.296 / log(1 + 1/(0.281 R)). R is in GeV^{-1}. */
static inline double bag_alphas(double r)
{
    return BAG_ALPHA_NUM / log(1.0 + 1.0 / (BAG_ALPHA_DEN * r));
}

/* Single-quark energy ω = sqrt(m^2 + x^2/R^2), in GeV. */
static inline double bag_kinematic(double m, double x, double r)
{
    return sqrt(m * m + x * x / (r * r));
}

static inline double bag_vacuum(double r, double B, double Z)
{
    return 4.0 * M_PI * r * r * r * B / 3.0 - Z / r;
}

/* Reduced quark magnetic moment μ = (R/6) (4 ω R + 2 m R - 3) / (2 ω R (ω R - 1) + m R). */
static inline double bag_momentum(double m, double x, double r)
{
    double wr = sqrt(m * m + x * x / (r * r)) * r;
    double mr = m * r;
    double num = (r / 6.0) * (4.0 * wr + 2.0 * mr - 3.0);
    double den = 2.0 * wr * (wr - 1.0) + mr;
    return (den != 0.0) ? num / den : NAN;
}

/* Ground-state frequency x(m, R). */
double bag_frequency(double m, double r);

/* Mass at fixed R. CMI uses the nprod (icolor,ispin) subspace of the once-built factor matrices. */
void bag_mass(const argsInput_t *input, const argsBAGModel_t *model, const matrix_t *color, const matrix_t *spin,
                int nprod, const int *icolor, const int *ispin, double radius, int state, bag_mit_t *mit);

/* Minimize mass over R in (1, 10) GeV^{-1}. */
void bag_variant(const argsInput_t *input, const argsBAGModel_t *model, const matrix_t *color, const matrix_t *spin,
                int nprod, const int *icolor, const int *ispin, int state, bag_mit_t *mit);

/* Load NONE/CMI, mix by st, and vary each state. order is a q/Q string, empty for NONE. */
void bag_compute(const argsInput_t *input, const argsBAGModel_t *model, bag_mit_t *state, int *nstate, char *order);

#endif
