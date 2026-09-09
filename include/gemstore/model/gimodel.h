/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef GEMSTORE_MODEL_GIMODEL
#define GEMSTORE_MODEL_GIMODEL

#include <gemstore/param/argset.h>

#include <math.h>

/* Default parameters for running strong couping constant */
extern const double GI_ALPHA_K[3], GI_GAMMA_K[3];

/* Get GI smearing parameter sigma_ij */
static inline double sigma_ij(double mi, double mj, double sigma0, double s);
static inline double sigma_ij(double mi, double mj, double sigma0, double s)
{
    double msum = mi + mj;
    double mprod = mi * mj;
    double frac1 = 4.0 * mprod / (msum * msum);
    double frac2 = 2.0 * mprod / msum;
    double term1 = sigma0 * sigma0 * (0.5 + 0.5 * pow(frac1, 4));
    double term2 = s * s * pow(frac2, 2);

    return sqrt(term1 + term2);
}

/* Get GI smearing parameters sigma_k_ij */
static inline void sigma_k_ij(double sigmaij, double sigmak[3]);
static inline void sigma_k_ij(double sigmaij, double sigmak[3])
{
    for (int k = 0; k < 3; k++) {
        sigmak[k] = GI_GAMMA_K[k] * sigmaij / sqrt(GI_GAMMA_K[k] * GI_GAMMA_K[k] + sigmaij * sigmaij);
    }
}

/* Pack GI static + dynamic parameters so they can be passed as integral ctx. */
typedef struct gi_pot_ctx {
    const argsGIModel_t *model;
    const argsGIModelDy_t *dyn;
} gi_pot_ctx_t;

/* Kinetic energy for GIScreen */
double GIVt(double p, void *ctx);

/* Spectator quark √(mi²+p²); baryon T_λ. GIVt is the two-body pair analogue. */
double GIVt_quark(double p, void *ctx);

/* GI smearing beta_ij for Vcoul */
double GIVbetaijcoul(double p, void *ctx);

/* GI smearing delta_ij for Vcont */
double GIVdeltaijcont(double p, void *ctx);

/* GI smearing delta_ii for Vsovi */
double GIVdeltaiisov(double p, void *ctx);

/* GI smearing delta_jj for Vsovj */
double GIVdeltajjsov(double p, void *ctx);

/* GI smearing delta_ij for Vsovij */
double GIVdeltaijsov(double p, void *ctx);

/* GI smearing delta_ii for Vsosi */
double GIVdeltaiisos(double p, void *ctx);

/* GI smearing delta_jj for Vsosj */
double GIVdeltajjsos(double p, void *ctx);

/* GI smearing delta_ij for Vtens */
double GIVdeltaijtens(double p, void *ctx);

/* Coulomb potential for GIScreen */
double GIVcoul(double r, void *ctx);

/* Confining potential for GIScreen */
double GIVconf(double r, void *ctx);

/* Contact potential for GIScreen */
double GIVcont(double r, void *ctx);

/* Spin-orbit coulping for GIScreen */
double GIVsovi(double r, void *ctx);

/* Spin-orbit coulping for GIScreen */
double GIVsovj(double r, void *ctx);

/* Spin-orbit coulping for GIScreen */
double GIVsovij(double r, void *ctx);

/* Thomas precession for GIScreen */
double GIVsosi(double r, void *ctx);

/* Thomas precession for GIScreen */
double GIVsosj(double r, void *ctx);

/* Tenser potential for GIScreen */
double GIVtens(double r, void *ctx);

#endif