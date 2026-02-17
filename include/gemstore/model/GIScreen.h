/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef GEMSTORE_MODEL_GISCREEN
#define GEMSTORE_MODEL_GISCREEN

#include <gemstore/model/model.h>

#include <math.h>

/* Default meson parameters for model GIScreen */
extern const argsModel_t argsGIScreen_meson;

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

/* Kinetic energy for GIScreen */
double GIScreen_T(double p, const argsModel_t *args_model, const argsModelDy_t *args_dynmc);

/* GI smearing beta_ij for Vcoul */
double GIScreen_betaij_coul(double p, const argsModel_t *args_model, const argsModelDy_t *args_dynmc);

/* GI smearing delta_ij for Vcont */
double GIScreen_deltaij_cont(double p, const argsModel_t *args_model, const argsModelDy_t *args_dynmc);

/* GI smearing delta_ii for Vsovi */
double GIScreen_deltaii_sov(double p, const argsModel_t *args_model, const argsModelDy_t *args_dynmc);

/* GI smearing delta_jj for Vsovj */
double GIScreen_deltajj_sov(double p, const argsModel_t *args_model, const argsModelDy_t *args_dynmc);

/* GI smearing delta_ij for Vsovij */
double GIScreen_deltaij_sov(double p, const argsModel_t *args_model, const argsModelDy_t *args_dynmc);

/* GI smearing delta_ii for Vsosi */
double GIScreen_deltaii_sos(double p, const argsModel_t *args_model, const argsModelDy_t *args_dynmc);

/* GI smearing delta_jj for Vsosj */
double GIScreen_deltajj_sos(double p, const argsModel_t *args_model, const argsModelDy_t *args_dynmc);

/* GI smearing delta_ij for Vtens */
double GIScreen_deltaij_tens(double p, const argsModel_t *args_model, const argsModelDy_t *args_dynmc);

/* Coulomb potential for GIScreen */
double GIScreen_Vcoul(double r, const argsModel_t *args_model, const argsModelDy_t *args_dynmc);

/* Confining potential for GIScreen */
double GIScreen_Vconf(double r, const argsModel_t *args_model, const argsModelDy_t *args_dynmc);

/* Contact potential for GIScreen */
double GIScreen_Vcont(double r, const argsModel_t *args_model, const argsModelDy_t *args_dynmc);

/* Spin-orbit coulping for GIScreen */
double GIScreen_Vsovi(double r, const argsModel_t *args_model, const argsModelDy_t *args_dynmc);

/* Spin-orbit coulping for GIScreen */
double GIScreen_Vsovj(double r, const argsModel_t *args_model, const argsModelDy_t *args_dynmc);

/* Spin-orbit coulping for GIScreen */
double GIScreen_Vsovij(double r, const argsModel_t *args_model, const argsModelDy_t *args_dynmc);

/* Thomas precession for GIScreen */
double GIScreen_Vsosi(double r, const argsModel_t *args_model, const argsModelDy_t *args_dynmc);

/* Thomas precession for GIScreen */
double GIScreen_Vsosj(double r, const argsModel_t *args_model, const argsModelDy_t *args_dynmc);

/* Tenser potential for GIScreen */
double GIScreen_Vtens(double r, const argsModel_t *args_model, const argsModelDy_t *args_dynmc);

#endif