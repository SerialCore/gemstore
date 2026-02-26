/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef GEMSTORE_MODEL_MODEL
#define GEMSTORE_MODEL_MODEL

#include <math.h>

typedef enum model_type {
    MODEL_NR_SCREEN,
    MODEL_GI_STRING,
    MODEL_GI_SCREEN,
    MODEL_GI_QUADRA
} model_type_t;

typedef enum system_type {
    SYSTEM_MESON,
    SYSTEM_BARYON,
    SYSTEM_MOLECULE
} system_type_t;

typedef struct argsModel {
    double mn;              /* mass of n quark */
    double ms;              /* mass of s quark */
    double mc;              /* mass of c quark */
    double mb;              /* mass of b quark */
    double mt;              /* mass of t quark */

    double alpha_s;         /* strong coupling constant */
    double b1;              /* string tension */
    double b2;              /* surface tension */
    double mu;              /* screen length */
    double c;               /* constant potential */
    double sigma;           /* short-range contribution */
    double sigma_0;         /* GI smearing parameter for sigma */
    double s;               /* GI smearing parameter for sigma */

    double epsilon_Coul;    /* GI smearing parameter for Coulumb */
    double epsilon_cont;    /* GI smearing parameter for contact */
    double epsilon_sov;     /* GI smearing parameter for spin-orbit */
    double epsilon_sos;     /* GI smearing parameter for Thomas */
    double epsilon_tens;    /* GI smearing parameter for tensor */
} argsModel_t;

typedef struct argsModelDy {
    model_type_t model;     /* model type */
    system_type_t system;   /* system type */
    double mi;              /* mass of particle i */
    double mj;              /* mass of particle j */
    double Cij;             /* color factor of pair ij */
    double OCent;           /* operator value of centor potential */
    double OSdS;            /* operator value of spin-spin coupling */
    double OLSi;            /* operator value of orbit-spini coupling */
    double OLSj;            /* operator value of orbit-spinj coupling */
    double OTens;           /* operator value of tensor potential */
    double Sigij;           /* GI smearing parameter sigma_ij */
    double Sigkij[3];       /* GI smearing parameters sigma_k_ij */
} argsModelDy_t;

/* Default meson parameters for model NRScreen */
extern const argsModel_t argsNRScreen_meson;

/* Default meson parameters for model GISstring */
extern const argsModel_t argsGIString_meson;

/* Default meson parameters for model GIScreen */
extern const argsModel_t argsGIScreen_meson;

/* Default meson parameters for model GISQuadra */
extern const argsModel_t argsGIQuadra_meson;

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

typedef double (*potential_t)(double x, const argsModel_t *args_model, const argsModelDy_t *args_dynmc);

/* Kinetic energy for NRScreen */
double NRVt(double p, const argsModel_t *args_model, const argsModelDy_t *args_dynmc);

/* Confining potential for NRScreen */
double NRVconf(double r, const argsModel_t *args_model, const argsModelDy_t *args_dynmc);

/* Contact potential for NRScreen */
double NRVcont(double r, const argsModel_t *args_model, const argsModelDy_t *args_dynmc);

/* Spin-orbit coulping for NRScreen */
double NRVsocm(double r, const argsModel_t *args_model, const argsModelDy_t *args_dynmc);

/* Thomas precession for NRScreen */
double NRVsotp(double r, const argsModel_t *args_model, const argsModelDy_t *args_dynmc);

/* Tenser potential for NRScreen */
double NRVtens(double r, const argsModel_t *args_model, const argsModelDy_t *args_dynmc);

/* Kinetic energy for GIScreen */
double GIVt(double p, const argsModel_t *args_model, const argsModelDy_t *args_dynmc);

/* GI smearing beta_ij for Vcoul */
double GIVbetaijcoul(double p, const argsModel_t *args_model, const argsModelDy_t *args_dynmc);

/* GI smearing delta_ij for Vcont */
double GIVdeltaijcont(double p, const argsModel_t *args_model, const argsModelDy_t *args_dynmc);

/* GI smearing delta_ii for Vsovi */
double GIVdeltaiisov(double p, const argsModel_t *args_model, const argsModelDy_t *args_dynmc);

/* GI smearing delta_jj for Vsovj */
double GIVdeltajjsov(double p, const argsModel_t *args_model, const argsModelDy_t *args_dynmc);

/* GI smearing delta_ij for Vsovij */
double GIVdeltaijsov(double p, const argsModel_t *args_model, const argsModelDy_t *args_dynmc);

/* GI smearing delta_ii for Vsosi */
double GIVdeltaiisos(double p, const argsModel_t *args_model, const argsModelDy_t *args_dynmc);

/* GI smearing delta_jj for Vsosj */
double GIVdeltajjsos(double p, const argsModel_t *args_model, const argsModelDy_t *args_dynmc);

/* GI smearing delta_ij for Vtens */
double GIVdeltaijtens(double p, const argsModel_t *args_model, const argsModelDy_t *args_dynmc);

/* Coulomb potential for GIScreen */
double GIVcoul(double r, const argsModel_t *args_model, const argsModelDy_t *args_dynmc);

/* Confining potential for GIScreen */
double GIVconf(double r, const argsModel_t *args_model, const argsModelDy_t *args_dynmc);

/* Contact potential for GIScreen */
double GIVcont(double r, const argsModel_t *args_model, const argsModelDy_t *args_dynmc);

/* Spin-orbit coulping for GIScreen */
double GIVsovi(double r, const argsModel_t *args_model, const argsModelDy_t *args_dynmc);

/* Spin-orbit coulping for GIScreen */
double GIVsovj(double r, const argsModel_t *args_model, const argsModelDy_t *args_dynmc);

/* Spin-orbit coulping for GIScreen */
double GIVsovij(double r, const argsModel_t *args_model, const argsModelDy_t *args_dynmc);

/* Thomas precession for GIScreen */
double GIVsosi(double r, const argsModel_t *args_model, const argsModelDy_t *args_dynmc);

/* Thomas precession for GIScreen */
double GIVsosj(double r, const argsModel_t *args_model, const argsModelDy_t *args_dynmc);

/* Tenser potential for GIScreen */
double GIVtens(double r, const argsModel_t *args_model, const argsModelDy_t *args_dynmc);

#endif