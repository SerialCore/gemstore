/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef GEMSTORE_MODEL_MODEL
#define GEMSTORE_MODEL_MODEL

typedef enum model_type {
    MODEL_NR_SCREEN,
    MODEL_GI_STRING,
    MODEL_GI_SCREEN,
    MODEL_GI_QUADRA_SCREEN
} model_type_t;

typedef struct argsModel {
    double mn;              /* mass of n */
    double ms;              /* mass of s */
    double mc;              /* mass of c */
    double mb;              /* mass of b */
    double mt;              /* mass of t */

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
    double epsilon_sonu;    /* GI smearing parameter for spin-orbit */
    double epsilon_sos;     /* GI smearing parameter for Thomas */
    double epsilon_tens;    /* GI smearing parameter for tensor */

    double alpha1;          /* GI smearing parameter for coupling constant */
    double alpha2;          /* GI smearing parameter for coupling constant */
    double alpha3;          /* GI smearing parameter for coupling constant */
    double gamma1;          /* GI smearing parameter for coupling constant */
    double gamma2;          /* GI smearing parameter for coupling constant */
    double gamma3;          /* GI smearing parameter for coupling constant */
} argsModel_t;

typedef struct argsSOC {
    double OCent;
    double OSdS;
    double OLS1;
    double OLS2;
    double OTens;
} argsSOC_t;

typedef struct argsFlavor {
    double m1;
    double m2;
} argsFlavor_t;

typedef double (*potential_t)(double x, const argsFlavor_t *args_flavor, const argsSOC_t *args_soc, const argsModel_t *args_model);

#endif