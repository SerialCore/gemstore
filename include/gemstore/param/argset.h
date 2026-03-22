/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef GEMSTORE_PARAM_ARGSET
#define GEMSTORE_PARAM_ARGSET

#include <gemstore/types.h>

typedef struct argsModel {
    model_type_t model;     /* model type */
    double mn;              /* mass of n quark */
    double ms;              /* mass of s quark */
    double mc;              /* mass of c quark */
    double mb;              /* mass of b quark */
} argsModel_t;

typedef struct argsGIModel {
    double mn;              /* mass of n quark */
    double ms;              /* mass of s quark */
    double mc;              /* mass of c quark */
    double mb;              /* mass of b quark */

    double b1;              /* string tension */
    double b2;              /* surface tension */
    double mu;              /* screen length */
    double c;               /* constant potential */
    double sigma_0;         /* GI smearing parameter for sigma */
    double s;               /* GI smearing parameter for sigma */

    double epsilon_cont;    /* GI smearing parameter for contact */
    double epsilon_sov;     /* GI smearing parameter for spin-orbit */
    double epsilon_sos;     /* GI smearing parameter for Thomas */
    double epsilon_tens;    /* GI smearing parameter for tensor */
} argsGIModel_t;

typedef struct argsGIModelDy {
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
} argsGIModelDy_t;

/* Default meson parameters for model GISstring */
extern const argsGIModel_t argsGIString_meson;

/* Default heavy meson parameters for model GIScreen */
extern const argsGIModel_t argsGIScreen_meson;

/* Default bbbar meson parameters for model GIScreen */
extern const argsGIModel_t argsGIScreen_bbbar;

/* Default ccbar meson parameters for model GIScreen */
extern const argsGIModel_t argsGIScreen_ccbar;

/* Default light meson parameters for model GIScreen */
extern const argsGIModel_t argsGIScreen_light;

/* Default heavy meson parameters for model GISQuadra */
extern const argsGIModel_t argsGIQuadra_meson;

/* Default bbbar meson parameters for model GISQuadra */
extern const argsGIModel_t argsGIQuadra_bbbar;

/* Default ccbar meson parameters for model GISQuadra */
extern const argsGIModel_t argsGIQuadra_ccbar;

/* Default light meson parameters for model GISQuadra */
extern const argsGIModel_t argsGIQuadra_light;

#endif