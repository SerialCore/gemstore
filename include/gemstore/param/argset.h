/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef GEMSTORE_PARAM_ARGSET
#define GEMSTORE_PARAM_ARGSET

#include <gemstore/input.h>
#include <gemstore/types.h>

typedef struct argsOrbit {
    int n;                  /* radial number & gaussian parameter */
    int l;                  /* orbital momentum */
    double scale;           /* scale factor, nu for GEM and beta for SHO */
    double param;           /* additional parameter */
} argsOrbit_t;

typedef struct argsModel {
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
    double b;               /* string tension */
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
    model_type_t model;
    system_type_t system;
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

typedef struct argsNRModel {
    double mn;              /* mass of n quark */
    double ms;              /* mass of s quark */
    double mc;              /* mass of c quark */
    double mb;              /* mass of b quark */
    double b;               /* string tension */
    double mu;              /* screen length */
    double c;               /* constant potential */
    double alpha_s;         /* strong coupling */
    double sigma;           /* contact smearing, GeV */
} argsNRModel_t;

typedef struct argsNRModelDy {
    model_type_t model;
    system_type_t system;
    double mi;              /* mass of particle i */
    double mj;              /* mass of particle j */
    double Cij;             /* color factor of pair ij */
    double OCent;           /* operator value of center potential */
    double OSdS;            /* operator value of spin-spin coupling */
    double OLSi;            /* operator value of orbit-spini coupling */
    double OLSj;            /* operator value of orbit-spinj coupling */
    double OTens;           /* operator value of tensor potential */
} argsNRModelDy_t;

typedef struct argsBAGModel {
    double mn;              /* mass of n quark */
    double ms;              /* mass of s quark */
    double mc;              /* mass of c quark */
    double mb;              /* mass of b quark */
    double Bcs;             /* binding energy c-s */
    double Bbs;             /* binding energy b-s */
    double Bcc;             /* binding energy c-c */
    double Bbc;             /* binding energy b-c */
    double Bbb;             /* binding energy b-b */
    double B;               /* bag constant */
    double Z;               /* ZPE constant */
} argsBAGModel_t;

/* Default meson parameters for model GIString */
extern const argsGIModel_t argsGIString_meson;

/* Defualt baryon parameters for model GIString */
extern const argsGIModel_t argsGIString_baryon;

/* Default heavy meson parameters for model GIScreen */
extern const argsGIModel_t argsGIScreen_meson;

/* Default bbbar meson parameters for model GIScreen */
extern const argsGIModel_t argsGIScreen_bbbar;

/* Default bcbar meson parameters for model GIScreen */
extern const argsGIModel_t argsGIScreen_bcbar;

/* Default bsbar meson parameters for model GIScreen */
extern const argsGIModel_t argsGIScreen_bsbar;

/* Default ccbar meson parameters for model GIScreen */
extern const argsGIModel_t argsGIScreen_ccbar;

/* Default csbar meson parameters for model GIScreen */
extern const argsGIModel_t argsGIScreen_csbar;

/* Default heavy meson parameters for model NRString */
extern const argsNRModel_t argsNRString_meson;

/* Default heavy meson parameters for model NRScreen */
extern const argsNRModel_t argsNRScreen_meson;

/* Default MIT bag model parameters */
extern const argsBAGModel_t argsBAGMIT_default;

/* Get GI model parameters from input */
argsGIModel_t argsGIModel_from(const argsInput_t *input);

/* Get NR model parameters from input */
argsNRModel_t argsNRModel_from(const argsInput_t *input);

/* Get MIT bag model parameters from input */
argsBAGModel_t argsBAGModel_from(const argsInput_t *input);

/* Get quark mass from any model */
static inline double getmq(int index, const void *args_model)
{
    argsModel_t *args = (argsModel_t *)args_model;
    switch (index) {
        case 0: return 0.0;
        case 1: return args->mn;
        case 2: return args->ms;
        case 3: return args->mc;
        case 4: return args->mb;
        default: return 0.0;
    }
}

#endif
