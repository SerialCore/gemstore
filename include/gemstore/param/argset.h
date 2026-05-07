/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef GEMSTORE_PARAM_ARGSET
#define GEMSTORE_PARAM_ARGSET

#include <gemstore/types.h>

typedef struct argsOrbit {
    int n;                  /* radial number & gaussian parameter */
    int l;                  /* orbital momentum */
    double scale;           /* scale factor, nu for GEM and beta for SHO */
    double param;           /* additional parameter, omega for CRG */
} argsOrbit_t;

typedef struct argsModel {
    model_type_t model;
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

typedef struct argsInput {
    task_type_t task;
    orbit_type_t orbit;
    model_type_t model;
    param_type_t param;
    system_type_t system;
    int print_pot;          /* if print potential */
    int print_wfn;          /* if print wavefunction */
    int f1;                 /* flavor 1 */
    int f2;                 /* flavor 2 */
    int f3;                 /* flavor 3 */
    int f4;                 /* flavor 4 */
    double S;               /* spin momentum S */
    double L;               /* orbit momentum L */
    double jl;              /* orbit momentum jl */
    double J;               /* total momentum J */
    int nmax;               /* Gaussian parameter */
    double rmax;            /* Gaussian parameter */
    double rmin;            /* Gaussian parameter */
    double beta;            /* harmonic oscillator parameter */
    double omega;           /* complex-range Gaussian parameter */
    char project[256];      /* project name */
    char param_file[256];   /* parameter file name */
} argsInput_t;

/* Default meson parameters for model GISstring */
extern const argsGIModel_t argsGIString_meson;

/* Default heavy meson parameters for model GIScreen */
extern const argsGIModel_t argsGIScreen_meson;

/* Default bbbar meson parameters for model GIScreen */
extern const argsGIModel_t argsGIScreen_bbbar;

/* Default ccbar meson parameters for model GIScreen */
extern const argsGIModel_t argsGIScreen_ccbar;

/* Get GI model parameters from input */
argsGIModel_t argsGIModel_from(const argsInput_t *input);

/* Get quark mass from GI model */
static inline double getmq(int index, const argsGIModel_t *args_model)
{
    switch (index) {
        case 1: return args_model->mn;
        case 2: return args_model->ms;
        case 3: return args_model->mc;
        case 4: return args_model->mb;
        default: return args_model->mn;
    }
}

#endif
