/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef GEMSTORE_INPUT
#define GEMSTORE_INPUT

#include <gemstore/types.h>

#define CMI_COLOR_MAX 8
#define CMI_SPIN_MAX 20
#define MOL3_PAIR_N 3

/* One pair in a three-body molecule: model1 → 1-2, model2 → 1-3, model3 → 2-3. */
typedef struct mol3_pair {
    model_type_t model;
    param_type_t param;
    char param_file[256];
} mol3_pair_t;

/* One spin channel. Unused couplings stay 0. */
typedef struct cmi_spin {
    double s12;
    double s34;
    double s123;
    double s45;
    double s456;
    double st;
    double st3;
} cmi_spin_t;

typedef struct argsInput {
    task_type_t task;
    basis_type_t basis;
    model_type_t model;
    param_type_t param;
    system_type_t system;
    fitting_type_t target;
    int print_pot;          /* if print potential */
    int print_wfn;          /* if print wavefunction */
    int f1;                 /* flavor 1 */
    int f2;                 /* flavor 2 */
    int f3;                 /* flavor 3 */
    int f4;                 /* flavor 4 */
    int f5;                 /* flavor 5 */
    int f6;                 /* flavor 6 */
    int nn;                 /* number of n quarks */
    int ns;                 /* number of s quarks */
    int nc;                 /* number of c quarks */
    int nb;                 /* number of b quarks */
    double S;               /* spin momentum S */
    double L;               /* orbit momentum L */
    double jl;              /* orbit momentum jl */
    double J;               /* total momentum J */
    int P;                  /* parity */
    double m1;              /* physical mass of hadron 1 */
    double m2;              /* physical mass of hadron 2 */
    double m3;              /* physical mass of hadron 3 */
    obe_leg_type_t type1;   /* Dirac / VECTOR / PSEUDOSCALAR of hadron 1 */
    obe_leg_type_t type2;   /* Dirac / VECTOR / PSEUDOSCALAR of hadron 2 */
    obe_leg_type_t type3;   /* Dirac / VECTOR / PSEUDOSCALAR of hadron 3 */
    double iso_I;           /* total isospin I. Not named I (complex.h). */
    double iso_t1;          /* particle isospin 1 */
    double iso_t2;          /* particle isospin 2 */
    double iso_t3;          /* particle isospin 3 */
    int sym12;              /* baryon 1<->2 exchange; MOLECUL3 identity uses type+mass */
    int Lmax;               /* baryon orbital truncation */
    int nmax;               /* Gaussian parameter */
    double rmax;            /* Gaussian parameter */
    double rmin;            /* Gaussian parameter */
    double beta;            /* harmonic oscillator parameter */
    double omega;           /* complex-range Gaussian parameter */
    double rho_b;           /* baryon density, fm^{-3}; BAG with QMCSTAD or QMCMOD */
    double temp;            /* temperature, GeV; BAG with QMCSTAD or QMCMOD */
    int ncolor;             /* number of CMI color channels */
    int nspin;              /* number of CMI spin channels */
    int color_id[CMI_COLOR_MAX];
    cmi_spin_t spin_ch[CMI_SPIN_MAX];
    char project[256];      /* project name */
    char param_file[256];   /* parameter file name */
    mol3_pair_t mol3_pairs[MOL3_PAIR_N];  /* MOLECUL3 pairs 12, 13, 23 */
} argsInput_t;

#endif
