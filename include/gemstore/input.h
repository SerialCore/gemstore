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
    int sym12;              /* symmetry under 1<->2 exchange */
    int Lmax;               /* baryon orbital truncation */
    int nmax;               /* Gaussian parameter */
    double rmax;            /* Gaussian parameter */
    double rmin;            /* Gaussian parameter */
    double beta;            /* harmonic oscillator parameter */
    double omega;           /* complex-range Gaussian parameter */
    int ncolor;             /* number of CMI color channels */
    int nspin;              /* number of CMI spin channels */
    int color_id[CMI_COLOR_MAX];
    cmi_spin_t spin_ch[CMI_SPIN_MAX];
    char project[256];      /* project name */
    char param_file[256];   /* parameter file name */
} argsInput_t;

#endif
