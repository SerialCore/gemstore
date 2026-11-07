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

typedef struct argsQMCModel {
    double mn;              /* mass of n quark */
    double ms;              /* mass of s quark */
    double mc;              /* mass of c quark */
    double mb;              /* mass of b quark */
    double B;               /* bag constant */
    double Z;               /* ZPE constant */
    double m_sigma;         /* mass of sigma meson */
    double m_omega;         /* mass of omega meson */
    double m_rho;           /* mass of rho meson */
    double g_sigma;         /* coupling of sigma meson */
    double g_omega;         /* coupling of omega meson */
    double g_rho;           /* coupling of rho meson */
    double lambda3;         /* scalar self-coupling; 0 if absent */
    double g_sigma_B;       /* coupling of sigma meson in bag constant */
    double delta;           /* exponent of bag constant */
} argsQMCModel_t;

typedef struct argsQMCModelDy {
    model_type_t model;
    double sigma;           /* scalar mean field, GeV */
    double omega;           /* isoscalar vector mean field, GeV */
    double rho3;            /* isovector vector mean field, GeV */
    double B;               /* bag constant at this field, GeV^4 */
    double mn_star;         /* in-medium light-quark mass, GeV */
} argsQMCModelDy_t;

typedef struct argsOBEModel {
    double lambda;          /* common monopole cutoff, GeV */
    double tensor_mass1;    /* M in f/(2M) on leg 1, not pole mass */
    double tensor_mass2;    /* M in f/(2M) on leg 2, not pole mass */
    double g_pi1, g_pi2;    /* g couplings for pion vertex */
    double g_sigma1, g_sigma2;
    double g_rho1, g_rho2, f_rho1, f_rho2;       /* f = κ g */
    double g_omega1, g_omega2, f_omega1, f_omega2;
    obe_leg_type_t type1;   /* DIRAC 1/2, VECTOR 1-, PSEUDOSCALAR 0- */
    obe_leg_type_t type2;
    int cc1;                /* meson C-parity; 0 no G flip, nonzero flips G=-1 exchanges for π,ω; Dirac ignored */
    int cc2;
    double eta;             /* MM C-parity on π and magnetic vector; write +1 or −1; 0 is read as +1; MB/BB ignored */
    double iso_pi;          /* ⟨τ1·τ2⟩ on π; I=0 → −3, I=1 → +1; 0 → +1 */
    double iso_rho;         /* ⟨τ1·τ2⟩ on ρ; same as iso_pi; σ/ω always 1 */
} argsOBEModel_t;

typedef struct argsOBEModelDy {
    double mi;              /* physical mass of hadron 1 */
    double mj;              /* physical mass of hadron 2 */
    double OCent;
    double OSdS;
    double OLSi;
    double OLSj;
    double OTens;
} argsOBEModelDy_t;

/* Local Yukawa, V = Σ (V_i e^{-μ_i r}/r) (W_i + B_i P_σ), P_σ=(1+σ1·σ2)/2.
 * V_i is in MeV·fm and μ_i in fm^{-1}. Radial callbacks convert to GeV. */
#define YKW_TERM_MAX 8

typedef struct argsYKWTerm {
    double V;               /* MeV·fm */
    double mu;              /* fm^{-1} */
    double wigner;          /* W, coefficient of 1 */
    double bartlett;        /* B, coefficient of P_σ */
} argsYKWTerm_t;

typedef struct argsYKWModel {
    int nterm;
    argsYKWTerm_t term[YKW_TERM_MAX];
} argsYKWModel_t;

typedef struct argsYKWModelDy {
    double mi;
    double mj;
    double OCent;
    double OSdS;
} argsYKWModelDy_t;

typedef struct argsEFTModel {
    double Ca;              /* contact amplitude in MeV */
    double Rc;              /* cut-off radius in fm */
} argsEFTModel_t;

typedef struct argsEFTModelDy {
    double mi;
    double mj;
    double OCent;
} argsEFTModelDy_t;

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

/* Free-space MIT bag plus physical ω, ρ and m_σ = 0.550 GeV.
 * Medium couplings are 0. QMCMOD cannot use this set. */
extern const argsQMCModel_t argsQMC_default;

extern const argsYKWModel_t argsYUKAWA_MT13;

/* Get GI model parameters from input */
argsGIModel_t argsGIModel_from(const argsInput_t *input);

/* Get NR model parameters from input */
argsNRModel_t argsNRModel_from(const argsInput_t *input);

/* Get MIT bag model parameters from input */
argsBAGModel_t argsBAGModel_from(const argsInput_t *input);

/* Get QMC parameters. QMCSTAD requires g_sigma_B = 0.
 * QMCMOD requires g_sigma_B ≠ 0 and delta > 0. */
argsQMCModel_t argsQMCModel_from(const argsInput_t *input);

/* Full OBE args for compute. MOLECUL2 ignores p and uses param_file plus type1/2.
 * MOLECUL3 uses mol3_pairs[p] (p = 0,1,2 → pairs 12, 13, 23) and the two legs of that pair. */
argsOBEModel_t argsOBEModel_from(const argsInput_t *input, int p);

/* Yukawa args. p is the MOLECUL3 pair index; MOLECUL2 ignores p. */
argsYKWModel_t argsYKWModel_from(const argsInput_t *input, int p);

/* EFT contact args. p is the MOLECUL3 pair index; MOLECUL2 ignores p. */
argsEFTModel_t argsEFTModel_from(const argsInput_t *input, int p);

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
