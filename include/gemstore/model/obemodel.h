/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef GEMSTORE_MODEL_OBEMODEL
#define GEMSTORE_MODEL_OBEMODEL

/* Static elastic OBE pair potential (MODEL_OBESTAT).
 *
 * Same limit as OBEPotential with StaticBaryon -> True (the library default):
 * Dirac mass -> infinity, Vls = 0, no q²/8m_B², Dirac magnetic g/(2m_B) drops.
 * Vector spin–spin and tensor use only f/(2M). Finite-mass Gordon (g+f) is
 * another potential and is not implemented here.
 *
 * V = Vc + Vss (s1·s2) + Vt S12, with s = σ on a Dirac leg and S on a
 * vector-meson leg. Units ħ = c = 1.
 *
 * The two legs are typed, not named as a hadron system:
 *   Dirac + Dirac                      OBE BaryonBaryon
 *   (0-, 0+ or 1-) + (0-, 0+ or 1-)    OBE MesonMeson
 *   Dirac + (0-, 0+ or 1-)             OBE MesonBaryon
 *
 * π, σ, ρ, ω (masses, G-parity, HeavySign) are fixed in obemodel.c.
 * Input is the two legs, Λ, M1, M2, η, isospin factors, and Liu couplings.
 * pion g is the table g; C converts with hidden mπ, fπ.
 * σ, ρ, ω g and f are OBE Coupling / TensorCoupling (electric g, magnetic f).
 * tensor_mass1/2 are Liu Table II M in f/(2M), one per leg (D* Ξ: 1.867 and
 * 1.318 GeV). They are not the D* pole mass, which does not enter V(r).
 * Elastic 0- pion exchange is identically zero. Elastic 0+ uses
 * OBE_LEG_PSEUDOSCALAR: OBELib Scalar and Pseudoscalar legs share the
 * elastic NR vertices.
 *
 * Not implemented, on purpose: P ↔ P* / MassGap, StaticBaryon -> False
 * (Gordon g+f, Vls, q²/8m_B²), AxialVector, dipole or Λ1 ≠ Λ2, extra
 * bosons (η, δ, a1). Those exist in OBEPotential; this C potential is
 * the static elastic πσρω operator.
 */

typedef enum obe_leg_type {
    OBE_LEG_DIRAC,
    OBE_LEG_PSEUDOSCALAR,
    OBE_LEG_VECTOR
} obe_leg_type_t;

typedef struct argsOBEStatModel {
    double lambda;          /* common monopole cutoff */
    double tensor_mass1;    /* M in f/(2M) on leg 1 */
    double tensor_mass2;    /* M in f/(2M) on leg 2 */
    double eta;             /* MM C-parity on π and magnetic vector; 0 is +1 */
    obe_leg_type_t type1;
    obe_leg_type_t type2;
    int cc1;                /* meson charge conjugation (G-parity ζ) */
    int cc2;
    double g_pi1, g_pi2;    /* Liu g; vertex conversion is internal */
    double g_sigma1, g_sigma2;
    double g_rho1, g_rho2, f_rho1, f_rho2;
    double g_omega1, g_omega2, f_omega1, f_omega2;
    double iso_pi;          /* isovector factor; 0 is read as 1 */
    double iso_rho;
} argsOBEStatModel_t;

typedef struct argsOBEStatModelDy {
    double OCent;
    double OSdS;
    double OLSi;
    double OLSj;
    double OTens;
} argsOBEStatModelDy_t;

typedef struct obe_pot_ctx {
    const argsOBEStatModel_t *model;
    const argsOBEStatModelDy_t *dyn;
} obe_pot_ctx_t;

/* Radial pieces without spin matrix elements: V = vc + vss (s1·s2) + vt S12. */
void obe_radial(const argsOBEStatModel_t *model, double r, double *vc, double *vss, double *vt);

double OBEVc(double r, void *ctx);
double OBEVss(double r, void *ctx);
double OBEVtens(double r, void *ctx);
double OBEVls1(double r, void *ctx);
double OBEVls2(double r, void *ctx);

#endif
