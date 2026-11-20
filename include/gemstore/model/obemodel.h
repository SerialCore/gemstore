/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef GEMSTORE_MODEL_OBEMODEL
#define GEMSTORE_MODEL_OBEMODEL

#include <gemstore/param/argset.h>

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
 * An antimeson is a meson leg with cc ≠ 0. That flips G = −1 exchanges
 * (π, ω) and leaves σ, ρ, and the isospin factor alone. cc does not
 * replace ⟨τ·T⟩. obe_tau_dot builds that factor from the two isospins
 * and the pair I: a doublet uses τ = 2T, an integer isospin uses T.
 * Two doublets give −3 or +1. Σc (t = 1) with D (t = 1/2) gives −2 at
 * I = 1/2 and +1 at I = 3/2. Two isoscalars give 0.
 *
 * Left out on purpose, and not a missing elastic channel: P ↔ P*
 * transitions and MassGap (μ² = m² − Δ²), and spin-3/2 baryons
 * (Σc*, Σc ↔ Σc*). A Dirac leg is spin 1/2. Also absent, as in the
 * static Liu operator: StaticBaryon -> False (Gordon g+f, Vls,
 * q²/8m_B²), AxialVector, dipole or Λ1 ≠ Λ2, and extra bosons
 * (η, δ, a1).
 *
 * Couplings and Dy live in argset.h (argsOBEStatModel / argsOBEStatModelDy).
 * Legs are obe_leg_type_t in types.h.
 */

typedef struct obe_pot_ctx {
    const argsOBEModel_t *model;
    const argsOBEModelDy_t *dyn;
} obe_pot_ctx_t;

/* Spin of an OBE leg: Dirac 1/2, vector 1, 0± 0. */
double obe_leg_spin(obe_leg_type_t type);

/* Intrinsic parity: Dirac +1, VECTOR (1−) and PSEUDOSCALAR (0−) −1.
 * Elastic 0+ currently shares PSEUDOSCALAR vertices and this parity. */
int obe_leg_parity(obe_leg_type_t type);

/* ⟨τi·τj⟩ for π and ρ. A t = 1/2 leg contributes τ = 2T; any other t
 * contributes the generator T. tij is the pair isospin. */
double obe_tau_dot(double ti, double tj, double tij);

/* Radial pieces without spin matrix elements: V = vc + vss (s1·s2) + vt S12. */
void obe_radial(const argsOBEModel_t *model, double r, double *vc, double *vss, double *vt);

double OBEVt(double p, void *ctx);
double OBEVc(double r, void *ctx);
double OBEVss(double r, void *ctx);
double OBEVtens(double r, void *ctx);
double OBEVls1(double r, void *ctx);
double OBEVls2(double r, void *ctx);

#endif
