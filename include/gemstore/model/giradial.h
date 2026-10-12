/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 *
 * Radial integrals of the baryon GI operators (recycle inteCenV).
 * Each function is ∫_0^∞ x^n exp(-b11 x^2) K(x) dx and is passed to
 * inteVcenPartA. For a coordinate operator x is the pair separation.
 * intep* and inteTi are reached through a momentum map, so x is that
 * momentum. ctx points at scdk_vargs_t.
 *
 * gimodel.c stores the same potentials as V(r) for the meson integral.
 * These kernels do not call it. GISTRING uses the closed erf integrals.
 * GISCREEN with mu > 0 evaluates confinement and the scalar spin-orbit
 * with the 50-point InteCenV rule from recycle, not integral_exp_rn.
 *
 * δ_ij and β_ij are raised to 0.5+ε. GIVdelta* also puts ε^2 in the
 * exponent; that term is absent here. Confinement uses C_ij = -2/3.
 * scdk_vargs_from_model copies the GI parameters and sets f = 0.5, eCoul = 0.
 */

#ifndef GEMSTORE_MODEL_GIRADIAL
#define GEMSTORE_MODEL_GIRADIAL

#include <gemstore/basis/basis.h>
#include <gemstore/param/argset.h>
#include <gemstore/types.h>

/* Masses, the three GI smearing nodes, and the potential parameters.
 * f is the string-smearing fraction kept at 1/2. eCoul is kept at 0,
 * so the Coulomb β_ij is the plain 1/2 power. */
typedef struct {
    double mn, ms, mc, mb;
    double alpha[3];
    double gamma[3];
    double b, c, sigma0, s, f, mu;
    double econt, etens, esov, esos, eCoul;
    model_type_t model;
} scdk_vargs_t;

void scdk_vargs_from_model(scdk_vargs_t *varg, const argsGIModel_t *model, model_type_t mtype);

/* Overlap, K = 1. */
double inteNfi(double b11, int n, basis_3b qf, basis_3b qi, void *ctx, int c);
/* Smeared one-gluon Coulomb, contact, and tensor. */
double inteVogeG(double b11, int n, basis_3b qf, basis_3b qi, void *ctx, int c);
double inteVcont(double b11, int n, basis_3b qf, basis_3b qi, void *ctx, int c);
double inteVtens(double b11, int n, basis_3b qf, basis_3b qi, void *ctx, int c);
/* Vector spin-orbit, one kernel for each spin-angular table. */
double inteVsovii(double b11, int n, basis_3b qf, basis_3b qi, void *ctx, int c);
double inteVsovjj(double b11, int n, basis_3b qf, basis_3b qi, void *ctx, int c);
double inteVsovij(double b11, int n, basis_3b qf, basis_3b qi, void *ctx, int c);
double inteVsovji(double b11, int n, basis_3b qf, basis_3b qi, void *ctx, int c);
/* Confinement, then the scalar spin-orbit it induces on quarks i and j. */
double inteVstring(double b11, int n, basis_3b qf, basis_3b qi, void *ctx, int c);
double inteVsosii(double b11, int n, basis_3b qf, basis_3b qi, void *ctx, int c);
double inteVsosjj(double b11, int n, basis_3b qf, basis_3b qi, void *ctx, int c);
/* β_ij and δ_ij. cbaryon sandwiches these with the momentum map. */
double intepogeG(double b11, int n, basis_3b qf, basis_3b qi, void *ctx, int c);
double intepcont(double b11, int n, basis_3b qf, basis_3b qi, void *ctx, int c);
double inteptens(double b11, int n, basis_3b qf, basis_3b qi, void *ctx, int c);
double intepsovii(double b11, int n, basis_3b qf, basis_3b qi, void *ctx, int c);
double intepsovjj(double b11, int n, basis_3b qf, basis_3b qi, void *ctx, int c);
double intepsovji(double b11, int n, basis_3b qf, basis_3b qi, void *ctx, int c);
double intepsovij(double b11, int n, basis_3b qf, basis_3b qi, void *ctx, int c);
double intepsosii(double b11, int n, basis_3b qf, basis_3b qi, void *ctx, int c);
double intepsosjj(double b11, int n, basis_3b qf, basis_3b qi, void *ctx, int c);
/* √(p²+m²) of one quark. c = 1, 2, 3 selects qf.m1, m2, m3. */
double inteTi(double b11, int n, basis_3b qf, basis_3b qi, void *ctx, int c);
/* Moment of x²: the same Gaussian integral at power n+2. */
double inteRMS(double b11, int n, basis_3b qf, basis_3b qi, void *ctx, int c);

#endif
