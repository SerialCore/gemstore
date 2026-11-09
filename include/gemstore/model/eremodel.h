/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef GEMSTORE_MODEL_EREMODEL
#define GEMSTORE_MODEL_EREMODEL

/* S-wave effective-range expansion. The potential is supplied by the caller.
 * Coulomb and three-body scattering are absent.
 *
 * Particle-physics ERE, Yuan–Liu Eq. (41) and Miyamoto Eq. (5.5):
 *   k cot δ0 = 1/a0 + (1/2) r0 k^2.
 * a0 > 0 is attraction that has not formed a bound state. a0 < 0 is repulsion.
 * The radial equation is Yuan–Liu Eq. (8) at l = 0. */

typedef double (*ere_pot_fn)(double r_fm, void *ctx);

typedef struct ere {
    double a0_fm;
    double r0_fm;
} ere_t;

/* v(r_fm, ctx) is the projected central potential in MeV.
 * mu_mev is the reduced mass. The solver uses ħc = 197.3 MeV·fm,
 * Numerov on r = 0…80 fm with step 0.01 fm, and a straight-line fit of
 * k cot δ0 over k = 0.05, 0.10, …, 0.50 MeV.
 * Returns 0 on success, -1 on a bad argument, -2 if the Numerov step is
 * singular, -3 if the outer wave does not determine a phase. */
int ere_swave(double mu_mev, ere_pot_fn v, void *ctx, ere_t *out);

/* S-wave phase shift in radians at momentum k_mev. Same grid and ħc as above. */
int ere_delta(double mu_mev, double k_mev, ere_pot_fn v, void *ctx, double *delta);

#endif
