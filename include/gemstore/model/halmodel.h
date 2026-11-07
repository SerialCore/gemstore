/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef GEMSTORE_MODEL_HALMODEL
#define GEMSTORE_MODEL_HALMODEL

#include <gemstore/param/argset.h>

/* Pack HAL QCD static + dynamic parameters for the integral callbacks. */
typedef struct hal_pot_ctx {
    const argsHALModel_t *model;
    const argsHALModelDy_t *dyn;
} hal_pot_ctx_t;

/* Kinetic energy m1 + m2 + p²/(2μ). */
double HALKin(double p, void *ctx);

/* V0(r), the spin-independent central piece. */
double HALVc(double r, void *ctx);

/* Vσ(r), the coefficient of σ1·σ2. The caller puts ⟨σ1·σ2⟩ in OSdS. */
double HALVss(double r, void *ctx);

/* VT(r), the coefficient of the Pauli tensor S12. The caller puts ⟨S12⟩ in OTens. */
double HALVten(double r, void *ctx);

/* S-wave spin-singlet combination V0 - 3 Vσ, Eq. (3.40). r is in fm. Returns GeV. */
double HALV1S0(const argsHALModel_t *model, double rfm);

#endif
