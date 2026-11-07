/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef GEMSTORE_MODEL_EFTMODEL
#define GEMSTORE_MODEL_EFTMODEL

#include <gemstore/param/argset.h>

/* Pack EFT contact static + dynamic parameters for the integral callbacks. */
typedef struct eft_pot_ctx {
    const argsEFTModel_t *model;
    const argsEFTModelDy_t *dyn;
} eft_pot_ctx_t;

/* Kinetic energy m1 + m2 + p²/(2μ). */
double EFTVt(double p, void *ctx);

/* Central contact V(r) = -Ca exp(-(r/Rc)^2), times the unit-operator matrix element.
 * Quadrature passes r in GeV^{-1}. The callback returns GeV. */
double EFTVc(double r, void *ctx);

#endif
