/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef GEMSTORE_MODEL_YKWMODEL
#define GEMSTORE_MODEL_YKWMODEL

#include <gemstore/param/argset.h>

/* Pack Yukawa static + dynamic parameters so they can be passed as integral ctx. */
typedef struct ykw_pot_ctx {
    const argsYKWModel_t *model;
    const argsYKWModelDy_t *dyn;
} ykw_pot_ctx_t;

/* Kinetic energy m1 + m2 + p²/(2μ) */
double YKWVt(double p, void *ctx);

/* Central Yukawa, Wigner plus half of Bartlett */
double YKWVc(double r, void *ctx);

/* Coefficient of σ1·σ2, half of Bartlett. Dirac legs multiply by kss. */
double YKWVss(double r, void *ctx);

#endif
