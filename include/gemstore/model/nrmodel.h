/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef GEMSTORE_MODEL_NRMODEL
#define GEMSTORE_MODEL_NRMODEL

#include <gemstore/param/argset.h>

/* Pack NR static + dynamic parameters so they can be passed as integral ctx. */
typedef struct nr_pot_ctx {
    const argsNRModel_t *model;
    const argsNRModelDy_t *dyn;
} nr_pot_ctx_t;

/* Kinetic energy Σ (m + p²/2m) */
double NRVt(double p, void *ctx);

/* Spectator quark m + p²/2m; baryon T_λ. NRVt is the two-body pair analogue. */
double NRVt_quark(double p, void *ctx);

/* Coulomb potential */
double NRVcoul(double r, void *ctx);

/* Confining potential, linear or screened */
double NRVconf(double r, void *ctx);

/* Colour contact interaction */
double NRVcont(double r, void *ctx);

/* Color-magnetic spin-orbit, diagonal and cross terms together */
double NRVsov(double r, void *ctx);

/* Thomas precession interaction */
double NRVsos(double r, void *ctx);

/* Colour tensor interaction */
double NRVtens(double r, void *ctx);

#endif
