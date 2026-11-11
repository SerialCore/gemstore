/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef GEMSTORE_MODEL_QPCMODEL
#define GEMSTORE_MODEL_QPCMODEL

#include <gemstore/param/argset.h>

/* Partial-wave amplitude of one 3P0 decay A -> B C.
 * Each meson is one SHO. p, beta, and the quark masses are in GeV. n starts at 1.
 * L is the outgoing orbital angular momentum. S is the total spin of B and C.
 * fl1 and fl2 weight diagrams 1-4-3-2 and 3-2-1-4.
 * The amplitude includes gamma and the Ackleh normalization. The constant sqrt(8)
 * in that factor stays. Energies and the charge multiplicity are not in this amplitude.
 */
double qpc_amplitude(const argsQPCModel_t *model, double p,
    const qpc_meson_t *a, const qpc_meson_t *b, const qpc_meson_t *c,
    int L, int S, double fl1, double fl2);

#endif
