/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef GEMSTORE_BASIS_COLOR
#define GEMSTORE_BASIS_COLOR

#include <gemstore/basis/intrin.h>

/* Get color wave functions of meson state */
intrin_wfn_t color_wfn_meson();

/* Get color wave functions of baryon state */
intrin_wfn_t color_wfn_baryon();

/* Get color wave functions of tetraquark state (1⊗1) */
intrin_wfn_t color_wfn_tetra1();

/* Get color wave functions of tetraquark state (8⊗8) */
intrin_wfn_t color_wfn_tetra8();

/* Get color wave functions of pentaquark state (1⊗1) */
intrin_wfn_t color_wfn_penta1();

/* Get color wave functions of pentaquark state (38⊗8) */
intrin_wfn_t color_wfn_penta38();

/* Get color wave functions of pentaquark state (68⊗8) */
intrin_wfn_t color_wfn_penta68();

/* Calculate the orthogonal degree of two color wave functions */
double OrthogonalColor(const intrin_wfn_t *wfn, const intrin_wfn_t *ref);

#endif