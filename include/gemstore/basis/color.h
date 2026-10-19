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
intrin_wfn_t color_wfn_tetra11();

/* Get color wave functions of tetraquark state (8⊗8) */
intrin_wfn_t color_wfn_tetra88();

/* Get color wave functions of tetraquark state (3⊗3) */
intrin_wfn_t color_wfn_tetra33();

/* Get color wave functions of tetraquark state (6⊗6) */
intrin_wfn_t color_wfn_tetra66();

/* Pentaquark [(12)^6 3]^8 (45̄)^8 */
intrin_wfn_t color_wfn_penta688();

/* Pentaquark [(12)^3̄ 3]^8 (45̄)^8 */
intrin_wfn_t color_wfn_penta388();

/* Pentaquark [(12)^3̄ 3]^1 (45̄)^1 */
intrin_wfn_t color_wfn_penta11();

/* Hexaquark [(12)^6 3]^10 [(4̄5̄)^6̄ 6̄]^10 */
intrin_wfn_t color_wfn_hexa1010();

/* Hexaquark [(12)^6 3]^8 [(4̄5̄)^6̄ 6̄]^8 */
intrin_wfn_t color_wfn_hexa6868();

/* Hexaquark [(12)^6 3]^8 [(4̄5̄)^3 6̄]^8 */
intrin_wfn_t color_wfn_hexa6838();

/* Hexaquark [(12)^3̄ 3]^8 [(4̄5̄)^6̄ 6̄]^8 */
intrin_wfn_t color_wfn_hexa3868();

/* Hexaquark [(12)^3̄ 3]^8 [(4̄5̄)^3 6̄]^8 */
intrin_wfn_t color_wfn_hexa3838();

/* Hexaquark [(12)^3̄ 3]^1 [(4̄5̄)^3 6̄]^1 */
intrin_wfn_t color_wfn_hexa11();

#endif