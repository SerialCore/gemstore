/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef GEMSTORE_BASIS_ISOSPIN
#define GEMSTORE_BASIS_ISOSPIN

#include <gemstore/basis/intrin.h>

/* Get isospin basis */
intrin_wfn_t isospin_basis(double st, const char config);

/* Get isospin wave function of meson state */
intrin_wfn_t isospin_wfn_meson(double st, double st3);

/* Get isospin wave function of diquark state */
intrin_wfn_t isospin_wfn_diquark(double it, double it3);

/* Get isospin wave function of baryon state */
intrin_wfn_t isospin_wfn_baryon(double s12, double st, double st3);

/* Get isospin wave function of tetraquark state */
intrin_wfn_t isospin_wfn_tetra(double s12, double s34, double st, double st3);

/* Get isospin wave function of pentaquark state */
intrin_wfn_t isospin_wfn_penta(double s12, double s123, double s45, double st, double st3);

/* Get isospin wave function of hexaquark state */
intrin_wfn_t isospin_wfn_hexa(double s12, double s123, double s45, double s456, double st, double st3);

#endif