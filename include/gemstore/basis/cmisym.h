/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef GEMSTORE_BASIS_CMISYM
#define GEMSTORE_BASIS_CMISYM

#include <gemstore/basis/intrin.h>
#include <gemstore/input.h>

/* Get single color wave function for given color id */
int cmi_color_wfn(int id, intrin_wfn_t *wfn);

/* Get single spin wave function for given spin configuration */
int cmi_spin_wfn(int nquark, const cmi_spin_t *channel, intrin_wfn_t *wfn);

/* CMI symmetry. 1 if color⊗spin is Pauli-allowed. Equal flavors are symmetric for now. */
int cmi_color_spin_allowed(int color_id, const cmi_spin_t *channel, int nf, const int flavor[]);

/* Get full color and spin wave functions from gemstore input */
int cmi_load_basis(const argsInput_t *input, intrin_wfn_t *color, intrin_wfn_t *spin);

#endif
