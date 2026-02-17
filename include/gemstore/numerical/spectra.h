/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef GEMSTORE_NUMERICAL_SPECTRA
#define GEMSTORE_NUMERICAL_SPECTRA

#include <gemstore/model/model.h>
#include <gemstore/numerical/matrix.h>

/* Calculate meson spectra in NRScreen model */
void spectra_meson_NRScreen(int f1, int f2, int S, int L, int J, int nmax, double rmax, double rmin, array_t *e_out, matrix_t *v_out, int v_len, argsModel_t *params);

/* Calculate meson spectra in GIScreen model */
void spectra_meson_GIScreen(int f1, int f2, int S, int L, int J, int nmax, double rmax, double rmin, array_t *e_out, matrix_t *v_out, int v_len, argsModel_t *params);

#endif
