/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef GEMSTORE_NUMERICAL_SPECTRA
#define GEMSTORE_NUMERICAL_SPECTRA

#include <gemstore/numerical/matrix.h>
#include <gemstore/numerical/model.h>

/* Calculate meson spectra in NRScreen model */
void spectra_meson_NR(int f1, int f2, int S, int L, int J, int nmax, double rmax, double rmin, 
    const argsModel_t *args_model, argsModelDy_t *args_dynmc,
    array_t *e_out, matrix_t *v_out, int v_len);

/* Calculate meson spectra in GIScreen model */
void spectra_meson_GI(int f1, int f2, int S, int L, int J, int nmax, double rmax, double rmin, 
    const argsModel_t *args_model, argsModelDy_t *args_dynmc,
    array_t *e_out, matrix_t *v_out, int v_len);

#endif
