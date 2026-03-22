/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef GEMSTORE_MODEL_SPECTRA
#define GEMSTORE_MODEL_SPECTRA

#include <gemstore/math/matrix.h>
#include <gemstore/param/argset.h>

/* Calculate meson spectra in GIScreen model */
void spectra_meson_GI(int f1, int f2, int S, int L, int J, int nmax, double rmax, double rmin, 
    const argsGIModel_t *args_model, argsGIModelDy_t *args_dynmc,
    array_t *e_out, matrix_t *v_out, int v_len);

#endif
