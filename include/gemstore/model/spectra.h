/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef GEMSTORE_MODEL_SPECTRA
#define GEMSTORE_MODEL_SPECTRA

#include <gemstore/math/matrix.h>
#include <gemstore/param/argset.h>

/* Calculate meson spectra in GIScreen model and return the eigenvalues and {v_len} of eigenvectors */
void spectra_meson_GEM(const argsInput_t *args_input, const argsGIModel_t *args_model, argsGIModelDy_t *args_dynmc,
    array_t *e_out, matrix_t *v_out, int v_len);

/* Calculate meson spectra using CRG (Complex-Range Gaussian / Hiyama's method) */
void spectra_meson_CRG(const argsInput_t *args_input, const argsGIModel_t *args_model, argsGIModelDy_t *args_dynmc,
    array_t *e_out, matrix_t *v_out, int v_len);

#endif
