/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef GEMSTORE_MODEL_MESONCRG
#define GEMSTORE_MODEL_MESONCRG

#include <gemstore/param/argset.h>
#include <gemstore/math/matrix.h>

/* Calculate meson spectra using CRG (Complex-Range Gaussian / Hiyama's method) */
void spectra_meson_CRG(const argsInput_t *args_input, const argsGIModel_t *args_model, argsGIModelDy_t *args_dynmc,
    array_t *e_out, matrix_t *v_out, int v_len);

/* Compute the root-mean-square radius of mesons with {len} of eigen vectors using CRG basis */
void radius_meson_CRG(const argsInput_t *input, const matrix_t *vector, array_t *radius, int len);

#endif