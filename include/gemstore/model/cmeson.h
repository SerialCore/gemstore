/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef GEMSTORE_MODEL_CMESON
#define GEMSTORE_MODEL_CMESON

#include <gemstore/param/argset.h>
#include <gemstore/math/matrix.h>

/* GI meson spectra on the GEM basis. Returns eigenvalues and v_len eigenvectors. */
void meson_gimodel_GEM(const argsInput_t *args_input, const argsGIModel_t *args_model, argsGIModelDy_t *args_dynmc,
    array_t *e_out, matrix_t *v_out, int v_len);

/* GI meson spectra on the SHO basis. */
void meson_gimodel_SHO(const argsInput_t *args_input, const argsGIModel_t *args_model, argsGIModelDy_t *args_dynmc,
    array_t *e_out, matrix_t *v_out, int v_len);

/* Non-relativistic meson spectra on the GEM basis. */
void meson_nrmodel_GEM(const argsInput_t *args_input, const argsNRModel_t *args_model, argsNRModelDy_t *args_dynmc,
    array_t *e_out, matrix_t *v_out, int v_len);

/* Non-relativistic meson spectra on the SHO basis. */
void meson_nrmodel_SHO(const argsInput_t *args_input, const argsNRModel_t *args_model, argsNRModelDy_t *args_dynmc,
    array_t *e_out, matrix_t *v_out, int v_len);

#endif
