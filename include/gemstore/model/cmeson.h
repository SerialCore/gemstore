/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef GEMSTORE_MODEL_CMESON
#define GEMSTORE_MODEL_CMESON

#include <gemstore/param/argset.h>
#include <gemstore/math/matrix.h>

/* GI meson spectra on the GEM basis.
 * GEM overlap is pruned with matrix_orthogonal, same cutoff as baryons.
 * e_out must be allocated with room for nmax; its length becomes the kept rank.
 * v_out, when not NULL, is nmax by nmax on entry and is trimmed to the stored rows.
 * n_out, when not NULL, is allocated here as the overlap in the original basis.
 * rms_out, when not NULL and v_out is set, is allocated here as the RMS radii. */
void meson_gimodel_GEM(const argsInput_t *args_input, const argsGIModel_t *args_model, argsGIModelDy_t *args_dynmc,
    array_t *e_out, matrix_t *v_out, matrix_t *n_out, array_t *rms_out);

/* GI meson spectra on the SHO basis.
 * n_out, when not NULL, is allocated here as the SHO overlap.
 * rms_out, when not NULL and v_out is set, is allocated here as the RMS radii. */
void meson_gimodel_SHO(const argsInput_t *args_input, const argsGIModel_t *args_model, argsGIModelDy_t *args_dynmc,
    array_t *e_out, matrix_t *v_out, matrix_t *n_out, array_t *rms_out);

/* Non-relativistic meson spectra on the GEM basis.
 * H is the sum of operator matrix elements, solved as H c = E N c.
 * n_out, when not NULL, is allocated here as that overlap.
 * rms_out, when not NULL and v_out is set, is allocated here as the RMS radii. */
void meson_nrmodel_GEM(const argsInput_t *args_input, const argsNRModel_t *args_model, argsNRModelDy_t *args_dynmc,
    array_t *e_out, matrix_t *v_out, matrix_t *n_out, array_t *rms_out);

/* Non-relativistic meson spectra on the SHO basis.
 * n_out, when not NULL, is allocated here as the SHO overlap.
 * rms_out, when not NULL and v_out is set, is allocated here as the RMS radii. */
void meson_nrmodel_SHO(const argsInput_t *args_input, const argsNRModel_t *args_model, argsNRModelDy_t *args_dynmc,
    array_t *e_out, matrix_t *v_out, matrix_t *n_out, array_t *rms_out);

#endif
