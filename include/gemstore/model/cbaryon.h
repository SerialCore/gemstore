/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef GEMSTORE_MODEL_CBARYON
#define GEMSTORE_MODEL_CBARYON

#include <gemstore/param/argset.h>
#include <gemstore/math/matrix.h>

/* Baryon SPECTRA on the three Jacobi GEM frames (c=1,2,3).
 * SCDK maps a bra/ket written in any frame onto the pair that a
 * potential depends on (e.g. V13 ← ρ of channel 2). */
void spectra_baryon_GEM(const argsInput_t *args_input, const argsGIModel_t *args_model, argsGIModelDy_t *args_dynmc,
    array_t *e_out, matrix_t *v_out, matrix_t *n_out,
    array_t *rms12, array_t *rms13, array_t *rms23);

#endif
