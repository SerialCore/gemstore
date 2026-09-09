/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef GEMSTORE_MODEL_CBARYON
#define GEMSTORE_MODEL_CBARYON

#include <gemstore/param/argset.h>
#include <gemstore/math/matrix.h>

/* Baryon SPECTRA in one Jacobi GEM frame (currently c=1 only).
 * Central V12/V13/V23 are included; off-channel L·S, tensor, and GI β(p),δ(p)
 * smearing are not. See doc/baryon-jacobi-gem.md. */
void spectra_baryon_GEM(const argsInput_t *args_input, const argsGIModel_t *args_model, argsGIModelDy_t *args_dynmc,
    array_t *e_out, matrix_t *v_out, matrix_t *n_out);

#endif
