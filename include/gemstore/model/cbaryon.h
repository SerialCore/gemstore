/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef GEMSTORE_MODEL_CBARYON
#define GEMSTORE_MODEL_CBARYON

#include <gemstore/math/matrix.h>
#include <gemstore/param/argset.h>

void spectra_baryon_GEM(const argsInput_t *input, array_t *e_out, matrix_t *v_out, matrix_t *n_out, int v_len,
    int *basis_len_out);

#endif
