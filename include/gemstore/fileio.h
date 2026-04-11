/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef GEMSTORE_FILEIO
#define GEMSTORE_FILEIO

#include <gemstore/math/matrix.h>
#include <gemstore/param/argset.h>

int write_meson_spectra(const argsInput_t *input, const array_t *mass, const array_t *radius, const matrix_t *vector, int len);

#endif
