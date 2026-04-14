/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef GEMSTORE_PRINT
#define GEMSTORE_PRINT

#include <gemstore/math/matrix.h>
#include <gemstore/param/argset.h>

/* Write {len} of meson spectra to a file, including mass, RMS radius, and eigenvectors */
int write_meson_spectra(const argsInput_t *input, const array_t *mass, const array_t *radius, const matrix_t *vector, int len);

#endif
