/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef GEMSTORE_MODEL_RADIUS
#define GEMSTORE_MODEL_RADIUS

#include <gemstore/param/argset.h>
#include <gemstore/math/matrix.h>

/* Compute the root-mean-square radius of mesons with {len} of eigen vectors */
void radius_meson_rms(const argsInput_t *input, const matrix_t *vector, array_t *radius, int len);

#endif