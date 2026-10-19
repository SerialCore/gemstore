/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef GEMSTORE_READ
#define GEMSTORE_READ

#include <gemstore/math/matrix.h>

/* Read meson spectra from a state JSON written by write_meson_spectra. */
void parse_meson_state(const char *filename, array_t *mass, array_t *radius, matrix_t *eigenvector);

#endif
