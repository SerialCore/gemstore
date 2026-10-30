/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef GEMSTORE_MODEL_WFNTRANS
#define GEMSTORE_MODEL_WFNTRANS

#include <gemstore/input.h>
#include <gemstore/math/matrix.h>

/* One GEM partial-wave block. coeff has nmax entries. r is in fm. */
double get_gem_channel_wfn(int L, int nmax, double rmax, double rmin, const double *coeff, double r);

/* One SHO partial-wave block. coeff has nmax entries. r is in fm. */
double get_sho_channel_wfn(int L, int nmax, double beta, const double *coeff, double r);

/* Compute effective SHO beta values that reproduce the target RMS radii. */
void get_effective_beta_sho(const argsInput_t *input, const array_t *radius, array_t *ebeta, int len);

#endif
