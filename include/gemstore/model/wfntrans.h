/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef GEMSTORE_MODEL_WFNTRANS
#define GEMSTORE_MODEL_WFNTRANS

#include <gemstore/input.h>
#include <gemstore/math/matrix.h>

/* Get the normalization factor for a given state vector. */
double get_normalized_factor(const argsInput_t *input, const double *vector);

/* Evaluate the overlap-normalized radial wave function at radius r in fm. */
double get_state_wfn_value(const argsInput_t *input, const double *vector, double normalized, double r);

/* Compute effective SHO beta values that reproduce the target RMS radii. */
void get_effective_beta_sho(const argsInput_t *input, const array_t *radius, array_t *ebeta, int len);

#endif
