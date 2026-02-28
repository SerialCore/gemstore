/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef GEMSTORE_MATH_CMI
#define GEMSTORE_MATH_CMI

#include <gemstore/basis/intrin.h>
#include <gemstore/math/matrix.h>

/* get matrix element of casimir operator sigma*sigma */
void operator_sigma2(const intrin_wfn_t swv[], int num_state, matrix_t *result);

/* get matrix element of casimir operator lambda*lambda */
void operator_lambda2(const intrin_wfn_t cwv[], int num_state, const char *config, matrix_t *result);

#endif