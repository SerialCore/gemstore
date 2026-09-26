/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 *
 * Shared three-body Jacobi GEM helpers (baryon now, molecules later):
 * pair identity, SCDK angular tables, overcomplete overlap reduction.
 */

#ifndef GEMSTORE_BASIS_THREEBODY
#define GEMSTORE_BASIS_THREEBODY

#include <gemstore/math/sumckdk.h>
#include <gemstore/math/matrix.h>

/* Jacobi channel c: 1=(12), 2=(31), 3=(23). */
int threebody_pair_identical(int id1, int id2, int id3, int c);

/* Recycle 1↔2 phase: f12 * (-1)^{1+sij+lρ}. */
double threebody_exchange_eta(int f12, double sij, int lrho);

void threebody_scdk_table_alloc(int *len_part, int len_list, sumckdk_scdk *****tab);
void threebody_scdk_table_free(int *len_part, int len_list, sumckdk_scdk *****tab);

/*
 * N-orthonormal rows of an overcomplete overlap.
 * Allocates *vt as (rank × n). Returns rank, or 0 if N has no kept modes.
 * Eigenvalues below rel_thresh * λ_max (and below 1e-12) are dropped.
 */
int threebody_overlap_basis(const matrix_t *N, double rel_thresh, matrix_t *vt);

#endif
