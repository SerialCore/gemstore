/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef GEMSTORE_MATH_MATRIX
#define GEMSTORE_MATH_MATRIX

typedef struct matrix {
	int row;      		/* count of rows */
	int col;      		/* count of cols */
	double **value; 	/* value of elements */
} matrix_t;

typedef struct array {
	int len;      		/* length of rows */
	double *value; 		/* value of elements */
} array_t;

/* Initialize a matrix with given dimensions, all elements 0 */
matrix_t matrix_init(int row, int col);

/* Initialize a matrix with random values between -1 and 1 */
matrix_t matrix_random(int row, int col);

/* (A + A^T)/2 in place */
void matrix_symmetrize(matrix_t *mat);

/* Copy src into dst (same shape) */
void matrix_copy(const matrix_t *sr, matrix_t *dst);

/* Calculate the inverse of a matrix */
void matrix_inverse(const matrix_t *mat, matrix_t *imat);

/* Transpose a matrix */
void matrix_transpose(const matrix_t *mat, matrix_t *tmat);

/* ⟨n|op|n⟩ with bra/ket = row n of vec (length vec->col) */
double matrix_expect(const matrix_t *vec, int n, const matrix_t *op);

/* N-orthonormal rows of an overcomplete overlap.
 * Allocates *vt as (rank × n). Returns rank, or 0 if N has no kept modes.
 * Eigenvalues below rel_thresh * λ_max (and below 1e-12) are dropped. */
int matrix_orthogonal(const matrix_t *N, double rel_thresh, matrix_t *vt);

/* Sum two matrices */
void matrix_sum(const matrix_t *matA, const matrix_t *matB, matrix_t *matC);

/* Product of two matrices */
void matrix_product(const matrix_t *matA, const matrix_t *matB, matrix_t *matC);

/* Product of three matrices: ABA^T */
void matrix_productT(const matrix_t *matA, const matrix_t *matB, matrix_t *matC);

/* Print the matrix */
void matrix_print(const matrix_t *mat);

/* Free the memory allocated for a matrix */
void matrix_free(matrix_t *mat);

/* Initialize an array of length len, all elements 0 */
array_t array_init(int len);

/* Print the array */
void array_print(const array_t *ary);

/* Free the memory allocated for an array */
void array_free(array_t *ary);

#endif