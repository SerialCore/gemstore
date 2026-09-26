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

/* Copy src into dst (same shape) */
void matrix_copy(matrix_t *dst, const matrix_t *src);

/* (A + A^T)/2 in place */
void matrix_symmetrize(matrix_t *mat);

/* ⟨n|op|n⟩ with bra/ket = row n of vec (length vec->col) */
double matrix_expect(const matrix_t *vec, int n, const matrix_t *op);

/* v ← p v p^T; tmp is workspace of the same shape as v */
void matrix_sandwich(matrix_t *v, const matrix_t *p, matrix_t *tmp);

/* Initialize a matrix with random values between -1 and 1 */
matrix_t matrix_random(int row, int col);

/* Calculate the inverse of a matrix */
void matrix_inverse(const matrix_t *mat, matrix_t *imat);

/* Calculate the inverse of a lower triangular matrix */
void matrix_inverse_lowertri(const matrix_t *mat, matrix_t *imat);

/* Perform Cholesky decomposition S = L * L^T with lower triangular matrix L */
void matrix_cholesky_decomp(const matrix_t *matS, matrix_t *matL);

/* Transpose a matrix */
void matrix_transpose(const matrix_t *mat, matrix_t *tmat);

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