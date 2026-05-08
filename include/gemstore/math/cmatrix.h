/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef GEMSTORE_MATH_CMATRIX
#define GEMSTORE_MATH_CMATRIX

#include <complex.h>

typedef struct cmatrix {
	int row;      		/* count of rows */
	int col;      		/* count of cols */
	complex **value; 	/* value of elements */
} cmatrix_t;

typedef struct carray {
	int len;      		/* length of rows */
	complex *value; 	/* value of elements */
} carray_t;

/* Initialize a matrix with given dimensions */
cmatrix_t cmatrix_init(int row, int col);

/* Initialize a matrix with random values between -1 and 1 */
cmatrix_t cmatrix_random(int row, int col);

/* Calculate the inverse of a matrix */
void cmatrix_inverse(const cmatrix_t *mat, cmatrix_t *imat);

/* Calculate the inverse of a lower triangular matrix */
void cmatrix_inverse_lowertri(const cmatrix_t *mat, cmatrix_t *imat);

/* Perform Cholesky decomposition S = L * L^T with lower triangular matrix L */
void cmatrix_cholesky_decomp(const cmatrix_t *matS, cmatrix_t *matL);

/* Transpose a matrix */
void cmatrix_transpose(const cmatrix_t *mat, cmatrix_t *tmat);

/* Sum two matrices */
void cmatrix_sum(const cmatrix_t *matA, const cmatrix_t *matB, cmatrix_t *matC);

/* Product of two matrices */
void cmatrix_product(const cmatrix_t *matA, const cmatrix_t *matB, cmatrix_t *matC);

/* Product of three matrices: ABA^T */
void cmatrix_productT(const cmatrix_t *matA, const cmatrix_t *matB, cmatrix_t *matC);

/* Print the matrix */
void cmatrix_print(const cmatrix_t *mat);

/* Free the memory allocated for a matrix */
void cmatrix_free(cmatrix_t *mat);

/* Initialize an array with given lengths */
carray_t carray_init(int len);

/* Print the array */
void carray_print(const carray_t *ary);

/* Free the memory allocated for an array */
void carray_free(carray_t *ary);

#endif