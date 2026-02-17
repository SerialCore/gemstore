/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef GEMSTORE_NUMERICAL_MATRIX
#define GEMSTORE_NUMERICAL_MATRIX

typedef struct matrix {
	int row;      		/* count of rows */
	int col;      		/* count of cols */
	double **value; 	/* value of elements */
} matrix_t;

typedef struct array {
	int len;      		/* length of rows */
	double *value; 		/* value of elements */
} array_t;

/* Initialize a matrix with given dimensions */
matrix_t matrix_init(int row, int col);

/* Initialize a matrix with random values */
matrix_t matrix_random(int row, int col);

/* Copy a matrix to a new one */
void matrix_copy(const matrix_t *mat, matrix_t *nmat);

/* Calculate the inverse of a matrix */
void matrix_inverse(const matrix_t *mat, matrix_t *imat);

/* Transpose a matrix */
void matrix_transpose(const matrix_t *mat, matrix_t *tmat);

/* Sum two matrices */
void matrix_sum(const matrix_t *matA, const matrix_t *matB, matrix_t *matC);

/* Product of two matrices */
void matrix_product(const matrix_t *matA, const matrix_t *matB, matrix_t *matC);

/* Product of three matrices: ABA^T */
void matrix_productT(const matrix_t *matA, const matrix_t *matB, matrix_t *matC);

/* Product of three matrices */
void matrix_product3(const matrix_t *matA, const matrix_t *matB, const matrix_t *matC, matrix_t *matD);

/* Calculate the trace of a square matrix */
double matrix_trace(const matrix_t *mat);

/* Calculate the norm of a matrix */
double matrix_norm(const matrix_t *mat);

/* Product a factor to the matrix */
void matrix_multi(double factor, matrix_t *matA);

/* Print the matrix */
void matrix_print(const matrix_t *mat);

/* Free the memory allocated for a matrix */
void matrix_free(matrix_t *mat);

/* Initialize an array with given lengths */
array_t array_init(int len);

/* Calculate the norm of an array */
double array_norm(const array_t *arr);

/* Print the array */
void array_print(const array_t *ary);

/* Free the memory allocated for an array */
void array_free(array_t *ary);

#endif