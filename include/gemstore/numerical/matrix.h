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
	int con;      		/* count of rows */
	double *value; 		/* value of elements */
} array_t;

/* Initialize a matrix with given dimensions */
matrix_t matrix_init(int row, int col);

/* Initialize a matrix with random values */
matrix_t matrix_random(int row, int col);

/* Transpose a matrix */
matrix_t matrix_transpose(const matrix_t *mat);

/* Sum two matrices */
void matrix_sum(const matrix_t *matA, const matrix_t *matB, matrix_t *matC);

/* Product of two matrices */
void matrix_product(const matrix_t *matA, const matrix_t *matB, matrix_t *matC);

/* Product of three matrices */
void matrix_product3(const matrix_t *matA, const matrix_t *matB, const matrix_t *matC, matrix_t *matD);

/* Calculate the diagonal sum of a square matrix */
double matrix_diagonal(const matrix_t *mat);

/* Push a new row to the matrix */
void matrix_push(matrix_t *mat);

/* Print the matrix */
void matrix_print(const matrix_t *mat);

/* Free the memory allocated for a matrix */
void matrix_free(matrix_t *mat);

/* Print the array */
void array_print(const array_t *ary);

/* Free the memory allocated for an array */
void array_free(array_t *ary);

#endif