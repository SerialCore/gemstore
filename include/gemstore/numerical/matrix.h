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

/* Initialize a matrix with given dimensions */
matrix_t matrix_init(int row, int col);

/* Initialize a matrix with random values */
matrix_t matrix_random(int row, int col);

/* Sum two matrices */
matrix_t matrix_sum(const matrix_t *a, const matrix_t *b);

/* Product of two matrices */
matrix_t matrix_product(const matrix_t *a, const matrix_t *b);

/* Transpose a matrix */
matrix_t matrix_transpose(const matrix_t *a);

/* Calculate the diagonal sum of a square matrix */
double matrix_diagonal(const matrix_t *a);

/* Push a new row to the matrix */
void matrix_push(matrix_t *matrix);

/* Print the matrix */
void matrix_print(const matrix_t *matrix);

/* Free the memory allocated for a matrix */
void matrix_free(matrix_t *matrix);

#endif