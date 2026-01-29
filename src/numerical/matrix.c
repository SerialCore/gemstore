/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gemstore/numerical/matrix.h>

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

matrix_t matrix_init(int row, int col)
{
	matrix_t matrix;

	double **value = (double**)malloc(row*sizeof(double*));
	if (value == NULL) {
		printf("error_initmatrix\n");
		exit(1);
	}
	for (int i = 0; i < row; i++) {
		value[i] = (double*)malloc(col*sizeof(double));
		if (value[i] == NULL) {
			printf("error_initmatrix\n");
			exit(1);
		}
	}
	matrix.value = value;
	matrix.row = row;
	matrix.col = col;

	return matrix;
}

matrix_t matrix_random(int row, int col)
{
	matrix_t matrix = matrix_init(row, col);

	srand(time(NULL));
	for (int i = 0; i < row; i++) {
		for (int j = 0; j < col; j++) {
			matrix.value[i][j] = (rand() / (double)RAND_MAX) * 2.0 - 1.0;
		}
	}

	return matrix;
}

matrix_t matrix_sum(const matrix_t *a, const matrix_t *b)
{
	if (a->row != b->row || a->col != b->col) {
		printf("error_matrix_sum: dimension mismatch\n");
		exit(1);
	}

	matrix_t c = matrix_init(a->row, a->col);
	for (int i = 0; i < a->row; i++) {
		for (int j = 0; j < a->col; j++) {
			c.value[i][j] = a->value[i][j] + b->value[i][j];
		}
	}

	return c;
}

matrix_t matrix_product(const matrix_t *a, const matrix_t *b)
{
	if (a->col != b->row) {
		printf("error_matrix_product: dimension mismatch\n");
		exit(1);
	}

	matrix_t c = matrix_init(a->row, b->col);
	for (int i = 0; i < a->row; i++) {
		for (int j = 0; j < b->col; j++) {
			c.value[i][j] = 0.0;
			for (int k = 0; k < a->col; k++) {
				c.value[i][j] += a->value[i][k] * b->value[k][j];
			}
		}
	}

	return c;
}

matrix_t matrix_transpose(const matrix_t *a)
{
	matrix_t b = matrix_init(a->col, a->row);
	for (int i = 0; i < a->row; i++) {
		for (int j = 0; j < a->col; j++) {
			b.value[j][i] = a->value[i][j];
		}
	}
	return b;
}

double matrix_diagonal(const matrix_t *a)
{
	if (a->row != a->col) {
		printf("error_matrix_diagonal: not a square matrix\n");
		exit(1);
	}

	double diag = 0.0;
	for (int i = 0; i < a->row; i++) {
		diag += a->value[i][i];
	}
	return diag;
}

void matrix_push(matrix_t *matrix)
{
	matrix->value = (double**)realloc(matrix->value, sizeof(double*)*(matrix->row+1));
	matrix->value[matrix->row] = (double*)malloc(sizeof(double)*matrix->col);
	matrix->row++;
}

void matrix_print(const matrix_t *matrix)
{
	double **value = matrix->value;
	int row = matrix->row, col = matrix->col;

	printf("     ");
	for (int j = 0; j < col; j++) {
		printf("   %4d    ", j+1);
	}
	printf("\n");
	for (int i = 0; i < row; i++) {
		printf("%3d: ", i+1);
		for (int j = 0; j < col; j++) {
			printf("%10.6f ", value[i][j]);
		}
		printf("\n");
	}
	printf("\n");
}

void matrix_free(matrix_t *matrix)
{
	double **value = matrix->value;
	int row = matrix->row;

	for (int i = 0; i < row; i++) {
		free(value[i]);
	}
	free(value);
	free(matrix);
}