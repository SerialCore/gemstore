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
	}
	for (int i = 0; i < row; i++) {
		value[i] = (double*)malloc(col*sizeof(double));
		if (value[i] == NULL) {
			printf("error_initmatrix\n");
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

void array_print(const double *a, int n)
{
	for (int i = 0; i < n; i++) {
		printf("%25.20f ", a[i]);
	}
	printf("\n");
}

void array_free(double *a)
{
	free(a);
}