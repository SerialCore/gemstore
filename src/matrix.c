/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gemstore/matrix.h>

#include <stdio.h>
#include <stdlib.h>

matrix_t *matrix_init(int row, int col)
{
	matrix_t *matrix = malloc(sizeof(matrix_t));
	double *p = malloc(row * col * sizeof(double));

	matrix->row = row;
	matrix->col = col;
	matrix->value = p;

	return matrix;
}

void matrix_push(matrix_t *matrix)
{
	double *p = matrix->value;
	int row = matrix->row, col = matrix->col;

	matrix->value = realloc(p, (row + 1) * col * sizeof(double));
	matrix->row++;
}

void matrix_print(matrix_t *matrix)
{
	double *p = matrix->value;
	int row = matrix->row, col = matrix->col;

	for (int j = 0; j < col; j++) {
		printf("   %4d    ", j + 1);
	}
	printf("\n");
	for (int i = 0; i < row; i++) {
		printf("%3d: ", i + 1);
		for (int j = 0; j < col; j++) {
			printf("%10.6f ", p[i*row + j]);
		}
		printf("\n");
	}
	printf("\n");
}

void matrix_free(matrix_t *matrix)
{
	free(matrix->value);
	free(matrix);
}