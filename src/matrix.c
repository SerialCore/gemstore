/*
 * Copyright (C) 2025, Wen-Xuan Zhang <serialcore@outlook.com>
 * Copyright (C) 2025, Si-Qiang Luo <luosq15@lzu.edu.cn>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <matrix.h>

#include <stdio.h>
#include <stdlib.h>

matrix_t *matrix_init(matrix_t *matrix, int n, int m)
{
	int i;
	double **p;
	p = (double**)malloc(n*sizeof(double*));
	if (p == NULL) {
		printf("error_initmatrix\n");
	}
	for (i = 0; i < n; i++) {
		p[i] = (double*)malloc(m*sizeof(double));
		if (p[i] == NULL) {
			printf("error_initmatrix\n");
		}
	}
	matrix->p = p;
	matrix->n = n;
	matrix->m = m;

	return matrix;
}

void matrix_push(matrix_t *matrix)
{
	matrix->p = (double**)realloc(matrix->p, sizeof(double*)*(matrix->n+1));
	matrix->p[matrix->n] = (double*)malloc(sizeof(double)*matrix->m);
	matrix->n++;
}

void matrix_print(matrix_t *matrix)
{
	double **p = matrix->p;
	int n = matrix->n, m = matrix->m;
	int i, j;
	printf("     ");
	for (j = 0; j < m; j++) {
		printf("   %4d    ", j+1);
	}
	printf("\n");
	for (i = 0; i < n; i++) {
		printf("%3d: ", i+1);
		for (j = 0; j < m; j++) {
			printf("%10.6f ", p[i][j]);
		}
		printf("\n");
	}
	printf("\n");
}

void matrix_free(matrix_t *matrix)
{
	double **p = matrix->p;
	int n = matrix->n;
	int i;
	for (i = 0; i < n; i++) {
		free(p[i]);
	}
	free(p);
	matrix->n = 0;
	matrix->m = 0;
}

void array_print(double *p, int n)
{
	for (int i = 0; i < n; i++) {
		printf("%25.20f ", p[i]);
	}
	printf("\n");
}

void array_free(double *p)
{
	free(p);
}