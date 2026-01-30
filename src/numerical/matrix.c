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
	matrix_t mat;

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
	mat.value = value;
	mat.row = row;
	mat.col = col;

	return mat;
}

matrix_t matrix_random(int row, int col)
{
	matrix_t mat = matrix_init(row, col);

	srand(time(NULL));
	for (int i = 0; i < row; i++) {
		for (int j = 0; j < col; j++) {
			mat.value[i][j] = (rand() / (double)RAND_MAX) * 2.0 - 1.0;
		}
	}

	return mat;
}

matrix_t matrix_transpose(const matrix_t *mat)
{
	matrix_t matT = matrix_init(mat->col, mat->row);
	for (int i = 0; i < mat->row; i++) {
		for (int j = 0; j < mat->col; j++) {
			matT.value[j][i] = mat->value[i][j];
		}
	}
	return matT;
}

void matrix_sum(const matrix_t *matA, const matrix_t *matB, matrix_t *matC)
{
	if (matA->row != matB->row || matA->col != matB->col
		|| matA->row != matC->row || matA->col != matC->col) {
		printf("error_matrix_sum: dimension mismatch\n");
		exit(1);
	}

	for (int i = 0; i < matA->row; i++) {
		for (int j = 0; j < matA->col; j++) {
			matC->value[i][j] = matA->value[i][j] + matB->value[i][j];
		}
	}
}

void matrix_product(const matrix_t *matA, const matrix_t *matB, matrix_t *matC)
{
	if (matA->col != matB->row || matA->row != matC->row || matB->col != matC->col) {
		printf("error_matrix_product: dimension mismatch\n");
		exit(1);
	}

	for (int i = 0; i < matA->row; i++) {
		for (int j = 0; j < matB->col; j++) {
			matC->value[i][j] = 0.0;
			for (int k = 0; k < matA->col; k++) {
				matC->value[i][j] += matA->value[i][k] * matB->value[k][j];
			}
		}
	}
}

void matrix_product3(const matrix_t *matA, const matrix_t *matB, const matrix_t *matC, matrix_t *matD)
{
	if (matA->col != matB->row || matB->col != matC->row
		|| matA->row != matD->row || matC->col != matD->col) {
		printf("error_matrix_product3: dimension mismatch\n");
		exit(1);
	}

	for (int i = 0; i < matA->row; i++) {
		for (int j = 0; j < matC->col; j++) {
			matD->value[i][j] = 0.0;
			for (int k = 0; k < matA->col; k++) {
				for (int l = 0; l < matB->col; l++) {
					matD->value[i][j] += matA->value[i][k] * matB->value[k][l] * matC->value[l][j];
				}
			}
		}
	}
}

double matrix_diagonal(const matrix_t *mat)
{
	if (mat->row != mat->col) {
		printf("error_matrix_diagonal: not a square matrix\n");
		exit(1);
	}

	double diag = 0.0;
	for (int i = 0; i < mat->row; i++) {
		diag += mat->value[i][i];
	}
	return diag;
}

void matrix_push(matrix_t *mat)
{
	mat->value = (double**)realloc(mat->value, sizeof(double*)*(mat->row+1));
	mat->value[mat->row] = (double*)malloc(sizeof(double)*mat->col);
	mat->row++;
}

void matrix_print(const matrix_t *mat)
{
	double **value = mat->value;
	int row = mat->row, col = mat->col;

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

void matrix_free(matrix_t *mat)
{
	double **value = mat->value;
	int row = mat->row;

	for (int i = 0; i < row; i++) {
		free(value[i]);
	}
	free(value);
}