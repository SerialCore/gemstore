/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gemstore/numerical/matrix.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

matrix_t matrix_init(int row, int col)
{
	matrix_t mat;

	double **value = (double**)malloc(row*sizeof(double*));
	for (int i = 0; i < row; i++) {
		value[i] = (double*)malloc(col*sizeof(double));
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

void matrix_copy(const matrix_t *mat, matrix_t *nmat)
{
	if (mat->row != nmat->row || mat->col != nmat->col) {
		printf("error_matrix_copy: dimension mismatch\n");
		return;
	}

	for (int i = 0; i < mat->row; i++) {
		for (int j = 0; j < mat->col; j++) {
			nmat->value[i][j] = mat->value[i][j];
		}
	}
}

void matrix_inverse(const matrix_t *mat, matrix_t *imat)
{
	if (mat->row != mat->col || imat->row != imat->col || mat->row != imat->row) {
		printf("error_matrix_inverse: dimension mismatch\n");
		return;
	}

	int n = mat->row;
    if (n <= 0) {
        printf("error_matrix_inverse: invalid size\n");
        return;
    }

 	/* construct [A | I] */
    double **aug = (double **)malloc((n) * sizeof(double *));
    for (int i = 0; i < n; i++) {
        aug[i] = (double *)malloc((2 * n) * sizeof(double));
    }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            aug[i][j] = mat->value[i][j];           /* matrix A */
            aug[i][j + n] = (i == j) ? 1.0 : 0.0;   /* matrix I */
        }
    }

    /* Gauss-Jordan elimination */
	double EPS = 1e-10;
    for (int p = 0; p < n; p++) {
		/* find the row with the largest pivot element */
		int max_row = p;
        for (int i = p + 1; i < n; i++) {
            if (fabs(aug[i][p]) > fabs(aug[max_row][p])) {
                max_row = i;
            }
        }

		/* strange matrix, return */
        if (fabs(aug[max_row][p]) < EPS) {
            printf("error_matrix_inverse: matrix is singular (or nearly singular)\n");
            for (int i = 0; i < n; i++) free(aug[i]);
            free(aug);
            return;
        }

		/* swap rows */
        if (max_row != p) {
            double *temp = aug[p];
            aug[p] = aug[max_row];
            aug[max_row] = temp;
        }

		/* pivot normalization */
        double pivot = aug[p][p];
        for (int j = 0; j < 2 * n; j++) {
            aug[p][j] /= pivot;
        }

		/* Gauss-Jordan elimination */
        for (int i = 0; i < n; i++) {
            if (i == p) continue;
            double factor = aug[i][p];
            for (int j = 0; j < 2 * n; j++) {
                aug[i][j] -= factor * aug[p][j];
            }
        }
    }

    /* extract the right half as the inverse matrix */
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            imat->value[i][j] = aug[i][j + n];
        }
    }

    for (int i = 0; i < n; i++) {
        free(aug[i]);
    }
    free(aug);
}

void matrix_transpose(const matrix_t *mat, matrix_t *tmat)
{
	if (mat->row != tmat->row || mat->col != tmat->col) {
		printf("error_matrix_transpose: dimension mismatch\n");
		return;
	}

	for (int i = 0; i < mat->row; i++) {
		for (int j = 0; j < mat->col; j++) {
			tmat->value[j][i] = mat->value[i][j];
		}
	}
}

void matrix_sum(const matrix_t *matA, const matrix_t *matB, matrix_t *matC)
{
	if (matA->row != matB->row || matA->col != matB->col
		|| matA->row != matC->row || matA->col != matC->col) {
		printf("error_matrix_sum: dimension mismatch\n");
		return;
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
		return;
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

void matrix_productT(const matrix_t *matA, const matrix_t *matB, matrix_t *matC)
{
	if (matA->col != matB->row || matB->col != matA->col
		|| matA->row != matC->row || matA->row != matC->col) {
		printf("error_matrix_productT: dimension mismatch\n");
		return;
	}

	for (int i = 0; i < matA->row; i++) {
		for (int j = 0; j < matC->col; j++) {
			matC->value[i][j] = 0.0;
			for (int k = 0; k < matA->col; k++) {
				for (int l = 0; l < matB->col; l++) {
					matC->value[i][j] += matA->value[i][k] * matB->value[k][l] * matA->value[j][l];
				}
			}
		}
	}
}

void matrix_product3(const matrix_t *matA, const matrix_t *matB, const matrix_t *matC, matrix_t *matD)
{
	if (matA->col != matB->row || matB->col != matC->row
		|| matA->row != matD->row || matC->col != matD->col) {
		printf("error_matrix_product3: dimension mismatch\n");
		return;
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

double matrix_trace(const matrix_t *mat)
{
	if (mat->row != mat->col) {
		printf("error_matrix_trace: not a square matrix\n");
		return 0.0;
	}

	double trace = 0.0;
	for (int i = 0; i < mat->row; i++) {
		trace += mat->value[i][i];
	}
	return trace;
}

double matrix_norm(const matrix_t *mat)
{
	double norm = 0.0;
	for (int i = 0; i < mat->row; i++) {
		for (int j = 0; j < mat->col; j++) {
			norm += mat->value[i][j] * mat->value[i][j];
		}
	}
	return sqrt(norm);
}

void matrix_multi(double factor, matrix_t *matA)
{
	for (int i = 0; i < matA->row; i++) {
		for (int j = 0; j < matA->col; j++) {
			matA->value[i][j] = matA->value[i][j] * factor;
		}
	}
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

array_t array_init(int len)
{
	array_t arr;

	double *value = (double*)malloc(len*sizeof(double));
	arr.value = value;
	arr.len = len;

	return arr;
}

double array_norm(const array_t *arr)
{
	double norm = 0.0;
	for (int i = 0; i < arr->len; i++) {
		norm += arr->value[i] * arr->value[i];
	}
	return sqrt(norm);
}

void array_print(const array_t *ary)
{
	for (int i = 0; i < ary->len; i++) {
		printf("%10.6f ", ary->value[i]);
	}
	printf("\n");
}

void array_free(array_t *ary)
{
	free(ary->value);
}