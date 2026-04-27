/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gemstore/math/cmatrix.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

cmatrix_t cmatrix_init(int row, int col)
{
	cmatrix_t mat;

	complex **value = (complex**)malloc(row*sizeof(complex*));
	for (int i = 0; i < row; i++) {
		value[i] = (complex*)malloc(col*sizeof(complex));
	}
	mat.value = value;
	mat.row = row;
	mat.col = col;

	return mat;
}

void cmatrix_inverse(const cmatrix_t *mat, cmatrix_t *imat)
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
    complex **aug = (complex **)malloc((n) * sizeof(complex *));
    for (int i = 0; i < n; i++) {
        aug[i] = (complex *)malloc((2 * n) * sizeof(complex));
    }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            aug[i][j] = mat->value[i][j];           /* matrix A */
            aug[i][j + n] = (i == j) ? 1.0 + 0.0*I : 0.0 + 0.0*I;   /* matrix I */
        }
    }

    /* Gauss-Jordan elimination */
	double EPS = 1e-10;
    for (int p = 0; p < n; p++) {
		/* find the row with the largest pivot element */
		int max_row = p;
        for (int i = p + 1; i < n; i++) {
            if (cabs(aug[i][p]) > cabs(aug[max_row][p])) {
                max_row = i;
            }
        }

		/* strange matrix, return */
        if (cabs(aug[max_row][p]) < EPS) {
            printf("error_matrix_inverse: matrix is singular (or nearly singular)\n");
            for (int i = 0; i < n; i++) free(aug[i]);
            free(aug);
            return;
        }

		/* swap rows */
        if (max_row != p) {
            complex *temp = aug[p];
            aug[p] = aug[max_row];
            aug[max_row] = temp;
        }

		/* pivot normalization */
        complex pivot = aug[p][p];
        for (int j = 0; j < 2 * n; j++) {
            aug[p][j] /= pivot;
        }

		/* Gauss-Jordan elimination */
        for (int i = 0; i < n; i++) {
            if (i == p) continue;
            complex factor = aug[i][p];
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

void cmatrix_inverse_lowertri(const cmatrix_t *mat, cmatrix_t *imat)
{
	if (mat->row != mat->col || imat->row != imat->col || mat->row != imat->row) {
		printf("error_matrix_inverse_lowertri: dimension mismatch\n");
		return;
	}

    int n = mat->row;
    complex **Linv = (complex**)malloc(n * sizeof(complex*));
    for (int i = 0; i < n; i++) Linv[i] = (complex*)malloc(n * sizeof(complex));

    for (int i = 0; i < n; i++) {
        /* initialize */
        for (int j = 0; j < n; j++) Linv[i][j] = 0.0 + 0.0*I;

		/* diagonal */
        Linv[i][i] = 1.0 / mat->value[i][i];

		/* lower triangle */
        for (int j = i - 1; j >= 0; j--) {
            complex sum = 0.0 + 0.0*I;
            for (int k = j; k < i; k++)
                sum += mat->value[i][k] * Linv[k][j];
            Linv[i][j] = -sum / mat->value[i][i];
        }

		/* set upper triangle to zero */
        for (int j = i + 1; j < n; j++)
            Linv[i][j] = 0.0 + 0.0*I;
    }

	/* write back to imat */
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            imat->value[i][j] = Linv[i][j];

    for (int i = 0; i < n; i++) free(Linv[i]);
    free(Linv);
}

void cmatrix_cholesky_decomp(const cmatrix_t *matS, cmatrix_t *matL)
{
	if (matS->row != matS->col || matL->row != matL->col || matS->row != matL->row) {
        printf("error_matrix_cholesky_decomp: dimension mismatch\n");
        return;
    }

    int n = matS->row;
    complex **L = (complex**)malloc(n * sizeof(complex*));
    for (int i = 0; i < n; i++) L[i] = (complex*)malloc(n * sizeof(complex));

    const double eps = 1e-12;

	/* Hermitianlize */
	for (int i = 0; i < n; i++) {
    	for (int j = 0; j < i; j++) {
        	complex avg = 0.5 * (matS->value[i][j] + conj(matS->value[j][i]));
        	matS->value[i][j] = avg;
        	matS->value[j][i] = conj(avg);
    	}
    	matS->value[i][i] = creal(matS->value[i][i]) + 0.0*I;
	}

    for (int i = 0; i < n; i++) {
        /* lower triangle */
        for (int j = 0; j < i; j++) {
            complex sum = 0.0 + 0.0*I;
            for (int k = 0; k < j; k++)
                sum += L[i][k] * conj(L[j][k]);
            L[i][j] = (matS->value[i][j] - sum) / L[j][j];
        }

        /* diagonal */
        complex sum = 0.0 + 0.0*I;
        for (int k = 0; k < i; k++)
            sum += L[i][k] * conj(L[i][k]);

        complex diag_val = matS->value[i][i] - sum;

		/* define tolerance  */
        double real_part = creal(diag_val);
        double imag_part = cimag(diag_val);
        if (real_part <= 0.0 || fabs(imag_part) > eps) {
            fprintf(stderr, "Cholesky failed at i=%d: real=%g, imag=%g\n", 
                    i, real_part, imag_part);
            exit(1);
        }

		/* only use real part of diagonal elements */
        L[i][i] = csqrt(real_part);

        /* upper triangle to zero */
        for (int j = i + 1; j < n; j++)
            L[i][j] = 0.0 + 0.0*I;
    }

    /* write back */
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            matL->value[i][j] = L[i][j];

    for (int i = 0; i < n; i++) free(L[i]);
    free(L);
}

void cmatrix_transpose(const cmatrix_t *mat, cmatrix_t *tmat)
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

void cmatrix_sum(const cmatrix_t *matA, const cmatrix_t *matB, cmatrix_t *matC)
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

void cmatrix_product(const cmatrix_t *matA, const cmatrix_t *matB, cmatrix_t *matC)
{
	if (matA->col != matB->row || matA->row != matC->row || matB->col != matC->col) {
		printf("error_matrix_product: dimension mismatch\n");
		return;
	}

	for (int i = 0; i < matA->row; i++) {
		for (int j = 0; j < matB->col; j++) {
			matC->value[i][j] = 0.0 + 0.0*I;
			for (int k = 0; k < matA->col; k++) {
				matC->value[i][j] += matA->value[i][k] * matB->value[k][j];
			}
		}
	}
}

void cmatrix_productT(const cmatrix_t *matA, const cmatrix_t *matB, cmatrix_t *matC)
{
	if (matA->col != matB->row || matB->col != matA->col
		|| matA->row != matC->row || matA->row != matC->col) {
		printf("error_matrix_productT: dimension mismatch\n");
		return;
	}

	for (int i = 0; i < matA->row; i++) {
		for (int j = 0; j < matC->col; j++) {
			matC->value[i][j] = 0.0 + 0.0*I;
			for (int k = 0; k < matA->col; k++) {
				for (int l = 0; l < matB->col; l++) {
					matC->value[i][j] += matA->value[i][k] * matB->value[k][l] * conj(matA->value[j][l]);
				}
			}
		}
	}
}

void cmatrix_print(const cmatrix_t *mat)
{
	/* Print a complex matrix with row/column indices and formatted complex elements
	 * Format: (real+imag*i) with 6 decimal places each
	 * Each complex number occupies approximately 24 characters per element
	 */
	complex **value = mat->value;
	int row = mat->row, col = mat->col;

	/* Print column headers */
	printf("     ");
	for (int j = 0; j < col; j++) {
		printf("        %3d        ", j + 1);
	}
	printf("\n");

	/* Print rows with row index and complex elements */
	for (int i = 0; i < row; i++) {
		printf("%3d: ", i + 1);
		for (int j = 0; j < col; j++) {
			printf("(%8.4f%+8.4fi) ", creal(value[i][j]), cimag(value[i][j]));
		}
		printf("\n");
	}
	printf("\n");
}

void cmatrix_free(cmatrix_t *mat)
{
	complex **value = mat->value;
	int row = mat->row;

	for (int i = 0; i < row; i++) {
		free(value[i]);
	}
	free(value);
}

carray_t carray_init(int len)
{
	carray_t arr;

	complex *value = (complex*)malloc(len*sizeof(complex));
	arr.value = value;
	arr.len = len;

	return arr;
}

void carray_print(const carray_t *ary)
{
	/* Print a complex array with formatted complex elements
	 * Format: (real+imag*i) with 6 decimal places each
	 * Each complex number occupies approximately 24 characters per element
	 */
	for (int i = 0; i < ary->len; i++) {
		printf("(%8.4f%+8.4fi) ", creal(ary->value[i]), cimag(ary->value[i]));
	}
	printf("\n");
}

void carray_free(carray_t *ary)
{
	free(ary->value);
}