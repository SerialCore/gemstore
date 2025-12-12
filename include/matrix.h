/*
 * Copyright (C) 2025 Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef MATRIX_H
#define MATRIX_H

typedef struct matrix {
	int n;      	/* count of rows */
	int m;      	/* count of cols */
	double **p; 	/* element values */
} matrix_t;

matrix_t *matrix_init(matrix_t *matrix, int n, int m);

void matrix_push(matrix_t *matrix);

void matrix_print(matrix_t *matrix);

void matrix_free(matrix_t *matrix);

void array_print(double *p, int n);

void array_free(double *p);

#endif