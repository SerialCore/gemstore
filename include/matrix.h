/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef MATRIX_H
#define MATRIX_H

typedef struct matrix {
	int row;      		/* count of rows */
	int col;      		/* count of cols */
	double *value; 		/* element values */
} matrix_t;

matrix_t *matrix_init(int row, int col);

void matrix_push(matrix_t *matrix);

void matrix_print(matrix_t *matrix);

void matrix_free(matrix_t *matrix);

#endif