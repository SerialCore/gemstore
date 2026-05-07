/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef GEMSTORE_PRINT
#define GEMSTORE_PRINT

#include <gemstore/math/matrix.h>
#include <gemstore/param/argset.h>

/* Print the GEMSTORE logo */
void print_logo();

/* Print the help message */
void print_help();

/* Print copyright and license information */
void print_copyright();

/* Print structured input parameters in professional form */
void print_input_parameters(const argsInput_t *input);

/* Print debug results for computed spectra with state analysis */
void print_meson_spectra(const array_t *eigenvalue, const array_t *rmsradius, const matrix_t *eigenvector, int nmax);

/* Write {len} of meson spectra to a file, including mass, RMS radius, and eigenvectors */
int write_meson_spectra(const argsInput_t *input, const array_t *mass, const array_t *radius, const matrix_t *vector, int len);

/* Write meson wavefunctions to a file */
int write_meson_wfn(const argsInput_t *input, const matrix_t *vector);

/* Write GI potential to a file */
int write_potential_GI(const argsInput_t *input, const argsGIModel_t *args_model, argsGIModelDy_t *args_dynmc);

#endif
