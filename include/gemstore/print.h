/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef GEMSTORE_PRINT
#define GEMSTORE_PRINT

#include <gemstore/input.h>
#include <gemstore/math/matrix.h>
#include <gemstore/model/bagmodel.h>

/* Print the GEMSTORE logo */
void print_logo();

/* Print the help message */
void print_help();

/* Print copyright and license information */
void print_copyright();

/* Print structured input parameters in professional form */
void print_input_parameters(const argsInput_t *input);

/* Print meson spectra. norm_check is c^T N c against overlap. */
void print_meson_spectra(const array_t *eigenvalue, const matrix_t *eigenvector, const matrix_t *overlap, 
    const array_t *rms12, int len);

/* Two-body molecule spectra. threshold is m1+m2; Ebind = Mass - threshold. */
void print_mol2_spectra(const array_t *eigenvalue, const matrix_t *eigenvector, const matrix_t *overlap,
    const array_t *rms12, double threshold, int len);

/* Print baryon spectra with overlap-metric normalization checks and pair RMS */
void print_baryon_spectra(const array_t *eigenvalue, const matrix_t *eigenvector, const matrix_t *overlap,
    const array_t *rms12, const array_t *rms13, const array_t *rms23, int len);

/* Print one variational bag result per color-spin state. order is NULL without CMI. */
void print_bag_spectra(const bag_mit_t *state, int nstate, const char *order);

#endif
