/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef GEMSTORE_WRITE
#define GEMSTORE_WRITE

#include <gemstore/input.h>
#include <gemstore/math/matrix.h>
#include <gemstore/model/bagmodel.h>
#include <gemstore/model/obemodel.h>
#include <gemstore/param/argset.h>

/* Write {len} of meson spectra to a file, including mass, RMS radius, and eigenvectors */
int write_meson_spectra(const argsInput_t *input, const array_t *mass, const array_t *radius, const matrix_t *vector, int len);

/* Write {len} of baryon spectra to a file, including mass, pair RMS, and eigenvectors */
int write_baryon_spectra(const argsInput_t *input, const array_t *mass, const matrix_t *vector,
    const array_t *rms12, const array_t *rms13, const array_t *rms23, int len);

/* Write meson wavefunctions to a file */
int write_meson_wfn(const argsInput_t *input, const matrix_t *vector);

/* Write GI potential to a file */
int write_meson_gi_pot(const argsInput_t *input, const argsGIModel_t *args_model, argsGIModelDy_t *args_dynmc);

/* Write non-relativistic potential to a file */
int write_meson_nr_pot(const argsInput_t *input, const argsNRModel_t *args_model, argsNRModelDy_t *args_dynmc);

/* Write baryon pair potentials V12, V13, V23 (central GI, colour −2/3) */
int write_baryon_gi_pot(const argsInput_t *input, const argsGIModel_t *args_model, argsGIModelDy_t *args_dynmc);

/* Write OBE radial pieces versus r in fm. Last column is V = Vc + Vss + Vt. */
int write_obe_pot(const char *path, const argsOBEStatModel_t *model);

/* Write bag spectra: mass, R, x, and color-spin eigenvector of each state. */
int write_bag_spectra(const argsInput_t *input, const bag_mit_t *state, int nstate, const char *order);

#endif