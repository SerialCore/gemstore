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
#include <gemstore/model/qmcmodel.h>
#include <gemstore/param/argset.h>

/* Write {len} of meson spectra to a file, including mass, RMS radius, and eigenvectors */
int write_meson_spectra(const argsInput_t *input, const array_t *mass, const array_t *radius, const matrix_t *vector, int len);

/* Write {len} of baryon spectra to a file, including mass, pair RMS, and eigenvectors */
int write_baryon_spectra(const argsInput_t *input, const array_t *mass, const matrix_t *vector,
    const array_t *rms12, const array_t *rms13, const array_t *rms23, int len);

/* Write three-body molecule spectra: mass, Ebind, pair RMS, eigenvectors. */
int write_mol3_spectra(const argsInput_t *input, const array_t *mass, const matrix_t *vector,
    const array_t *rms12, const array_t *rms13, const array_t *rms23, double threshold, int len);

/* Write bag spectra. matter is NULL for a vacuum bag and set for QMC. */
int write_bag_spectra(const argsInput_t *input, const bag_mit_t *state, int nstate, const char *order, const qmc_matter_t *matter);

/* Write meson wavefunctions to a file. overlap is the solver N. */
int write_meson_wfn(const argsInput_t *input, const matrix_t *vector, const matrix_t *overlap);

/* Write two-body molecule wavefunctions: one R_a(r) column per partial wave. overlap is the solver N. */
int write_mol2_wfn(const argsInput_t *input, const matrix_t *vector, const matrix_t *overlap);

/* Write GI potential to a file */
int write_meson_gi_pot(const argsInput_t *input, const argsGIModel_t *args_model, argsGIModelDy_t *args_dynmc);

/* Write non-relativistic potential to a file */
int write_meson_nr_pot(const argsInput_t *input, const argsNRModel_t *args_model, argsNRModelDy_t *args_dynmc);

/* Write baryon pair potentials V12, V13, V23 (central GI, colour −2/3) */
int write_baryon_gi_pot(const argsInput_t *input, const argsGIModel_t *args_model, argsGIModelDy_t *args_dynmc);

/* Write OBE radial pieces versus r in fm. Last column is V = Vc + Vss + Vt. */
int write_mol2_obe_pot(const argsInput_t *input, const argsOBEModel_t *args_model);

/* Write Yukawa Vc, Vss versus r in fm. Vt is zero. V = Vc + Vss. */
int write_mol2_ykw_pot(const argsInput_t *input, const argsYKWModel_t *args_model);

/* Write the EFT contact. Vss and Vt are zero. V = -Ca exp(-(r/Rc)^2). */
int write_mol2_eft_pot(const argsInput_t *input, const argsEFTModel_t *args_model);

/* Write three pair OBE curves into {project}.pot.dat. */
int write_mol3_pot(const argsInput_t *input);

#endif