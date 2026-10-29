/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef GEMSTORE_MODEL_CMOL2
#define GEMSTORE_MODEL_CMOL2

#include <gemstore/input.h>
#include <gemstore/math/matrix.h>
#include <gemstore/param/argset.h>

#define MOL2_CHMAX 6

/* S-wave plus D-wave (L and L+2) with every total spin that couples to J.
 * Returns the number of channels. L[] and S[] may be NULL. */
int mol2_build_channels(const argsInput_t *input, int *L, double *S, int maxch);
int mol2_channel_count(const argsInput_t *input);

/* Two-body molecular GEM, MODEL_OBESTAT.
 * Block basis: nch partial waves times nmax Gaussians. Tensor S12 mixes
 * the L and L+2 blocks (Liu ²S_{1/2}–⁴D_{1/2} or ⁴S–²D–⁴D at J=3/2).
 * H = T + Vc + Vss + Vt, T = m1+m2+p²/2μ, solved as H c = E N c.
 * e_out / v_out must hold nch*nmax. Eigenvalues are the molecule mass. */
void mol2_obestat_GEM(const argsInput_t *args_input, const argsOBEModel_t *args_model, argsOBEModelDy_t *args_dynmc, 
    array_t *e_out, matrix_t *v_out, matrix_t *n_out, array_t *rms_out);

#endif
