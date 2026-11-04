/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef GEMSTORE_MODEL_CMOL3
#define GEMSTORE_MODEL_CMOL3

#include <gemstore/input.h>
#include <gemstore/math/matrix.h>

/* Three-body molecular GEM, SYSTEM_MOLECUL3 + MODEL_HYBRID.
 * Each pair is OBESTAT or YUKAWA (or another radial kernel registered
 * in mol3_load_pairs). Jacobi charts c=1,2,3 with physical masses and
 * hadron spins 0, 1/2, 1.
 * H = Σ m_i + Σ p_i²/(2m_i) + Σ_pairs (Vc + Vss s_i·s_j + Vt S12), Vls=0.
 * system.I is required: pair tij is a basis QN recoupled to total I;
 * pair iso_pi/iso_rho are ignored. t1,t2,t3 default to 1/2.
 * sym12=0: three independent charts. sym12=-1 on three identical Dirac
 * legs: each Pauli-allowed component is the rearrangement sum of charts
 * c=1,2,3 at the same (l, sij, tij), with (−1)^{lρ+sij+tij}=−1.
 * Overcomplete N is pruned, then H' u = E u. Eigenvalues are the mass. */
void mol3_hybrid_GEM(const argsInput_t *args_input,
    array_t *e_out, matrix_t *v_out, matrix_t *n_out,
    array_t *rms12, array_t *rms13, array_t *rms23);

#endif
