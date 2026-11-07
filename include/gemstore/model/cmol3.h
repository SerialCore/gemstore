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
 * Each pair is OBESTAT, YUKAWA, or EFTCONT (or another radial kernel
 * registered in mol3_load_pairs). Jacobi charts c=1,2,3 with physical masses and
 * hadron spins 0, 1/2, 1.
 * H = Σ m_i + Σ p_i²/(2m_i) + Σ_pairs (Vc + Vss s_i·s_j + Vt S12), Vls=0.
 * system.I is required: pair tij is a basis QN recoupled to total I;
 * pair iso_pi/iso_rho are ignored. t1,t2,t3 default to 1/2.
 * π/ρ use c_i c_j Ti·Tj with c=2 on a doublet (τ) and c=1 on integer isospin.
 * Identical legs: same type and mass. Dirac is a fermion, VECTOR and
 * PSEUDOSCALAR are bosons. A Jacobi ρ-pair that is identical is kept
 * only when (−1)^{lρ+(2s−sij)+(2t−tij)} equals −1 (fermion) or +1
 * (boson). Three identical legs: one S3 component, rearrangement sum of
 * charts c=1,2,3 at the same (l, sij, tij). Two or three distinguishable
 * legs: independent charts; Bose/Pauli only on the identical ρ-pair.
 * sym12 is not used for this selection (kept in JSON for baryons).
 * Overcomplete N is pruned, then H' u = E u. Eigenvalues are the mass. */
void mol3_hybrid_GEM(const argsInput_t *args_input,
    array_t *e_out, matrix_t *v_out, matrix_t *n_out,
    array_t *rms12, array_t *rms13, array_t *rms23);

#endif
