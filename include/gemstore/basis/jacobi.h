/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef GEMSTORE_BASIS_JACOBI
#define GEMSTORE_BASIS_JACOBI

/* Jacobi channel c: 1 = pair(1,2) spectator 3,
 *                   2 = pair(3,1) spectator 2,
 *                   3 = pair(2,3) spectator 1.
 *
 * Unweighted coordinates:
 *   ρ_c = r_i - r_j
 *   λ_c = (m_i r_i + m_j r_j)/(m_i+m_j) - r_k
 *
 * Linear map of a basis channel `from` onto an operator channel `to`:
 *   ρ_from = α ρ_to + β λ_to
 *   λ_from = γ ρ_to + δ λ_to
 */

void jacobi_pair_mass(double m1, double m2, double m3, int c,
    double *mi, double *mj, double *mk);

void jacobi_r_map(double m1, double m2, double m3, int from, int to,
    double *alpha, double *beta, double *gamma, double *delta);

double jacobi_mu_rho(double mi, double mj);

double jacobi_mu_lam(double mi, double mj, double mk);

/* Angle of the mass-weighted rotation that takes channel `from` into `to`. */
double jacobi_angle(double m1, double m2, double m3, int from, int to);

/* Raynal–Revai coefficient ⟨lρ' lλ' L | lρ lλ L⟩_φ. */
double raynal_revai(int lrho, int llam, int L, int lrhop, int llamp, double phi);

/* Overlap of two s-wave Gaussians exp(−νρ ρ² − νλ λ²) in different Jacobi frames.
 * Angular factors for ℓ≠0 are not included; use solidharm_central_me with V=1. */
double jacobi_gaussian_overlap(double m1, double m2, double m3,
    int ca, double nurho_a, double nulam_a,
    int cb, double nurho_b, double nulam_b);

/* Complete-the-square widths after mapping both Gaussians onto operator channel `pair`.
 * b11 is the remaining relative-coordinate Gaussian; aRR is the spectator width. */
int jacobi_gaussian_reduce(
    double m1, double m2, double m3,
    int from_a, double nurho_a, double nulam_a,
    int from_b, double nurho_b, double nulam_b,
    int pair,
    double *b11, double *aRR);

/* Same quadratic reduction, plus maps after the spectator shift R' = R + κ r
 * with κ = a_{rR}/(2 a_{RR}):
 *   ρ_a = al_a r + be_a R',  λ_a = ga_a r + de_a R'   (and likewise for b).
 * Exponentials then factor as exp(−b11 r² − aRR R'²). */
typedef struct {
    double b11;
    double aRR;
    double al_a, be_a, ga_a, de_a;
    double al_b, be_b, ga_b, de_b;
} jacobi_shift_t;

int jacobi_gaussian_shift(
    double m1, double m2, double m3,
    int from_a, double nurho_a, double nulam_a,
    int from_b, double nurho_b, double nulam_b,
    int pair,
    jacobi_shift_t *out);

#endif
