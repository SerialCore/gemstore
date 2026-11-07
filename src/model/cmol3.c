/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 *
 * Three-body molecular GEM. SCDK kinematics as in cbaryon. Each pair
 * supplies its own radial kernel (OBESTAT, YUKAWA, or EFTCONT) through
 * integral_exp_rn. Kinetics are non-relativistic.
 *
 * H = Σ m_i + Σ p_i²/(2m_i) + Σ_pairs (Vc + Vss s_i·s_j + Vt S12).
 * Pair tij is a basis QN recoupled to required system.I.
 */

#include <gemstore/model/cmol3.h>
#include <gemstore/model/obemodel.h>
#include <gemstore/model/ykwmodel.h>
#include <gemstore/model/eftmodel.h>
#include <gemstore/model/giradial.h>

#include <gemstore/basis/basis.h>
#include <gemstore/basis/channel.h>
#include <gemstore/basis/orbit.h>

#include <math.h>

#include <gemstore/math/kinemtc.h>
#include <gemstore/math/spinang.h>
#include <gemstore/math/ckdk.h>
#include <gemstore/math/soc.h>
#include <gemstore/math/eigen.h>
#include <gemstore/math/matrix.h>
#include <gemstore/math/integral.h>

#include <gemstore/param/argset.h>
#include <gemstore/thread.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define MOL3_NOP 9
#define MOL3_ISO_K 4

/* One isospin slot of one pair. ctx is owned by mol3_vargs and stays put. */
typedef struct mol3_kernel {
    potential_t vc;
    potential_t vss;
    potential_t vt;
    void *ctx;
} mol3_kernel_t;

typedef struct mol3_vargs {
    int n_tp[3];
    double tp_val[3][MOL3_ISO_K];
    int has_vt[3];
    mol3_kernel_t ker[3][MOL3_ISO_K];
    argsOBEModel_t obe_m[3][MOL3_ISO_K];
    argsOBEModelDy_t obe_d[3][MOL3_ISO_K];
    obe_pot_ctx_t obe_c[3][MOL3_ISO_K];
    argsYKWModel_t ykw_m[3];
    argsYKWModelDy_t ykw_d[3];
    ykw_pot_ctx_t ykw_c[3];
    argsEFTModel_t eft_m[3];
    argsEFTModelDy_t eft_d[3];
    eft_pot_ctx_t eft_c[3];
} mol3_vargs_t;

typedef struct mol3_iv {
    void *ctx;
    potential_t pot;
} mol3_iv_t;

typedef struct mol3_job {
    bases_3b qnlist_spfy;
    bases_3b qnlist_full;
    matrix_t **mlsj;
    ckdk_scdk ****scdk[MOL3_NOP];
    mol3_vargs_t varg;
    double **iso_ov;
    double ***iso_U;
    matrix_t Nfi;
    matrix_t T[3];
    matrix_t Vc[3];
    matrix_t Vss[3];
    matrix_t Vt[3];
    matrix_t rmsr[3];
} mol3_job_t;

static int mol3_triangle(double a, double b, double c)
{
    return (a + b + 1e-10 >= c) && (a + c + 1e-10 >= b) && (b + c + 1e-10 >= a);
}

/* ---------- pair OBE kernels ---------- */

/* s_i · s_j with hadron spins on operator pair c. vcont hard-codes s=1/2. */
static void vss_hadron(ckdk_scdk *scdk, double coef, double coei, double msai, double msaj, double msak,
    double msbi, double msbj, double msbk, int l1, int l2, int l3, int l4, int ml1, int ml2, int ml3, int ml4,
    basis_3b qa, basis_3b qb, pthread_mutex_t *mutex, int lock, int c)
{
    double si, sj, sk, coe = 0.0;

    (void)qb;
    channel_3b_ijk(qa.s1, qa.s2, qa.s3, &si, &sj, &sk, c);
    (void)sk;
    if (si < 1e-12 || sj < 1e-12) {
        return;
    }
    if (msak == msbk) {
        for (int M = -1; M <= 1; M++) {
            coe -= sqrt(3.0) * sqrt(si * (si + 1.0)) * sqrt(sj * (sj + 1.0))
                * clebsch_gordan(1, M, 1, -M, 0, 0)
                * clebsch_gordan(si, msbi, 1, M, si, msai)
                * clebsch_gordan(sj, msbj, 1, -M, sj, msaj);
        }
    }
    if (coe != 0.0) {
        ckdk_scdk_calc(scdk, 4.0 * M_PI * coe * coef * coei,
            l1, l2, l3, l4, 0, 0, ml1, ml2, ml3, ml4, 0, 0, mutex, lock);
    }
}

static void vtens_hadron(ckdk_scdk *scdk, double coef, double coei, double msai, double msaj, double msak,
    double msbi, double msbj, double msbk, int l1, int l2, int l3, int l4, int ml1, int ml2, int ml3, int ml4,
    basis_3b qa, basis_3b qb, pthread_mutex_t *mutex, int lock, int c)
{
    double si, sj, sk;

    (void)qb;
    channel_3b_ijk(qa.s1, qa.s2, qa.s3, &si, &sj, &sk, c);
    (void)sk;
    if (si < 1e-12 || sj < 1e-12 || msak != msbk) {
        return;
    }
    for (int M = -2; M <= 2; M++) {
        for (int mu = -1; mu <= 1; mu++) {
            double coe = clebsch_gordan(2, -M, 2, M, 0, 0)
                * clebsch_gordan(1, mu, 1, -M - mu, 2, -M)
                * sqrt(si * (si + 1.0)) * sqrt(sj * (sj + 1.0))
                * clebsch_gordan(si, msbi, 1, mu, si, msai)
                * clebsch_gordan(sj, msbj, 1, -M - mu, sj, msaj);
            if (coe != 0.0) {
                ckdk_scdk_calc(scdk, sqrt(4.0 * M_PI) * coe * coef * coei,
                    l1, l2, l3, l4, 2, 0, ml1, ml2, ml3, ml4, M, 0, mutex, lock);
            }
        }
    }
}

/* ops 0–2 central, 3–5 spin–spin, 6–8 tensor; pair = op % 3 + 1. */
static const vtype mol3_vfn[MOL3_NOP] = {
    vcent, vcent, vcent,
    vss_hadron, vss_hadron, vss_hadron,
    vtens_hadron, vtens_hadron, vtens_hadron
};
static const int mol3_fg[MOL3_NOP] = {
    1, 1, 1, 1, 1, 1, 2, 2, 2
};
static const int mol3_pair[MOL3_NOP] = {
    1, 2, 3, 1, 2, 3, 1, 2, 3
};

/* ⟨p²⟩ / (2 m_c); tpi_cent puts x in momentum space, inteRMS raises the power by 2. */
static double mol3_inte_tnr(double b11, int n, basis_3b qf, basis_3b qi, void *ctx, int c)
{
    double mi = (c == 1) ? qf.m1 : (c == 2) ? qf.m2 : (c == 3) ? qf.m3 : 0.0;

    (void)ctx;
    if (mi <= 0.0) {
        return 0.0;
    }
    return inteRMS(b11, n, qf, qi, ctx, c) / (2.0 * mi);
}

static double mol3_inte_v(double b11, int n, basis_3b qf, basis_3b qi, void *ctx, int c)
{
    mol3_iv_t *iv = (mol3_iv_t *)ctx;

    (void)qf;
    (void)qi;
    if (c < 1 || c > 3) {
        return 0.0;
    }
    return integral_exp_rn(iv->pot, b11, n, iv->ctx);
}

/* ---------- isospin recoupling ---------- */

/* Isovector factor at the vertex: a doublet uses τ = 2T, an integer isospin
 * (Σc is t=1) uses the generator T. Two doublets give τ1·τ2 (−3 or +1).
 * A doublet–triplet pair gives T·τ (−2 for I=1/2, +1 for I=3/2). */
static double mol3_iso_vertex(double t)
{
    return (fabs(t - 0.5) < 1e-8) ? 2.0 : 1.0;
}

static double mol3_tau_dot(double ti, double tj, double tij)
{
    double tt = tij * (tij + 1.0) - ti * (ti + 1.0) - tj * (tj + 1.0);

    return 0.5 * mol3_iso_vertex(ti) * mol3_iso_vertex(tj) * tt;
}

static double mol3_iso_scheme_amp(double t1, double t2, double t3, int c,
    double tij, double T, double MT, double m1, double m2, double m3)
{
    double ti, tj, tk, mi, mj, mk;

    channel_3b_ijk(t1, t2, t3, &ti, &tj, &tk, c);
    channel_3b_ijk(m1, m2, m3, &mi, &mj, &mk, c);
    return clebsch_gordan(ti, mi, tj, mj, tij, mi + mj)
        * clebsch_gordan(tij, mi + mj, tk, mk, T, MT);
}

static double mol3_iso_bracket(double t1, double t2, double t3, double T,
    int c_a, double tij_a, int c_b, double tij_b)
{
    double acc = 0.0;

    if (T < -1e-12) {
        return 0.0;
    }
    for (double m1 = -t1; m1 <= t1 + 1e-12; m1 += 1.0) {
        for (double m2 = -t2; m2 <= t2 + 1e-12; m2 += 1.0) {
            for (double m3 = -t3; m3 <= t3 + 1e-12; m3 += 1.0) {
                if (fabs(m1 + m2 + m3 - T) > 1e-9) {
                    continue;
                }
                acc += mol3_iso_scheme_amp(t1, t2, t3, c_a, tij_a, T, T, m1, m2, m3)
                    * mol3_iso_scheme_amp(t1, t2, t3, c_b, tij_b, T, T, m1, m2, m3);
            }
        }
    }
    return acc;
}

static void mol3_iso_cache_free(mol3_job_t *job, int nsp)
{
    if (job->iso_ov == NULL) {
        return;
    }
    for (int i = 0; i < nsp; i++) {
        free(job->iso_ov[i]);
        for (int p = 0; p < 3; p++) {
            free(job->iso_U[i][p]);
        }
        free(job->iso_U[i]);
    }
    free(job->iso_ov);
    free(job->iso_U);
    job->iso_ov = NULL;
    job->iso_U = NULL;
}

static void mol3_iso_cache_build(mol3_job_t *job)
{
    int nsp = job->qnlist_spfy.len_list;
    basis_3b q0 = job->qnlist_spfy.qnum[0][0];
    double t1 = q0.t1, t2 = q0.t2, t3 = q0.t3, T = q0.T;
    double u00, u01, same, orth, nrm;

    job->iso_ov = (double **)malloc(sizeof(double *) * (size_t)nsp);
    job->iso_U = (double ***)malloc(sizeof(double **) * (size_t)nsp);
    for (int i = 0; i < nsp; i++) {
        basis_3b qi = job->qnlist_spfy.qnum[i][0];
        job->iso_ov[i] = (double *)malloc(sizeof(double) * (size_t)nsp);
        job->iso_U[i] = (double **)malloc(sizeof(double *) * 3);
        for (int p = 0; p < 3; p++) {
            job->iso_U[i][p] = (double *)calloc((size_t)MOL3_ISO_K, sizeof(double));
            for (int k = 0; k < job->varg.n_tp[p]; k++) {
                job->iso_U[i][p][k] = mol3_iso_bracket(t1, t2, t3, T,
                    qi.c, qi.tij, p + 1, job->varg.tp_val[p][k]);
            }
        }
        for (int j = 0; j < nsp; j++) {
            basis_3b qj = job->qnlist_spfy.qnum[j][0];
            job->iso_ov[i][j] = mol3_iso_bracket(t1, t2, t3, T,
                qi.c, qi.tij, qj.c, qj.tij);
        }
    }

    u00 = mol3_iso_bracket(0.5, 0.5, 0.5, 0.5, 1, 0.0, 3, 0.0);
    u01 = mol3_iso_bracket(0.5, 0.5, 0.5, 0.5, 1, 0.0, 3, 1.0);
    same = mol3_iso_bracket(0.5, 0.5, 0.5, 0.5, 1, 0.0, 1, 0.0);
    orth = mol3_iso_bracket(0.5, 0.5, 0.5, 0.5, 1, 0.0, 1, 1.0);
    nrm = u00 * u00 + u01 * u01;
    printf("  isospin recoupling check (t=1/2, I=1/2): ⟨c1,t12=0|c1,0⟩=%.6f  ⟨0|1⟩=%.6f\n",
        same, orth);
    printf("  ⟨c1,t12=0|c3,t23=0⟩=%.6f  ⟨0|1⟩=%.6f  |U|^2=%.6f\n", u00, u01, nrm);
    if (fabs(same - 1.0) > 1e-8 || fabs(orth) > 1e-8 || fabs(nrm - 1.0) > 1e-8) {
        fprintf(stderr, "Error: isospin recoupling failed the t=1/2 unitarity check.\n");
        exit(1);
    }
    fflush(stdout);
}

/* ---------- GEM basis ---------- */

static obe_leg_type_t mol3_leg_type(const argsInput_t *input, int a)
{
    return (a == 1) ? input->type1 : (a == 2) ? input->type2 : input->type3;
}

static double mol3_leg_mass(const argsInput_t *input, int a)
{
    return (a == 1) ? input->m1 : (a == 2) ? input->m2 : input->m3;
}

static int mol3_identical(const argsInput_t *input, int a, int b)
{
    return mol3_leg_type(input, a) == mol3_leg_type(input, b)
        && fabs(mol3_leg_mass(input, a) - mol3_leg_mass(input, b)) < 1e-9;
}

static int mol3_chart_pair_identical(int id12, int id31, int id23, int c)
{
    switch (c) {
        case 1: return id12;
        case 2: return id31;
        case 3: return id23;
        default: return 0;
    }
}

/* (−1)^{ℓ+(2s−S)+(2t−T)} = −1 fermion, +1 boson. */
static int mol3_pair_sym_ok(int ell, double s, double t, double S, double T, int fermion)
{
    int two_s = (int)(2.0 * s + 0.5);
    int two_t = (int)(2.0 * t + 0.5);
    int iS = (int)(S + 0.5);
    int iT = (int)(T + 0.5);
    int expo = ell + (two_s - iS) + (two_t - iT);
    int odd = (expo % 2 != 0);
    return fermion ? odd : !odd;
}

static void mol3_basis_build(const argsInput_t *input, bases_3b *spfy, bases_3b *full)
{
    double m1 = input->m1, m2 = input->m2, m3 = input->m3;
    double s1 = obe_leg_spin(input->type1);
    double s2 = obe_leg_spin(input->type2);
    double s3 = obe_leg_spin(input->type3);
    double J = input->J, totI = input->iso_I;
    double t1 = input->iso_t1, t2 = input->iso_t2, t3 = input->iso_t3;
    int P = input->P, Lmax = input->Lmax;
    int Pint = obe_leg_parity(input->type1) * obe_leg_parity(input->type2)
        * obe_leg_parity(input->type3);
    int id12 = mol3_identical(input, 1, 2);
    int id31 = mol3_identical(input, 3, 1);
    int id23 = mol3_identical(input, 2, 3);
    int s3_ident = id12 && id23;
    int any_ident = id12 || id23 || id31;
    int ident_fermion = (id12 && input->type1 == OBE_LEG_DIRAC)
        || (id23 && input->type2 == OBE_LEG_DIRAC)
        || (id31 && input->type3 == OBE_LEG_DIRAC);

    bases_3b_init(spfy);
    bases_3b_init(full);

    if (P != 1 && P != -1) {
        fprintf(stderr, "Error: MOLECUL3 parity P must be +1 or -1.\n");
        exit(1);
    }
    if (Lmax < 0) {
        fprintf(stderr, "Error: MOLECUL3 Lmax must be non-negative.\n");
        exit(1);
    }
    if (totI < -1e-12) {
        fprintf(stderr, "Error: MOLECUL3 total I must be non-negative.\n");
        exit(1);
    }

    for (int lrho = 0; lrho <= Lmax; lrho++) {
        for (int llam = 0; llam <= Lmax - lrho; llam++) {
            int spatial = ((lrho + llam) % 2 == 0) ? 1 : -1;
            int cmax;
            if (Pint * spatial != P) {
                continue;
            }
            cmax = s3_ident ? 1 : 3;
            for (int L = abs(lrho - llam); L <= lrho + llam; L++) {
                for (int c = 1; c <= cmax; c++) {
                    double si, sj, sk, ti, tj, tk;
                    int pair_fermion;
                    channel_3b_ijk(s1, s2, s3, &si, &sj, &sk, c);
                    channel_3b_ijk(t1, t2, t3, &ti, &tj, &tk, c);
                    if (c == 1) {
                        pair_fermion = (input->type1 == OBE_LEG_DIRAC);
                    } else if (c == 2) {
                        pair_fermion = (input->type3 == OBE_LEG_DIRAC);
                    } else {
                        pair_fermion = (input->type2 == OBE_LEG_DIRAC);
                    }
                    for (double sij = fabs(si - sj); sij <= si + sj + 1e-9; sij += 1.0) {
                        for (double jl = fabs((double)L - sij); jl <= (double)L + sij + 1e-9; jl += 1.0) {
                            if (!mol3_triangle(jl, sk, J)) {
                                continue;
                            }
                            for (double tij = fabs(ti - tj); tij <= ti + tj + 1e-9; tij += 1.0) {
                                if (!mol3_triangle(tij, tk, totI)) {
                                    continue;
                                }
                                if (mol3_chart_pair_identical(id12, id31, id23, c)
                                    && !mol3_pair_sym_ok(lrho, si, ti, sij, tij, pair_fermion)) {
                                    continue;
                                }
                                bases_3b_push(spfy, 1, -1, -1, 1.0, m1, m2, m3, s1, s2, s3,
                                    t1, t2, t3, tij, totI, c, lrho, llam, L, sij, jl, J,
                                    0, 0, 0.0, 0.0);
                                /* Three identical hadrons: add the two cyclic rearrangement
                                 * charts to this component. A cluster can then sit on any
                                 * pair. Pair Fermi/Bose is already on (lρ, sij, tij), and a
                                 * 3-cycle is even, so the chart weight is +1. */
                                if (s3_ident) {
                                    bases_3b_push(spfy, 0, -1, -1, 1.0, m1, m2, m3, s1, s2, s3,
                                        t1, t2, t3, tij, totI, 2, lrho, llam, L, sij, jl, J,
                                        0, 0, 0.0, 0.0);
                                    bases_3b_push(spfy, 0, -1, -1, 1.0, m1, m2, m3, s1, s2, s3,
                                        t1, t2, t3, tij, totI, 3, lrho, llam, L, sij, jl, J,
                                        0, 0, 0.0, 0.0);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    bases_3b_push_full(spfy, full, input->rmin, input->rmax, input->nmax);
    /* Unweighted Jacobi: ρ=r1−r2, λ=(r1+r2)/2−r3. Equal-mass S3-symmetric
     * Gaussians need νλ/νρ = μλ/μρ = 4/3 so that ⟨r12²⟩=⟨r13²⟩=⟨r23²⟩. */
    if (s3_ident && fabs(m1 - m2) < 1e-9 && fabs(m2 - m3) < 1e-9) {
        for (int i = 0; i < full->len_list; i++) {
            for (int j = 0; j < full->len_part[i]; j++) {
                full->qnum[i][j].nulam *= 4.0 / 3.0;
            }
        }
        printf("  equal-mass S3: νλ=(4/3)νρ so nρ=nλ is permutation-symmetric\n");
    }
    if (s3_ident) {
        printf("Molecule3 GEM basis: S3  rearrangement c=1+2+3, %s, P_int=%d, lρ+lλ≤%d\n",
            ident_fermion ? "fermion Pauli" : "boson Bose", Pint, Lmax);
    } else if (any_ident) {
        printf("Molecule3 GEM basis: independent c=1,2,3, identical-pair %s, P_int=%d, lρ+lλ≤%d\n",
            ident_fermion ? "Pauli" : "Bose", Pint, Lmax);
    } else {
        printf("Molecule3 GEM basis: distinguishable  independent c=1,2,3, P_int=%d, lρ+lλ≤%d\n",
            Pint, Lmax);
    }
    printf("  total I=%.1f  t=(%.1f,%.1f,%.1f); pair tij recoupled, pair iso_pi ignored\n",
        totI, t1, t2, t3);
    if (any_ident) {
        printf("  identical ρ-pair (−1)^{lρ+(2s−sij)+(2t−tij)}=%s\n",
            ident_fermion ? "−1" : "+1");
    }
    printf("  angular channels: %d  (n_full=%d)\n", spfy->len_list, full->len_list);
    fflush(stdout);
}

static void mlsj_push(matrix_t *m, double coe, double ms1, double ms2, double ms3, double mrho, double mlam)
{
    int n = m->row;
    m->value = (double **)realloc(m->value, sizeof(double *) * (size_t)(n + 1));
    m->value[n] = (double *)malloc(sizeof(double) * 6);
    m->value[n][0] = coe;
    m->value[n][1] = ms1;
    m->value[n][2] = ms2;
    m->value[n][3] = ms3;
    m->value[n][4] = mrho;
    m->value[n][5] = mlam;
    m->row = n + 1;
    m->col = 6;
}

static void mol3_mlsj_jl(mol3_job_t *job)
{
    bases_3b *qnlist = &job->qnlist_spfy;

    job->mlsj = (matrix_t **)malloc(sizeof(matrix_t *) * (size_t)qnlist->len_list);
    for (int i = 0; i < qnlist->len_list; i++) {
        job->mlsj[i] = (matrix_t *)malloc(sizeof(matrix_t) * (size_t)qnlist->len_part[i]);
        for (int j = 0; j < qnlist->len_part[i]; j++) {
            basis_3b q = qnlist->qnum[i][j];
            double si, sj, sk;
            job->mlsj[i][j].value = (double **)malloc(sizeof(double *) * 0);
            job->mlsj[i][j].row = 0;
            job->mlsj[i][j].col = 6;
            channel_3b_ijk(q.s1, q.s2, q.s3, &si, &sj, &sk, q.c);
            for (double msi = -si; msi <= si + 1e-12; msi += 1.0) {
                for (double msj = -sj; msj <= sj + 1e-12; msj += 1.0) {
                    for (double msk = -sk; msk <= sk + 1e-12; msk += 1.0) {
                        for (int mrho = -q.lrho; mrho <= q.lrho; mrho++) {
                            for (int mlam = -q.llam; mlam <= q.llam; mlam++) {
                                double cgf = q.coe
                                    * clebsch_gordan(si, msi, sj, msj, q.sij, msi + msj)
                                    * clebsch_gordan(q.lrho, mrho, q.llam, mlam, q.L, mrho + mlam)
                                    * clebsch_gordan(q.sij, msi + msj, q.L, mrho + mlam, q.jl,
                                        msi + msj + mrho + mlam)
                                    * clebsch_gordan(q.jl, msi + msj + mrho + mlam, sk, msk, q.J, q.J);
                                if (cgf != 0.0) {
                                    double ms1p, ms2p, ms3p;
                                    channel_3b_123(&ms1p, &ms2p, &ms3p, msi, msj, msk, q.c);
                                    mlsj_push(&job->mlsj[i][j], cgf, ms1p, ms2p, ms3p,
                                        (double)mrho, (double)mlam);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

/* ---------- pair models and matrix fill ---------- */

static void mol3_pair_mass(const argsInput_t *input, int p, double *mi, double *mj)
{
    if (p == 0) {
        *mi = input->m1;
        *mj = input->m2;
    }
    else if (p == 1) {
        *mi = input->m1;
        *mj = input->m3;
    }
    else {
        *mi = input->m2;
        *mj = input->m3;
    }
}

static double mol3_pair_kss(const argsInput_t *input, int p)
{
    obe_leg_type_t t1, t2;

    if (p == 0) {
        t1 = input->type1;
        t2 = input->type2;
    }
    else if (p == 1) {
        t1 = input->type1;
        t2 = input->type3;
    }
    else {
        t1 = input->type2;
        t2 = input->type3;
    }
    return (t1 == OBE_LEG_DIRAC ? 2.0 : 1.0) * (t2 == OBE_LEG_DIRAC ? 2.0 : 1.0);
}

/* OBESTAT: one radial kernel per pair isospin, because π/ρ carry ⟨τ·τ⟩.
 * YUKAWA: spin dependence sits in Vss, so every tij of a pair shares one kernel.
 * EFTCONT: the file is already one projected channel, and V is central,
 * so every tij of that pair shares one kernel. */
static void mol3_load_pairs(const argsInput_t *input, mol3_vargs_t *varg)
{
    for (int p = 0; p < 3; p++) {
        double mi, mj, kss, ti, tj, tk;
        int c = p + 1;
        model_type_t kind = input->mol3_pairs[p].model;

        mol3_pair_mass(input, p, &mi, &mj);
        kss = mol3_pair_kss(input, p);
        varg->has_vt[p] = 0;
        varg->n_tp[p] = 0;

        if (kind == MODEL_YUKAWA) {
            varg->ykw_m[p] = argsYKWModel_from(input, p);
            varg->ykw_d[p].mi = mi;
            varg->ykw_d[p].mj = mj;
            varg->ykw_d[p].OCent = 1.0;
            varg->ykw_d[p].OSdS = kss;
            varg->ykw_c[p].model = &varg->ykw_m[p];
            varg->ykw_c[p].dyn = &varg->ykw_d[p];
            varg->has_vt[p] = 0;
        }
        else if (kind == MODEL_EFTCONT) {
            varg->eft_m[p] = argsEFTModel_from(input, p);
            varg->eft_d[p].mi = mi;
            varg->eft_d[p].mj = mj;
            varg->eft_d[p].OCent = 1.0;
            varg->eft_c[p].model = &varg->eft_m[p];
            varg->eft_c[p].dyn = &varg->eft_d[p];
            varg->has_vt[p] = 0;
        }
        else if (kind != MODEL_OBESTAT) {
            fprintf(stderr, "Error: MOLECUL3 pair %d model %s is not OBESTAT, YUKAWA, or EFTCONT.\n",
                p + 1, model_type_str[kind]);
            exit(1);
        }

        channel_3b_ijk(input->iso_t1, input->iso_t2, input->iso_t3, &ti, &tj, &tk, c);
        for (double tij = fabs(ti - tj); tij <= ti + tj + 1e-9; tij += 1.0) {
            int k;
            double iso;
            if (!mol3_triangle(tij, tk, input->iso_I)) {
                continue;
            }
            if (varg->n_tp[p] >= MOL3_ISO_K) {
                fprintf(stderr, "Error: MOLECUL3 pair isospin slots exceeded.\n");
                exit(1);
            }
            k = varg->n_tp[p];
            iso = mol3_tau_dot(ti, tj, tij);
            varg->tp_val[p][k] = tij;
            if (kind == MODEL_YUKAWA) {
                varg->ker[p][k].vc = YKWVc;
                varg->ker[p][k].vss = YKWVss;
                varg->ker[p][k].vt = NULL;
                varg->ker[p][k].ctx = &varg->ykw_c[p];
            }
            else if (kind == MODEL_EFTCONT) {
                varg->ker[p][k].vc = EFTVc;
                varg->ker[p][k].vss = NULL;
                varg->ker[p][k].vt = NULL;
                varg->ker[p][k].ctx = &varg->eft_c[p];
            }
            else {
                varg->obe_m[p][k] = argsOBEModel_from(input, p);
                varg->obe_m[p][k].iso_pi = iso;
                varg->obe_m[p][k].iso_rho = iso;
                varg->obe_d[p][k].mi = mi;
                varg->obe_d[p][k].mj = mj;
                varg->obe_d[p][k].OCent = 1.0;
                varg->obe_d[p][k].OSdS = kss;
                varg->obe_d[p][k].OLSi = 0.0;
                varg->obe_d[p][k].OLSj = 0.0;
                varg->obe_d[p][k].OTens = kss;
                varg->obe_c[p][k].model = &varg->obe_m[p][k];
                varg->obe_c[p][k].dyn = &varg->obe_d[p][k];
                varg->ker[p][k].vc = OBEVc;
                varg->ker[p][k].vss = OBEVss;
                varg->ker[p][k].vt = OBEVtens;
                varg->ker[p][k].ctx = &varg->obe_c[p][k];
                varg->has_vt[p] = 1;
            }
            varg->n_tp[p]++;
        }
        printf("  pair %d recoupled tij:", p + 1);
        for (int k = 0; k < varg->n_tp[p]; k++) {
            printf(" %.1f(⟨τ·τ⟩=%.3f)", varg->tp_val[p][k],
                mol3_tau_dot(ti, tj, varg->tp_val[p][k]));
        }
        printf("\n");
    }
    fflush(stdout);
}

static void mol3_scdk_table_alloc(int *len_part, int len_list, ckdk_scdk *****tab)
{
    *tab = (ckdk_scdk ****)malloc(sizeof(ckdk_scdk ***) * (size_t)len_list);
    for (int i = 0; i < len_list; i++) {
        (*tab)[i] = (ckdk_scdk ***)malloc(sizeof(ckdk_scdk **) * (size_t)len_part[i]);
        for (int j = 0; j < len_part[i]; j++) {
            (*tab)[i][j] = (ckdk_scdk **)malloc(sizeof(ckdk_scdk *) * (size_t)len_list);
            for (int k = 0; k < len_list; k++) {
                (*tab)[i][j][k] = (ckdk_scdk *)malloc(sizeof(ckdk_scdk) * (size_t)len_part[k]);
            }
        }
    }
}

static void mol3_scdk_table_free(int *len_part, int len_list, ckdk_scdk *****tab)
{
    for (int i = 0; i < len_list; i++) {
        for (int j = 0; j < len_part[i]; j++) {
            for (int k = 0; k < len_list; k++) {
                for (int l = 0; l < len_part[k]; l++) {
                    ckdk_scdk_free(&((*tab)[i][j][k][l]));
                }
                free((*tab)[i][j][k]);
            }
            free((*tab)[i][j]);
        }
        free((*tab)[i]);
    }
    free(*tab);
    *tab = NULL;
}

static void calc_scdk_mt(void *args)
{
    argsThread_t *mtargs = (argsThread_t *)args;
    mol3_job_t *job = (mol3_job_t *)mtargs->args;
    int ith = mtargs->index;
    int nf = ith / MOL3_NOP;
    int op = ith % MOL3_NOP;
    int *len_part = job->qnlist_spfy.len_part;
    int len_list = job->qnlist_spfy.len_list;

    for (int nfp = 0; nfp < len_part[nf]; nfp++) {
        for (int ni = 0; ni < len_list; ni++) {
            for (int nip = 0; nip < len_part[ni]; nip++) {
                ckdk_scdk_vtype(&job->scdk[op][nf][nfp][ni][nip],
                    mol3_vfn[op], mol3_fg[op], job->qnlist_spfy, job->mlsj,
                    nf, nfp, ni, nip, mtargs->mutex, mtargs->lock, mol3_pair[op]);
            }
        }
    }
}

static void getmfi(void *args)
{
    argsThread_t *mtarg = (argsThread_t *)args;
    mol3_job_t *job = (mol3_job_t *)mtarg->args;
    int nf = mtarg->index;

    for (int ni = 0; ni < job->qnlist_full.len_list; ni++) {
        double nfi = 0.0;
        double T[3] = {0}, Vc[3] = {0}, Vss[3] = {0}, Vt[3] = {0}, rmsr[3] = {0};

        /* Isospin is per rearrangement component. Charts of one basis
         * vector are different coupling schemes, so a single factor from
         * part 0 would drop the cross-chart norm and the pair projectors. */
        for (int nfp = 0; nfp < job->qnlist_full.len_part[nf]; nfp++) {
            for (int nip = 0; nip < job->qnlist_full.len_part[ni]; nip++) {
                int mapf1 = job->qnlist_full.qnum[nf][nfp].map1;
                int mapf2 = job->qnlist_full.qnum[nf][nfp].map2;
                int mapi1 = job->qnlist_full.qnum[ni][nip].map1;
                int mapi2 = job->qnlist_full.qnum[ni][nip].map2;
                basis_3b qf = job->qnlist_full.qnum[nf][nfp];
                basis_3b qi = job->qnlist_full.qnum[ni][nip];
                ckdk_scdk *sc0 = &job->scdk[0][mapf1][mapf2][mapi1][mapi2];
                double iso = mol3_iso_bracket(qf.t1, qf.t2, qf.t3, qf.T,
                    qf.c, qf.tij, qi.c, qi.tij);

                if (fabs(iso) > 1e-16) {
                    nfi += iso * inteVcenPartA(tir_cent, *sc0, qf, qi, inteNfi, &job->varg, 1);
                }
                for (int p = 0; p < 3; p++) {
                    int c = p + 1;
                    ckdk_scdk *scp = &job->scdk[p][mapf1][mapf2][mapi1][mapi2];
                    if (fabs(iso) > 1e-16) {
                        T[p] += iso * inteVcenPartA(tpi_cent, *scp, qf, qi, mol3_inte_tnr, &job->varg, c);
                        rmsr[p] += iso * inteVcenPartA(tir_cent, *scp, qf, qi, inteRMS, &job->varg, c);
                    }
                    for (int k = 0; k < job->varg.n_tp[p]; k++) {
                        double uf = mol3_iso_bracket(qf.t1, qf.t2, qf.t3, qf.T,
                            qf.c, qf.tij, c, job->varg.tp_val[p][k]);
                        double ui = mol3_iso_bracket(qi.t1, qi.t2, qi.t3, qi.T,
                            qi.c, qi.tij, c, job->varg.tp_val[p][k]);
                        double w = uf * ui;
                        mol3_kernel_t *ker = &job->varg.ker[p][k];
                        mol3_iv_t iv;
                        if (fabs(w) < 1e-16) {
                            continue;
                        }
                        iv.ctx = ker->ctx;
                        iv.pot = ker->vc;
                        Vc[p] += w * inteVcenPartA(tir_cent, *scp, qf, qi, mol3_inte_v, &iv, c);
                        if (ker->vss != NULL) {
                            iv.pot = ker->vss;
                            Vss[p] += w * inteVcenPartA(tir_cent,
                                job->scdk[3 + p][mapf1][mapf2][mapi1][mapi2],
                                qf, qi, mol3_inte_v, &iv, c);
                        }
                        if (job->varg.has_vt[p] && ker->vt != NULL) {
                            iv.pot = ker->vt;
                            Vt[p] += w * inteVcenPartA(t1r_tens,
                                job->scdk[6 + p][mapf1][mapf2][mapi1][mapi2],
                                qf, qi, mol3_inte_v, &iv, c);
                        }
                    }
                }
            }
        }

        job->Nfi.value[nf][ni] = nfi;
        for (int p = 0; p < 3; p++) {
            job->T[p].value[nf][ni] = T[p];
            job->Vc[p].value[nf][ni] = Vc[p];
            job->Vss[p].value[nf][ni] = Vss[p];
            job->Vt[p].value[nf][ni] = Vt[p];
            job->rmsr[p].value[nf][ni] = rmsr[p];
        }
    }
}

static void mol3_addH(matrix_t *H, const matrix_t *vt, const matrix_t *M, matrix_t *tmp, double scale)
{
    int i, j;

    matrix_productT(vt, M, tmp);
    if (scale != 1.0) {
        for (i = 0; i < tmp->row; i++) {
            for (j = 0; j < tmp->col; j++) {
                tmp->value[i][j] *= scale;
            }
        }
    }
    matrix_sum(H, tmp, H);
}

/* ---------- driver ---------- */

void mol3_hybrid_GEM(const argsInput_t *args_input,
    array_t *e_out, matrix_t *v_out, matrix_t *n_out,
    array_t *rms12, array_t *rms13, array_t *rms23)
{
    mol3_job_t job = {0};
    matrix_t *mats[16];
    matrix_t vt, temp, Hfi, ut;
    double mtot;
    int nbas, nsp, nmats, nkeep;

    mol3_load_pairs(args_input, &job.varg);
    mol3_basis_build(args_input, &job.qnlist_spfy, &job.qnlist_full);

    nbas = job.qnlist_full.len_list;
    if (nbas <= 0) {
        fprintf(stderr, "Error: empty MOLECUL3 basis. Check J, P, Lmax, sym12, I, and leg types.\n");
        exit(1);
    }
    mol3_iso_cache_build(&job);
    mol3_mlsj_jl(&job);

    nsp = job.qnlist_spfy.len_list;
    for (int op = 0; op < MOL3_NOP; op++) {
        mol3_scdk_table_alloc(job.qnlist_spfy.len_part, nsp, &job.scdk[op]);
    }
    thread_load(calc_scdk_mt, &job, MOL3_NOP * nsp, getNumCores());

    nmats = 0;
    mats[nmats++] = &job.Nfi;
    for (int p = 0; p < 3; p++) {
        mats[nmats++] = &job.T[p];
        mats[nmats++] = &job.Vc[p];
        mats[nmats++] = &job.Vss[p];
        mats[nmats++] = &job.Vt[p];
        mats[nmats++] = &job.rmsr[p];
    }
    for (int i = 0; i < nmats; i++) {
        *mats[i] = matrix_init(nbas, nbas);
    }

    thread_load(getmfi, &job, nbas, getNumCores());
    for (int i = 0; i < nmats; i++) {
        matrix_symmetrize(mats[i]);
    }

    nkeep = matrix_orthogonal(&job.Nfi, 1e-8, &vt);
    if (nkeep <= 0) {
        fprintf(stderr, "Error: MOLECUL3 overlap N has no positive eigenvalues.\n");
        exit(1);
    }
    printf("  overlap rank: %d / %d\n", nkeep, nbas);
    fflush(stdout);

    temp = matrix_init(nkeep, nkeep);
    Hfi = matrix_init(nkeep, nkeep);
    mtot = args_input->m1 + args_input->m2 + args_input->m3;
    mol3_addH(&Hfi, &vt, &job.Nfi, &temp, mtot);
    for (int p = 0; p < 3; p++) {
        mol3_addH(&Hfi, &vt, &job.T[p], &temp, 1.0);
        mol3_addH(&Hfi, &vt, &job.Vc[p], &temp, 1.0);
        mol3_addH(&Hfi, &vt, &job.Vss[p], &temp, 1.0);
        mol3_addH(&Hfi, &vt, &job.Vt[p], &temp, 1.0);
    }

    *e_out = array_init(nkeep);
    ut = (matrix_t){0};
    if (v_out != NULL) {
        ut = matrix_init(nkeep, nkeep);
    }
#ifdef LAPACKE
    lapack_standard(Hfi.value, nkeep, e_out->value, (v_out != NULL) ? ut.value : NULL,
        (v_out != NULL) ? nkeep : 0);
#else
    eigen_standard(Hfi.value, nkeep, e_out->value, (v_out != NULL) ? ut.value : NULL,
        (v_out != NULL) ? nkeep : 0);
#endif

    if (v_out != NULL) {
        *v_out = matrix_init(nkeep, nbas);
        matrix_product(&ut, &vt, v_out);
        matrix_free(&ut);
        if (rms12 != NULL) *rms12 = array_init(nkeep);
        if (rms13 != NULL) *rms13 = array_init(nkeep);
        if (rms23 != NULL) *rms23 = array_init(nkeep);
        for (int n = 0; n < nkeep; n++) {
            if (rms12 != NULL) {
                rms12->value[n] = sqrt(fabs(matrix_expect(v_out, n, &job.rmsr[0]))) / GEMSTORE_FM;
            }
            if (rms13 != NULL) {
                rms13->value[n] = sqrt(fabs(matrix_expect(v_out, n, &job.rmsr[1]))) / GEMSTORE_FM;
            }
            if (rms23 != NULL) {
                rms23->value[n] = sqrt(fabs(matrix_expect(v_out, n, &job.rmsr[2]))) / GEMSTORE_FM;
            }
        }
    }
    if (n_out != NULL) {
        *n_out = matrix_init(nbas, nbas);
        matrix_copy(&job.Nfi, n_out);
    }

    if (v_out != NULL && args_input->print_wfn) {
        char path[280];
        sprintf(path, "%s.basis.dat", args_input->project);
        FILE *bf = fopen(path, "w");
        if (bf) {
            fprintf(bf, "# i  c  nrho nlam lrho llam L  sij  jl  tij  T  coe\n");
            for (int i = 0; i < nbas; i++) {
                basis_3b q = job.qnlist_full.qnum[i][0];
                fprintf(bf, "%d  %d  %d %d  %d %d %d  %.1f  %.1f  %.1f  %.1f  %.6f\n",
                    i, q.c, q.nrho, q.nlam, q.lrho, q.llam, q.L, q.sij, q.jl, q.tij, q.T, q.coe);
            }
            fclose(bf);
        }
        for (int n = 0; n < nkeep; n++) {
            sprintf(path, "%s.wfn.%d.dat", args_input->project, n + 1);
            FILE *wf = fopen(path, "w");
            if (!wf) {
                continue;
            }
            fprintf(wf, "# i  coefficient   (see %s.basis.dat)\n", args_input->project);
            for (int i = 0; i < nbas; i++) {
                fprintf(wf, "%d  %.16e\n", i, v_out->value[n][i]);
            }
            fclose(wf);
        }
    }

    mol3_iso_cache_free(&job, nsp);
    for (int op = 0; op < MOL3_NOP; op++) {
        mol3_scdk_table_free(job.qnlist_spfy.len_part, nsp, &job.scdk[op]);
    }
    for (int i = 0; i < nsp; i++) {
        for (int j = 0; j < job.qnlist_spfy.len_part[i]; j++) {
            matrix_free(&job.mlsj[i][j]);
        }
        free(job.mlsj[i]);
    }
    free(job.mlsj);
    bases_3b_free(&job.qnlist_spfy);
    bases_3b_free(&job.qnlist_full);
    matrix_free(&temp);
    matrix_free(&vt);
    matrix_free(&Hfi);
    for (int i = 0; i < nmats; i++) {
        matrix_free(mats[i]);
    }
}
