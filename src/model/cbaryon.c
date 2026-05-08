/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gemstore/model/cbaryon.h>

#include <gemstore/basis/basis.h>
#include <gemstore/basis/color.h>
#include <gemstore/basis/orbit.h>
#include <gemstore/basis/spin.h>
#include <gemstore/math/cmi.h>
#include <gemstore/math/eigen.h>
#include <gemstore/math/matrix.h>
#include <gemstore/model/gimodel.h>
#include <gemstore/print.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define BARYON_RADIAL_ORDER 16
#define BARYON_ANGLE_ORDER 16

typedef struct {
    double x;
    double y;
    double z;
} vec3_t;

typedef struct {
    int num_spec;
    int **index_map;
    intrin_wfn_t *spin_wf;
    matrix_t spin_overlap;
    matrix_t spin_sigma;
} baryon_spin_cache_t;

static vec3_t vec3_make(double x, double y, double z)
{
    vec3_t v = {x, y, z};
    return v;
}

static vec3_t vec3_add(vec3_t a, vec3_t b)
{
    return vec3_make(a.x + b.x, a.y + b.y, a.z + b.z);
}

static vec3_t vec3_sub(vec3_t a, vec3_t b)
{
    return vec3_make(a.x - b.x, a.y - b.y, a.z - b.z);
}

static vec3_t vec3_scale(vec3_t a, double s)
{
    return vec3_make(a.x * s, a.y * s, a.z * s);
}

static double vec3_dot(vec3_t a, vec3_t b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

static double vec3_norm2(vec3_t a)
{
    return vec3_dot(a, a);
}

static double vec3_norm(vec3_t a)
{
    return sqrt(vec3_norm2(a));
}

static void gauss_legendre_rule(int n, double *nodes, double *weights)
{
    const double eps = 1e-14;
    const int m = (n + 1) / 2;

    for (int i = 0; i < m; i++) {
        double z = cos(M_PI * (i + 0.75) / (n + 0.5));

        for (;;) {
            double p1 = 1.0;
            double p2 = 0.0;
            for (int j = 1; j <= n; j++) {
                double p3 = p2;
                p2 = p1;
                p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3) / j;
            }
            double pp = n * (z * p1 - p2) / (z * z - 1.0);
            double z1 = z;
            z = z1 - p1 / pp;
            if (fabs(z - z1) < eps) {
                nodes[i] = -z;
                nodes[n - 1 - i] = z;
                weights[i] = 2.0 / ((1.0 - z * z) * pp * pp);
                weights[n - 1 - i] = weights[i];
                break;
            }
        }
    }
}

static void jacobi_indices(int c, int *i, int *j, int *k)
{
    if (c == 1) {
        *i = 0;
        *j = 1;
        *k = 2;
    }
    else if (c == 2) {
        *i = 0;
        *j = 2;
        *k = 1;
    }
    else {
        *i = 1;
        *j = 2;
        *k = 0;
    }
}

static void positions_from_jacobi(int c, const double mass[3], vec3_t rho, vec3_t lambda, vec3_t r[3])
{
    int i, j, k;
    jacobi_indices(c, &i, &j, &k);

    const double mij = mass[i] + mass[j];
    const double mtot = mass[0] + mass[1] + mass[2];

    r[i] = vec3_add(vec3_scale(lambda, -mass[k] / mtot), vec3_scale(rho, mass[j] / mij));
    r[j] = vec3_add(vec3_scale(lambda, -mass[k] / mtot), vec3_scale(rho, -mass[i] / mij));
    r[k] = vec3_scale(lambda, mij / mtot);
}

static void jacobi_from_positions(int c, const double mass[3], const vec3_t r[3], vec3_t *rho, vec3_t *lambda)
{
    int i, j, k;
    jacobi_indices(c, &i, &j, &k);

    const double mij = mass[i] + mass[j];

    *rho = vec3_sub(r[i], r[j]);
    *lambda = vec3_sub(r[k], vec3_scale(vec3_add(vec3_scale(r[i], mass[i]), vec3_scale(r[j], mass[j])), 1.0 / mij));
}

static void momenta_from_jacobi(int c, const double mass[3], vec3_t prho, vec3_t plam, vec3_t p[3])
{
    int i, j, k;
    jacobi_indices(c, &i, &j, &k);

    const double mij = mass[i] + mass[j];

    p[i] = vec3_add(prho, vec3_scale(plam, -mass[i] / mij));
    p[j] = vec3_add(vec3_scale(prho, -1.0), vec3_scale(plam, -mass[j] / mij));
    p[k] = plam;
}

static void jacobi_from_momenta(int c, const double mass[3], const vec3_t p[3], vec3_t *prho, vec3_t *plam)
{
    int i, j, k;
    jacobi_indices(c, &i, &j, &k);

    const double mij = mass[i] + mass[j];

    *prho = vec3_scale(vec3_sub(vec3_scale(p[i], mass[j]), vec3_scale(p[j], mass[i])), 1.0 / mij);
    *plam = p[k];
}

static double gaussian_coord_norm(double nu)
{
    return GRnlr(1.0, 1, 0, nu);
}

static double gaussian_momentum_norm(double nu)
{
    return creal(GRnlp(1.0, 1, 0, nu));
}

static double baryon_orbital_coord(const basis_qnum *q, const vec3_t r[3])
{
    vec3_t rho;
    vec3_t lambda;
    const double mass[3] = {q->m1, q->m2, q->m3};

    jacobi_from_positions(q->c, mass, r, &rho, &lambda);
    return gaussian_coord_norm(q->nurho) * gaussian_coord_norm(q->nulam)
        * exp(-q->nurho * vec3_norm2(rho) - q->nulam * vec3_norm2(lambda));
}

static double baryon_orbital_momentum(const basis_qnum *q, const vec3_t p[3])
{
    vec3_t prho;
    vec3_t plam;
    const double mass[3] = {q->m1, q->m2, q->m3};

    jacobi_from_momenta(q->c, mass, p, &prho, &plam);
    return gaussian_momentum_norm(q->nurho) * gaussian_momentum_norm(q->nulam)
        * exp(-vec3_norm2(prho) / (4.0 * q->nurho) - vec3_norm2(plam) / (4.0 * q->nulam));
}

static intrin_wfn_t spin_wfn_baryon_pair(int c, double sij, double st, double st3)
{
    intrin_wfn_t wf = intrin_wfn_init(3);
    cg_table_t couple = {0};
    int order[3];

    if (c == 1) {
        order[0] = 0;
        order[1] = 1;
        order[2] = 2;
    }
    else if (c == 2) {
        order[0] = 0;
        order[1] = 2;
        order[2] = 1;
    }
    else {
        order[0] = 1;
        order[1] = 2;
        order[2] = 0;
    }

    if (sij == 0.0) {
        intrin_wfn_t pair = spin_wfn_meson(0.0, 0.0);
        intrin_wfn_t third = spin_basis(st3);
        intrin_wfn_t source = intrin_wfn_init(3);

        intrin_wfn_product(&pair, &third, 1.0, &source);
        for (int i = 0; i < source.num_terms; i++) {
            char config[4] = {'0', '0', '0', '\0'};
            config[order[0]] = source.configs[i][0];
            config[order[1]] = source.configs[i][1];
            config[order[2]] = source.configs[i][2];
            intrin_wfn_push(&wf, source.coeffs[i], config);
        }

        intrin_wfn_free(&source);
        intrin_wfn_free(&pair);
        intrin_wfn_free(&third);
    }
    else {
        couple = CallCGTable(sij, 0.5, st, st3);
        for (int i = 0; i < couple.num; i++) {
            intrin_wfn_t pair = spin_wfn_meson(sij, couple.tuples[i].ms1);
            intrin_wfn_t third = spin_basis(couple.tuples[i].ms2);
            intrin_wfn_t source = intrin_wfn_init(3);

            intrin_wfn_product(&pair, &third, couple.tuples[i].cg, &source);
            for (int j = 0; j < source.num_terms; j++) {
                char config[4] = {'0', '0', '0', '\0'};
                config[order[0]] = source.configs[j][0];
                config[order[1]] = source.configs[j][1];
                config[order[2]] = source.configs[j][2];
                intrin_wfn_push(&wf, source.coeffs[j], config);
            }

            intrin_wfn_free(&source);
            intrin_wfn_free(&pair);
            intrin_wfn_free(&third);
        }
        cg_table_free(&couple);
    }

    intrin_wfn_trim(&wf);
    return wf;
}

static int pair_column(int i, int j)
{
    return (i * (2 * 3 - i - 1)) / 2 + (j - i - 1);
}

static double baryon_pair_sdots(const baryon_spin_cache_t *cache, int spec_f, int spec_i, int pair)
{
    return 0.25 * cache->spin_sigma.value[spec_f * cache->num_spec + spec_i][pair];
}

static double baryon_pair_potential(double r, double mi, double mj, double spin_overlap, double sdot,
    model_type_t model_type, const argsGIModel_t *args_model)
{
    argsGIModelDy_t dyn = {0};
    double sigmaij;

    dyn.model = model_type;
    dyn.system = SYSTEM_BARYON;
    dyn.mi = mi;
    dyn.mj = mj;
    dyn.Cij = -2.0 / 3.0;
    dyn.OCent = spin_overlap;
    dyn.OSdS = sdot;
    sigmaij = sigma_ij(mi, mj, args_model->sigma_0, args_model->s);
    dyn.Sigij = sigmaij;
    sigma_k_ij(sigmaij, dyn.Sigkij);

    return GIVconf(r, args_model, &dyn)
        + GIVcoul(r, args_model, &dyn)
        + GIVcont(r, args_model, &dyn);
}

static double baryon_coordinate_cutoff(const basis_qnum *qf, const basis_qnum *qi)
{
    double min_nu = qf->nurho;
    if (qf->nulam < min_nu) min_nu = qf->nulam;
    if (qi->nurho < min_nu) min_nu = qi->nurho;
    if (qi->nulam < min_nu) min_nu = qi->nulam;
    return 8.0 / sqrt(min_nu);
}

static double baryon_momentum_cutoff(const basis_qnum *qf, const basis_qnum *qi)
{
    double max_nu = qf->nurho;
    if (qf->nulam > max_nu) max_nu = qf->nulam;
    if (qi->nurho > max_nu) max_nu = qi->nurho;
    if (qi->nulam > max_nu) max_nu = qi->nulam;
    return 8.0 * sqrt(max_nu);
}

static double integrate_baryon_overlap(const basis_qnum *qf, const basis_qnum *qi,
    const double *radial_nodes, const double *radial_weights,
    const double *angle_nodes, const double *angle_weights)
{
    const double mass[3] = {qf->m1, qf->m2, qf->m3};
    const double cutoff = baryon_coordinate_cutoff(qf, qi);
    double sum = 0.0;

    for (int ir = 0; ir < BARYON_RADIAL_ORDER; ir++) {
        double rho = 0.5 * cutoff * (radial_nodes[ir] + 1.0);
        double wrho = 0.5 * cutoff * radial_weights[ir];
        for (int il = 0; il < BARYON_RADIAL_ORDER; il++) {
            double lambda = 0.5 * cutoff * (radial_nodes[il] + 1.0);
            double wlambda = 0.5 * cutoff * radial_weights[il];
            for (int ix = 0; ix < BARYON_ANGLE_ORDER; ix++) {
                double costh = angle_nodes[ix];
                double sinth = sqrt(fmax(0.0, 1.0 - costh * costh));
                vec3_t rho_vec = vec3_make(0.0, 0.0, rho);
                vec3_t lambda_vec = vec3_make(lambda * sinth, 0.0, lambda * costh);
                vec3_t r[3];
                double measure = 8.0 * M_PI * M_PI * rho * rho * lambda * lambda;

                positions_from_jacobi(qf->c, mass, rho_vec, lambda_vec, r);
                sum += wrho * wlambda * angle_weights[ix] * measure
                    * baryon_orbital_coord(qf, r)
                    * baryon_orbital_coord(qi, r);
            }
        }
    }

    return sum;
}

static double integrate_baryon_kinetic(const basis_qnum *qf, const basis_qnum *qi,
    const double *radial_nodes, const double *radial_weights,
    const double *angle_nodes, const double *angle_weights)
{
    const double mass[3] = {qf->m1, qf->m2, qf->m3};
    const double cutoff = baryon_momentum_cutoff(qf, qi);
    double sum = 0.0;

    for (int ir = 0; ir < BARYON_RADIAL_ORDER; ir++) {
        double prho = 0.5 * cutoff * (radial_nodes[ir] + 1.0);
        double wrho = 0.5 * cutoff * radial_weights[ir];
        for (int il = 0; il < BARYON_RADIAL_ORDER; il++) {
            double plam = 0.5 * cutoff * (radial_nodes[il] + 1.0);
            double wlambda = 0.5 * cutoff * radial_weights[il];
            for (int ix = 0; ix < BARYON_ANGLE_ORDER; ix++) {
                double costh = angle_nodes[ix];
                double sinth = sqrt(fmax(0.0, 1.0 - costh * costh));
                vec3_t prho_vec = vec3_make(0.0, 0.0, prho);
                vec3_t plam_vec = vec3_make(plam * sinth, 0.0, plam * costh);
                vec3_t p[3];
                double kinetic;
                double measure = 8.0 * M_PI * M_PI * prho * prho * plam * plam;

                momenta_from_jacobi(qf->c, mass, prho_vec, plam_vec, p);
                kinetic = vec3_norm2(p[0]) / (2.0 * mass[0])
                    + vec3_norm2(p[1]) / (2.0 * mass[1])
                    + vec3_norm2(p[2]) / (2.0 * mass[2]);
                sum += wrho * wlambda * angle_weights[ix] * measure
                    * baryon_orbital_momentum(qf, p)
                    * baryon_orbital_momentum(qi, p)
                    * kinetic;
            }
        }
    }

    return sum;
}

static double integrate_baryon_potential(const basis_qnum *qf, const basis_qnum *qi,
    const baryon_spin_cache_t *spin_cache, model_type_t model_type, const argsGIModel_t *args_model,
    const double *radial_nodes, const double *radial_weights,
    const double *angle_nodes, const double *angle_weights)
{
    const double mass[3] = {qf->m1, qf->m2, qf->m3};
    const double cutoff = baryon_coordinate_cutoff(qf, qi);
    const int spec_f = spin_cache->index_map[qf->map1][qf->map2];
    const int spec_i = spin_cache->index_map[qi->map1][qi->map2];
    const double spin_overlap = spin_cache->spin_overlap.value[spec_f][spec_i];
    double sum = 0.0;

    for (int ir = 0; ir < BARYON_RADIAL_ORDER; ir++) {
        double rho = 0.5 * cutoff * (radial_nodes[ir] + 1.0);
        double wrho = 0.5 * cutoff * radial_weights[ir];
        for (int il = 0; il < BARYON_RADIAL_ORDER; il++) {
            double lambda = 0.5 * cutoff * (radial_nodes[il] + 1.0);
            double wlambda = 0.5 * cutoff * radial_weights[il];
            for (int ix = 0; ix < BARYON_ANGLE_ORDER; ix++) {
                double costh = angle_nodes[ix];
                double sinth = sqrt(fmax(0.0, 1.0 - costh * costh));
                vec3_t rho_vec = vec3_make(0.0, 0.0, rho);
                vec3_t lambda_vec = vec3_make(lambda * sinth, 0.0, lambda * costh);
                vec3_t r[3];
                double measure = 8.0 * M_PI * M_PI * rho * rho * lambda * lambda;
                double v = 0.0;

                positions_from_jacobi(qf->c, mass, rho_vec, lambda_vec, r);
                v += baryon_pair_potential(vec3_norm(vec3_sub(r[0], r[1])), mass[0], mass[1], spin_overlap,
                    baryon_pair_sdots(spin_cache, spec_f, spec_i, pair_column(0, 1)), model_type, args_model);
                v += baryon_pair_potential(vec3_norm(vec3_sub(r[0], r[2])), mass[0], mass[2], spin_overlap,
                    baryon_pair_sdots(spin_cache, spec_f, spec_i, pair_column(0, 2)), model_type, args_model);
                v += baryon_pair_potential(vec3_norm(vec3_sub(r[1], r[2])), mass[1], mass[2], spin_overlap,
                    baryon_pair_sdots(spin_cache, spec_f, spec_i, pair_column(1, 2)), model_type, args_model);

                sum += wrho * wlambda * angle_weights[ix] * measure
                    * baryon_orbital_coord(qf, r)
                    * baryon_orbital_coord(qi, r)
                    * v;
            }
        }
    }

    return sum;
}

static void build_baryon_specific_basis(const argsInput_t *input, const argsGIModel_t *args_model,
    basis_list *qnlist_spfy, basis_list *qnlist_full)
{
    const double m1 = getmq(input->f1, args_model);
    const double m2 = getmq(input->f2, args_model);
    const double m3 = getmq(input->f3, args_model);
    const double si = 0.5;
    const double sj = 0.5;
    const double sk = 0.5;
    const double t1 = 0.0;
    const double t2 = 0.0;
    const double t3 = 0.0;
    const double tij = 0.0;
    const double T = 0.0;

    basis_list_init(qnlist_spfy);
    basis_list_init(qnlist_full);

    for (int lrho = 0; lrho <= input->Lmax; lrho++) {
        for (int llam = 0; llam <= input->Lmax; llam++) {
            for (int L = abs(lrho - llam); L <= lrho + llam; L++) {
                for (double sij = 0.0; sij <= 1.0; sij += 1.0) {
                    for (double S = fabs(sij - sk); S <= fabs(sij + sk) + 1e-12; S += 1.0) {
                        if (lrho + llam > input->Lmax) continue;
                        if (input->P != (int)pow(-1.0, lrho + llam)) continue;
                        if (!(fabs(S - L) <= input->J && input->J <= S + L)) continue;

                        if (input->f12 == (int)pow(-1.0, 1.0 + sij + lrho)) {
                            basis_list_push(qnlist_spfy, 1, -1, -1, 1.0,
                                m1, m2, m3,
                                si, sj, sk,
                                t1, t2, t3, tij, T,
                                1, lrho, llam, L, sij, S, input->J,
                                0, 0, 0.0, 0.0);
                        }

                        basis_list_push(qnlist_spfy, 1, -1, -1, 1.0,
                            m1, m2, m3,
                            si, sj, sk,
                            t1, t2, t3, tij, T,
                            2, lrho, llam, L, sij, S, input->J,
                            0, 0, 0.0, 0.0);

                        basis_list_push(qnlist_spfy, 0, -1, -1, input->f12 * pow(-1.0, 1.0 + sij + lrho),
                            m1, m2, m3,
                            si, sj, sk,
                            t1, t2, t3, tij, T,
                            3, lrho, llam, L, sij, S, input->J,
                            0, 0, 0.0, 0.0);
                    }
                }
            }
        }
    }

    basis_list_push_full(qnlist_spfy, qnlist_full, input->rmin, input->rmax, input->nmax);
}

static void baryon_spin_cache_init(baryon_spin_cache_t *cache, const basis_list *qnlist_spfy)
{
    cache->num_spec = 0;
    cache->index_map = (int **)malloc(sizeof(int *) * qnlist_spfy->len_list);
    for (int i = 0; i < qnlist_spfy->len_list; i++) {
        cache->index_map[i] = (int *)malloc(sizeof(int) * qnlist_spfy->len_part[i]);
        for (int j = 0; j < qnlist_spfy->len_part[i]; j++) {
            cache->index_map[i][j] = cache->num_spec++;
        }
    }

    cache->spin_wf = (intrin_wfn_t *)malloc(sizeof(intrin_wfn_t) * cache->num_spec);
    for (int i = 0; i < qnlist_spfy->len_list; i++) {
        for (int j = 0; j < qnlist_spfy->len_part[i]; j++) {
            int idx = cache->index_map[i][j];
            const basis_qnum *q = &qnlist_spfy->qnum[i][j];
            cache->spin_wf[idx] = spin_wfn_baryon_pair(q->c, q->sij, q->jl, q->J);
        }
    }

    cache->spin_overlap = matrix_init(cache->num_spec, cache->num_spec);
    cache->spin_sigma = matrix_init(cache->num_spec * cache->num_spec, 3);
    for (int i = 0; i < cache->num_spec; i++) {
        for (int j = 0; j < cache->num_spec; j++) {
            cache->spin_overlap.value[i][j] = intrin_wfn_overlap(&cache->spin_wf[i], &cache->spin_wf[j]);
        }
    }
    for (int i = 0; i < cache->spin_sigma.row; i++) {
        for (int j = 0; j < cache->spin_sigma.col; j++) {
            cache->spin_sigma.value[i][j] = 0.0;
        }
    }
    operator_sigma2(cache->spin_wf, cache->num_spec, &cache->spin_sigma);
}

static void baryon_spin_cache_free(baryon_spin_cache_t *cache, const basis_list *qnlist_spfy)
{
    for (int i = 0; i < cache->num_spec; i++) {
        intrin_wfn_free(&cache->spin_wf[i]);
    }
    free(cache->spin_wf);
    for (int i = 0; i < qnlist_spfy->len_list; i++) {
        free(cache->index_map[i]);
    }
    free(cache->index_map);
    matrix_free(&cache->spin_overlap);
    matrix_free(&cache->spin_sigma);
}

void spectra_baryon_GEM(const argsInput_t *input, array_t *e_out, matrix_t *v_out, matrix_t *n_out, int v_len,
    int *basis_len_out)
{
    argsGIModel_t args_model = argsGIModel_from(input);
    basis_list qnlist_spfy;
    basis_list qnlist_full;
    baryon_spin_cache_t spin_cache;
    double radial_nodes[BARYON_RADIAL_ORDER];
    double radial_weights[BARYON_RADIAL_ORDER];
    double angle_nodes[BARYON_ANGLE_ORDER];
    double angle_weights[BARYON_ANGLE_ORDER];

    if (input->Lmax != 0) {
        fprintf(stderr, "Error: native baryon GEM currently supports only Lmax = 0.\n");
        exit(1);
    }

    build_baryon_specific_basis(input, &args_model, &qnlist_spfy, &qnlist_full);
    if (qnlist_full.len_list <= 0) {
        fprintf(stderr, "Error: baryon basis generation produced no states.\n");
        basis_list_free(&qnlist_spfy);
        basis_list_free(&qnlist_full);
        exit(1);
    }

    if (basis_len_out != NULL) {
        *basis_len_out = qnlist_full.len_list;
    }

    gauss_legendre_rule(BARYON_RADIAL_ORDER, radial_nodes, radial_weights);
    gauss_legendre_rule(BARYON_ANGLE_ORDER, angle_nodes, angle_weights);
    baryon_spin_cache_init(&spin_cache, &qnlist_spfy);

    matrix_t Hfi = matrix_init(qnlist_full.len_list, qnlist_full.len_list);
    matrix_t Nfi = matrix_init(qnlist_full.len_list, qnlist_full.len_list);
    for (int i = 0; i < qnlist_full.len_list; i++) {
        for (int j = 0; j < qnlist_full.len_list; j++) {
            Hfi.value[i][j] = 0.0;
            Nfi.value[i][j] = 0.0;
        }
    }

    for (int nf = 0; nf < qnlist_full.len_list; nf++) {
        for (int ni = nf; ni < qnlist_full.len_list; ni++) {
            double hij = 0.0;
            double nij = 0.0;

            for (int nfp = 0; nfp < qnlist_full.len_part[nf]; nfp++) {
                const basis_qnum *qf = &qnlist_full.qnum[nf][nfp];
                const int spec_f = spin_cache.index_map[qf->map1][qf->map2];
                for (int nip = 0; nip < qnlist_full.len_part[ni]; nip++) {
                    const basis_qnum *qi = &qnlist_full.qnum[ni][nip];
                    const int spec_i = spin_cache.index_map[qi->map1][qi->map2];
                    const double coeff = qf->coe * qi->coe;
                    const double spin_overlap = spin_cache.spin_overlap.value[spec_f][spec_i];
                    double overlap_orb;
                    double kinetic_orb;
                    double potential_orb;

                    overlap_orb = integrate_baryon_overlap(qf, qi, radial_nodes, radial_weights, angle_nodes, angle_weights);
                    kinetic_orb = integrate_baryon_kinetic(qf, qi, radial_nodes, radial_weights, angle_nodes, angle_weights);
                    potential_orb = integrate_baryon_potential(qf, qi, &spin_cache, input->model,
                        &args_model, radial_nodes, radial_weights, angle_nodes, angle_weights);

                    nij += coeff * spin_overlap * overlap_orb;
                    hij += coeff * (spin_overlap * kinetic_orb + potential_orb);
                }
            }

            Nfi.value[nf][ni] = nij;
            Nfi.value[ni][nf] = nij;
            Hfi.value[nf][ni] = hij;
            Hfi.value[ni][nf] = hij;
        }
    }

    array_t overlap_eval = array_init(qnlist_full.len_list);
    matrix_t overlap_vec = matrix_init(qnlist_full.len_list, qnlist_full.len_list);
    eigen_standard(Nfi.value, qnlist_full.len_list, overlap_eval.value, overlap_vec.value, qnlist_full.len_list);

    double max_overlap = overlap_eval.value[qnlist_full.len_list - 1];
    double overlap_cut = fmax(max_overlap * 1e-10, 1e-12);
    int keep = 0;
    for (int i = 0; i < qnlist_full.len_list; i++) {
        if (overlap_eval.value[i] > overlap_cut) {
            keep++;
        }
    }
    if (keep <= 0) {
        fprintf(stderr, "Error: baryon overlap matrix has no stable subspace.\n");
        exit(1);
    }

    matrix_t vt = matrix_init(keep, qnlist_full.len_list);
    int row = 0;
    for (int i = 0; i < qnlist_full.len_list; i++) {
        if (overlap_eval.value[i] <= overlap_cut) {
            continue;
        }
        double norm = sqrt(overlap_eval.value[i]);
        for (int j = 0; j < qnlist_full.len_list; j++) {
            vt.value[row][j] = overlap_vec.value[i][j] / norm;
        }
        row++;
    }

    matrix_t Horth = matrix_init(keep, keep);
    matrix_productT(&vt, &Hfi, &Horth);

    if (v_len < keep || e_out == NULL || v_out == NULL || n_out == NULL) {
        fprintf(stderr, "Error: baryon output buffers are invalid or too small.\n");
        exit(1);
    }

    array_t eigenvalue = array_init(keep);
    matrix_t ut = matrix_init(keep, keep);
    matrix_t eigenvector = matrix_init(keep, qnlist_full.len_list);

    eigen_standard(Horth.value, keep, eigenvalue.value, ut.value, keep);
    matrix_product(&ut, &vt, &eigenvector);

    e_out->len = keep;
    for (int i = 0; i < keep; i++) {
        e_out->value[i] = eigenvalue.value[i];
    }

    for (int i = 0; i < qnlist_full.len_list; i++) {
        for (int j = 0; j < qnlist_full.len_list; j++) {
            n_out->value[i][j] = Nfi.value[i][j];
        }
    }

    for (int i = 0; i < keep; i++) {
        for (int j = 0; j < qnlist_full.len_list; j++) {
            v_out->value[i][j] = eigenvector.value[i][j];
        }
    }
    for (int i = keep; i < v_out->row; i++) {
        for (int j = 0; j < v_out->col; j++) {
            v_out->value[i][j] = 0.0;
        }
    }

    array_free(&overlap_eval);
    array_free(&eigenvalue);
    matrix_free(&overlap_vec);
    matrix_free(&vt);
    matrix_free(&Horth);
    matrix_free(&ut);
    matrix_free(&eigenvector);
    matrix_free(&Hfi);
    matrix_free(&Nfi);
    baryon_spin_cache_free(&spin_cache, &qnlist_spfy);
    basis_list_free(&qnlist_spfy);
    basis_list_free(&qnlist_full);
}
