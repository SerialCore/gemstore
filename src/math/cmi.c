/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gemstore/math/cmi.h>
#include <gemstore/math/matrix.h>

#include <gemstore/basis/intrin.h>

#include <ctype.h>
#include <complex.h>

/* Helper to get spin vector */
static inline void get_swv(char c, double v[2]);
static inline void get_swv(char c, double v[2])
{
    if (c == '1') {
        v[0] = 1.0;
        v[1] = 0.0;
    } else if (c == '0') {
        v[0] = 0.0;
        v[1] = 1.0;
    } else {
        // Error handling
        return;
    }
}

/* Helper to get color vector */
static inline void get_cwv(char c, double v[3]);
static inline void get_cwv(char c, double v[3])
{
    char lc = tolower(c);
    if (lc == 'r') {
        v[0] = 1.0; v[1] = 0.0; v[2] = 0.0;
    } else if (lc == 'g') {
        v[0] = 0.0; v[1] = 1.0; v[2] = 0.0;
    } else if (lc == 'b') {
        v[0] = 0.0; v[1] = 0.0; v[2] = 1.0;
    } else {
        // Error handling
        return;
    }
}

/* Trace for spin-sigma-spin */
static inline complex trace_spin(const double p[2], const complex m[4], const double v[2]);
static inline complex trace_spin(const double p[2], const complex m[4], const double v[2])
{
    complex sum = 0.0 + 0.0 * I;
    for (int a = 0; a < 2; a++) {
        complex mv = 0.0 + 0.0 * I;
        for (int b = 0; b < 2; b++) {
            mv += m[a * 2 + b] * v[b];
        }
        sum += p[a] * mv;
    }
    return sum;
}

/* Trace for color-lambda-color */
static inline complex trace_color(const double p[3], const complex m[9], const double v[3]);
static inline complex trace_color(const double p[3], const complex m[9], const double v[3])
{
    complex sum = 0.0 + 0.0 * I;
    for (int a = 0; a < 3; a++) {
        complex mv = 0.0 + 0.0 * I;
        for (int b = 0; b < 3; b++) {
            mv += m[a * 3 + b] * v[b];
        }
        sum += p[a] * mv;
    }
    return sum;
}

/* Compute σ_i · σ_j matrix element between single quark states */
static double sigma_dot_sigma(char fi, char ii, char fj, char ij);
static double sigma_dot_sigma(char fi, char ii, char fj, char ij)
{
    double sum = 0.0;
    double fvi[2] = {0.0, 0.0};
    double ivi[2] = {0.0, 0.0};
    double fvj[2] = {0.0, 0.0};
    double ivj[2] = {0.0, 0.0};
    complex sigma[3][4] = {
        {
            0.0 + 0.0 * I, 1.0 + 0.0 * I,
            1.0 + 0.0 * I, 0.0 + 0.0 * I
        },
        {
            0.0 + 0.0 * I, 0.0 - 1.0 * I,
            0.0 + 1.0 * I, 0.0 + 0.0 * I
        },
        {
            1.0 + 0.0 * I, 0.0 + 0.0 * I,
            0.0 + 0.0 * I, -1.0 + 0.0 * I
        }
    };

    get_swv(fi, fvi);
    get_swv(ii, ivi);
    get_swv(fj, fvj);
    get_swv(ij, ivj);
    for (int x = 0; x < 3; x++) {
        complex tr1 = trace_spin(fvi, sigma[x], ivi);
        complex tr2 = trace_spin(fvj, sigma[x], ivj);
        sum += creal(tr1 * tr2);
    }

    return sum;
}

/* Compute λ_i · λ_j matrix element between single quark states */
static double lambda_dot_lambda(char fi, char ii, char fj, char ij);
static double lambda_dot_lambda(char fi, char ii, char fj, char ij)
{
    double sum = 0.0;
    double fvi[3] = {0.0, 0.0, 0.0};
    double ivi[3] = {0.0, 0.0, 0.0};
    double fvj[3] = {0.0, 0.0, 0.0};
    double ivj[3] = {0.0, 0.0, 0.0};
    double complex lambda[8][9] = {
        {
            0.0 + 0.0 * I, 1.0 + 0.0 * I, 0.0 + 0.0 * I,
            1.0 + 0.0 * I, 0.0 + 0.0 * I, 0.0 + 0.0 * I,
            0.0 + 0.0 * I, 0.0 + 0.0 * I, 0.0 + 0.0 * I
        },
        {
            0.0 + 0.0 * I, 0.0 - 1.0 * I, 0.0 + 0.0 * I,
            0.0 + 1.0 * I, 0.0 + 0.0 * I, 0.0 + 0.0 * I,
            0.0 + 0.0 * I, 0.0 + 0.0 * I, 0.0 + 0.0 * I
        },
        {
            1.0 + 0.0 * I, 0.0 + 0.0 * I, 0.0 + 0.0 * I,
            0.0 + 0.0 * I, -1.0 + 0.0 * I, 0.0 + 0.0 * I,
            0.0 + 0.0 * I, 0.0 + 0.0 * I, 0.0 + 0.0 * I
        },
        {
            0.0 + 0.0 * I, 0.0 + 0.0 * I, 1.0 + 0.0 * I,
            0.0 + 0.0 * I, 0.0 + 0.0 * I, 0.0 + 0.0 * I,
            1.0 + 0.0 * I, 0.0 + 0.0 * I, 0.0 + 0.0 * I
        },
        {
            0.0 + 0.0 * I, 0.0 + 0.0 * I, 0.0 - 1.0 * I,
            0.0 + 0.0 * I, 0.0 + 0.0 * I, 0.0 + 0.0 * I,
            0.0 + 1.0 * I, 0.0 + 0.0 * I, 0.0 + 0.0 * I
        },
        {
            0.0 + 0.0 * I, 0.0 + 0.0 * I, 0.0 + 0.0 * I,
            0.0 + 0.0 * I, 0.0 + 0.0 * I, 1.0 + 0.0 * I,
            0.0 + 0.0 * I, 1.0 + 0.0 * I, 0.0 + 0.0 * I
        },
        {
            0.0 + 0.0 * I, 0.0 + 0.0 * I, 0.0 + 0.0 * I,
            0.0 + 0.0 * I, 0.0 + 0.0 * I, 0.0 - 1.0 * I,
            0.0 + 0.0 * I, 0.0 + 1.0 * I, 0.0 + 0.0 * I
        },
        {
            0.57735026918962576451 + 0.0 * I, 0.0 + 0.0 * I, 0.0 + 0.0 * I,
            0.0 + 0.0 * I, 0.57735026918962576451 + 0.0 * I, 0.0 + 0.0 * I,
            0.0 + 0.0 * I, 0.0 + 0.0 * I, -1.1547005383792515290 + 0.0 * I
        }
    }; 

    get_cwv(fi, fvi);
    get_cwv(ii, ivi);
    get_cwv(fj, fvj);
    get_cwv(ij, ivj);
    for (int x = 0; x < 8; x++) {
        complex tr1 = trace_color(fvi, lambda[x], ivi);
        complex tr2 = trace_color(fvj, lambda[x], ivj);
        sum += creal(tr1 * tr2);
    }

    return sum;
}

/* match configs except positions 1 and 2 */
static inline int configs_match_except_positions(const char *config1, const char *config2, int num_configs, int skip1, int skip2);
static inline int configs_match_except_positions(const char *config1, const char *config2, int num_configs, int skip1, int skip2)
{
    for (int k = 0; k < num_configs; k++) {
        /* skip positions */
        if (k == skip1 || k == skip2) {
            continue;
        }
        /* match the rest positions */
        if (config1[k] != config2[k]) {
            return 0;
        }
    }
    return 1;
}

void operator_sigma2(const intrin_wfn_t swv[], int num_state, matrix_t *result)
{
    int num_configs = swv[0].num_configs;       /* assume all state has the same number of configs */

    /* traversing for hadron state */
    for (int nf = 0; nf < num_state; nf++) {
        for (int ni = 0; ni < num_state; ni++) {
            /* traversing for configs in each term */
            for (int i = 0; i < num_configs - 1; i++) {
                for (int j = i + 1; j < num_configs; j++) {
                    double ss_sum = 0.0;
                    /* traversing for terms in each state */
                    for (int cf = 0; cf < swv[nf].num_terms; cf++) {
                        for (int ci = 0; ci < swv[ni].num_terms; ci++) {
                            /* match configs except positions i and j */
                            if (!configs_match_except_positions(swv[nf].configs[cf], swv[ni].configs[ci], num_configs, i, j)) {
                                continue;
                            }
                            char fi = swv[nf].configs[cf][i];
                            char ii = swv[ni].configs[ci][i];
                            char fj = swv[nf].configs[cf][j];
                            char ij = swv[ni].configs[ci][j];
                            ss_sum += swv[nf].coeffs[cf] * swv[ni].coeffs[ci] * sigma_dot_sigma(fi, ii, fj, ij);
                        }
                    }
                    result->value[nf * num_state + ni][(i * (2 * num_configs - i - 1)) / 2 + (j - i - 1)] = ss_sum;
                }
            }
        }
    }
}

void operator_lambda2(const intrin_wfn_t cwv[], int num_state, const char *config, matrix_t *result)
{
    int num_configs = cwv[0].num_configs;       /* assume all state has the same number of configs */

    /* traversing for hadron state */
    for (int nf = 0; nf < num_state; nf++) {
        for (int ni = 0; ni < num_state; ni++) {
            /* traversing for configs in each term */
            for (int i = 0; i < num_configs - 1; i++) {
                for (int j = i + 1; j < num_configs; j++) {
                    double ss_sum = 0.0;
                    /* traversing for terms in each state */
                    for (int cf = 0; cf < cwv[nf].num_terms; cf++) {
                        for (int ci = 0; ci < cwv[ni].num_terms; ci++) {
                            /* match configs except positions i and j */
                            if (!configs_match_except_positions(cwv[nf].configs[cf], cwv[ni].configs[ci], num_configs, i, j)) {
                                continue;
                            }
                            char fi = cwv[nf].configs[cf][i];
                            char ii = cwv[ni].configs[ci][i];
                            char fj = cwv[nf].configs[cf][j];
                            char ij = cwv[ni].configs[ci][j];
                            if (config[i] == 'q' && config[j] == 'q') {
                                ss_sum += cwv[nf].coeffs[cf] * cwv[ni].coeffs[ci] * lambda_dot_lambda(fi, ii, fj, ij);
                            } else if (config[i] == 'q' && config[j] == 'Q') {
                                ss_sum -= cwv[nf].coeffs[cf] * cwv[ni].coeffs[ci] * lambda_dot_lambda(fi, ii, ij, fj);
                            } else if (config[i] == 'Q' && config[j] == 'q') {
                                ss_sum -= cwv[nf].coeffs[cf] * cwv[ni].coeffs[ci] * lambda_dot_lambda(ii, fi, fj, ij);
                            } else if (config[i] == 'Q' && config[j] == 'Q') {
                                ss_sum += cwv[nf].coeffs[cf] * cwv[ni].coeffs[ci] * lambda_dot_lambda(ii, fi, ij, fj);
                            }
                        }
                    }
                    result->value[nf * num_state + ni][(i * (2 * num_configs - i - 1)) / 2 + (j - i - 1)] = ss_sum;
                }
            }
        }
    }
}