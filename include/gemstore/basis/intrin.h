/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef GEMSTORE_BASIS_INTRIN
#define GEMSTORE_BASIS_INTRIN

typedef struct intrin_wfn
{
    /* number of terms in the wavefunction */
    int num_terms;
    /* number of configurations in each term */
    int num_configs;
    /* coefficient of each term, like 1, 1/2, and 1/sqrt(2) */
    double *coeffs;
    /* configuration of each term, like "10", "rgb", and "uds" */
    char **configs;
} intrin_wfn_t;

/* Initialize intrinsic wave function */
intrin_wfn_t intrin_wfn_init(int num_configs);

/* Product of two intrinsic wave functions */
void intrin_wfn_product(const intrin_wfn_t *wfnA, const intrin_wfn_t *wfnB, intrin_wfn_t *wfnC, double factor);

/* Trim intrinsic wave function */
void intrin_wfn_trim(intrin_wfn_t *wfn);

/* Push a term to the intrinsic wave function */
void intrin_wfn_push(intrin_wfn_t *wfn, double coeff, const char *config);

/* Print intrinsic wave function */
void intrin_wfn_print(const intrin_wfn_t *wfn);

/* Free intrinsic wave function */
void intrin_wfn_free(intrin_wfn_t *wfn);

typedef struct cg_tuple
{
    double ms1;             /* the third component of s1 */
    double ms2;             /* the third component of s2 */
    double cg;              /* Clebsch-Gordan coefficient */
} cg_tuple_t;

typedef struct cg_table
{
    int num;                /* number of entries */
    cg_tuple_t *tuples;     /* Clebsch-Gordan table */
} cg_table_t;

/* Call Clebsch-Gordan table */
cg_table_t CallCGTable(double s1, double s2, double st, double st3);

/* Free Clebsch-Gordan table */
void cg_table_free(cg_table_t *t);

#endif