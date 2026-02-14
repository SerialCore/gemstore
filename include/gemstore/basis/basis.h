/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef GEMSTORE_BASIS_BASIS
#define GEMSTORE_BASIS_BASIS

typedef enum system_type {
    SYSTEM_MESON,
    SYSTEM_BARYON
} system_type_t;

typedef struct basis_base {
    system_type_t type;
    double J;           /* total momentum J */
    double P;           /* parity P */
} basis_base_t;

typedef struct basis_meson {
    basis_base_t base;
    double m1;          /* mass of quark 1 */
    double m2;          /* mass of quark 2 */
    double m3;          /* mass of quark 3 */
    double s1;          /* spin of quark 1 */
    double s2;          /* spin of quark 2 */
    double s3;          /* spin of quark 3 */
    double t1;          /* isospin of quark 1 */
    double t2;          /* isospin of quark 2 */
    double t3;          /* isospin of quark 3 */
    double tij;         /* isospin of quark ij */
    double T;           /* total isospin */
    int c;              /* index of pair: 1=12, 2=13, 3=23 */
    int lrho;           /* orbital momentum of rho */
    int llam;           /* orbital momentum of lambda */
    int L;              /* total orbital momentum */
    double sij;         /* spin of quark ij */
    double jl;          /* jj coupling */
    double J;           /* total J */
    int nrho;           /* Gaussian quantum number of rho */
    int nlam;           /* Gaussian quantum number of lambda */
    double nurho;       /* Gaussian parameter of rho */
    double nulam;       /* Gaussian parameter of lambda */
} basis_meson_t;

typedef struct basis_baryon {
    basis_base_t base;
    int map1;           /* index map of Jacobi coordinate */
    int map2;           /* index map of Jacobi coordinate */
    double coe;         /* coefficient of basis */
    double m1;          /* mass of quark 1 */
    double m2;          /* mass of quark 2 */
    double m3;          /* mass of quark 3 */
    double s1;          /* spin of quark 1 */
    double s2;          /* spin of quark 2 */
    double s3;          /* spin of quark 3 */
    double t1;          /* isospin of quark 1 */
    double t2;          /* isospin of quark 2 */
    double t3;          /* isospin of quark 3 */
    double tij;         /* isospin of quark ij */
    double T;           /* total isospin */
    int c;              /* index of pair: 1=12, 2=13, 3=23 */
    int lrho;           /* orbital momentum of rho */
    int llam;           /* orbital momentum of lambda */
    int L;              /* total orbital momentum */
    double sij;         /* spin of quark ij */
    double jl;          /* jj coupling */
    double J;           /* total J */
    int nrho;           /* Gaussian quantum number of rho */
    int nlam;           /* Gaussian quantum number of lambda */
    double nurho;       /* Gaussian parameter of rho */
    double nulam;       /* Gaussian parameter of lambda */
} basis_baryon_t;

typedef struct {
    basis_base_t **basis;
    int len_list;
    int *len_part;
} basis_list;

/*
void basis_list_init(basis_list *qnlist);

void basis_list_logs(basis_list qnlist);

void basis_list_push(basis_list *qnlist, int newqnQ, int map1, int map2, double coe, double m1, double m2, double m3, double s1, double s2, double s3, double t1, double t2, double t3, double tij, double T, int c, int lrho, int llam, int L, double sij, double jl, double J, int nrho, int nlam, double nurho, double nulam);

void basis_list_push_full(basis_list *qnlist_spfy, basis_list *qnlist_full, double rmin, double rmax, int nmax);
*/

#endif