/*
 * Copyright (C) 2025 Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef BASIS_H
#define BASIS_H

typedef struct {
    int map1;
    int map2;
    double coe;
    double m1;
    double m2;
    double m3;
    double s1;
    double s2;
    double s3;
    double t1;
    double t2;
    double t3;
    double tij;
    double T;
    int c;
    int lrho;
    int llam;
    int L;
    double sij;
    double jl;
    double J;
    int nrho;
    int nlam;
    double nurho;
    double nulam;
} basis_qnum;

typedef struct {
    basis_qnum **qnum;
    int *len_part;
    int len_list;
} basis_list;

void basis_list_init(basis_list *qnlist);

void basis_list_logs(basis_list qnlist);

void basis_list_push(basis_list *qnlist, int newqnQ, int map1, int map2, double coe, double m1, double m2, double m3, double s1, double s2, double s3, double t1, double t2, double t3, double tij, double T, int c, int lrho, int llam, int L, double sij, double jl, double J, int nrho, int nlam, double nurho, double nulam);

double getnu(double xmin, double xmax, int N, int n);

void basis_list_push_full(basis_list *qnlist_spfy, basis_list *qnlist_full, double rmin, double rmax, int nmax);

#endif