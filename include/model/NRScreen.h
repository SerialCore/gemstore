/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef NRSCREEN_H
#define NRSCREEN_H

typedef struct argsNRScreen {
    double mn;              /* mass of n */
    double ms;              /* mass of s */
    double mc;              /* mass of c */
    double mb;              /* mass of b */
    double mt;              /* mass of t */

    double alpha_s;         /* strong coupling constant */
    double b;               /* string tension */
    double mu;              /* screen length */
    double c;               /* constant */
    double sigma;           /* short-range */
} argsNRScreen_t;

argsNRScreen_t argsNRScreen_meson = {
    0.606, 0.780, 1.984, 5.368, 172.57, 
    0.3930, 0.2312, 0.069, -1.1711, 1.842
};

typedef double (*Vr)(double);

typedef double (*Vp)(double);

double NRScreen_T(double p, double m1, double m2);

double NRScreen_Vconf(double r, double alpha_s, double b, double mu, double c);

double NRScreen_Vcont(double r, double alpha_s, double sigma, double m1, double m2, double sds);

double NRScreen_Vsocm(double r, double alpha_s, double m1, double m2, double lds1, double lds2);

double NRScreen_Vsotp(double r, double alpha_s, double b, double mu, double m1, double m2, double lds1, double lds2);

double NRScreen_Vtens(double r, double alpha_s, double m1, double m2, double tens);

#endif