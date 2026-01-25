/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef GIQUADRASCREEN_H
#define GIQUADRASCREEN_H

typedef struct argsGIQuadraScreen {
    double mn;              /* mass of n */
    double ms;              /* mass of s */
    double mc;              /* mass of c */
    double mb;              /* mass of b */
    double mt;              /* mass of t */

    double b1;              /* string tension */
    double b2;              /* surface tension */
    double mu;              /* screen length */
    double c;               /* constant */
    double sigma_0;         /* GI smearing parameter for sigma_ij */
    double s;               /* GI smearing parameter for sigma_ij */
    
    double epsilon_Coul;    /* GI smearing parameter for Coulumb */
    double epsilon_cont;    /* GI smearing parameter for contact */
    double epsilon_sonu;    /* GI smearing parameter for spin-orbit */
    double epsilon_sos;     /* GI smearing parameter for Thomas */
    double epsilon_tens;    /* GI smearing parameter for tensor */
} argsGIQuadraScreen_t;

argsGIQuadraScreen_t argsGIQuadraScreen_meson = {
    0.220, 0.419, 1.628, 4.977, 172.57, 
    0.16, 0.02, 0.15, -0.253, 1.8, 1.55, 
    0.0, -0.168, -0.035, 0.055, 0.025
};

#endif