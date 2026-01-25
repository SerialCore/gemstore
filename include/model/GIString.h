/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef GISTRING_H
#define GISTRING_H

typedef struct argsGIString {
    double md;              /* mass of d */
    double mu;              /* mass of u */
    double ms;              /* mass of s */
    double mc;              /* mass of c */
    double mb;              /* mass of b */
    double mt;              /* mass of t */

    double b;               /* string tension */
    double c;               /* constant */
    double sigma_0;         /* GI smearing parameter for sigma_ij */
    double s;               /* GI smearing parameter for sigma_ij */
    
    double epsilon_Coul;    /* GI smearing parameter for Coulumb */
    double epsilon_cont;    /* GI smearing parameter for contact */
    double epsilon_sonu;    /* GI smearing parameter for spin-orbit */
    double epsilon_sos;     /* GI smearing parameter for Thomas */
    double epsilon_tens;    /* GI smearing parameter for tensor */
} argsGIString_t;

argsGIString_t argsGIString_meson = {
    0.220, 0.220, 0.419, 1.628, 4.977, 172.57, 
    0.18, -0.253, 1.8, 1.55, 
    0.0, -0.168, -0.035, 0.055, 0.025
};

#endif