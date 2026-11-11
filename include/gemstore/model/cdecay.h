/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef GEMSTORE_MODEL_CDECAY
#define GEMSTORE_MODEL_CDECAY

#include <gemstore/input.h>

#define QPC_WAVE_MAX 32

/* One partial wave. */
typedef struct qpc_wave {
    int L;       /* outgoing orbital angular momentum */
    int S;       /* total spin of the two final mesons */
    double M;    /* partial-wave amplitude */
} qpc_wave_t;

/* One channel. Momenta and energies are in GeV. The width is in MeV. */
typedef struct qpc_result {
    int open;                      /* 0 below threshold; 1 when the channel is open */
    int nwave;                     /* number of entries stored in wave[] */
    double p;                      /* center-of-mass momentum */
    double EB;                     /* energy of daughter B */
    double EC;                     /* energy of daughter C */
    double width;                  /* 2 pi (EB EC / MA) p F sum |M|^2, in MeV */
    qpc_wave_t wave[QPC_WAVE_MAX]; /* partial waves; the first nwave entries are used */
} qpc_result_t;

/* Open-flavor QPC width of one meson into two mesons, one SHO each.
 * Calls qpc_amplitude for each allowed partial wave, sum up |M|^2.
 * Width = 2 pi (EB EC / MA) p F sum |M|^2, converted from GeV to MeV.
 */
void decay_qpc_meson(const argsInput_t *input, qpc_result_t *result);

#endif
