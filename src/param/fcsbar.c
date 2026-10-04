/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gemstore/param/fitting.h>

static const fit_state_t k_states[] = {
    // Ds (s c-bar) 9
    {"Ds(1^1S0)", 2, 3, 1, 0, 0, 0, 1968.4, 1},     /* Ds(1^1S0) */
    {"Ds*(1^3S1)", 2, 3, 1, 1, 0, 1, 2112.2, 1},    /* Ds*(1^3S1) */
    {"Ds0*(2317)", 2, 3, 1, 1, 1, 0, 2317.8, 1},    /* Ds0*(2317), candidate 1^3P0 */
    {"Ds1(2460)", 2, 3, 1, 1, 1, 1, 2459.5, 1},     /* Ds1(2460), mixed 1^1P1/1^3P1, stored here as dominant ^3P1-like state */
    {"Ds1(2536)", 2, 3, 1, 0, 1, 1, 2535.1, 1},     /* Ds1(2536), 1^1P1 */
    {"Ds2*(2573)", 2, 3, 1, 1, 1, 2, 2569.1, 1},    /* Ds2*(2573), 1^3P2 */
    {"Ds0(2590)", 2, 3, 2, 0, 0, 0, 2591.0, 9},     /* Ds0(2590), 2^1S0 (LHCb, JP=0-) */
    {"Ds3*(2860)", 2, 3, 1, 1, 2, 3, 2860.5, 7},    /* Ds3*(2860), 1^3D3 */
};

static const fit_param_t k_params[] = {
    {"mn", 0.4560806209112, 0.01, 0.1, 0.5, 1},
    {"ms", 0.6173440068792, 0.01, 0.3, 0.7, 1},
    {"mc", 1.805387067165, 0.01, 1.5, 2.0, 1},
    {"mb", 5.151269542382, 0.01, 4.5, 5.5, 1},
    {"b", 0.2522010221331, 0.01, 0.1, 0.3, 0},
    {"mu", 0.16255, 0.01, 0.1, 0.2, 1},
    {"c", -0.6482214863381, 0.01, -2.0, 0.0, 1},
    {"sigma_0", 1.770545357386, 0.01, 1.0, 3.0, 1},
    {"s", 1.146340880694, 0.01, 1.0, 3.0, 1},
    {"epsilon_cont", -0.3199525941084, 0.01, -0.5, 0.0, 0},
    {"epsilon_sov", -0.3428467195716, 0.01, -1.0, 1.0, 0},
    {"epsilon_sos", 0.7905135472165, 0.01, -1.0, 1.0, 0},
    {"epsilon_tens", -0.5000469699086, 0.01, -1.0, 1.0, 0},
};

const fit_target_t fit_target_giscreen_csbar = {
    MODEL_GISCREEN,
    k_states,
    (int)(sizeof k_states / sizeof k_states[0]),
    k_params,
    (int)(sizeof k_params / sizeof k_params[0])
};
