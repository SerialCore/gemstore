/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef GEMSTORE_MODEL_MODEL
#define GEMSTORE_MODEL_MODEL

typedef enum model_enum {
    NRScreen,
    GIString,
    GIScreen,
    GIQuadraScreen
} model_enum_t;

typedef double (*potential_r)(double r, void* args);
typedef double (*potential_p)(double p, void* args);

#endif