/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef GEMSTORE_TYPES
#define GEMSTORE_TYPES

#include <complex.h>

typedef enum orbit_type {
    ORBIT_GEM,
    ORBIT_SHO
} orbit_type_t;

typedef enum system_type {
    SYSTEM_MESON,
    SYSTEM_BARYON,
    SYSTEM_MOLECULE
} system_type_t;

typedef enum model_type {
    MODEL_GI_STRING,
    MODEL_GI_SCREEN,
    MODEL_GI_QUADRA
} model_type_t;

typedef struct argsModel {
    model_type_t model;     /* model type */
    double mn;              /* mass of n quark */
    double ms;              /* mass of s quark */
    double mc;              /* mass of c quark */
    double mb;              /* mass of b quark */
    double mt;              /* mass of t quark */
} argsModel_t;

#endif