/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef GEMSTORE_TYPES
#define GEMSTORE_TYPES

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

#endif