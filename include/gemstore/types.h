/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef GEMSTORE_TYPES
#define GEMSTORE_TYPES

typedef enum orbit_type {
    ORBIT_GEM,
    ORBIT_CSM,
    ORBIT_CRG,
    ORBIT_SHO
} orbit_type_t;

/* used for determining potential and dispatching computing task */
typedef enum system_type {
    SYSTEM_MESON,
    SYSTEM_BARYON,
    SYSTEM_MOLECULE
} system_type_t;

/* used for determining potential and dispatching fitting task */
typedef enum model_type {
    MODEL_GISTRING,
    MODEL_GISCREEN
} model_type_t;

/* used for determining parameter type */
typedef enum param_type {
    PARAM_GISTRING_MESON,
    PARAM_GISCREEN_MESON,
    PARAM_GISCREEN_BBBAR,
    PARAM_GISCREEN_CCBAR,
    PARAM_GISCREEN_LIGHT
} param_type_t;

/* used for dispatching task */
typedef enum task_type {
    TASK_SPECTRA,
    TASK_DECAY3P0,
    TASK_COUPLCHN,
    TASK_SCATTER
} task_type_t;

#endif
