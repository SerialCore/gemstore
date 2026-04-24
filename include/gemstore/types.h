/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef GEMSTORE_TYPES
#define GEMSTORE_TYPES

/* used for dispatching task */
typedef enum task_type {
    TASK_SPECTRA,
    TASK_DECAY3P0,
    TASK_COUPLCHN,
    TASK_SCATTER
} task_type_t;

extern const char *task_type_str[];

/* used for determining orbit type */
typedef enum orbit_type {
    ORBIT_SHO,
    ORBIT_GEM,
    ORBIT_CRG,
    ORBIT_CSM
} orbit_type_t;

extern const char *orbit_type_str[];

/* used for determining potential and dispatching fitting task */
typedef enum model_type {
    MODEL_GISTRING,
    MODEL_GISCREEN
} model_type_t;

extern const char *model_type_str[];

/* used for determining parameter type */
typedef enum param_type {
    PARAM_GISTRING_MESON,
    PARAM_GISCREEN_MESON,
    PARAM_GISCREEN_BBBAR,
    PARAM_GISCREEN_CCBAR,
    PARAM_GISTRING_CUSTOM,
    PARAM_GISCREEN_CUSTOM
} param_type_t;

extern const char *param_type_str[];

/* used for determining potential and dispatching computing task */
typedef enum system_type {
    SYSTEM_MESON,
    SYSTEM_BARYON,
    SYSTEM_MOLECULE
} system_type_t;

extern const char *system_type_str[];

#endif
