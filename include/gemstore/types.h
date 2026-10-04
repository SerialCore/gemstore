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
    TASK_SCATTER,
    TASK_FITTING
} task_type_t;

extern const char *task_type_str[];

/* used for determining potential and dispatching fitting task */
typedef enum model_type {
    MODEL_GISTRING,
    MODEL_GISCREEN,
    MODEL_NRSTRING,
    MODEL_NRSCREEN
} model_type_t;

extern const char *model_type_str[];

/* used for determining potential and dispatching computing task */
typedef enum system_type {
    SYSTEM_MESON,
    SYSTEM_BARYON,
    SYSTEM_MOLECULE
} system_type_t;

extern const char *system_type_str[];

/* used for determining orbit type */
typedef enum orbit_type {
    ORBIT_SHO,
    ORBIT_GEM,
    ORBIT_NONE
} orbit_type_t;

extern const char *orbit_type_str[];

/* used for determining fitting task */
typedef enum fitting_type {
    FITTING_GISCREEN_MESON,
    FITTING_GISCREEN_BBBAR,
    FITTING_GISCREEN_BCBAR,
    FITTING_GISCREEN_BSBAR,
    FITTING_GISCREEN_CCBAR,
    FITTING_GISCREEN_CSBAR
} fitting_type_t;

extern const char *fitting_type_str[];

/* used for determining parameter type */
typedef enum param_type {
    PARAM_GISTRING_MESON,
    PARAM_GISTRING_BARYON,
    PARAM_GISCREEN_MESON,
    PARAM_GISCREEN_BBBAR,
    PARAM_GISCREEN_BCBAR,
    PARAM_GISCREEN_BSBAR,
    PARAM_GISCREEN_CCBAR,
    PARAM_GISCREEN_CSBAR,
    PARAM_GISTRING_CUSTOM,
    PARAM_GISCREEN_CUSTOM,
    PARAM_NRSTRING_MESON,
    PARAM_NRSCREEN_MESON,
    PARAM_NRSTRING_CUSTOM,
    PARAM_NRSCREEN_CUSTOM
} param_type_t;

extern const char *param_type_str[];

#endif
