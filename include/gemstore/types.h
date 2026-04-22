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
    MODEL_GI_STRING,
    MODEL_GI_SCREEN
} model_type_t;

/* used for dispatching task */
typedef enum task_type {
    TASK_SPECTRA,
    TASK_DECAY3P0,
    TASK_COUPLCHN,
    TASK_SCATTER
} task_type_t;

/* used for parsing input */
typedef enum input_section {
    SECTION_NONE,
    SECTION_GLOBAL,
    SECTION_SYSTEM,
    SECTION_PARAMS,
    SECTION_QUANTUM,
    SECTION_GAUSS,
} input_section_t;

#endif