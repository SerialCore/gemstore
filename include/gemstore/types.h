/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef GEMSTORE_TYPES
#define GEMSTORE_TYPES

/* used for determining potential and dispatching computing task */
typedef enum system_type {
    SYSTEM_MESON,
    SYSTEM_BARYON,
    SYSTEM_MOLECULE
} system_type_t;

/* used for determining potential and dispatching fitting task */
typedef enum model_type {
    MODEL_GI_STRING,
    MODEL_GI_SCREEN,
    MODEL_GI_QUADRA
} model_type_t;

/* used for dispatching task */
typedef enum task_type {
    TASK_SPECTRA,
    TASK_RADIUS,
    TASK_DECAY3P0,
    TASK_COUPLCHN,
    TASK_SCATTER,
    TASK_FITTING,
    TASK_DEBUG,
    TASK_PRINT
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