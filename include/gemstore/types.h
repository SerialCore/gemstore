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
    TASK_FITTING
} task_type_t;

extern const char *task_type_str[];

/* used for determining potential and dispatching fitting task */
typedef enum model_type {
    MODEL_GISTRING,
    MODEL_GISCREEN,
    MODEL_NRSTRING,
    MODEL_NRSCREEN,
    MODEL_BAGMIT,
    MODEL_OBESTAT
} model_type_t;

extern const char *model_type_str[];

/* used for determining potential and dispatching computing task */
typedef enum system_type {
    SYSTEM_BAG,
    SYSTEM_MESON,
    SYSTEM_BARYON,
    SYSTEM_MOLECUL2,
    SYSTEM_MOLECUL3
} system_type_t;

extern const char *system_type_str[];

/* used for determining basis type */
typedef enum basis_type {
    BASIS_SHO,
    BASIS_GEM,
    BASIS_CMI,
    BASIS_NONE
} basis_type_t;

extern const char *basis_type_str[];

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
    PARAM_NRSCREEN_CUSTOM,
    PARAM_BAGMIT_DEFAULT,
    PARAM_BAGMIT_CUSTOM,
} param_type_t;

extern const char *param_type_str[];

/* used for selecting a CMI color channel */
typedef enum cmicolor_type {
    CMI_COLOR_MESON,
    CMI_COLOR_BARYON,
    CMI_COLOR_TETRA11,
    CMI_COLOR_TETRA33,
    CMI_COLOR_TETRA66,
    CMI_COLOR_TETRA88,
    CMI_COLOR_PENTA688,
    CMI_COLOR_PENTA388,
    CMI_COLOR_PENTA11,
    CMI_COLOR_HEXA1010,
    CMI_COLOR_HEXA6868,
    CMI_COLOR_HEXA6838,
    CMI_COLOR_HEXA3868,
    CMI_COLOR_HEXA3838,
    CMI_COLOR_HEXA11
} cmicolor_type_t;

extern const char *cmicolor_type_str[];

#endif
