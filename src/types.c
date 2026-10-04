/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gemstore/types.h>

const char *task_type_str[] = {
    "SPECTRA",
    "DECAY3P0",
    "COUPLCHN",
    "SCATTER",
    "FITTING"
};

const char *model_type_str[] = {
    "GISTRING",
    "GISCREEN",
    "NRSTRING",
    "NRSCREEN"
};

const char *system_type_str[] = {
    "MESON",
    "BARYON",
    "MOLECULE"
};

const char *orbit_type_str[] = {
    "SHO",
    "GEM",
    "NONE"
};

const char *fitting_type_str[] = {
    "GISCREEN_MESON",
    "GISCREEN_BBBAR",
    "GISCREEN_BCBAR",
    "GISCREEN_BSBAR",
    "GISCREEN_CCBAR",
    "GISCREEN_CSBAR"
};

const char *param_type_str[] = {
    "GISTRING_MESON",
    "GISTRING_BARYON",
    "GISCREEN_MESON",
    "GISCREEN_BBBAR",
    "GISCREEN_BCBAR",
    "GISCREEN_BSBAR",
    "GISCREEN_CCBAR",
    "GISCREEN_CSBAR",
    "GISTRING_CUSTOM",
    "GISCREEN_CUSTOM",
    "NRSTRING_MESON",
    "NRSCREEN_MESON",
    "NRSTRING_CUSTOM",
    "NRSCREEN_CUSTOM"
};