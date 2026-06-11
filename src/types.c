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
    "SCATTER"
};

const char *orbit_type_str[] = {
    "SHO",
    "GEM",
    "CRG"
};

const char *model_type_str[] = {
    "GISTRING",
    "GISCREEN"
};

const char *param_type_str[] = {
    "GISTRING_MESON",
    "GISCREEN_MESON",
    "GISCREEN_BBBAR",
    "GISCREEN_BCBAR",
    "GISCREEN_BSBAR",
    "GISCREEN_CCBAR",
    "GISCREEN_CSBAR",
    "GISTRING_CUSTOM",
    "GISCREEN_CUSTOM"
};

const char *system_type_str[] = {
    "MESON",
    "BARYON",
    "MOLECULE"
};