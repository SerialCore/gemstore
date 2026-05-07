/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef GEMSTORE_PARSE
#define GEMSTORE_PARSE

#include <gemstore/math/matrix.h>
#include <gemstore/param/argset.h>

/* parse input file */
void parse_input_file(const char *filename, argsInput_t *input);

/* parse GISTRING parameters */
void parse_param_GISTRING(const char *filename, argsGIModel_t *args_model);

/* parse GISCREEN parameters */
void parse_param_GISCREEN(const char *filename, argsGIModel_t *args_model);

/* parse meson state file (mass, rms_radius, eigenvector) */
void parse_meson_state(const char *filename, array_t *mass, array_t *radius, matrix_t *eigenvector);

#endif
