/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gemstore/entry.h>
#include <gemstore/debug.h>
#include <gemstore/print.h>
#include <gemstore/parse.h>
#include <gemstore/model/compute.h>
#include <gemstore/param/argset.h>
#include <gemstore/param/fitting.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void entry_compute(const char* arg)
{
    argsInput_t input = {0};
    parse_input_file(arg, &input);
    print_input_parameters(&input);

    if (input.task == TASK_SPECTRA) {
        if (input.system == SYSTEM_MESON) {
            if (input.orbit != ORBIT_GEM && input.orbit != ORBIT_CRG && input.orbit != ORBIT_SHO) {
                fprintf(stderr, "SPECTRA meson compute supports: GEM, CRG, SHO basis\n");
                exit(1);
            }
            compute_spectra_meson(&input);
        }
        else if (input.system == SYSTEM_BARYON) {
            if (input.orbit != ORBIT_GEM) {
                fprintf(stderr, "SPECTRA baryon compute supports: GEM basis\n");
                exit(1);
            }
            compute_spectra_baryon(&input);
        }
        else {
            fprintf(stderr, "SPECTRA currently supports only MESON and BARYON systems\n");
            exit(1);
        }
    }
}

void entry_fitting(const char* arg)
{
    /* copy the content after ‘_’ */
    size_t prefix_len = strcspn(arg, "_");
    const char *suffix = arg + prefix_len;
    suffix++;

    /* copy the content before ‘_’ */
    char prefix[10];
    strncpy(prefix, arg, prefix_len);

    if (strcmp(prefix, "GIScreen") == 0) call_minuit2_GIScreen(suffix);
    else {fprintf(stderr, "Unknown fitting model: %s\n", arg); exit(1);}
}

void entry_debug(const char* arg)
{
    if (strcmp(arg, "su3_product") == 0) debug_su3_product();
    else if (strcmp(arg, "soc_operator") == 0) debug_soc_operator();
    else if (strcmp(arg, "casimir_operator") == 0) debug_casimir_operator();
    else if (strcmp(arg, "color_wfn") == 0) debug_color_wfn();
    else if (strcmp(arg, "spin_wfn") == 0) debug_spin_wfn();
    else if (strcmp(arg, "isospin_wfn") == 0) debug_isospin_wfn();
    else if (strcmp(arg, "orbit_wfn") == 0) debug_orbit_wfn();
    else if (strcmp(arg, "eigen_system") == 0) debug_eigen_system();
    else if (strcmp(arg, "eigen_system_complex") == 0) debug_eigen_system_complex();
    else {fprintf(stderr, "Unknown debug unit: %s\n", arg); exit(1);}
}
