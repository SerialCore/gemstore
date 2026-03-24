/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gemstore/entry.h>
#include <gemstore/debug.h>
#include <gemstore/fileio.h>
#include <gemstore/model/compute.h>
#include <gemstore/param/argset.h>
#include <gemstore/param/fitting.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <ctype.h>

static void trim_line(char *s)
{
    if (!s || !*s) return;
    
    /* Fully handle # and // comments - truncate at first comment */
    char *comment = strchr(s, '#');
    if (!comment) comment = strstr(s, "//");
    if (comment) *comment = '\0';
    
    /* Trim all leading and trailing whitespace */
    char *end = s + strlen(s) - 1;
    while (end >= s && isspace((unsigned char)*end)) *end-- = '\0';
    
    char *start = s;
    while (isspace((unsigned char)*start)) start++;
    if (start != s) memmove(s, start, strlen(start) + 1);
    
    /* Extra pass to remove any remaining spaces around key (for strtok safety) */
    for (char *p = s; *p; p++) {
        if (isspace((unsigned char)*p)) {
            memmove(p, p+1, strlen(p));
            p--;
        }
    }
}

static int parse_line(char *line, input_section_t sec, argsInput_t *input)
{
    trim_line(line);
    if (strlen(line) == 0) return 1;
    char key[64] = {0}, val[128] = {0};
    char *token = strtok(line, "=");
    if (token) {
        strncpy(key, token, 63);
        trim_line(key);
        token = strtok(NULL, "=");
        if (token) {
            strncpy(val, token, 127);
            trim_line(val);
        }
    }

    if (sec == SECTION_GLOBAL) {
        if (strcmp(key, "project") == 0) {
            strncpy(input->project, val, 255);
        } else if (strcmp(key, "task") == 0) {
            if (strcmp(val, "SPECTRA") == 0) input->task = TASK_SPECTRA;
            else if (strcmp(val, "RADIUS") == 0) input->task = TASK_RADIUS;
            else if (strcmp(val, "DECAY3P0") == 0) input->task = TASK_DECAY3P0;
            else if (strcmp(val, "COUPLCHN") == 0) input->task = TASK_COUPLCHN;
            else if (strcmp(val, "SCATTER") == 0) input->task = TASK_SCATTER;
            else {fprintf(stderr, "Unknown task: %s\n", val); exit(1);}
        }
    } else if (sec == SECTION_SYSTEM) {
        if (strcmp(key, "model") == 0) {
            if (strcmp(val, "GI_STRING") == 0) input->model = MODEL_GI_STRING;
            else if (strcmp(val, "GI_SCREEN") == 0) input->model = MODEL_GI_SCREEN;
            else if (strcmp(val, "GI_QUADRA") == 0) input->model = MODEL_GI_QUADRA;
            else {fprintf(stderr, "Unknown model: %s\n", val); exit(1);}
        } else if (strcmp(key, "system") == 0) {
            if (strcmp(val, "MESON") == 0) input->system = SYSTEM_MESON;
            else if (strcmp(val, "BARYON") == 0) input->system = SYSTEM_BARYON;
            else if (strcmp(val, "MOLECULE") == 0) input->system = SYSTEM_MOLECULE;
            else {fprintf(stderr, "Unknown system: %s\n", val); exit(1);}
        }
    } else if (sec == SECTION_PARAMS) {
        if (strcmp(key, "params") == 0) {
            if (strcmp(val, "GIString_meson") == 0) input->params = argsGIString_meson;
            else if (strcmp(val, "GIScreen_meson") == 0) input->params = argsGIScreen_meson;
            else if (strcmp(val, "GIScreen_bbbar") == 0) input->params = argsGIScreen_bbbar;
            else if (strcmp(val, "GIScreen_ccbar") == 0) input->params = argsGIScreen_ccbar;
            else if (strcmp(val, "GIScreen_light") == 0) input->params = argsGIScreen_light;
            else if (strcmp(val, "GIQuadra_meson") == 0) input->params = argsGIQuadra_meson;
            else if (strcmp(val, "GIQuadra_bbbar") == 0) input->params = argsGIQuadra_bbbar;
            else if (strcmp(val, "GIQuadra_ccbar") == 0) input->params = argsGIQuadra_ccbar;
            else if (strcmp(val, "GIQuadra_light") == 0) input->params = argsGIQuadra_light;
            else {fprintf(stderr, "Unknown paramset: %s\n", val); exit(1);}
        }
        else if (strcmp(key, "mn") == 0) input->params.mn = atof(val);
        else if (strcmp(key, "ms") == 0) input->params.ms = atof(val);
        else if (strcmp(key, "mc") == 0) input->params.mc = atof(val);
        else if (strcmp(key, "mb") == 0) input->params.mb = atof(val);
        else if (strcmp(key, "b1") == 0) input->params.b1 = atof(val);
        else if (strcmp(key, "b2") == 0) input->params.b2 = atof(val);
        else if (strcmp(key, "mu") == 0) input->params.mu = atof(val);
        else if (strcmp(key, "c") == 0) input->params.c = atof(val);
        else if (strcmp(key, "sigma_0") == 0) input->params.sigma_0 = atof(val);
        else if (strcmp(key, "s") == 0) input->params.s = atof(val);
        else if (strcmp(key, "epsilon_cont") == 0) input->params.epsilon_cont = atof(val);
        else if (strcmp(key, "epsilon_sov") == 0) input->params.epsilon_sov = atof(val);
        else if (strcmp(key, "epsilon_sos") == 0) input->params.epsilon_sos = atof(val);
        else if (strcmp(key, "epsilon_tens") == 0) input->params.epsilon_tens = atof(val);
        else {fprintf(stderr, "Unknown param: %s\n", key); exit(1);}
    } else if (sec == SECTION_QUANTUM) {
        if (strcmp(key, "f1") == 0) input->f1 = atoi(val);
        else if (strcmp(key, "f2") == 0) input->f2 = atoi(val);
        else if (strcmp(key, "f3") == 0) input->f3 = atoi(val);
        else if (strcmp(key, "f4") == 0) input->f4 = atoi(val);
        else if (strcmp(key, "S") == 0) input->S = atof(val);
        else if (strcmp(key, "L") == 0) input->L = atof(val);
        else if (strcmp(key, "jl") == 0) input->jl = atof(val);
        else if (strcmp(key, "J") == 0) input->J = atof(val);
        else {fprintf(stderr, "Unknown quantum number: %s\n", key); exit(1);}
    } else if (sec == SECTION_GAUSS) {
        if (strcmp(key, "nmax") == 0) input->nmax = atoi(val);
        else if (strcmp(key, "rmax") == 0) input->rmax = atof(val);
        else if (strcmp(key, "rmin") == 0) input->rmin = atof(val);
        else {fprintf(stderr, "Unknown gauss parameter: %s\n", key); exit(1);}
    }

    return 1;
}

static int parse_input_file(const char *filename, argsInput_t *input)
{;
    FILE *f = fopen(filename, "r");
    if (!f) {
        fprintf(stderr, "Cannot open %s\n", filename);
        exit(1);
    }
    char line[512];

    input_section_t curr_sec = SECTION_NONE;
    while (fgets(line, sizeof(line), f)) {
        if (strstr(line, "&GLOBAL")) curr_sec = SECTION_GLOBAL;
        else if (strstr(line, "&SYSTEM")) curr_sec = SECTION_SYSTEM;
        else if (strstr(line, "&PARAMS")) curr_sec = SECTION_PARAMS;
        else if (strstr(line, "&QUANTUM")) curr_sec = SECTION_QUANTUM;
        else if (strstr(line, "&GAUSS")) curr_sec = SECTION_GAUSS;
        else if (strstr(line, "&END") || strstr(line, "/")) curr_sec = SECTION_NONE;
        else if (curr_sec != SECTION_NONE) parse_line(line, curr_sec, input);
        else {fprintf(stderr, "Unknown input section: %s\n", line); exit(1);}
    }
    fclose(f);
    
    return 1;
}

void entry_compute(const char* arg)
{
    argsInput_t input = {0};
    parse_input_file(arg, &input);

    if (input.task == TASK_SPECTRA) compute_spectra_meson(&input);
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
    else if (strcmp(prefix, "GIQuadra") == 0) call_minuit2_GIQuadra(suffix);
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
}

void entry_print(const char* arg)
{
    //
}