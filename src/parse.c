/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gemstore/parse.h>
#include <gemstore/types.h>
#include <gemstore/param/argset.h>

#include "cJSON.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

typedef struct input_flags {
    bool project;
    bool task;
    bool model;
    bool system;
    bool basis;
    bool f1;
    bool f2;
    bool f3;
    bool f4;
    bool S;
    bool L;
    bool jl;
    bool J;
    bool nmax;
    bool rmax;
    bool rmin;
    bool omega;
    bool beta;
    bool theta;
} input_flags_t;

static void check_flag(bool condition, const char *message)
{
    if (!condition) {
        fprintf(stderr, "%s\n", message);
        exit(1);
    }
}

static char *read_input_file(const char *filename)
{
    FILE *f = fopen(filename, "rb");
    char *buffer;
    long size;

    if (!f) {
        fprintf(stderr, "Cannot open %s\n", filename);
        exit(1);
    }

    if (fseek(f, 0, SEEK_END) != 0) {
        fclose(f);
        fprintf(stderr, "Cannot seek %s\n", filename);
        exit(1);
    }

    size = ftell(f);
    if (size < 0) {
        fclose(f);
        fprintf(stderr, "Cannot determine size of %s\n", filename);
        exit(1);
    }

    if (fseek(f, 0, SEEK_SET) != 0) {
        fclose(f);
        fprintf(stderr, "Cannot rewind %s\n", filename);
        exit(1);
    }

    buffer = (char *)malloc((size_t)size + 1);
    if (!buffer) {
        fclose(f);
        fprintf(stderr, "Out of memory while reading %s\n", filename);
        exit(1);
    }

    if (fread(buffer, 1, (size_t)size, f) != (size_t)size) {
        free(buffer);
        fclose(f);
        fprintf(stderr, "Cannot read %s\n", filename);
        exit(1);
    }

    buffer[size] = '\0';
    fclose(f);
    return buffer;
}

static cJSON *read_object_item(const cJSON *object, const char *key)
{
    cJSON *item = cJSON_GetObjectItemCaseSensitive(object, key);
    if (!cJSON_IsObject(item)) {
        fprintf(stderr, "Missing or invalid object: %s\n", key);
        exit(1);
    }
    return item;
}

static cJSON *read_string_item(const cJSON *object, const char *key)
{
    cJSON *item = cJSON_GetObjectItemCaseSensitive(object, key);
    if (!cJSON_IsString(item) || item->valuestring == NULL) {
        fprintf(stderr, "Missing or invalid string: %s\n", key);
        exit(1);
    }
    return item;
}

static cJSON *read_number_item(const cJSON *object, const char *key)
{
    cJSON *item = cJSON_GetObjectItemCaseSensitive(object, key);
    if (!cJSON_IsNumber(item)) {
        fprintf(stderr, "Missing or invalid number: %s\n", key);
        exit(1);
    }
    return item;
}

static void parse_task_string(const char *value, argsInput_t *input, input_flags_t *flags)
{
    if (strcmp(value, "SPECTRA") == 0) input->task = TASK_SPECTRA;
    else if (strcmp(value, "DECAY3P0") == 0) input->task = TASK_DECAY3P0;
    else if (strcmp(value, "COUPLCHN") == 0) input->task = TASK_COUPLCHN;
    else if (strcmp(value, "SCATTER") == 0) input->task = TASK_SCATTER;
    else {
        fprintf(stderr, "Unknown task: %s\n", value);
        exit(1);
    }

    flags->task = true;
}

static void parse_system_section(const cJSON *root, argsInput_t *input, input_flags_t *flags)
{
    cJSON *system_json = read_object_item(root, "system");
    const char *type = read_string_item(system_json, "type")->valuestring;

    if (strcmp(type, "MESON") == 0) input->system = SYSTEM_MESON;
    else {
        fprintf(stderr, "Unsupported system type: %s\n", type);
        exit(1);
    }
    flags->system = true;

    input->f1 = read_number_item(system_json, "f1")->valueint;
    flags->f1 = true;
    input->f2 = read_number_item(system_json, "f2")->valueint;
    flags->f2 = true;
    input->S = read_number_item(system_json, "S")->valuedouble;
    flags->S = true;
    input->L = read_number_item(system_json, "L")->valuedouble;
    flags->L = true;
    input->J = read_number_item(system_json, "J")->valuedouble;
    flags->J = true;
}

static void parse_model_section(const cJSON *root, argsInput_t *input, input_flags_t *flags)
{
    cJSON *model_json = read_object_item(root, "model");
    const char *type = read_string_item(model_json, "type")->valuestring;
    cJSON *item;

    if (strcmp(type, "GISTRING") == 0) input->model = MODEL_GISTRING;
    else if (strcmp(type, "GISCREEN") == 0) input->model = MODEL_GISCREEN;
    else {
        fprintf(stderr, "Unknown model type: %s\n", type);
        exit(1);
    }
    flags->model = true;

    item = cJSON_GetObjectItemCaseSensitive(model_json, "param");
    if (cJSON_IsString(item) && item->valuestring != NULL) {
        if (strcmp(item->valuestring, "GISTRING_MESON") == 0) input->params = argsGIString_meson;
        else if (strcmp(item->valuestring, "GISCREEN_MESON") == 0) input->params = argsGIScreen_meson;
        else if (strcmp(item->valuestring, "GISCREEN_BBBAR") == 0) input->params = argsGIScreen_bbbar;
        else if (strcmp(item->valuestring, "GISCREEN_CCBAR") == 0) input->params = argsGIScreen_ccbar;
        else {
            fprintf(stderr, "Unknown model preset: %s\n", item->valuestring);
            exit(1);
        }
    }
}

static void parse_basis_section(const cJSON *root, argsInput_t *input, input_flags_t *flags)
{
    cJSON *basis_json = read_object_item(root, "basis");
    const char *type = read_string_item(basis_json, "type")->valuestring;

    if (strcmp(type, "GEM") == 0) input->orbit = ORBIT_GEM;
    else if (strcmp(type, "CRG") == 0) input->orbit = ORBIT_CRG;
    else if (strcmp(type, "CSM") == 0) input->orbit = ORBIT_CSM;
    else if (strcmp(type, "SHO") == 0) input->orbit = ORBIT_SHO;
    else {
        fprintf(stderr, "Unknown basis type: %s\n", type);
        exit(1);
    }
    flags->basis = true;

    if (input->orbit == ORBIT_GEM || input->orbit == ORBIT_CRG || input->orbit == ORBIT_CSM) {
        input->nmax = read_number_item(basis_json, "nmax")->valueint;
        flags->nmax = true;
        input->rmax = read_number_item(basis_json, "rmax")->valuedouble;
        flags->rmax = true;
        input->rmin = read_number_item(basis_json, "rmin")->valuedouble;
        flags->rmin = true;
    }

    if (input->orbit == ORBIT_CRG) {
        input->omega = read_number_item(basis_json, "omega")->valuedouble;
        flags->omega = true;
    }

    if (input->orbit == ORBIT_CSM) {
        input->theta = read_number_item(basis_json, "theta")->valuedouble;
        flags->theta = true;
    }

    if (input->orbit == ORBIT_SHO) {
        input->beta = read_number_item(basis_json, "beta")->valuedouble;
        flags->beta = true;
    }
}

static void validate_input(const argsInput_t *input, const input_flags_t *flags)
{
    check_flag(flags->project, "Missing required global field: project");
    check_flag(flags->task, "Missing required global field: task");
    check_flag(flags->model, "Missing required system field: model");
    check_flag(flags->system, "Missing required system field: system");
    check_flag(flags->basis, "Missing required basis field: basis");

    if (input->system == SYSTEM_MESON) {
        check_flag(flags->f1, "Missing required meson quantum number: f1");
        check_flag(flags->f2, "Missing required meson quantum number: f2");
        check_flag(flags->S, "Missing required meson quantum number: S");
        check_flag(flags->L, "Missing required meson quantum number: L");
        check_flag(flags->J, "Missing required meson quantum number: J");
    }
    else {
        fprintf(stderr, "Unsupported system type: %d\n", input->system);
        exit(1);
    }

    if (input->orbit == ORBIT_GEM) {
        check_flag(flags->nmax, "Missing required GEM basis parameter: nmax");
        check_flag(flags->rmax, "Missing required GEM basis parameter: rmax");
        check_flag(flags->rmin, "Missing required GEM basis parameter: rmin");
        return;
    }

    if (input->orbit == ORBIT_CRG) {
        check_flag(flags->nmax, "Missing required CRG basis parameter: nmax");
        check_flag(flags->rmax, "Missing required CRG basis parameter: rmax");
        check_flag(flags->rmin, "Missing required CRG basis parameter: rmin");
        check_flag(flags->omega, "Missing required CRG basis parameter: omega");
        return;
    }

    if (input->orbit == ORBIT_CSM) {
        check_flag(flags->nmax, "Missing required CSM basis parameter: nmax");
        check_flag(flags->rmax, "Missing required CSM basis parameter: rmax");
        check_flag(flags->rmin, "Missing required CSM basis parameter: rmin");
        check_flag(flags->theta, "Missing required CSM basis parameter: theta");
        return;
    }

    if (input->orbit == ORBIT_SHO) {
        check_flag(flags->beta, "Missing required SHO basis parameter: beta");
        return;
    }

    fprintf(stderr, "Unsupported basis type: %d\n", input->orbit);
    exit(1);
}

int parse_input_file(const char *filename, argsInput_t *input)
{
    char *json_text = read_input_file(filename);
    const char *parse_error = NULL;
    cJSON *root = cJSON_Parse(json_text);
    input_flags_t flags = {0};

    if (!root) {
        parse_error = cJSON_GetErrorPtr();
        fprintf(stderr, "Invalid JSON input in %s", filename);
        if (parse_error) fprintf(stderr, " near: %.40s", parse_error);
        fprintf(stderr, "\n");
        free(json_text);
        exit(1);
    }

    strncpy(input->project, read_string_item(root, "project")->valuestring, 255);
    input->project[255] = '\0';
    flags.project = true;

    parse_task_string(read_string_item(root, "task")->valuestring, input, &flags);
    parse_system_section(root, input, &flags);
    parse_model_section(root, input, &flags);
    parse_basis_section(root, input, &flags);

    validate_input(input, &flags);

    cJSON_Delete(root);
    free(json_text);

    return 1;
}
