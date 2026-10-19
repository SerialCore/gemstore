/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gemstore/parse.h>
#include <gemstore/input.h>
#include <gemstore/types.h>

#include "cJSON.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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

static cJSON *read_array_item(const cJSON *object, const char *key)
{
    cJSON *item = cJSON_GetObjectItemCaseSensitive(object, key);
    if (!cJSON_IsArray(item)) {
        fprintf(stderr, "Missing or invalid array: %s\n", key);
        exit(1);
    }
    return item;
}

static double read_number_or(const cJSON *object, const char *key, double fallback)
{
    cJSON *item = cJSON_GetObjectItemCaseSensitive(object, key);

    if (item == NULL) return fallback;
    if (!cJSON_IsNumber(item)) {
        fprintf(stderr, "Invalid number: %s\n", key);
        exit(1);
    }
    return item->valuedouble;
}

static void parse_task_string(const char *value, argsInput_t *input)
{
    if (strcmp(value, "SPECTRA") == 0) input->task = TASK_SPECTRA;
    else if (strcmp(value, "DECAY3P0") == 0) input->task = TASK_DECAY3P0;
    else if (strcmp(value, "COUPLCHN") == 0) input->task = TASK_COUPLCHN;
    else if (strcmp(value, "SCATTER") == 0) input->task = TASK_SCATTER;
    else if (strcmp(value, "FITTING") == 0) input->task = TASK_FITTING;
    else {
        fprintf(stderr, "Unknown task: %s\n", value);
        exit(1);
    }
}

static void parse_model_section(const cJSON *root, argsInput_t *input)
{
    cJSON *model_json = read_object_item(root, "model");
    const char *type = read_string_item(model_json, "type")->valuestring;
    const char *param;

    if (strcmp(type, "GISTRING") == 0) input->model = MODEL_GISTRING;
    else if (strcmp(type, "GISCREEN") == 0) input->model = MODEL_GISCREEN;
    else if (strcmp(type, "NRSTRING") == 0) input->model = MODEL_NRSTRING;
    else if (strcmp(type, "NRSCREEN") == 0) input->model = MODEL_NRSCREEN;
    else if (strcmp(type, "BAGMIT") == 0) input->model = MODEL_BAGMIT;
    else {
        fprintf(stderr, "Unknown model type: %s\n", type);
        exit(1);
    }

    if (input->task == TASK_FITTING) {
        return;
    }

    param = read_string_item(model_json, "param")->valuestring;

    if (input->model == MODEL_GISTRING) {
        if (strcmp(param, "GISTRING_MESON") == 0) input->param = PARAM_GISTRING_MESON;
        else if (strcmp(param, "GISTRING_BARYON") == 0) input->param = PARAM_GISTRING_BARYON;
        else if (strcmp(param, "GISTRING_CUSTOM") == 0) input->param = PARAM_GISTRING_CUSTOM;
        else {
            fprintf(stderr, "Unknown GISTRING parameter set: %s\n", param);
            exit(1);
        }
    }
    else if (input->model == MODEL_GISCREEN) {
        if (strcmp(param, "GISCREEN_MESON") == 0) input->param = PARAM_GISCREEN_MESON;
        else if (strcmp(param, "GISCREEN_BBBAR") == 0) input->param = PARAM_GISCREEN_BBBAR;
        else if (strcmp(param, "GISCREEN_BCBAR") == 0) input->param = PARAM_GISCREEN_BCBAR;
        else if (strcmp(param, "GISCREEN_BSBAR") == 0) input->param = PARAM_GISCREEN_BSBAR;
        else if (strcmp(param, "GISCREEN_CCBAR") == 0) input->param = PARAM_GISCREEN_CCBAR;
        else if (strcmp(param, "GISCREEN_CSBAR") == 0) input->param = PARAM_GISCREEN_CSBAR;
        else if (strcmp(param, "GISCREEN_CUSTOM") == 0) input->param = PARAM_GISCREEN_CUSTOM;
        else {
            fprintf(stderr, "Unknown GISCREEN parameter set: %s\n", param);
            exit(1);
        }
    }
    else if (input->model == MODEL_NRSTRING) {
        if (strcmp(param, "NRSTRING_MESON") == 0) input->param = PARAM_NRSTRING_MESON;
        else if (strcmp(param, "NRSTRING_CUSTOM") == 0) input->param = PARAM_NRSTRING_CUSTOM;
        else {
            fprintf(stderr, "Unknown NRSTRING parameter set: %s\n", param);
            exit(1);
        }
    }
    else if (input->model == MODEL_NRSCREEN) {
        if (strcmp(param, "NRSCREEN_MESON") == 0) input->param = PARAM_NRSCREEN_MESON;
        else if (strcmp(param, "NRSCREEN_CUSTOM") == 0) input->param = PARAM_NRSCREEN_CUSTOM;
        else {
            fprintf(stderr, "Unknown NRSCREEN parameter set: %s\n", param);
            exit(1);
        }
    }
    else if (input->model == MODEL_BAGMIT) {
        if (strcmp(param, "BAGMIT_DEFAULT") == 0) input->param = PARAM_BAGMIT_DEFAULT;
        else if (strcmp(param, "BAGMIT_CUSTOM") == 0) input->param = PARAM_BAGMIT_CUSTOM;
        else {
            fprintf(stderr, "Unknown BAGMIT parameter set: %s\n", param);
            exit(1);
        }
    }

    if (input->param == PARAM_GISTRING_CUSTOM || input->param == PARAM_GISCREEN_CUSTOM
        || input->param == PARAM_NRSTRING_CUSTOM || input->param == PARAM_NRSCREEN_CUSTOM
        || input->param == PARAM_BAGMIT_CUSTOM) {
        const char *file = read_string_item(model_json, "file")->valuestring;
        if (strlen(file) > 0) {
            strncpy(input->param_file, file, 256);
            input->param_file[255] = '\0';
        }
        else {
            fprintf(stderr, "Missing parameter file for custom parameter set\n");
            exit(1);
        }
    }
}

static void parse_system_section(const cJSON *root, argsInput_t *input)
{
    cJSON *system_json = read_object_item(root, "system");
    const char *type = read_string_item(system_json, "type")->valuestring;

    if (strcmp(type, "BAG") == 0) input->system = SYSTEM_BAG;
    else if (strcmp(type, "MESON") == 0) input->system = SYSTEM_MESON;
    else if (strcmp(type, "BARYON") == 0) input->system = SYSTEM_BARYON;
    else {
        fprintf(stderr, "Unsupported system type: %s\n", type);
        exit(1);
    }

    if (input->task == TASK_FITTING) {
        return;
    }

    if (input->system == SYSTEM_BAG) {
        input->nn = read_number_item(system_json, "nn")->valueint;
        input->ns = read_number_item(system_json, "ns")->valueint;
        input->nc = read_number_item(system_json, "nc")->valueint;
        input->nb = read_number_item(system_json, "nb")->valueint;
        input->f1 = (int)read_number_or(system_json, "f1", 0.0);
        input->f2 = (int)read_number_or(system_json, "f2", 0.0);
        input->f3 = (int)read_number_or(system_json, "f3", 0.0);
        input->f4 = (int)read_number_or(system_json, "f4", 0.0);
        input->f5 = (int)read_number_or(system_json, "f5", 0.0);
        input->f6 = (int)read_number_or(system_json, "f6", 0.0);
    }
    else if (input->system == SYSTEM_MESON) {
        input->f1 = read_number_item(system_json, "f1")->valueint;
        input->f2 = read_number_item(system_json, "f2")->valueint;
        input->S = read_number_item(system_json, "S")->valuedouble;
        input->L = read_number_item(system_json, "L")->valuedouble;
        input->J = read_number_item(system_json, "J")->valuedouble;
    }
    else if (input->system == SYSTEM_BARYON) {
        input->f1 = read_number_item(system_json, "f1")->valueint;
        input->f2 = read_number_item(system_json, "f2")->valueint;
        input->f3 = read_number_item(system_json, "f3")->valueint;
        input->J = read_number_item(system_json, "J")->valuedouble;
        input->P = read_number_item(system_json, "P")->valueint;
        input->sym12 = read_number_item(system_json, "sym12")->valueint;
        input->Lmax = read_number_item(system_json, "Lmax")->valueint;
        input->jl = read_number_item(system_json, "jl")->valueint;
    }
}

static void parse_cmi_basis(const cJSON *basis_json, argsInput_t *input)
{
    cJSON *color_json = read_array_item(basis_json, "color");
    cJSON *spin_json = read_array_item(basis_json, "spin");
    cJSON *item;
    int n;

    n = cJSON_GetArraySize(color_json);
    if (n < 1 || n > CMI_COLOR_MAX) {
        fprintf(stderr, "CMI color basis must contain 1 to %d channels\n", CMI_COLOR_MAX);
        exit(1);
    }
    input->ncolor = n;
    for (int i = 0; i < n; i++) {
        const char *name;
        int id;

        item = cJSON_GetArrayItem(color_json, i);
        if (!cJSON_IsString(item) || item->valuestring == NULL) {
            fprintf(stderr, "CMI color channel %d is not a string\n", i + 1);
            exit(1);
        }
        name = item->valuestring;
        if (strcmp(name, "meson") == 0) id = CMI_COLOR_MESON;
        else if (strcmp(name, "baryon") == 0) id = CMI_COLOR_BARYON;
        else if (strcmp(name, "tetra11") == 0) id = CMI_COLOR_TETRA11;
        else if (strcmp(name, "tetra33") == 0) id = CMI_COLOR_TETRA33;
        else if (strcmp(name, "tetra66") == 0) id = CMI_COLOR_TETRA66;
        else if (strcmp(name, "tetra88") == 0) id = CMI_COLOR_TETRA88;
        else if (strcmp(name, "penta688") == 0) id = CMI_COLOR_PENTA688;
        else if (strcmp(name, "penta388") == 0) id = CMI_COLOR_PENTA388;
        else if (strcmp(name, "penta11") == 0) id = CMI_COLOR_PENTA11;
        else if (strcmp(name, "hexa1010") == 0) id = CMI_COLOR_HEXA1010;
        else if (strcmp(name, "hexa6868") == 0) id = CMI_COLOR_HEXA6868;
        else if (strcmp(name, "hexa6838") == 0) id = CMI_COLOR_HEXA6838;
        else if (strcmp(name, "hexa3868") == 0) id = CMI_COLOR_HEXA3868;
        else if (strcmp(name, "hexa3838") == 0) id = CMI_COLOR_HEXA3838;
        else if (strcmp(name, "hexa11") == 0) id = CMI_COLOR_HEXA11;
        else {
            fprintf(stderr, "Unknown CMI color basis: %s\n", name);
            exit(1);
        }
        input->color_id[i] = id;
    }

    n = cJSON_GetArraySize(spin_json);
    if (n < 1 || n > CMI_SPIN_MAX) {
        fprintf(stderr, "CMI spin basis must contain 1 to %d channels\n", CMI_SPIN_MAX);
        exit(1);
    }
    input->nspin = n;
    for (int i = 0; i < n; i++) {
        item = cJSON_GetArrayItem(spin_json, i);
        if (!cJSON_IsObject(item)) {
            fprintf(stderr, "CMI spin channel %d is not an object\n", i + 1);
            exit(1);
        }
        input->spin_ch[i].s12 = read_number_or(item, "s12", 0.0);
        input->spin_ch[i].s34 = read_number_or(item, "s34", 0.0);
        input->spin_ch[i].s123 = read_number_or(item, "s123", 0.0);
        input->spin_ch[i].s45 = read_number_or(item, "s45", 0.0);
        input->spin_ch[i].s456 = read_number_or(item, "s456", 0.0);
        input->spin_ch[i].st = read_number_or(item, "st", 0.0);
        input->spin_ch[i].st3 = read_number_or(item, "st3", 0.0);
    }
}

static void parse_basis_section(const cJSON *root, argsInput_t *input)
{
    cJSON *basis_json = read_object_item(root, "basis");
    const char *type = read_string_item(basis_json, "type")->valuestring;

    if (strcmp(type, "GEM") == 0) input->basis = BASIS_GEM;
    else if (strcmp(type, "SHO") == 0) input->basis = BASIS_SHO;
    else if (strcmp(type, "CMI") == 0) input->basis = BASIS_CMI;
    else if (strcmp(type, "NONE") == 0) input->basis = BASIS_NONE;
    else {
        fprintf(stderr, "Unknown basis type: %s\n", type);
        exit(1);
    }

    if (input->basis == BASIS_GEM) {
        input->nmax = read_number_item(basis_json, "nmax")->valueint;
        input->rmax = read_number_item(basis_json, "rmax")->valuedouble;
        input->rmin = read_number_item(basis_json, "rmin")->valuedouble;
    }
    else if (input->basis == BASIS_SHO) {
        input->nmax = read_number_item(basis_json, "nmax")->valueint;
        input->beta = read_number_item(basis_json, "beta")->valuedouble;
    }
    else if (input->basis == BASIS_CMI) {
        parse_cmi_basis(basis_json, input);
    }
    else if (input->basis == BASIS_NONE) {
        input->ncolor = 0;
        input->nspin = 0;
    }
}

static void parse_print_section(const cJSON *root, argsInput_t *input)
{
    cJSON *print_json = read_object_item(root, "print");
    const char *pot_str = read_string_item(print_json, "pot")->valuestring;
    const char *wfn_str = read_string_item(print_json, "wfn")->valuestring;

    input->print_pot = (strcmp(pot_str, "true") == 0) ? 1 : 0;
    input->print_wfn = (strcmp(wfn_str, "true") == 0) ? 1 : 0;
}

static void parse_fit_section(const cJSON *root, argsInput_t *input)
{
    cJSON *fit_json = read_object_item(root, "fit");
    const char *target = read_string_item(fit_json, "target")->valuestring;

    if (strcmp(target, "GISCREEN_MESON") == 0) input->target = FITTING_GISCREEN_MESON;
    else if (strcmp(target, "GISCREEN_BBBAR") == 0) input->target = FITTING_GISCREEN_BBBAR;
    else if (strcmp(target, "GISCREEN_BCBAR") == 0) input->target = FITTING_GISCREEN_BCBAR;
    else if (strcmp(target, "GISCREEN_BSBAR") == 0) input->target = FITTING_GISCREEN_BSBAR;
    else if (strcmp(target, "GISCREEN_CCBAR") == 0) input->target = FITTING_GISCREEN_CCBAR;
    else if (strcmp(target, "GISCREEN_CSBAR") == 0) input->target = FITTING_GISCREEN_CSBAR;
    else {
        fprintf(stderr, "Unknown fit target: %s\n", target);
        exit(1);
    }
}

void parse_input_file(const char *filename, argsInput_t *input)
{
    char *json_text = read_input_file(filename);
    const char *parse_error = NULL;
    cJSON *root = cJSON_Parse(json_text);

    if (!root) {
        parse_error = cJSON_GetErrorPtr();
        fprintf(stderr, "Invalid JSON input in %s", filename);
        if (parse_error) fprintf(stderr, " near: %.40s", parse_error);
        fprintf(stderr, "\n");
        free(json_text);
        exit(1);
    }

    strncpy(input->project, read_string_item(root, "project")->valuestring, 256);
    input->project[255] = '\0';

    parse_task_string(read_string_item(root, "task")->valuestring, input);
    parse_model_section(root, input);
    parse_system_section(root, input);
    parse_basis_section(root, input);
    parse_print_section(root, input);
    if (input->task == TASK_FITTING) {
        parse_fit_section(root, input);
    }

    cJSON_Delete(root);
    free(json_text);
}

void parse_param_GISTRING(const char *filename, argsGIModel_t *args_model)
{
    char *json_text = read_input_file(filename);
    const char *parse_error = NULL;
    cJSON *root = cJSON_Parse(json_text);
    cJSON *param_json;

    if (!root) {
        parse_error = cJSON_GetErrorPtr();
        fprintf(stderr, "Invalid GISTRING parameter JSON in %s", filename);
        if (parse_error) fprintf(stderr, " near: %.40s", parse_error);
        fprintf(stderr, "\n");
        free(json_text);
        exit(1);
    }

    param_json = read_object_item(root, "param");
    args_model->mn = read_number_item(param_json, "mn")->valuedouble;
    args_model->ms = read_number_item(param_json, "ms")->valuedouble;
    args_model->mc = read_number_item(param_json, "mc")->valuedouble;
    args_model->mb = read_number_item(param_json, "mb")->valuedouble;
    args_model->b = read_number_item(param_json, "b")->valuedouble;
    args_model->c = read_number_item(param_json, "c")->valuedouble;
    args_model->sigma_0 = read_number_item(param_json, "sigma_0")->valuedouble;
    args_model->s = read_number_item(param_json, "s")->valuedouble;
    args_model->epsilon_cont = read_number_item(param_json, "epsilon_cont")->valuedouble;
    args_model->epsilon_sov = read_number_item(param_json, "epsilon_sov")->valuedouble;
    args_model->epsilon_sos = read_number_item(param_json, "epsilon_sos")->valuedouble;
    args_model->epsilon_tens = read_number_item(param_json, "epsilon_tens")->valuedouble;
    args_model->mu = 0.0;

    cJSON_Delete(root);
    free(json_text);
}

void parse_param_GISCREEN(const char *filename, argsGIModel_t *args_model)
{
    char *json_text = read_input_file(filename);
    const char *parse_error = NULL;
    cJSON *root = cJSON_Parse(json_text);
    cJSON *param_json;

    if (!root) {
        parse_error = cJSON_GetErrorPtr();
        fprintf(stderr, "Invalid GISCREEN parameter JSON in %s", filename);
        if (parse_error) fprintf(stderr, " near: %.40s", parse_error);
        fprintf(stderr, "\n");
        free(json_text);
        exit(1);
    }

    param_json = read_object_item(root, "param");
    args_model->mn = read_number_item(param_json, "mn")->valuedouble;
    args_model->ms = read_number_item(param_json, "ms")->valuedouble;
    args_model->mc = read_number_item(param_json, "mc")->valuedouble;
    args_model->mb = read_number_item(param_json, "mb")->valuedouble;
    args_model->b = read_number_item(param_json, "b")->valuedouble;
    args_model->mu = read_number_item(param_json, "mu")->valuedouble;
    args_model->c = read_number_item(param_json, "c")->valuedouble;
    args_model->sigma_0 = read_number_item(param_json, "sigma_0")->valuedouble;
    args_model->s = read_number_item(param_json, "s")->valuedouble;
    args_model->epsilon_cont = read_number_item(param_json, "epsilon_cont")->valuedouble;
    args_model->epsilon_sov = read_number_item(param_json, "epsilon_sov")->valuedouble;
    args_model->epsilon_sos = read_number_item(param_json, "epsilon_sos")->valuedouble;
    args_model->epsilon_tens = read_number_item(param_json, "epsilon_tens")->valuedouble;

    cJSON_Delete(root);
    free(json_text);
}

void parse_param_NRSTRING(const char *filename, argsNRModel_t *args_model)
{
    char *json_text = read_input_file(filename);
    const char *parse_error = NULL;
    cJSON *root = cJSON_Parse(json_text);
    cJSON *param_json;

    if (!root) {
        parse_error = cJSON_GetErrorPtr();
        fprintf(stderr, "Invalid NRSTRING parameter JSON in %s", filename);
        if (parse_error) fprintf(stderr, " near: %.40s", parse_error);
        fprintf(stderr, "\n");
        free(json_text);
        exit(1);
    }

    param_json = read_object_item(root, "param");
    args_model->mn = read_number_item(param_json, "mn")->valuedouble;
    args_model->ms = read_number_item(param_json, "ms")->valuedouble;
    args_model->mc = read_number_item(param_json, "mc")->valuedouble;
    args_model->mb = read_number_item(param_json, "mb")->valuedouble;
    args_model->b = read_number_item(param_json, "b")->valuedouble;
    args_model->c = read_number_item(param_json, "c")->valuedouble;
    args_model->alpha_s = read_number_item(param_json, "alpha_s")->valuedouble;
    args_model->sigma = read_number_item(param_json, "sigma")->valuedouble;
    args_model->mu = 0.0;

    cJSON_Delete(root);
    free(json_text);
}

void parse_param_NRSCREEN(const char *filename, argsNRModel_t *args_model)
{
    char *json_text = read_input_file(filename);
    const char *parse_error = NULL;
    cJSON *root = cJSON_Parse(json_text);
    cJSON *param_json;

    if (!root) {
        parse_error = cJSON_GetErrorPtr();
        fprintf(stderr, "Invalid NRSCREEN parameter JSON in %s", filename);
        if (parse_error) fprintf(stderr, " near: %.40s", parse_error);
        fprintf(stderr, "\n");
        free(json_text);
        exit(1);
    }

    param_json = read_object_item(root, "param");
    args_model->mn = read_number_item(param_json, "mn")->valuedouble;
    args_model->ms = read_number_item(param_json, "ms")->valuedouble;
    args_model->mc = read_number_item(param_json, "mc")->valuedouble;
    args_model->mb = read_number_item(param_json, "mb")->valuedouble;
    args_model->b = read_number_item(param_json, "b")->valuedouble;
    args_model->mu = read_number_item(param_json, "mu")->valuedouble;
    args_model->c = read_number_item(param_json, "c")->valuedouble;
    args_model->alpha_s = read_number_item(param_json, "alpha_s")->valuedouble;
    args_model->sigma = read_number_item(param_json, "sigma")->valuedouble;

    cJSON_Delete(root);
    free(json_text);
}

void parse_param_BAGMIT(const char *filename, argsBAGModel_t *args_model)
{
    char *json_text = read_input_file(filename);
    const char *parse_error = NULL;
    cJSON *root = cJSON_Parse(json_text);
    cJSON *param_json;

    if (!root) {
        parse_error = cJSON_GetErrorPtr();
        fprintf(stderr, "Invalid MITBAG parameter JSON in %s", filename);
        if (parse_error) fprintf(stderr, " near: %.40s", parse_error);
        fprintf(stderr, "\n");
        free(json_text);
        exit(1);
    }

    param_json = read_object_item(root, "param");
    args_model->mn = read_number_item(param_json, "mn")->valuedouble;
    args_model->ms = read_number_item(param_json, "ms")->valuedouble;
    args_model->mc = read_number_item(param_json, "mc")->valuedouble;
    args_model->mb = read_number_item(param_json, "mb")->valuedouble;
    args_model->Bcs = read_number_item(param_json, "Bcs")->valuedouble;
    args_model->Bbs = read_number_item(param_json, "Bbs")->valuedouble;
    args_model->Bcc = read_number_item(param_json, "Bcc")->valuedouble;
    args_model->Bbc = read_number_item(param_json, "Bbc")->valuedouble;
    args_model->Bbb = read_number_item(param_json, "Bbb")->valuedouble;
    args_model->B = read_number_item(param_json, "B")->valuedouble;
    args_model->Z = read_number_item(param_json, "Z")->valuedouble;

    cJSON_Delete(root);
    free(json_text);
}
