/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gemstore/read.h>

#include "cJSON.h"

#include <stdio.h>
#include <stdlib.h>

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

void parse_meson_state(const char *filename, array_t *mass, array_t *radius, matrix_t *eigenvector)
{
    if (filename == NULL || mass == NULL || radius == NULL || eigenvector == NULL) {
        fprintf(stderr, "Error: Invalid input parameters to parse_meson_state()\n");
        exit(1);
    }

    if (mass->value == NULL || radius->value == NULL || eigenvector->value == NULL) {
        fprintf(stderr, "Error: Arrays must be pre-allocated (non-NULL values)\n");
        exit(1);
    }

    char *json_text = read_input_file(filename);
    const char *parse_error = NULL;
    cJSON *root = cJSON_Parse(json_text);

    if (!root) {
        parse_error = cJSON_GetErrorPtr();
        fprintf(stderr, "Invalid JSON in state file %s", filename);
        if (parse_error) fprintf(stderr, " near: %.40s", parse_error);
        fprintf(stderr, "\n");
        free(json_text);
        exit(1);
    }

    /* Get states array using read_array_item */
    cJSON *states_json = read_array_item(root, "states");
    int num_states = cJSON_GetArraySize(states_json);
    if (num_states <= 0) {
        fprintf(stderr, "Error: Empty states array in %s\n", filename);
        cJSON_Delete(root);
        free(json_text);
        exit(1);
    }

    /* Verify array dimensions match */
    if (num_states != mass->len || num_states != radius->len || num_states != eigenvector->row) {
        fprintf(stderr, "Error: State count mismatch (expected %d, got %d states in file)\n",
                mass->len, num_states);
        cJSON_Delete(root);
        free(json_text);
        exit(1);
    }

    /* Get eigenvector length from first state */
    cJSON *first_state = cJSON_GetArrayItem(states_json, 0);
    cJSON *first_eigenvector = read_array_item(first_state, "eigenvector");
    int eigenvector_len = cJSON_GetArraySize(first_eigenvector);
    if (eigenvector_len != eigenvector->col) {
        fprintf(stderr, "Error: Eigenvector length mismatch (expected %d, got %d)\n",
                eigenvector->col, eigenvector_len);
        cJSON_Delete(root);
        free(json_text);
        exit(1);
    }

    /* Parse states data */
    for (int n = 0; n < num_states; n++) {
        cJSON *state = cJSON_GetArrayItem(states_json, n);
        if (!cJSON_IsObject(state)) {
            fprintf(stderr, "Error: Invalid state object at index %d\n", n);
            cJSON_Delete(root);
            free(json_text);
            exit(1);
        }

        /* Parse mass */
        cJSON *mass_item = read_number_item(state, "mass");
        mass->value[n] = mass_item->valuedouble;

        /* Parse rms_radius */
        cJSON *radius_item = read_number_item(state, "rms_radius");
        radius->value[n] = radius_item->valuedouble;

        /* Parse eigenvector using read_array_item */
        cJSON *eigen_array = read_array_item(state, "eigenvector");
        int current_eigen_len = cJSON_GetArraySize(eigen_array);
        if (current_eigen_len != eigenvector_len) {
            fprintf(stderr, "Error: Eigenvector length mismatch in state %d (expected %d, got %d)\n",
                    n, eigenvector_len, current_eigen_len);
            cJSON_Delete(root);
            free(json_text);
            exit(1);
        }

        /* Extract eigenvector coefficients */
        for (int d = 0; d < eigenvector_len; d++) {
            cJSON *coeff = cJSON_GetArrayItem(eigen_array, d);
            if (!cJSON_IsNumber(coeff)) {
                fprintf(stderr, "Error: Invalid eigenvector coefficient at state %d, index %d\n", n, d);
                cJSON_Delete(root);
                free(json_text);
                exit(1);
            }
            eigenvector->value[n][d] = coeff->valuedouble;
        }
    }

    cJSON_Delete(root);
    free(json_text);
}
