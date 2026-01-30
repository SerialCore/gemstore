/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gemstore/basis/intrin.h>
#include <gemstore/basis/soc.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

intrin_wfn_t intrin_wfn_init(int num_configs)
{
    intrin_wfn_t wf;
    wf.num_terms = 0;
    wf.num_configs = num_configs;
    wf.coeffs = NULL;
    wf.configs = NULL;
    return wf;
}

void intrin_wfn_product(const intrin_wfn_t *wfnA, const intrin_wfn_t *wfnB, intrin_wfn_t *wfnC, double factor)
{
    for (int i = 0; i < wfnA->num_terms; i++) {
        for (int j = 0; j < wfnB->num_terms; j++) {
            double coeff =  factor * wfnA->coeffs[i] * wfnB->coeffs[j];
            char config[16];
            sprintf(config, "%s%s", wfnA->configs[i], wfnB->configs[j]);
            intrin_wfn_push(wfnC, coeff, config);
        }
    }
}

void intrin_wfn_trim(intrin_wfn_t *wfn)
{
    intrin_wfn_t trimmed = intrin_wfn_init(wfn->num_configs);

    for (int i = 0; i < wfn->num_terms; i++)
    {
        double coeff = wfn->coeffs[i];
        const char *config = wfn->configs[i];

        if (coeff == 0.0) {
            continue;
        }

        int found = 0;
        for (int j = 0; j < trimmed.num_terms; j++) {
            if (strcmp(trimmed.configs[j], config) == 0) {
                trimmed.coeffs[j] += coeff;
                found = 1;
                break;
            }
        }
        if (!found) {
            intrin_wfn_push(&trimmed, coeff, config);
        }
    }

    intrin_wfn_free(wfn);
    *wfn = trimmed;
}

void intrin_wfn_push(intrin_wfn_t *wfn, double coeff, const char *config)
{
    if (strlen(config) != (size_t)wfn->num_configs)
    {
        // Error: invalid config length
        return;
    }

    wfn->num_terms++;
    wfn->coeffs = (double *)realloc(wfn->coeffs, wfn->num_terms * sizeof(double));
    wfn->configs = (char **)realloc(wfn->configs, wfn->num_terms * sizeof(char *));
    wfn->coeffs[wfn->num_terms - 1] = coeff;
    wfn->configs[wfn->num_terms - 1] = strdup(config);
}

void intrin_wfn_print(const intrin_wfn_t *wfn)
{
    printf("{");
    for (int i = 0; i < wfn->num_terms; i++) {
        if (i > 0) printf(", ");
        printf("{%.6f, \"%s\"}", wfn->coeffs[i], wfn->configs[i]);
    }
    printf("}\n");
}

void intrin_wfn_free(intrin_wfn_t *wfn)
{
    for (int i = 0; i < wfn->num_terms; i++)
    {
        free(wfn->configs[i]);
    }
    free(wfn->coeffs);
    free(wfn->configs);
    wfn->num_terms = 0;
    wfn->num_configs = 0;
    wfn->coeffs = NULL;
    wfn->configs = NULL;
}

cg_table_t CallCGTable(double s1, double s2, double st, double st3)
{
    cg_table_t couple;
    couple.num = 0;
    couple.tuples = NULL;

    for (double ms1 = -s1; ms1 <= s1 + 1e-10; ms1 += 1.0) {
        for (double ms2 = -s2; ms2 <= s2 + 1e-10; ms2 += 1.0) {
            if (ms1 + ms2 - st3 == 0.0) {
                double cg = clebsch_gordan(s1, ms1, s2, ms2, st, st3);
                if (cg != 0.0) {
                    couple.num++;
                    couple.tuples = (cg_tuple_t *)realloc(couple.tuples, couple.num * sizeof(cg_tuple_t));
                    couple.tuples[couple.num - 1].ms1 = ms1;
                    couple.tuples[couple.num - 1].ms2 = ms2;
                    couple.tuples[couple.num - 1].cg = cg;
                }
            }
        }
    }

    return couple;
}

void cg_table_free(cg_table_t *t)
{
    free(t->tuples);
    t->num = 0;
    t->tuples = NULL;
}