/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gemstore/math/interplt.h>
#include <gemstore/math/matrix.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct {
    int n;          /* order of the data point */
    double value;   /* original value (input as data, will be modified internally) */
    int anomaly;    /* anomaly flag (1=anomaly, 0=normal) */
} Point;

void interpolate_quadratic(array_t *data)
{
if (data == NULL || data->len < 4 || data->value == NULL) {
        return;
    }

    int N = data->len;
    
    /* temporary Point array */
    Point *points = (Point *)malloc(N * sizeof(Point));
    if (points == NULL) return;

    /* initialize */
    for (int i = 0; i < N; i++) {
        points[i].n = i + 1;
        points[i].value = data->value[i];
        points[i].anomaly = 0;
    }

    /* first step: detect anomalies (i=1 to N-2) */
    for (int i = 1; i < N - 1; i++) {
        double delta_prev = points[i].value - points[i - 1].value;
        double delta_curr = (i + 1 < N) ? (points[i + 1].value - points[i].value) : 0.0;

        /* violate monotonic increase or non-increasing difference → mark as anomaly */
        if (points[i].value <= points[i - 1].value || 
            (i + 1 < N && delta_curr <= delta_prev)) {
            points[i].anomaly = 1;
        }
    }

    /* second step: fix anomalies and write back to data->value */
    for (int i = 1; i < N - 1; i++) {
        if (!points[i].anomaly) {
            data->value[i] = points[i].value;   // normal state keeps original value
            continue;
        }

        /* find the two previous clean points: j1 < j2 < i */
        int j2 = i - 1;
        while (j2 >= 0 && points[j2].anomaly) j2--;
        int j1 = j2 - 1;
        while (j1 >= 0 && points[j1].anomaly) j1--;

        /* find the next clean point: k > i */
        int k = i + 1;
        while (k < N && points[k].anomaly) k++;

        double fixed_value;

        if (j1 >= 0 && j2 >= 0 && k < N) {
            /* quadratic Lagrange interpolation */
            double x1 = points[j1].n, y1 = points[j1].value;
            double x2 = points[j2].n, y2 = points[j2].value;
            double x3 = points[k].n,  y3 = points[k].value;
            double x  = points[i].n;

            double l1 = ((x - x2)*(x - x3)) / ((x1 - x2)*(x1 - x3));
            double l2 = ((x - x1)*(x - x3)) / ((x2 - x1)*(x2 - x3));
            double l3 = ((x - x1)*(x - x2)) / ((x3 - x1)*(x3 - x2));

            fixed_value = y1 * l1 + y2 * l2 + y3 * l3;
        }
        else if (j2 >= 0 && k < N) {
            /* fallback to linear interpolation */
            double x1 = points[j2].n, y1 = points[j2].value;
            double x2 = points[k].n,  y2 = points[k].value;
            double x  = points[i].n;

            fixed_value = y1 + (y2 - y1) * (x - x1) / (x2 - x1);
        }
        else {
            /* extreme case: copy the nearest clean point */
            fixed_value = (j2 >= 0) ? points[j2].value : points[k].value;
        }

        /* write back the result */
        data->value[i] = fixed_value;
        points[i].value = fixed_value;   // update temporary array to avoid using fixed value later
    }

    free(points);
}