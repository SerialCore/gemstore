/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gemstore/math/interplt.h>

#include <stdlib.h>

typedef struct {
    int n;          /* order of the data point */
    double value;   /* original value (input as data, will be modified internally) */
    int anomaly;    /* anomaly flag (1=anomaly, 0=normal) */
} Point;

static double interpolate_quadratic(const Point *points, const int *indices, double x)
{
    double result = 0.0;

    for (int i = 0; i < 3; i++) {
        double xi = points[indices[i]].n;
        double yi = points[indices[i]].value;
        double basis = 1.0;

        for (int j = 0; j < 3; j++) {
            if (i == j) {
                continue;
            }

            basis *= (x - points[indices[j]].n) / (xi - points[indices[j]].n);
        }

        result += yi * basis;
    }

    return result;
}

static double interpolate_cubic(const Point *points, const int *indices, double x)
{
    double result = 0.0;

    for (int i = 0; i < 4; i++) {
        double xi = points[indices[i]].n;
        double yi = points[indices[i]].value;
        double basis = 1.0;

        for (int j = 0; j < 4; j++) {
            if (i == j) {
                continue;
            }

            basis *= (x - points[indices[j]].n) / (xi - points[indices[j]].n);
        }

        result += yi * basis;
    }

    return result;
}

static double interpolate_linear(const Point *points, const int *indices, double x)
{
    double x1 = points[indices[0]].n;
    double y1 = points[indices[0]].value;
    double x2 = points[indices[1]].n;
    double y2 = points[indices[1]].value;

    return y1 + (y2 - y1) * (x - x1) / (x2 - x1);
}

void interpolate_fix_divergence(double *data, int n)
{
    if (data == NULL || n < 4) {
        return;
    }

    Point *points = (Point *)malloc(n * sizeof(Point));
    if (points == NULL) {
        return;
    }

    for (int i = 0; i < n; i++) {
        points[i].n = i + 1;
        points[i].value = data[i];
        points[i].anomaly = 0;
    }

    for (int i = 1; i < n - 1; i++) {
        if (points[i].value <= points[i - 1].value) {
            points[i].anomaly = 1;
        }
    }

    for (int i = 1; i < n - 1; i++) {
        if (!points[i].anomaly) {
            data[i] = points[i].value;
            continue;
        }

        int indices[4];
        int count = 0;

        for (int j = i - 1; j >= 0 && count < 4; j--) {
            if (!points[j].anomaly) {
                indices[count++] = j;
            }
        }

        for (int j = i + 1; j < n && count < 4; j++) {
            if (!points[j].anomaly) {
                indices[count++] = j;
            }
        }

        if (count >= 4) {
            double fixed_value = interpolate_cubic(points, indices, points[i].n);
            data[i] = fixed_value;
            points[i].value = fixed_value;
        }
        else if (count == 3) {
            double fixed_value = interpolate_quadratic(points, indices, points[i].n);
            data[i] = fixed_value;
            points[i].value = fixed_value;
        }
        else if (count == 2) {
            double fixed_value = interpolate_linear(points, indices, points[i].n);
            data[i] = fixed_value;
            points[i].value = fixed_value;
        }
        else if (count == 1) {
            data[i] = points[indices[0]].value;
            points[i].value = data[i];
        }
    }

    free(points);
}

void interpolate_fix_convergence(double *data, int n)
{
    if (data == NULL || n < 4) {
        return;
    }

    Point *points = (Point *)malloc(n * sizeof(Point));
    if (points == NULL) {
        return;
    }

    for (int i = 0; i < n; i++) {
        points[i].n = i + 1;
        points[i].value = data[i];
        points[i].anomaly = 0;
    }

    for (int i = 1; i < n - 1; i++) {
        if (points[i].value <= points[i - 1].value) {
            points[i].anomaly = 1;
        }
    }

    for (int i = 1; i < n - 1; i++) {
        if (!points[i].anomaly) {
            data[i] = points[i].value;
            continue;
        }

        int indices[4];
        int count = 0;

        for (int j = i - 1; j >= 0 && count < 4; j--) {
            if (!points[j].anomaly) {
                indices[count++] = j;
            }
        }

        for (int j = i + 1; j < n && count < 4; j++) {
            if (!points[j].anomaly) {
                indices[count++] = j;
            }
        }

        if (count >= 4) {
            double fixed_value = interpolate_cubic(points, indices, points[i].n);
            data[i] = fixed_value;
            points[i].value = fixed_value;
        }
        else if (count == 3) {
            double fixed_value = interpolate_quadratic(points, indices, points[i].n);
            data[i] = fixed_value;
            points[i].value = fixed_value;
        }
        else if (count == 2) {
            double fixed_value = interpolate_linear(points, indices, points[i].n);
            data[i] = fixed_value;
            points[i].value = fixed_value;
        }
        else if (count == 1) {
            data[i] = points[indices[0]].value;
            points[i].value = data[i];
        }
    }

    free(points);
}
