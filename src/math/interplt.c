/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gemstore/math/interplt.h>

#include <math.h>
#include <stdlib.h>

typedef struct {
    int n;          /* order of the data point */
    double value;   /* original value (input as data, will be modified internally) */
    int anomaly;    /* anomaly flag (1=anomaly, 0=normal) */
} Point;

static double max_double(double a, double b)
{
    return (a > b) ? a : b;
}

static double min_double(double a, double b)
{
    return (a < b) ? a : b;
}

static double clamp_double(double x, double lo, double hi)
{
    if (x < lo) {
        return lo;
    }
    if (x > hi) {
        return hi;
    }
    return x;
}

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

static double interpolate_linear(const Point *points, const int *indices, double x)
{
    double x1 = points[indices[0]].n;
    double y1 = points[indices[0]].value;
    double x2 = points[indices[1]].n;
    double y2 = points[indices[1]].value;

    return y1 + (y2 - y1) * (x - x1) / (x2 - x1);
}

void interpolate_divergence(double *data, int n)
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

    for (int pass = 0; pass < n; pass++) {
        int changed = 0;

        for (int i = 0; i < n; i++) {
            points[i].anomaly = 0;
        }

        for (int i = 2; i < n; i++) {
            double slope_prev = points[i - 1].value - points[i - 2].value;
            double slope_next = points[i].value - points[i - 1].value;
            double slope_ratio = slope_next / slope_prev;

            if (slope_next <= 0.0 || slope_ratio < 0.75) {
                points[i].anomaly = 1;
            }
        }

        for (int i = 2; i < n; i++) {
            if (!points[i].anomaly) {
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

            double fixed_value;
            if (count >= 3) {
                fixed_value = interpolate_quadratic(points, indices, points[i].n);
            }
            else if (count == 2) {
                fixed_value = interpolate_linear(points, indices, points[i].n);
            }
            else if (count == 1) {
                fixed_value = points[indices[0]].value;
            }
            else {
                continue;
            }

            if (fabs(fixed_value - points[i].value) > 1e-12) {
                points[i].value = fixed_value;
                data[i] = fixed_value;
                changed = 1;
            }
        }

        if (!changed) {
            break;
        }
    }

    free(points);
}

void extrapolate_divergence(double *data, int n)
{
    if (data == NULL || n < 4) {
        return;
    }

    for (int i = 3; i < n; i++) {
        double s_prev2 = data[i - 2] - data[i - 3];
        double s_prev1 = data[i - 1] - data[i - 2];
        double s_curr = data[i] - data[i - 1];
        double trend = 0.5 * (s_prev2 + s_prev1);

        if (s_curr <= 0.0) {
            int pivot = i;
            int anchor = pivot - 1;
            if (anchor < 2) {
                anchor = 2;
            }

            double s0 = data[anchor] - data[anchor - 1];
            double s1 = data[anchor - 1] - data[anchor - 2];
            double slope = max_double(s0, s1);
            double accel = max_double(s0 - s1, 0.0);

            if (slope <= 0.0) {
                slope = max_double(data[anchor] * 0.05, 1e-6);
            }

            if (accel <= 0.0) {
                accel = max_double(slope * 0.12, 1e-6);
            }

            double accel_growth = clamp_double(1.0 + accel / max_double(slope, 1e-6), 1.02, 1.18);

            for (int j = pivot; j < n; j++) {
                double min_step = max_double(slope, max_double(data[j - 1] * 1e-4, 1e-6));
                double candidate = data[j - 1] + min_step;

                if (data[j] > data[j - 1]) {
                    double ceiling = data[j - 1] + 1.35 * min_step;
                    candidate = max_double(candidate, min_double(data[j], ceiling));
                }

                data[j] = candidate;
                accel = max_double(accel * accel_growth, slope * 0.08);
                slope += accel;
            }

            return;
        }

        if (trend > 0.0 && s_curr < 0.55 * trend) {
            int pivot = i;
            int anchor = pivot - 1;
            if (anchor < 2) {
                anchor = 2;
            }

            double s0 = data[anchor] - data[anchor - 1];
            double s1 = data[anchor - 1] - data[anchor - 2];
            double slope = max_double(s0, s1);
            double accel = max_double(s0 - s1, 0.0);

            if (slope <= 0.0) {
                slope = max_double(data[anchor] * 0.05, 1e-6);
            }

            if (accel <= 0.0) {
                accel = max_double(slope * 0.12, 1e-6);
            }

            double accel_growth = clamp_double(1.0 + accel / max_double(slope, 1e-6), 1.02, 1.18);

            for (int j = pivot; j < n; j++) {
                double min_step = max_double(slope, max_double(data[j - 1] * 1e-4, 1e-6));
                double candidate = data[j - 1] + min_step;

                if (data[j] > data[j - 1]) {
                    double ceiling = data[j - 1] + 1.35 * min_step;
                    candidate = max_double(candidate, min_double(data[j], ceiling));
                }

                data[j] = candidate;
                accel = max_double(accel * accel_growth, slope * 0.08);
                slope += accel;
            }

            return;
        }
    }
}
