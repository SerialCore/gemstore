/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gemstore/entry.h>
#include <gemstore/fitting.h>

void call_gem(double *a)
{
    a[0] = 0.1;
    a[1] = 0.2;
    a[2] = 0.3;
}

void chi2_test()
{
    double a[] = {1.1, 1.2, 1.3};
    test(a);
}