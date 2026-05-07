/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef GEMSTORE_ENTRY
#define GEMSTORE_ENTRY

#include <gemstore/types.h>
#include <gemstore/param/argset.h>

/* dispatch the compute tasks */
void entry_compute(const char* arg);

/* dispatch the fittiing tasks */
void entry_fitting(const char* arg);

/* dispatch the debug tasks */
void entry_debug(const char* arg);

#endif