/*
 * Copyright (C) 2025 Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef PARALLEL_H
#define PARALLEL_H

#include <unistd.h>
#include <pthread.h>

typedef void *(*mt_fun)(void *);

typedef struct {
	int ith;
	int len;
	int *cal;
	pthread_mutex_t *mutex;
	int lock;
	mt_fun f;
	void *p;
}mt_args;

int getNumProcessors();

void *mt_run(void *args);

void mt_load(int len, mt_fun f, void *p, int lth);

#endif