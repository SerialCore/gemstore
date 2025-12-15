/*
 * Copyright (C) 2025 Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef THREAD_H
#define THREAD_H

#include <unistd.h>
#include <pthread.h>

typedef void *(*thread_function)(void *);

typedef struct thread {
	int ith;
	int len;
	int *cal;
	pthread_mutex_t *mutex;
	int lock;
	thread_function f;
	void *p;
} thread_t;

int getNumCores();

void *thread_run(void *args);

void thread_load(int len, thread_function f, void *p, int lth);

#endif