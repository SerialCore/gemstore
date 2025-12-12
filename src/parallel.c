/*
 * Copyright (C) 2025, Wen-Xuan Zhang <serialcore@outlook.com>
 * Copyright (C) 2025, Si-Qiang Luo <luosq15@lzu.edu.cn>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <parallel.h>

#include <stdlib.h>
#include <unistd.h>
#include <pthread.h>

int getNumProcessors()
{
    return sysconf(_SC_NPROCESSORS_ONLN);
}

void *mt_run(void *args)
{
	mt_args *arg = (mt_args*)args;
	int len = arg->len;
	int *cal = arg->cal;
	pthread_mutex_t *mutex = arg->mutex;
	mt_fun f = arg->f;
	int i, pos;
	
    for (i = 0; i < len; i++) {
		pthread_mutex_lock(mutex);
		if (*cal < len) {
			pos = *cal;
			*cal = *cal + 1;
		}
        else {
			pos = -1;
		}
		pthread_mutex_unlock(mutex);
		if (pos >= 0) {
			f(&arg[pos]);
		}
	}

	return NULL;
}

void mt_load(int len, mt_fun f, void *p, int lth)
{
	pthread_t *tid = (pthread_t*)malloc(sizeof(pthread_t)*lth);
	int cal = 0;
	pthread_mutex_t mutex;
	mt_args *arg = (mt_args*)malloc(sizeof(mt_args)*len);
	int i;
	
	pthread_mutex_init(&mutex, NULL);
	
	for (i = 0; i < len; i++) {
		arg[i].ith = i;
		arg[i].len = len;
		arg[i].cal = &cal;
		arg[i].mutex = &mutex;
		arg[i].lock = 0;
		arg[i].p = p;
		arg[i].f = f;
	}
	for (i = 0; i < lth; i++) {
		pthread_create(&tid[i], NULL, mt_run, arg);
	}
	for (i = 0; i < lth; i++) {
		pthread_join(tid[i], NULL);
	}

	pthread_mutex_destroy(&mutex);
	
	free(tid);
	free(arg);
}
