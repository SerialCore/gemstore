/*
 * Copyright (C) 2025 Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <task.h>

#include <stdio.h>

void print_help();

int main(int argc, char **argv)
{
	printf("Lambda_s:\n");
	task_baryon(1, 1, 2, -1, 0.5, +1, 2);
    task_baryon(1, 1, 2, -1, 1.5, +1, 2);
    task_baryon(1, 1, 2, -1, 2.5, +1, 2);
    task_baryon(1, 1, 2, -1, 0.5, -1, 1);
    task_baryon(1, 1, 2, -1, 1.5, -1, 1);
    task_baryon(1, 1, 2, -1, 2.5, -1, 1);
	printf("\n");
	
	printf("Sigma_s: \n");
	task_baryon(1, 1, 2, +1, 0.5, +1, 2);
    task_baryon(1, 1, 2, +1, 1.5, +1, 2);
    task_baryon(1, 1, 2, +1, 2.5, +1, 2);
    task_baryon(1, 1, 2, +1, 0.5, -1, 1);
    task_baryon(1, 1, 2, +1, 1.5, -1, 1);
    task_baryon(1, 1, 2, +1, 2.5, -1, 1);
	printf("\n");

	printf("Xi_ss:   \n");
	task_baryon(2, 2, 1, +1, 0.5, +1, 2);
    task_baryon(2, 2, 1, +1, 1.5, +1, 2);
    task_baryon(2, 2, 1, +1, 2.5, +1, 2);
    task_baryon(2, 2, 1, +1, 0.5, -1, 1);
    task_baryon(2, 2, 1, +1, 1.5, -1, 1);
    task_baryon(2, 2, 1, +1, 2.5, -1, 1);
	printf("\n");

	printf("Lambda_c:\n");
	task_baryon(1, 1, 3, -1, 0.5, +1, 2);
    task_baryon(1, 1, 3, -1, 1.5, +1, 2);
    task_baryon(1, 1, 3, -1, 2.5, +1, 2);
    task_baryon(1, 1, 3, -1, 0.5, -1, 1);
    task_baryon(1, 1, 3, -1, 1.5, -1, 1);
    task_baryon(1, 1, 3, -1, 2.5, -1, 1);
	printf("\n");

	printf("Sigma_c: \n");
	task_baryon(1, 1, 3, +1, 0.5, +1, 2);
    task_baryon(1, 1, 3, +1, 1.5, +1, 2);
    task_baryon(1, 1, 3, +1, 2.5, +1, 2);
    task_baryon(1, 1, 3, +1, 0.5, -1, 1);
    task_baryon(1, 1, 3, +1, 1.5, -1, 1);
    task_baryon(1, 1, 3, +1, 2.5, -1, 1);
	printf("\n");

	printf("Omega_c: \n");
	task_baryon(2, 2, 3, +1, 0.5, +1, 2);
    task_baryon(2, 2, 3, +1, 1.5, +1, 2);
    task_baryon(2, 2, 3, +1, 2.5, +1, 2);
    task_baryon(2, 2, 3, +1, 0.5, -1, 1);
    task_baryon(2, 2, 3, +1, 1.5, -1, 1);
    task_baryon(2, 2, 3, +1, 2.5, -1, 1);
	printf("\n");

	printf("Lambda_b:\n");
	task_baryon(1, 1, 4, -1, 0.5, +1, 2);
    task_baryon(1, 1, 4, -1, 1.5, +1, 2);
    task_baryon(1, 1, 4, -1, 2.5, +1, 2);
    task_baryon(1, 1, 4, -1, 0.5, -1, 1);
    task_baryon(1, 1, 4, -1, 1.5, -1, 1);
    task_baryon(1, 1, 4, -1, 2.5, -1, 1);
	printf("\n");
	
	printf("Sigma_b: \n");
	task_baryon(1, 1, 4, +1, 0.5, +1, 2);
    task_baryon(1, 1, 4, +1, 1.5, +1, 2);
    task_baryon(1, 1, 4, +1, 2.5, +1, 2);
    task_baryon(1, 1, 4, +1, 0.5, -1, 1);
    task_baryon(1, 1, 4, +1, 1.5, -1, 1);
    task_baryon(1, 1, 4, +1, 2.5, -1, 1);
	printf("\n");
	
	printf("Omega_b: \n");
	task_baryon(2, 2, 4, +1, 0.5, +1, 2);
    task_baryon(2, 2, 4, +1, 1.5, +1, 2);
    task_baryon(2, 2, 4, +1, 2.5, +1, 2);
    task_baryon(2, 2, 4, +1, 0.5, -1, 1);
    task_baryon(2, 2, 4, +1, 1.5, -1, 1);
    task_baryon(2, 2, 4, +1, 2.5, -1, 1);
	printf("\n");

	printf("Xi_cc:   \n");
	task_baryon(3, 3, 1, +1, 0.5, +1, 2);
    task_baryon(3, 3, 1, +1, 1.5, +1, 2);
    task_baryon(3, 3, 1, +1, 2.5, +1, 2);
    task_baryon(3, 3, 1, +1, 0.5, -1, 1);
    task_baryon(3, 3, 1, +1, 1.5, -1, 1);
    task_baryon(3, 3, 1, +1, 2.5, -1, 1);
	printf("\n");
	
	printf("Xi_bb:   \n");
	task_baryon(4, 4, 1, +1, 0.5, +1, 2);
    task_baryon(4, 4, 1, +1, 1.5, +1, 2);
    task_baryon(4, 4, 1, +1, 2.5, +1, 2);
    task_baryon(4, 4, 1, +1, 0.5, -1, 1);
    task_baryon(4, 4, 1, +1, 1.5, -1, 1);
    task_baryon(4, 4, 1, +1, 2.5, -1, 1);
	printf("\n");

	printf("Omega_cc:\n");
	task_baryon(3, 3, 2, +1, 0.5, +1, 2);
    task_baryon(3, 3, 2, +1, 1.5, +1, 2);
    task_baryon(3, 3, 2, +1, 2.5, +1, 2);
    task_baryon(3, 3, 2, +1, 0.5, -1, 1);
    task_baryon(3, 3, 2, +1, 1.5, -1, 1);
    task_baryon(3, 3, 2, +1, 2.5, -1, 1);
	printf("\n");
	
	printf("Omega_bb:\n");
	task_baryon(4, 4, 2, +1, 0.5, +1, 2);
    task_baryon(4, 4, 2, +1, 1.5, +1, 2);
    task_baryon(4, 4, 2, +1, 2.5, +1, 2);
    task_baryon(4, 4, 2, +1, 0.5, -1, 1);
    task_baryon(4, 4, 2, +1, 1.5, -1, 1);
    task_baryon(4, 4, 2, +1, 2.5, -1, 1);
	printf("\n");

    return 0;
}

void print_help()
{
    printf("This is the help page.\n");
}