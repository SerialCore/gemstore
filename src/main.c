/*
 * Copyright (C) 2025 Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <debug.h>

#include <stdio.h>

void print_help();

int main(int argc, char **argv)
{
	printf("Lambda_s:\n");debug02(1,1,2,-1);printf("\n");
	printf("Sigma_s: \n");debug02(1,1,2,+1);printf("\n");
	printf("Xi_ss:   \n");debug02(2,2,1,+1);printf("\n");

	printf("Lambda_c:\n");debug02(1,1,3,-1);printf("\n");
	printf("Sigma_c: \n");debug02(1,1,3,+1);printf("\n");
	printf("Omega_c: \n");debug02(2,2,3,+1);printf("\n");

	printf("Lambda_b:\n");debug02(1,1,4,-1);printf("\n");
	printf("Sigma_b: \n");debug02(1,1,4,+1);printf("\n");
	printf("Omega_b: \n");debug02(2,2,4,+1);printf("\n");

	printf("Xi_cc:   \n");debug02(3,3,1,+1);printf("\n");
	printf("Xi_bb:   \n");debug02(4,4,1,+1);printf("\n");

	printf("Omega_cc:\n");debug02(3,3,2,+1);printf("\n");
	printf("Omega_bb:\n");debug02(4,4,2,+1);printf("\n");

    return 0;
}

void print_help()
{
    printf("This is the help page.\n");
}