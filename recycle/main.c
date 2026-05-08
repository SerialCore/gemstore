#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#ifdef _WIN32   
#include <windows.h>
#include "pthread.h"
#pragma comment(lib, "pthreadVC2.lib")
#else
#include <unistd.h>
#include <pthread.h>
#endif


#include "jc.h"
#include "mt.h"
#include "eig.h"
#include "matrixarray.h"
#include "sumckdk.h"
#include "basis.h"
#include "vtype.h"
#include "inteCenV.h"
#include "mfi.h"
#include "res.h"
#include "debug.h"

int main()
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
