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
//	debug02(1,1,3);
	debug02(1,2,3);
//	debug02(2,2,3);
//	debug02(1,1,4);
//	debug02(1,2,4);
//	debug02(2,2,4);
	return 0;
}
