typedef struct{
	int n;
	int m;
	double **p;
}matrix;

void initmatrix(matrix *mat,int n,int m)
{
	int i,j;
	mat->p=(double**)malloc(n*sizeof(double*));
	if(NULL==mat->p)
	{
		printf("error_initmatrix\n");
	}
	for(i=0;i<n;i++)
	{
		mat->p[i]=(double*)malloc(m*sizeof(double));
		if(NULL==mat->p[i])
		{
			printf("error_initmatrix\n");
		}
		else
		{
			for(j=0;j<m;j++)
			{
				mat->p[i][j]=0;
			}
		}
	}
	mat->n=n;
	mat->m=m;
}

void pushmatrix(matrix *mat)
{
	mat->p=(double**)realloc(mat->p,sizeof(double*)*(mat->n+1));
	mat->p[mat->n]=(double*)malloc(sizeof(double)*mat->m);
	mat->n++;
}

void printmatrix(matrix mat)
{
	double **p=mat.p;
	int n=mat.n,m=mat.m;
	int i,j;
	printf("     ");
	for(j=0;j<m;j++)
	{
		printf("   %4d    ",j+1);
	}
	printf("\n");
	for(i=0;i<n;i++)
	{
		printf("%3d: ",i+1);
		for(j=0;j<m;j++)
		{
			printf("%10.6f ",p[i][j]);
		}
		printf("\n");
	}
	printf("\n");
}

void freematrix(matrix *mat)
{
	int n=mat->n;
	int i;
	for(i=0;i<n;i++)
	{
		free(mat->p[i]);
	}
	free(mat->p);
	mat->n=0;
	mat->m=0;
}

void printarrayd1(double *p,int n)
{
	int i;
	for(i=0;i<n;i++)
	{
		printf("%25.20f ",p[i]);
	}
	printf("\n");
}

void printarrayd2(double **p,int n,int m)
{
	int i,j;
	printf("     ");
	for(j=0;j<m;j++)
	{
		printf("   %4d    ",j+1);
	}
	printf("\n");
	for(i=0;i<n;i++)
	{
		printf("%3d: ",i+1);
		for(j=0;j<m;j++)
		{
			printf("%10.6f ",p[i][j]);
		}
		printf("\n");
	}
	printf("\n");
}
void printarraye1(double *p,int n)
{
	int i;
	for(i=0;i<n;i++)
	{
		printf("%25.16E ",p[i]);
	}
	printf("\n");
}

void printarraya1(double *p,int n)
{
	int i;
	for(i=0;i<n;i++)
	{
		printf("%20A ",p[i]);
	}
	printf("\n");
}


void printarrayh1(double *p,int n)
{
	int i,j,k;
	unsigned char *c=(unsigned char*)p;
	k=0;
	for(i=0;i<n;i++)
	{
		for(j=0;j<sizeof(double);j++)
		{
			k++;
			printf("%02X",c[k]);
		}
		printf(" ");
	}
	printf("\n");
}

void freearrayd1(double **p)
{
	free(*p);
}


void freearrayd2(double ***p,int n)
{
/*
	int i;
	for(i=0;i<n;i++)
	{
		free(*(p[i]));
	}
	free(*p);*/
}
