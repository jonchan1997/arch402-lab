#include "stdio.h"
#include "stdlib.h"
#include "stdbool.h"
#include "string.h"
#include "c_timer.h"
#include "cblas.h" 
#include "math.h"
#include "x86intrin.h"

#define UNROLL (4)

void unoptimized_dgemm(int n, double* A, double* B, double* C)
{
	int i,j,k;
	for(i = 0; i < n; ++i)
	{
		for(j = 0; j< n; ++j)
		{
			double cij = C[i+j*n]; /* cij = C[i][j] */
			for(k = 0; k < n; k++)
			{
				cij += A[i+k*n] * B[k+j*n]; /* cij += A[i][k]*B[k][j] */
				C[i+j*n] = cij; /* C[i][j] = cij */
			}
		}
	}
}
void compare_matrix(int n, double* A1, double* A2){
	int i, j;
	double d1, d2;
	for(i = 0; i< n; i++)
	{
		for(j = 0; j<n; j++)
		{
			d1 = *(A1 + i*n + j);
			d2 = *(A2 + i*n + j);
			if(fabs(d2-d1)/(d2) > 1e-6)
			{
				printf("ERROR: %f(my)<>%f(dgemm)\n \n", d1, d2);
				exit(1);
			}
		}
	}
	printf("Correct result! :-) \n \n");
	
}
void init_matrix(int n, double* A){
	int i,j;
	
	for(i=0; i<n; i++)
	{
		
		for(j=0; j<n; j++)
		{
			
			*(A + i*n + j) = rand() / (RAND_MAX * 1.0);
		}
	}
	
}
void sse_dgemm (int n, double* A, double* B, double* C)
{
	int i,j,k;
	for ( i = 0; i < n; i+=4 )
	{
		for ( j = 0; j < n; j++ )
		{
			__m256d c0 = _mm256_load_pd(C+i+j*n); /* c0 = C[i][j] */
			for( k = 0; k < n; k++ )
			{
				c0 = _mm256_add_pd(c0, /* c0 += A[i][k]*B[k][j] */_mm256_mul_pd(_mm256_load_pd(A+i+k*n),  _mm256_broadcast_sd(B+k+j*n)));
				_mm256_store_pd(C+i+j*n, c0); /* C[i][j] = c0 */
			}
			 
		}
	}
}
void sse_dgemm_with_unrolling (int n, double* A, double* B, double* C)
{
	int i,j,x,k;
	for (i = 0; i < n; i+=UNROLL*4 )
	{
		for (j = 0; j < n; j++ )
		{
			__m256d c[4];
			for (x = 0; x < UNROLL; x++ )
			{
				c[x]=_mm256_load_pd(C+i+x*4+j*n);
			}
			for(k = 0; k < n; k++ )
			{
				__m256d b = _mm256_broadcast_sd(B+k+j*n);
				for (x = 0; x < UNROLL; x++)
				{
					c[x]= _mm256_add_pd(c[x],_mm256_mul_pd(_mm256_load_pd(A+n*k+x*4+i), b));
				}
			}
			for (x = 0; x < UNROLL; x++ )
			{
				_mm256_store_pd(C+i+x*4+j*n, c[x]);
			}
		}
	}
}

int main (int arc, char **argv)
{
	srand(time(NULL));
	bool Quit = false;
	int number = 0;
	while(!Quit)
	{
		printf("Please Type the Number of your Selection then Press ENTER: \n \n");
		printf("1. Test unoptimized_dgemm \n");
		printf("2. Test sse_dgemm \n");
		printf("3. Test sse_dgemm_with_unrolling \n");
		printf("4. Generate Gflop Data \n");
		printf("5. Quit \n\n");
		printf("~");
		scanf("%d", &number);
		if( number == 1)
		{
			printf("what is the value of N? \n");
			printf("~");
			scanf("%d", &number);
			int n = number;
			double *A, *B, *C, *D;
			posix_memalign((void**)&A, 32, n*n*sizeof(double));
			posix_memalign((void**)&B, 32, n*n*sizeof(double));
			posix_memalign((void**)&C, 32, n*n*sizeof(double));
			posix_memalign((void**)&D, 32, n*n*sizeof(double));
			init_matrix(n, A);
			init_matrix(n, B);
			unoptimized_dgemm(n, A, B, C);
			cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, A, n, B, n, 0.0, D, n); 
			compare_matrix(n, C, D);
			continue;
		}
		else if( number == 2)
		{
			
			int n;
			for(n = 64; n <= 1025; n*=2 )
			{
				printf("N: %d \n", n);
				double *A, *B, *C, *D;
				posix_memalign((void**)&A, 32, n*n*sizeof(double));
				posix_memalign((void**)&B, 32, n*n*sizeof(double));
				posix_memalign((void**)&C, 32, n*n*sizeof(double));
				posix_memalign((void**)&D, 32, n*n*sizeof(double));
				init_matrix(n, A);
				init_matrix(n, B);
				sse_dgemm(n, A, B, C);
				cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, A, n, B, n, 0.0, D, n);
				compare_matrix(n, C, D);
			}
			continue;
		}
		else if( number == 3)
		{
			
			int n;
			for(n = 64; n <= 1025; n*=2 )
			{
				printf("N: %d \n", n);
				double *A, *B, *C, *D;
				posix_memalign((void**)&A, 32, n*n*sizeof(double));
				posix_memalign((void**)&B, 32, n*n*sizeof(double));
				posix_memalign((void**)&C, 32, n*n*sizeof(double));
				posix_memalign((void**)&D, 32, n*n*sizeof(double));
				init_matrix(n, A);
				init_matrix(n, B);
				sse_dgemm_with_unrolling(n, A, B, C);
				cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, A, n, B, n, 0.0, D, n);
				compare_matrix(n, C, D);
			}
			continue;
		}
		else if( number == 4)
		{
			int n;
			for(n = 64; n <= 1025; n*=2)
			{
				printf("N: %d \n", n);
				double *A, *B, *C, *D, *E;
				posix_memalign((void**)&A, 32, n*n*sizeof(double));
				posix_memalign((void**)&B, 32, n*n*sizeof(double));
				posix_memalign((void**)&C, 32, n*n*sizeof(double));
				posix_memalign((void**)&D, 32, n*n*sizeof(double));
				posix_memalign((void**)&E, 32, n*n*sizeof(double));
				init_matrix(n, A);
				init_matrix(n, B);
				//unoptimized
				double N = (double) n;
				double t0 = get_cur_time();
				unoptimized_dgemm(n, A, B, C);
				double t1 = get_cur_time();
				double exTime = t1 - t0;
				double unoptimized_Gflops =(double)(2.0*N*N*N)/(exTime*1000000000.0);
				t0 = get_cur_time();
				sse_dgemm(n, A, B, D);
				t1 = get_cur_time();
				exTime = t1 - t0;
				double sse_Gflops =(double)(2.0*N*N*N)/(exTime*1000000000.0);
				t0 = get_cur_time();
				sse_dgemm_with_unrolling(n, A, B, E);
				t1 = get_cur_time();
				exTime = t1 - t0;
				double sse_with_rolling_Gflops =(double)(2.0*N*N*N)/(exTime*1000000000.0);
				printf("unoptimized_dgemm: N = %d Gflops = %lf \n", n, unoptimized_Gflops);
				printf("sse_dgemm: N = %d Gflops = %lf \n", n, sse_Gflops);
				printf("sse_with_rolling_dgemm: N = %d Gflops = %lf \n \n", n, sse_with_rolling_Gflops);
			}
			continue;
		}
		else if( number == 5)
		{
			Quit = true;
			break;
		}
		else
		{
			printf(" \n Invalaid Input, Try Again! \n \n \n");
			continue;
		}
	}	
	return 0;
}

