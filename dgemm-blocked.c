/* 
    Please include compiler name below (you may also include any other modules you would like to be loaded)

COMPILER= gnu

    Please include All compiler flags and libraries as you want them run. You can simply copy this over from the Makefile's first few lines
 
CC = cc
OPT = -O3
CFLAGS = -Wall -std=gnu99 $(OPT)
MKLROOT = /opt/intel/composer_xe_2013.1.117/mkl
LDLIBS = -lrt -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_sequential.a $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm

*/

#include<time.h>
#include<stdio.h>
#include<stdlib.h>
#include<limits.h>

const char* dgemm_desc = "Simple blocked dgemm.";

#define min(a,b) (((a)<(b))?(a):(b))

/* This auxiliary subroutine performs a smaller dgemm operation
 *  C := C + A * B
 * where C is M-by-N, A is M-by-K, and B is K-by-N. */
static void do_block (int lda, int M, int N, int K, double* A, double* B, double* C)
{
  /* For each row i of A */
  for (int i = 0; i < M; ++i)
    /* For each column j of B */ 
    for (int j = 0; j < N; ++j) 
    {
      /* Compute C(i,j) */
      double cij = C[i+j*lda];
      for (int k = 0; k < K; ++k)
	cij += A[i+k*lda] * B[k+j*lda];
      C[i+j*lda] = cij;
    }
}


/*
void square_dgemm(int lda, double* A, double* B, double* C, const int ROW_BLOCK_SIZE, const int COL_BLOCK_SIZE) {
	// for-each block row of A 
	for(int i=0; i<lda; i+=ROW_BLOCK_SIZE) {
		int M = min(ROW_BLOCK_SIZE, lda-i);
		// for-each block column of B 
		for(int j=0; j<lda; j+=COL_BLOCK_SIZE) {
			int N = min(COL_BLOCK_SIZE, lda-j);
			// for-each block of C
			for(int k=0; k<lda; k+=ROW_BLOCK_SIZE) {
				int K = min(ROW_BLOCK_SIZE, lda-k);
				
				//printf("do it, %d, %d, %d\n", M, N, K);
				//I double checked this line
				do_block(lda, M, N, K, A + i + k*lda, B + k + j*lda, C + i + j*lda);

			}
		}
	}
}

*/

void square_dgemm (int lda, double* A, double* B, double* C, int ROW_BLOCK_SIZE, int COL_BLOCK_SIZE)
{
  /* For each block-row of A */ 
  for (int i = 0; i < lda; i += ROW_BLOCK_SIZE)
    /* For each block-column of B */
    for (int j = 0; j < lda; j += COL_BLOCK_SIZE)
      /* Accumulate block dgemms into block of C */
      for (int k = 0; k < lda; k += ROW_BLOCK_SIZE)
      {
	/* Correct block dimensions if block "goes off edge of" the matrix */
	int M = min (ROW_BLOCK_SIZE, lda-i);
	int N = min (COL_BLOCK_SIZE, lda-j);
	int K = min (ROW_BLOCK_SIZE, lda-k);

	/* Perform individual block dgemm */
	do_block(lda, M, N, K, A + i + k*lda, B + k + j*lda, C + i + j*lda);
      }
}

void test_square_dgemm(int lda, double* A, double* B, double* C) {
	int best_row = 0, best_col = 0, best_time = INT_MAX;

	for(int ROW_BLOCK_SIZE = 1; ROW_BLOCK_SIZE <= lda; ROW_BLOCK_SIZE++) {
		for(int COL_BLOCK_SIZE = 1; COL_BLOCK_SIZE <= lda; COL_BLOCK_SIZE++) {
			clock_t start = clock();
			square_dgemm(lda, A, B, C, ROW_BLOCK_SIZE, COL_BLOCK_SIZE);
			clock_t diff = clock() - start;
			int msec = diff * 1000 / CLOCKS_PER_SEC;
			//printf("Row = %d, Col = %d, time = %d\n", ROW_BLOCK_SIZE, COL_BLOCK_SIZE, msec);
			
			if(msec < best_time) {
				best_time = msec;
				best_row = ROW_BLOCK_SIZE;
				best_col = COL_BLOCK_SIZE;
				printf("row: %d, col: %d, time: %d\n", best_row, best_col, best_time);
			}
		}
	}

	printf("best row: %d, best col: %d\n", best_row, best_col);

}

void fill_array_with_random_doubles(double* arr, int len) {
	for(int i = 0; i < len; i++) {
		printf("%d\n", i);
		arr[i] = rand() * 0.01;
	}
}

int main() {
	srand(120);

	int lda = 100;
	int len = lda * lda;
	double* A = (double*)calloc(len, sizeof(double));
	double* B = (double*)calloc(len, sizeof(double));
	double* C = (double*)calloc(len, sizeof(double));

	
	fill_array_with_random_doubles(A, len);
	printf("done 1\n");
	fill_array_with_random_doubles(B, len);
	printf("done 2\n");

	test_square_dgemm(lda, A, B, C);

	/*
	free(A);
	free(B);
	free(C);
	*/
}


void dgemm_ (int lda, double* A, double* B, double* C)
{
  int ROW_BLOCK_SIZE = 25;
  int COL_BLOCK_SIZE = 25;
  /* For each block-row of A */ 
  for (int i = 0; i < lda; i += ROW_BLOCK_SIZE)
    /* For each block-column of B */
    for (int j = 0; j < lda; j += COL_BLOCK_SIZE)
      /* Accumulate block dgemms into block of C */
      for (int k = 0; k < lda; k += ROW_BLOCK_SIZE)
      {
	/* Correct block dimensions if block "goes off edge of" the matrix */
	int M = min (ROW_BLOCK_SIZE, lda-i);
	int N = min (COL_BLOCK_SIZE, lda-j);
	int K = min (ROW_BLOCK_SIZE, lda-k);

	/* Perform individual block dgemm */
	do_block(lda, M, N, K, A + i + k*lda, B + k + j*lda, C + i + j*lda);
      }
}
