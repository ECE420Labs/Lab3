/*
Perform Gauss-Jordan elimanation on a square matrix
-----
Compiling:
"Lab3IO.c" should be included and "-fopenmp" tag is needed, like
> gcc lab3_parallel_sol.c Lab3IO.c -o main -fopenmp
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Lab3IO.h"
#include "timer.h"
#include <omp.h>

void Gauss_elim(int);

double **Au;  //pointer to the augmented matrix
int size;
double* X;
int* ind;

int main (int argc, char* argv[]) {
    double start, end;

    // Get thread count from input
    int thread_count_;
    if (argc < 2) {
        printf("Please indicate the number of threads!\n");
        return 1;
    }
    thread_count_ = strtol(argv[1], NULL, 10);

    Lab3LoadInput(&Au, &size);

    X = CreateVec(size);
    ind = malloc(size * sizeof(int));
    int i;
    for (i = 0; i < size; ++i)
    ind[i] = i;

    GET_TIME(start);
    if (size == 1)
    X[0] = Au[0][1] / Au[0][0];
    else{
        {
            /*Gaussian elimination*/
            Gauss_elim(thread_count_);
        }

    }
    GET_TIME(end);

    Lab3SaveOutput(X, size, end-start);
    printf("%f\n", end-start);
    DestroyMat(Au, size);
    return 0;
}

void Gauss_elim(int nt) {
    int k;
    int j=0;
    int i=0;
    double temp=0;
    #pragma omp parallel num_threads(nt)
    {
        /* Gaussian elimination */
        #pragma omp for ordered private(temp, i, j)
        for (k = 0; k < size - 1; ++k){
            temp = 0;
            j = 0;
            #pragma omp ordered
            /*Pivoting*/
            {
                for (i = k; i < size; ++i) {// Find row with largest kth element
                    if (temp < Au[ind[i]][k] * Au[ind[i]][k]){ // square value to make it positive
                        temp = Au[ind[i]][k] * Au[ind[i]][k];
                        j = i; // j is the row ind with the largest element for column k
                    }
                }
                if (j != k)/*swap*/{
                    i = ind[j];
                    ind[j] = ind[k];
                    ind[k] = i;
                }
            }
            /*calculating*/
            for (i = k + 1; i < size; ++i){
                temp = Au[ind[i]][k] / Au[ind[k]][k];
                for (j = k; j < size + 1; ++j) {
                    Au[ind[i]][j] -= Au[ind[k]][j] * temp;
                }
            }
        }
        /* Jordan elimination */
        #pragma omp for ordered private(temp, i)
        for (k = size - 1; k > 0; --k){
            #pragma omp ordered
            for (i = k - 1; i >= 0; --i ){
                temp = Au[ind[i]][k] / Au[ind[k]][k];
                Au[ind[i]][k] -= temp * Au[ind[k]][k];
                Au[ind[i]][size] -= temp * Au[ind[k]][size]; // output vector
            }
        }
        /* Solve */
        #pragma omp for
        for (k=0; k< size; ++k) {
            X[k] = Au[ind[k]][size] / Au[ind[k]][k];
        }
    }
}
