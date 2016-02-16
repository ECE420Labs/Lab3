/*
Perform Gauss-Jordan elimanation on a square matrix

-----
Compiling:
    "Lab3IO.c" should be included and "-lm" tag is needed, like
    > gcc serialtester.c Lab3IO.c -o serialtester -lm
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Lab3IO.h"
#include "timer.h"
#include <omp.h>

int **Au;  //pointer to the augmented matrix
int size;
double* X;
int* index;

int main (int argc, char* argv[]) {
    double start, end;

    // Get thread count from input
    int thread_count_;
    if (argc < 2) {
        printf("Please indicate the number of threads!\n");
        return 1;
    }
    thread_count_ = strtol(argv[1], NULL, 10);

    Lab2_loadinput(&Au, &size);

    X = CreateVec(size);
    index = malloc(size * sizeof(int));
    for (i = 0; i < size; ++i)
        index[i] = i;

    GET_TIME(start);
    if (size == 1)
        X[0] = Au[0][1] / Au[0][0];
    else{
        /*Gaussian elimination*/
        #pragma omp parallel num_threads(thread_count_)
        Gauss_elim();
        Jordan_elim(); // serial
        solve();
    }
    GET_TIME(end);

    Lab2_saveoutput(Au, size, end-start);

    DestroyMat(AU, size);
    return 0;
}

void Gauss_elim(void) {
    int my_rank = omp_get_thread_num();
    int thread_count = omp_get_num_threads();

    for (k = my_rank*size/thread_count; k < (my_rank + 1)*size/thread_count; ++k){
        /*Pivoting*/
        temp = 0;
        for (i = k, j = 0; i < size; ++i) // Find row with largest kth element
            if (temp < Au[index[i]][k] * Au[index[i]][k]){ // square value to make it positive
                temp = Au[index[i]][k] * Au[index[i]][k];
                j = i; // j is the row index with the largest element for column k
            }
        if (j != k)/*swap*/{
            i = index[j];
            index[j] = index[k];
            index[k] = i;
        }
        /*calculating*/
        #pragma omp critical
        for (i = k + 1; i < size; ++i){  // serialize for now
            temp = Au[index[i]][k] / Au[index[k]][k];
            for (j = k; j < size + 1; ++j)
                Au[index[i]][j] -= Au[index[k]][j] * temp;
        }
    }
}

void Jordan_elim() {
    /*Jordan elimination*/
    #pragma omp critical
    for (k = size - 1; k > 0; --k){
        for (i = k - 1; i >= 0; --i ){
            temp = Au[index[i]][k] / Au[index[k]][k];
            Au[index[i]][k] -= temp * Au[index[k]][k];
            Au[index[i]][size] -= temp * Au[index[k]][size]; // output vector
        }
    }
}

void solve() {
    /*solution*/
    int my_rank = omp_get_thread_num();
    int thread_count = omp_get_num_threads();

    for (k = my_rank*size/thread_count; k < (my_rank + 1)*size/thread_count; ++k){
        X[k] = Au[index[k]][size] / Au[index[k]][k];
}
