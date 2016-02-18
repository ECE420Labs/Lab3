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

void Gauss_elim(int);
void Jordan_elim();
void solve(int);

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
        int k;
        int j;
        int i;
        double temp;

        for (k = 0; k < size - 1; ++k){
            /*Pivoting*/
            temp = 0;
            j = 0;
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
        #pragma omp parallel num_threads(thread_count_)
{
        /*Gaussian elimination*/
        Gauss_elim(thread_count_);
        /*printf("%f ", Au[0][0]);printf("%f ", Au[0][1]);printf("%f\n", Au[0][2]);
        printf("%f ", Au[1][0]);printf("%f ", Au[1][1]);printf("%f\n", Au[1][2]);
        printf("%f ", Au[2][0]);printf("%f ", Au[2][1]);printf("%f\n", Au[2][2]);*/
        //Jordan_elim(); // serial

        //solve(thread_count_);
        //printf("%f ", X[0]);printf("%f ", X[1]);printf("%f\n", X[2]);
        }

    }
    GET_TIME(end);

    Lab3SaveOutput(X, size, end-start);

    DestroyMat(Au, size);
    return 0;
}

void Gauss_elim(int nt) {
    int k;
    int j;
    int i;
    double temp;

    for (k = 0; k < size - 1; ++k){

        /*calculating*/
        #pragma omp for private(j) private(temp)
        for (i = k + 1; i < size; ++i){
            temp = Au[ind[i]][k] / Au[ind[k]][k];
            for (j = k; j < size + 1; ++j) {
                Au[ind[i]][j] -= Au[ind[k]][j] * temp;
            }
        }
    }

    /*Jordan elimination*/
    for (k = size - 1; k > 0; --k){
        for (i = k - 1; i >= 0; --i ){
            temp = Au[ind[i]][k] / Au[ind[k]][k];
            Au[ind[i]][k] -= temp * Au[ind[k]][k];
            Au[ind[i]][size] -= temp * Au[ind[k]][size]; // output vector
        }
    }

    /*solution*/
    #pragma omp for
    for (k=0; k< size; ++k) {
        X[k] = Au[ind[k]][size] / Au[ind[k]][k];
    }
}

void Jordan_elim() {

}

void solve(int nt) {


}
