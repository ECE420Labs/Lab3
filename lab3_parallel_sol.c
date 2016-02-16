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

int thread_count_;
int **Au;  //pointer to the augmented matrix
int size;

int main (int argc, char* argv[]) {
    if (argc < 2) {
        printf("Please indicate the number of threads!\n");
        return 1;
    }
    thread_count_ = strtol(argv[1], NULL, 10);

    Lab2_loadinput(&Au, &size);

    if (size == 1)
        X[0] = Au[0][1] / Au[0][0];
    else{
        /*Gaussian elimination*/
        for (k = 0; k < size - 1; ++k){
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
            for (i = k + 1; i < size; ++i){  // serialize for now
                temp = Au[index[i]][k] / Au[index[k]][k];
                for (j = k; j < size + 1; ++j)
                    Au[index[i]][j] -= Au[index[k]][j] * temp;
            }
        }
        /*Jordan elimination*/
        for (k = size - 1; k > 0; --k){
            for (i = k - 1; i >= 0; --i ){
                temp = Au[index[i]][k] / Au[index[k]][k];
                Au[index[i]][k] -= temp * Au[index[k]][k];
                Au[index[i]][size] -= temp * Au[index[k]][size]; // output vector
            }
        }
        /*solution*/
        for (k=0; k< size; ++k)
            X[k] = Au[index[k]][size] / Au[index[k]][k];
    }
}
