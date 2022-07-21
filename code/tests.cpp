#include <iostream>
#include <vector>
#include <cstring>
#include "RCGrid.h"
#include "NCMesh.h"
#include "schemes.h"

double* allocateDoubleArray(int n, int m) {

    double* a = (double*) calloc(n*m, sizeof(double));
    if(!a) {
        printf("Error: could not allocate enough memory\n");
        return NULL;
    }

    for(int i = 0; i < n*m; i++)
        a[i] = i;

    return a;
}

void allocateDoubleArray(double* &a, int n, int m) {

    a = (double*) calloc(n*m, sizeof(double));
    if(!a) {
        printf("Error: could not allocate enough memory for u\n");
        return;
    }

    for(int i = 0; i < n*m; i++)
        a[i] = i;

}

int main(int argc, char* argv[]) {

    // double* a;
    // int n = 5;
    // int m = 3;

    // allocateDoubleArray(a, n, m);
    //
    // if(!a)
    //     printf("Error\n");
    //
    // for(int j = 0; j < m; j++) {
    //     for(int i = 0; i < n; i++) {
    //         printf("%7.0f", a[i+j*m]);
    //     }
    //     printf("\n");
    // }

    // double* b = allocateDoubleArray(n, m);
    // if(!b)
    //     printf("Error");
    //
    // for(int j = m-1; j >= 0; j--) {
    //     for(int i = 0; i < n; i++) {
    //         printf("%7.0f", b[i+j*m]);
    //     }
    //     printf("\n");
    // }

    int n = 5;

    double* a = (double*) calloc(n, sizeof(double));
    double* b = (double*) calloc(n, sizeof(double));

    for(int i = 0; i < n; i++)
        a[i] = i;

    for(int i = 0; i < n; i++)
        b[i] = -i;

    printf("\n\n");
    for(int i = 0; i < n; i++)
        printf("%10.0f%5s%10.0f\n", a[i], "", b[i]);

    for(int i = 0; i < n; i++)
        b[i] = a[i];

    printf("\n\n");
    for(int i = 0; i < n; i++)
        printf("%10.0f%5s%10.0f\n", a[i], "", b[i]);

    for(int i = 0; i < n; i++)
        a[i] = i*i;

    printf("\n\n");
    for(int i = 0; i < n; i++)
        printf("%10.0f%5s%10.0f\n", a[i], "", b[i]);



}
