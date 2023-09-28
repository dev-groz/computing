#include <stdio.h>
#include <stdlib.h>

#include "printing.h"
#include "lagrange.h"

void print_table(double x[], double y[], int n){
    for (int i = 0; i < n + 1; i++)
    {
        printf("\t%f", x[i]);
    }
    printf("\n");
    for (int i = 0; i < n + 1; i++)
    {
        printf("\t%f", y[i]);
    }
}

void print_by_rows(double* x, double* y, int n){
    for (int i = 0; i < n + 1; i++)
    {
        printf("%f; %f\n", x[i], y[i]);
    }
}
