#include <stdio.h>
#include <stdlib.h>

#include "printing.h"
#include "lagrange.h"
#include "newton.h"

double mysin(double x, int n){
    double result = 0;
    double elem = x;

    result += elem;

    for (int i = 1; i <= n; i++){
        elem = elem * (-1) * (x * x) / ((2 * i + 1) * (2 * i));
        result += elem;
    }

    return result;
}

double sin_of_square(double x){
    return mysin(x * x, 20);
}

int main(){
    double a = -3.0;
    double b = 3.0;
    int n = 20;
    double step = (b - a) / n;

    double* x = (double*)malloc(sizeof(double) * (n + 1));
    double* y = (double*)malloc(sizeof(double) * (n + 1));

    int new_n = 500;

    double* new_x = (double*)malloc(sizeof(double) * (new_n + 1));

    for (int i = 0; i < n + 1; i++)
    {
        x[i] = a + step * i;
        y[i] = sin_of_square(x[i]);
    }

    double new_step = (b - a) / new_n;

    for (int i = 0; i < new_n + 1; i++)
    {
        new_x[i] = a + new_step * i;
    }


    // FILE* file = fopen("output_lag.dat", "w");
    FILE* file = fopen("output_new.dat", "w");

    if (file == NULL){
        exit(1);
    }

    for (int i = 0; i < new_n + 1; i++)
    {
        // double result = polynom_lagrange(x, y, n, new_x[i]);
        double result = polynom_newton(x, y, n, new_x[i]);
        fprintf(file, "%f\t %f\n", new_x[i], result);
    }

    fclose(file);


    // for (int i = 0; i < new_n + 1; i++)
    // {
    // //     double result = polynom_lagrange(x, y, n, new_x[i]);
    //     double result = polynom_newton(x, y, n, new_x[i]);
    //     fprintf(file1, "%f\t %f\n", new_x[i], result);
    // }

    // printf("Newton: \n");    

    //printf("\nTaylor values: \n");

    // print_table(x, y, n);
    //print_by_rows(x, y, n);
    //output_in_file(x, y, n, "output.dat");

    free(x);
    free(y);
    free(new_x);
}