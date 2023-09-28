#include <stdio.h>
#include <stdlib.h>

#include "printing.h"
#include "lagrange.h"

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


double divided_difference(double* x, double* y, int j, int k, int n){

    if(k == 0){
        return y[j];
    }

    return (divided_difference(x, y, j + 1, k - 1, n) - divided_difference(x, y, j, k - 1, n)) / (x[j + k] - x[j]);
}

double polynom_newton(double* x, double* y, int n, double x_point){
    double result = 0;

    for (int i = 0; i < n; i++)
    {
        double product = 1;
        for (int j = 0; j < i; j++)
        {
            product *= x_point - x[j];
        }
        product *= divided_difference(x, y, 0, i, n);

        result += product;
    }

    return result;   
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

    // printf("Polynom: \n");


    FILE* file = fopen("output_lag.dat", "w");

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