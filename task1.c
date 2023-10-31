#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "printing.h"
#include "lagrange.h"
#include "newton.h"

#define MAX_INT 2147483647
#define PI 3.14159265358979323846

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

int factorial(int n){

    int result = 1;

    for (int i = 1; i <= n; i++){
        result *= i;
    }

    if (result < 0)
        return MAX_INT;

    return result;
}

double sin_of_square(double x){
    return mysin(x * x, 20);
}

double derivative(double xi, int n){
    if (n == 1){
        return 2 * xi;
    }

    if (n == 2){
        return 2 + 2 * 2 * xi * xi;
    }

    return n * (n - 1) * derivative(xi, n - 2) + pow(2, n) * pow(xi, n);
}

double find_error(double* x, int n, double xi, double point){
    double result = derivative(xi, n + 1);

    double fact = factorial(n + 1);

    result /= fact;

    double omega = 1;
    for (int i = 0; i < n; i++)
    {
        omega *= point - x[i];
    }
    
    omega = omega > 0 ? omega : -omega;

    result *= omega;

    return result;
}

void output_in_file_polynom_lag(double* x, double* y, double* new_x, int n, int new_n, char* filename){
    FILE* file = fopen(filename, "w");

    if (file == NULL){
        exit(1);
    }

    for (int i = 0; i < new_n + 1; i++)
    {
        double result = polynom_lagrange(x, y, n, new_x[i]);
        fprintf(file, "%f\t %f\n", new_x[i], result);
    }

    fclose(file);
}

void output_in_file_polynom_new(double* x, double* y, double* new_x, int n, int new_n, char* filename){
    FILE* file = fopen(filename, "w");

    if (file == NULL){
        exit(1);
    }

    for (int i = 0; i < new_n + 1; i++)
    {
        double result = polynom_newton(x, y, n + 1, new_x[i]);
        fprintf(file, "%f\t %f\n", new_x[i], result);
    }

    fclose(file);
}

void output_real_error(double* x, double* y, int n, double* new_x, int new_n){

    FILE* file = fopen("outputs\\real_error.dat", "w");

    for (int i = 0; i < new_n + 1; i++)
    {
        double result = sin(new_x[i] * new_x[i]) - polynom_lagrange(x, y, n, new_x[i]);
        result = result > 0 ? result : -result;
        fprintf(file, "%f\t %f\n", new_x[i], result);
    }
    fclose(file);
}

void output_approximate_error(double* x, double* new_x, int new_n){
    FILE* file = fopen("outputs\\approximate_error.dat", "w");

    for (int i = 0; i < new_n; i++)
    {
        fprintf(file, "%f\t%f\n", new_x[i], find_error(x, 3, 3, new_x[i]));
    }

    fclose(file);
}

int main(){
    double a = -3.0;
    double b = 3.0;
    int n = 7;
    double step = (b - a) / n;

    double* x = (double*)malloc(sizeof(double) * (n + 1));
    double* y = (double*)malloc(sizeof(double) * (n + 1));

    int new_n = 300;

    double* new_x = (double*)malloc(sizeof(double) * (new_n + 1));

    for (int i = 0; i < n + 1; i++)
    {
        //equidistant values
        x[i] = a + step * i;
        //Chebyshev values
        //x[i] = 0.5 * (b + a) + 0.5 * (b - a) * cos((2 * i + 1) * PI /(2 * (n + 1))); 
        // y[i] = sin_of_square(x[i]);
        y[i] = pow(x[i], 10);
    }

    double new_step = (b - a) / new_n;

    for (int i = 0; i < new_n + 1; i++)
    {
        new_x[i] = a + new_step * i;
    }   

    output_in_file_polynom_lag(x, y, new_x, n, new_n, "outputs\\output_lag.dat");
    output_in_file_polynom_new(x, y, new_x, n, new_n, "outputs\\output_new.dat");

    // print_table(x, y, n);
    //print_by_rows(x, y, n);
    //output_in_file(x, y, n + 1, "outputs\\output.dat");

    //output_approximate_error(x, new_x, new_n);

    // output_real_error(x,y,n, new_x, new_n);

    free(x);
    free(y);
    free(new_x);
}
