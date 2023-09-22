#include <stdio.h>
#include <stdlib.h>

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


void print_table(double x[], double y[], int n){
    for (int i = 0; i < n + 1; i++)
    {
        printf("\t%.3f", x[i]);
    }
    printf("\n");
    for (int i = 0; i < n + 1; i++)
    {
        printf("\t%.3f", y[i]);
    }
}

void print_by_rows(double* x, double* y, int n){
    for (int i = 0; i < n + 1; i++)
    {
        printf("%.3f; %.3f\n", x[i], y[i]);
    }
}

double polynom_lagrange(double* x, double* y, int n, double x_point){
    
    double result = 0;

    for (int i = 0; i < n + 1; i++){
        double product = 1;

        for (int j = 0; j < n + 1; j++)
        {
            if (i == j) continue;
            double numerator = (x_point - x[j]);
            double denominator = (x[i] - x[j]);

            product *= (numerator / denominator);
        }
        product *= y[i];
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

    for (int i = 0; i < n + 1; i++)
    {
        x[i] = a + step * i;
        y[i] = sin_of_square(x[i]);

    }

    printf("Lagrange: \n");

    for (int i = 0; i < n + 1; i++)
    {
        double result = polynom_lagrange(x, y, n, x[i]);
        printf("%f; %f\n", x[i], result);
    }
    
    
    printf("\nTaylor values: \n");

    // print_table(x, y, n);
    print_by_rows(x, y, n);
}