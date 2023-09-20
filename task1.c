#include <stdio.h>

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

void print_by_rows(double x[], double y[], int n){
    for (int i = 0; i < n + 1; i++)
    {
        printf("%.3f; %.3f\n ", x[i], y[i]);
    }
}

void polynom_lagrange(double x[], double y[], int n){
    
}

int main(){
    double a = -3.0;
    double b = 3.0;
    int n = 20;
    double step = (b - a) / n;

    double x[100];
    double y[100];

    for (int i = 0; i <= n; i++)
    {
        x[i] = a + step * i;
        y[i] = sin_of_square(x[i]);

    }

    //print_table(x, y, n);
    print_by_rows(x, y, n);
}