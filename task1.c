#include <stdio.h>

// int fact(int n){
//     double result = 1;
//     for (int i = 1; i <= n; i++)
//     {
//         result *= i;
//     }
    
//     return result;
// }

// double mypow(double x, int n){
//     double result = 1;
//     for (int i = 0; i < n; i++)
//     {
//         result *= x;
//     }
    
//     return result;
// }

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

    double x[100];
    double y[100];

    for (int i = 0; i <= n; i++)
    {
        x[i] = a + step * i;
        y[i] = sin_of_square(x[i]);

    }

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