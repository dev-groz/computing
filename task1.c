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

    for (int i = 0; i <= n; i++)
    {
        double pos = a + step * i;
        printf("%f; %f\n", pos, sin_of_square(pos));

    }
}