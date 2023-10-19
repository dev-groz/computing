#include <stdio.h>

#include "lagrange.h"

double polynom_lagrange(double* x, double* y, int n, double x_point){
    
    double result = 0;

    for (int i = 0; i < n + 1; i++){
        double product = 1;

        double full_num = 1;
        double full_den = 1;

        for (int j = 0; j < n + 1; j++)
        {
            if (i == j) continue;
            double numerator = (x_point - x[j]);
            double denominator = (x[i] - x[j]);

            full_num *= numerator;
            full_den *= denominator;
            // product *= (numerator / denominator);
        }
        product *= full_num/full_den;
        product *= y[i];
        result += product;
    }

    return result;
}