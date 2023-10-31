#include <stdio.h>
#include <stdlib.h>

double divided_difference(double* x, double* y, int j, int k, int n){

    if(k == 0){
        return y[j];
    }

    return (divided_difference(x, y, j + 1, k - 1, n) - divided_difference(x, y, j, k - 1, n)) / (x[j + k] - x[j]);
}

double** divided_difference_table(double* x, double* y, int n){
    double **table = (double**) malloc(sizeof(double*) * n);

    for (int i = 0; i < n; i++)
    {
        table[i] = (double*) malloc(sizeof(double) * (n + 1));
    }

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n+1; j++)
        {
            table[i][j] = 0;
        }
    }       
    
    for (int i = 0; i < n; i++)
    {
        table[i][0] = x[i];
        table[i][1] = y[i];
    }

    for (int j = 2; j < n; j++)
    {
        for (int i = 0; i < n - j + 1; i++)
        {
            table[i][j] = (table[i + 1][j - 1] - table[i][j - 1]) / (table[i + j - 1][0] - table[i][0]);
        }
    }

    // for (int i = 0; i < n; i++)
    // {
    //     for (int j = 0; j < n + 1; j++)
    //     {
    //         printf("%.3f ", table[i][j]);
    //     }
    //     printf("\n");
    // }
    
    return table;
}


// double polynom_newton(double* x, double* y, int n, double x_point){
//     double result = 0;

//     for (int i = 0; i < n; i++)
//     {
//         double product = 1;
//         for (int j = 0; j < i; j++)
//         {
//             product *= x_point - x[j];
//         }
//         product *= divided_difference(x, y, 0, i, n);

//         result += product;
//     }

//     return result;   
// }

double polynom_newton(double* x, double* y, int n, double x_point){
    double result = 0;
    double ** div_table = divided_difference_table(x, y, n);

    for (int i = 0; i < n; i++)
    {
        double product = 1;
        for (int j = 0; j < i; j++)
        {
            product *= x_point - x[j];
        }
        //product *= divided_difference(x, y, 0, i, n);
        product *= div_table[0][i + 1];

        result += product;
    }

    for (int i = 0; i < n; i++)
    {
        free(div_table[i]);
    }
    free(div_table);

    return result;   
}
