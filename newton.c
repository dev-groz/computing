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
