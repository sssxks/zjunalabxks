/**
2 100 1
---
1.68
 */
#include <stdio.h>
#include <math.h>

double f0(double x, double l, double t)
{
    return sqrt(1.0 + l * l * t * t * cos(t * x) * cos(t * x));
}

double Integral(double a, double b, double (*f)(double x, double y, double z),
                double eps, double l, double t);

int main()
{
    double a = 0.0, b, eps = 0.005, l, t;

    scanf("%lf %lf %lf", &l, &b, &t);
    printf("%.2f\n", Integral(a, b, f0, eps, l, t));

    return 0;
}

#include <stdlib.h>

/**
 * @brief Calculate the integral of a function using the Romberg method.
 * 
 * @param a The lower limit of the integral.
 * @param b The upper limit of the integral.
 * @param f The function to be integrated.(double x, double y, double z), only x is used.
 * @param eps The precision of the result.
 * @param l a constant of the function.
 * @param t a constant of the function.
 */
double Integral(double a, double b, double (*f)(double x, double y, double z), double eps, double l, double t)
{
    double *r[2];
    r[0] = (double *)malloc(sizeof(double) * 1);
    r[1] = (double *)malloc(sizeof(double) * 2);

    double h = b - a;
    r[0][0] = (f(a, l, t) + f(b, l, t)) * h / 2.0;

    int count = 0;
    int pow_i = 1;
    int i = 1;
    do
    {
        // use trapzoidal rule
        double sum = 0.0;
        for (int j = 1; j <= pow_i; j++)
        {
            sum += f(a + (j - 0.5) * h, l, t);
        }
        r[1][0] = r[0][0] / 2.0 + sum * h / 2.0;

        // use Richardson extrapolation
        for (int j = 1; j <= i; j++)
        {
            r[1][j] = r[1][j - 1] + (r[1][j - 1] - r[0][j - 1]) / (pow(4, j) - 1);
        }

        // Check for desired precision
        if (fabs(r[1][i] - r[0][i - 1]) < eps)
        {
            double result = r[1][i];
            free(r[0]);
            free(r[1]);
            return result;
        }

        pow_i *= 2;
        i++;
        h /= 2.0;

        free(r[0]);
        r[0] = r[1];
        r[1] = (double *)malloc(sizeof(double) * (i + 1));
        count++;
    } while (count <= 10);

    double result = r[0][i - 1];
    free(r[0]);
    free(r[1]);
    return result;
}
