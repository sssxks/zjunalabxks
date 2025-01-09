/**
2 1 1.5 -1
0 2
---
0.5000
 */
#include <stdio.h>
#include <math.h>

#define ZERO 1e-13 /* X is considered to be 0 if |X|<ZERO */
#define MAXN 11    /* Max Polynomial Degree + 1 */

double Polynomial_Root(int n, double c[], double a, double b, double EPS);

int main()
{
    int n;
    double c[MAXN], a, b;
    double EPS = 0.00005;
    int i;

    scanf("%d", &n);
    for (i = n; i >= 0; i--)
        scanf("%lf", &c[i]);
    scanf("%lf %lf", &a, &b);
    printf("%.4f\n", Polynomial_Root(n, c, a, b, EPS));

    return 0;
}

/* Your function will be put here */


double f(double x, double c[], int n);
double f_prime(double x, double c[], int n);
double f_double_prime(double x, double c[], int n);
double Polynomial_Root(int n, double c[], double a, double b, double EPS);


double f(double x, double c[], int n)
{
    double ret = c[(int)n];
    for (int i = n - 1; i >= 0; i--)
    {
        ret = ret * x + c[i];
    }
    return ret;
}

double f_prime(double x, double c[], int n)
{
    double ret = c[(int)n] * n;
    for (int i = n - 1; i >= 1; i--)
    {
        ret = ret * x + c[i] * i;
    }
    return ret;
}

double f_double_prime(double x, double c[], int n)
{
    double ret = c[(int)n] * n * (n - 1);
    for (int i = n - 1; i >= 2; i--)
    {
        ret = ret * x + c[i] * i * (i - 1);
    }
    return ret;
}

double Polynomial_Root(int n, double c[], double a, double b, double EPS)
{
    const int N = 128;
    const int MAX_ITER = 1024;

    if (a > b)
    {
        double temp = a;
        a = b;
        b = temp;
    }

    for (double i = 0; i < N; i++)
    {
        double x = a + (b - a) * (double)i / (double)N;
        for (int iter = 0; iter < MAX_ITER; iter++)
        {
            double x_last = x;

            double f_val = f(x, c, n);
            double f_prime_val = f_prime(x, c, n);
            double f_double_prime_val = f_double_prime(x, c, n);
            x = x - (f_prime_val * f_val) / (f_prime_val * f_prime_val - f_val * f_double_prime_val);
            if (x < a || x > b)
            {
                break;
            }

            if (fabs(x - x_last) < EPS)
            {
                return x;
            }
        }
    }
}