/**
12 23.64 17.39 12.77 16.62 10.67 14.85 12.68 26.90 28.30 15.59 37.99 23.18
---
9.20 5.58 3.24 7.00 1.99 6.36 2.25 10.01 11.52 0.50 17.65 4.88
 */
#include <stdio.h>

#define Max_size 10000 /* max number of dishes */

void Price(int n, double p[]);

int main()
{
    int n, i;
    double p[Max_size];

    scanf("%d", &n);
    for (i = 0; i < n; i++)
        scanf("%lf", &p[i]);
    Price(n, p);
    for (i = 0; i < n; i++)
        printf("%.2f ", p[i]);
    printf("\n");

    return 0;
}


#include <math.h>
#define bound pow(2, 127)
#define ZERO 1e-9 /* X is considered to be 0 if |X|<ZERO */

int Jacobi(int n, double a[][Max_size], double b[], double x[], double TOL, int MAXN)
{
    double local_a[Max_size][Max_size];
    double local_b[Max_size];
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
            local_a[i][j] = a[i][j];
        local_b[i] = b[i];
    }

    for (int i = 0; i < n; i++)
    {
        int max_row = i;
        double max_val = fabs(local_a[i][i]);
        for (int k = i + 1; k < n; k++)
        {
            if (fabs(local_a[k][i]) > max_val)
            {
                max_val = fabs(local_a[k][i]);
                max_row = k;
            }
        }
        if (max_val < ZERO)
        {
            max_val = 0.0;
            max_row = -1;
            for (int k = 0; k < i; k++)
            {
                if (fabs(local_a[k][i]) > max_val)
                {
                    max_val = fabs(local_a[k][i]);
                    max_row = k;
                }
            }
            if (max_val < ZERO)
            {
                return -1; // Zero column
            }
            else
            {
                for (int j = 0; j < n; j++)
                    local_a[i][j] += local_a[max_row][j];
                local_b[i] += local_b[max_row];
            }
        }
        else if (max_row != i)
        {
            for (int j = 0; j < n; j++)
            {
                double temp = local_a[i][j];
                local_a[i][j] = local_a[max_row][j];
                local_a[max_row][j] = temp;
            }
            double temp_b = local_b[i];
            local_b[i] = local_b[max_row];
            local_b[max_row] = temp_b;
        }
    }

    double x_old[Max_size];
    for (int i = 0; i < n; i++)
        x_old[i] = x[i];

    for (int k = 1; k <= MAXN; k++)
    {
        double x_new[Max_size];
        for (int i = 0; i < n; i++)
        {
            double sum = 0.0;
            for (int j = 0; j < n; j++)
            {
                if (j != i)
                    sum += local_a[i][j] * x_old[j];
            }
            x_new[i] = (local_b[i] - sum) / local_a[i][i];
        }

        double max_diff = 0.0;
        for (int i = 0; i < n; i++)
        {
            double diff = fabs(x_new[i] - x_old[i]);
            if (diff > max_diff)
                max_diff = diff;
        }
        if (max_diff < TOL)
        {
            for (int i = 0; i < n; i++)
                x[i] = x_new[i];
            return k;
        }

        int divergence = 0;
        for (int i = 0; i < n; i++)
        {
            if (x_new[i] < -bound || x_new[i] > bound)
            {
                divergence = 1;
                break;
            }
        }
        if (divergence)
            return -2;

        for (int i = 0; i < n; i++)
            x_old[i] = x_new[i];
    }

    return 0;
}

void Price(int n, double p[])
{
    double a[Max_size][Max_size] = {0};
    double b[Max_size] = {0};
    double x[Max_size] = {0};

    a[0][n - 1] = 0.5;
    a[n - 1][0] = 0.5;
    for (int i = 0; i < n; i++)
    {
        a[i][i] = 2;
        if (i > 0)
        {
            a[i][i - 1] = 0.5;
        }
        if (i < n - 1)
        {
            a[i][i + 1] = 0.5;
        }
    }

    for (int i = 0; i < n; i++)
    {
        b[i] = p[i];
    }

    // Call Jacobi function
    int iterations = Jacobi(n, a, b, x, 1e-8, 1000);
    if (iterations < 0)
    {
        fprintf(stderr, "Jacobi method failed to converge.\n");
    }

    for (int i = 0; i < n; i++)
    {
        p[i] = x[i];
    }
}
