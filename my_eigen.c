/**
3
1 2 3
2 3 4
3 4 5
0.0000000001 1000
1
-0.6 1 1 1
---
-0.62347538
  1.00000000   0.17206558  -0.65586885
===
2
0 1
1 0
0.0000000001 10
2
1.0 1 1
100.0 1 0
---
  1.00000000 is an eigenvalue.
Maximum number of iterations exceeded.
 */
#include <stdio.h>

#define MAX_SIZE 10

int EigenV(int n, double a[][MAX_SIZE], double *lambda, double v[], double TOL, int MAXN);

int main()
{
    int n, MAXN, m, i, j, k;
    double a[MAX_SIZE][MAX_SIZE], v[MAX_SIZE];
    double lambda, TOL;

    scanf("%d", &n);
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            scanf("%lf", &a[i][j]);
    scanf("%lf %d", &TOL, &MAXN);
    scanf("%d", &m);
    for (i = 0; i < m; i++)
    {
        scanf("%lf", &lambda);
        for (j = 0; j < n; j++)
            scanf("%lf", &v[j]);
        switch (EigenV(n, a, &lambda, v, TOL, MAXN))
        {
        case -1:
            printf("%12.8f is an eigenvalue.\n", lambda);
            break;
        case 0:
            printf("Maximum number of iterations exceeded.\n");
            break;
        case 1:
            printf("%12.8f\n", lambda);
            for (k = 0; k < n; k++)
                printf("%12.8f ", v[k]);
            printf("\n");
            break;
        }
    }

    return 0;
}

/* Your function will be put here */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

const double EPS = 1e-12;

double infinity_norm(double v[], int n)
{
    double max_abs = 0.0;
    for (int i = 0; i < n; i++)
    {
        if (fabs(v[i]) > max_abs)
            max_abs = fabs(v[i]);
    }
    return max_abs;
}

void matrix_vector_product(double a[][MAX_SIZE], double v[], double result[], int n)
{
    for (int i = 0; i < n; i++)
    {
        result[i] = 0.0;
        for (int j = 0; j < n; j++)
        {
            result[i] += a[i][j] * v[j];
        }
    }
}

double dot_product(double u[], double v[], int n)
{
    double sum = 0.0;
    for (int i = 0; i < n; i++)
    {
        sum += u[i] * v[i];
    }
    return sum;
}

void form_M(double a[][MAX_SIZE], double m[][MAX_SIZE], double lambda, int n)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            m[i][j] = a[i][j] - (i == j ? lambda : 0.0);
        }
    }
}

int lu_decompose(double m[][MAX_SIZE], int n, int pivot[])
{
    for (int i = 0; i < n; i++)
        pivot[i] = i;

    for (int i = 0; i < n; i++)
    {
        double max = fabs(m[i][i]);
        int max_row = i;
        for (int k = i + 1; k < n; k++)
        {
            if (fabs(m[k][i]) > max)
            {
                max = fabs(m[k][i]);
                max_row = k;
            }
        }

        if (max < EPS)
            return -1;

        if (max_row != i)
        {
            int temp = pivot[i];
            pivot[i] = pivot[max_row];
            pivot[max_row] = temp;

            for (int k = 0; k < n; k++)
            {
                double tmp = m[i][k];
                m[i][k] = m[max_row][k];
                m[max_row][k] = tmp;
            }
        }

        for (int k = i + 1; k < n; k++)
        {
            m[k][i] /= m[i][i];
            for (int j = i + 1; j < n; j++)
            {
                m[k][j] -= m[k][i] * m[i][j];
            }
        }
    }
    return 0;
}

void lu_solve(double m[][MAX_SIZE], int n, int pivot[], double b[], double x[])
{
    double y[MAX_SIZE];
    for (int i = 0; i < n; i++)
    {
        y[i] = b[pivot[i]];
        for (int j = 0; j < i; j++)
        {
            y[i] -= m[i][j] * y[j];
        }
    }

    for (int i = n - 1; i >= 0; i--)
    {
        x[i] = y[i];
        for (int j = i + 1; j < n; j++)
        {
            x[i] -= m[i][j] * x[j];
        }
        x[i] /= m[i][i];
    }
}

/**
 * Computes the eigenvalue and eigenvector of a given matrix.
 *
 * @param n The order of the matrix.
 * @param a The matrix for which the eigenvalue and eigenvector are to be computed. It is a 2D array of size n x n.
 * @param lambda an initial approximation p of the eigenvalue and is supposed to be returned as a more accurate eigenvalue
 * @param v initial approximation of the eigenvector and is supposed to be returned as a more accurate eigenvector with unit norm
 * @param TOL The tolerance level for the convergence of the algorithm.
 * @param MAXN The maximum number of iterations allowed for the algorithm to converge.
 * @return An integer indicating the success or failure of the computation. Typically, 0 for success and non-zero for failure.
 */
int EigenV(int n, double a[][MAX_SIZE], double *lambda, double v[], double TOL, int MAXN)
{
    for (int iteration = 0; iteration < MAXN; iteration++)
    {
        double M[MAX_SIZE][MAX_SIZE];
        form_M(a, M, *lambda, n);

        // solve Mw = v
        int pivot[n];
        double w[MAX_SIZE];
        if (lu_decompose(M, n, pivot))
        {
            return -1;
        }
        lu_solve(M, n, pivot, v, w);

        double norm = infinity_norm(w, n);
        if (norm == 0)
        {
            return 0;
        }
        double diff[MAX_SIZE];
        for (int i = 0; i < n; i++)
        {
            double new_v = w[i] / norm;
            diff[i] = new_v - v[i];
            v[i] = new_v;
        }

        double A_times_v[MAX_SIZE];
        matrix_vector_product(a, v, A_times_v, n);
        *lambda = dot_product(v, A_times_v, n);

        if (infinity_norm(diff, n) < TOL)
        {
            return 1;
        }
    }

    return 0;
}