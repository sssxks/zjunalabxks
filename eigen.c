#include <math.h>

/**
 * @brief Computes the infinity norm (maximum absolute value) of a vector.
 * 
 * @param v The vector whose infinity norm is to be computed.
 * @param n The number of elements in the vector.
 * @return The infinity norm of the vector.
 */
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

/**
 * @brief Finds the index of the maximum absolute value in a vector.
 * 
 * @param v The vector to search.
 * @param N The number of elements in the vector.
 * @return The index of the maximum absolute value.
 */
int max_index(double* v, int N)
{
    double max = infinity_norm(v, N);
    int index = 0;
    for (int i = 0; i < N; i++)
    {
        if (fabs(v[i]) == max)
        {
            index = i;
            break;
        }
    }
    return index;
}

/**
 * @brief Performs LU decomposition and solves the system of linear equations.
 * 
 * @param n The number of equations.
 * @param a The coefficient matrix.
 * @param x The solution vector.
 * @param b The right-hand side vector.
 * @return 0 if successful, -1 if the matrix is singular.
 */
int lu(int n, double a[][MAX_SIZE], double x[], double b[])
{
    double l[MAX_SIZE][MAX_SIZE];
    double u[MAX_SIZE][MAX_SIZE];
    double y[MAX_SIZE] = {0};
    for (int i = 0; i < n; i++)
    {
        l[i][i] = 1;
        for (int j = i; j < n; j++)
        {
            double sum = 0;
            for (int k = 0; k < i; k++)
            {
                sum += l[i][k] * u[k][j];
            }
            u[i][j] = a[i][j] - sum;
        }
        for (int j = i + 1; j < n; j++)
        {
            double sum = 0;
            for (int k = 0; k < i; k++)
            {
                sum += l[j][k] * u[k][i];
            }
            if (u[i][i] == 0) return -1;
            l[j][i] = (a[j][i] - sum) / u[i][i];
        }
    }
    y[0] = b[0];
    for (int i = 1; i < n; i++)
    {
        double sum = 0;
        for (int j = 0; j < i; j++)
        {
            sum += l[i][j] * y[j];
        }
        y[i] = b[i] - sum;
    }
    for (int i = n - 1; i >= 0; i--)
    {
        if (u[i][i] == 0) return -1;
        double sum = 0;
        for (int j = i+1; j < n; j++)
        {
            sum += u[i][j] * x[j];
        }
        x[i] = (y[i] - sum) / u[i][i];
    }
    return 0;
}

/**
 * @brief Computes the eigenvalue and eigenvector of a matrix using the inverse power method.
 * 
 * @param n The size of the matrix.
 * @param a The matrix whose eigenvalue and eigenvector are to be computed.
 * @param lambda The initial guess for the eigenvalue, updated with the computed eigenvalue.
 * @param x0 The initial guess for the eigenvector, updated with the computed eigenvector.
 * @param TOL The tolerance for convergence.
 * @param MAXN The maximum number of iterations.
 * @return 1 if successful, 0 if the method did not converge, -1 if the matrix is singular.
 */
int EigenV(int n, double a[][MAX_SIZE], double* lambda, double x0[], double TOL, int MAXN)
{
    double y[MAX_SIZE] = { 0 };
    double a2[MAX_SIZE][MAX_SIZE] = { 0 };
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            a2[i][j] =(i==j)? a[i][i] - *lambda:a[i][j];

    int k = 1;
    int p = max_index(x0, n);
    double max = infinity_norm(x0, n);
    for (int i = 0; i < n; i++) x0[i] /= max;

    while (k <= MAXN)
    {
        if (lu(n, a2,  y, x0) == -1)
            return -1;
        double u = y[p];
        p = max_index(y, n);
        
        double error[MAX_SIZE];
        for(int i=0;i<n;i++)
            error[i] = x0[i]-y[i]/y[p];
        
        for (int i = 0; i < n; i++) x0[i] = y[i] / y[p];

        if (infinity_norm(error, n) < TOL)
        {
            *lambda = 1.0 / u + *lambda;
            return 1;
        }
        k++;
    }
    return 0;
}
