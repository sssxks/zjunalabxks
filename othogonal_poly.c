/**
 * Orthogonal Polynomial Approximation (OPA) implementation
 */
#include <stdio.h>
#include <math.h>

#define MAX_m 200
#define MAX_n 5

double f1(double x)
{
    return sin(x);
}

double f2(double x)
{
    return exp(x);
}

int OPA(double (*f)(double t), int m, double x[], double w[], double c[], double *eps);

void print_results(int n, double c[], double eps)
{
    int i;

    printf("%d\n", n);
    for (i = 0; i <= n; i++)
        printf("%12.4e ", c[i]);
    printf("\n");
    printf("error = %9.2e\n", eps);
    printf("\n");
}

int main()
{
    int m, i, n;
    double x[MAX_m], w[MAX_m], c[MAX_n + 1], eps;

    m = 90;
    for (i = 0; i < m; i++)
    {
        x[i] = 3.1415926535897932 * (double)(i + 1) / 180.0;
        w[i] = 1.0;
    }
    eps = 0.001;
    n = OPA(f1, m, x, w, c, &eps);
    print_results(n, c, eps);

    m = 200;
    for (i = 0; i < m; i++)
    {
        x[i] = 0.01 * (double)i;
        w[i] = 1.0;
    }
    eps = 0.001;
    n = OPA(f2, m, x, w, c, &eps);
    print_results(n, c, eps);

    return 0;
}

// above is test code, do not modify it

#include <stdio.h>
#include <math.h>

/**
 * Computes the inner product of two functions.
 *
 * @param f A pointer to a function that takes a double and returns a double. ignored when flag is 0
 * @param m The number of points in the arrays x, w, p, and q.
 * @param x An array of double values representing the points at which the functions are evaluated.
 * @param w An array of double values representing the weights for the inner product computation.
 * @param p An array of double values representing the values of the first function at the points in x. ignored when flag is 1
 * @param q An array of double values representing the values of the second function at the points in x.
 * @param flag A flag to determine which function to compute the inner product with. 0 for p and q, 1 for f and q.
 * @return The computed inner product as a double.
 */
double inner_product(double (*f)(double), int m, double x[], double w[], double p[], double q[], int flag)
{
    double sum = 0.0;
    if (flag == 0)
    {
        for (int i = 0; i < m; i++)
        {
            double p_val = 0.0, q_val = 0.0;
            for (int k = 0; k <= MAX_n; k++)
            {
                p_val += p[k] * pow(x[i], k);
                q_val += q[k] * pow(x[i], k);
            }
            sum += w[i] * p_val * q_val;
        }
    }
    else
    {
        for (int i = 0; i < m; i++)
        {
            double q_val = 0.0;
            for (int k = 0; k <= MAX_n; k++)
            {
                q_val += q[k] * pow(x[i], k);
            }
            sum += w[i] * f(x[i]) * q_val;
        }
    }
    return sum;
}

/**
 * Multiplies a polynomial by the variable x.
 *
 * This function takes a polynomial represented by the array `p` and multiplies it by the variable x,
 * storing the result in the array `result`.
 *
 * @param p The input polynomial coefficients array. The array should be of size n, where n is the degree of the polynomial plus one.
 * @param result The output array where the resulting polynomial coefficients will be stored. This array should be of size n+1.
 */
void polynomial_multiply_x(double p[], double result[])
{
    for (int i = MAX_n; i > 0; i--)
    {
        result[i] = p[i - 1];
    }
    result[0] = 0.0;
}

/**
 * @brief Subtracts two polynomials represented by arrays.
 *
 * This function takes two arrays representing polynomials `p` and `q`,
 * and computes the polynomial subtraction `p - q`, storing the result
 * in the `result` array.
 *
 * @param p The first polynomial array.
 * @param q The second polynomial array.
 * @param result The array where the result of the subtraction will be stored.
 *
 * @note The arrays `p`, `q`, and `result` must have a size of at least `MAX_n + 1`.
 */
void polynomial_subtract(double p[], double q[], double result[])
{
    for (int i = 0; i <= MAX_n; i++)
    {
        result[i] = p[i] - q[i];
    }
}

/**
 * @brief Multiplies each coefficient of a polynomial by a constant.
 *
 * This function takes a polynomial represented by an array of coefficients `p[]`,
 * multiplies each coefficient by a given constant, and stores the result in the `result[]` array.
 *
 * @param p[] The array of coefficients representing the polynomial.
 * @param constant The constant value to multiply each coefficient by.
 * @param result[] The array where the resulting coefficients will be stored.
 *
 * @note The size of the arrays `p[]` and `result[]` should be at least `MAX_n + 1`.
 */
void polynomial_multiply_constant(double p[], double constant, double result[])
{
    for (int i = 0; i <= MAX_n; i++)
    {
        result[i] = p[i] * constant;
    }
}

/**
 * @brief Computes the error between a given function and its polynomial approximation.
 *
 * This function calculates the weighted sum of squared errors between the values of a given function `f`
 * and its polynomial approximation at specified points.
 *
 * @param f A pointer to the function to be approximated.
 * @param m The number of points at which the error is computed.
 * @param x An array of points at which the function and polynomial are evaluated.
 * @param w An array of weights corresponding to each point.
 * @param c An array of coefficients for the polynomial approximation.
 * @return The computed error as a double.
 */
double compute_error(double (*f)(double), int m, double x[], double w[], double c[])
{
    double err = 0.0;
    for (int i = 0; i < m; i++)
    {
        double p_val = 0.0;
        for (int k = 0; k <= MAX_n; k++)
        {
            p_val += c[k] * pow(x[i], k);
        }
        err += w[i] * pow(f(x[i]) - p_val, 2);
    }
    return err;
}

/**
 * @brief Computes the orthogonal polynomial approximation (OPA) of a given function.
 *
 * This function computes the orthogonal polynomial approximation of a given function `f` using
 * Legendre polynomials. The approximation is computed up to a maximum degree `MAX_n`.
 *
 * @param f The function to approximate.
 * @param m The number of sample points.
 * @param x An array of sample points.
 * @param w An array of weights corresponding to the sample points.
 * @param c An array to store the coefficients of the orthogonal polynomial approximation.
 * @param eps A pointer to the tolerance value. On return, it contains the final error.
 * @return The degree of the orthogonal polynomial approximation.
 */
int OPA(double (*f)(double), int m, double x[], double w[], double c[], double *eps)
{
    double phi[MAX_n + 1][MAX_n + 1] = {0};
    double a[MAX_n + 1] = {0};
    double B[MAX_n + 1] = {0};
    double C[MAX_n + 1] = {0};

    // Initialize phi_0(x) = 1
    phi[0][0] = 1.0;

    // Compute a_0
    a[0] = inner_product(f, m, x, w, phi[0], phi[0], 1) / inner_product(f, m, x, w, phi[0], phi[0], 0);

    // Initialize P_n(x) with a_0 * phi_0(x)
    for (int i = 0; i <= MAX_n; i++)
    {
        c[i] = a[0] * phi[0][i];
    }

    // Compute the initial error
    // If the error is already within tolerance, return n = 0
    double err = compute_error(f, m, x, w, c);
    if (err <= *eps)
    {
        *eps = err;
        return 0;
    }

    // Iterate to find the best n
    for (int n = 1; n <= MAX_n; n++)
    {
        double xMultiplied[MAX_n + 1] = {0};
        polynomial_multiply_x(phi[n - 1], xMultiplied);

        // Compute B_n
        double numerator = inner_product(f, m, x, w, xMultiplied, phi[n - 1], 0);
        double denominator = inner_product(f, m, x, w, phi[n - 1], phi[n - 1], 0);
        B[n] = numerator / denominator;

        // Compute C_n
        if (n > 1)
        {
            numerator = inner_product(f, m, x, w, phi[n - 1], phi[n - 1], 0);
            denominator = inner_product(f, m, x, w, phi[n - 2], phi[n - 2], 0);
            C[n] = numerator / denominator;
        }

        // Compute phi_n(x)
        double BMultiplied[MAX_n + 1] = {0};
        polynomial_multiply_constant(phi[n - 1], B[n], BMultiplied);

        double LastTerm[MAX_n + 1] = {0};
        polynomial_subtract(xMultiplied, BMultiplied, LastTerm);

        if (n > 1)
        {
            double SecondToLastTerm[MAX_n + 1] = {0};
            polynomial_multiply_constant(phi[n - 2], C[n], SecondToLastTerm);

            polynomial_subtract(LastTerm, SecondToLastTerm, phi[n]);
        }
        else
        {
            for (int i = 0; i <= MAX_n; i++)
            {
                phi[n][i] = LastTerm[i];
            }
        }

        // Compute a_n
        a[n] = inner_product(f, m, x, w, phi[n], phi[n], 1) / inner_product(f, m, x, w, phi[n], phi[n], 0);

        // Update P_n(x)
        for (int i = 0; i <= MAX_n; i++)
        {
            c[i] += a[n] * phi[n][i];
        }

        // Compute the new error
        // If the error is within tolerance, return n
        err = compute_error(f, m, x, w, c);
        if (err <= *eps)
        {
            *eps = err;
            return n;
        }
    }

    // If Max_n is reached, return Max_n
    *eps = err;
    return MAX_n;
}