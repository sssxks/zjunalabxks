#include <stdlib.h>
#include <math.h>
#include <complex.h>

#define PI 3.14159265358979323846

void solve_circulant(complex double *c, complex double *b, complex double *x, int n);
void dft(complex double *x, int n, complex double *X);
void idft(complex double *X, int n, complex double *x);

void Price(int n, double p[])
{
    // Allocate memory for padded input
    complex double *p_complex = calloc((size_t)n, sizeof(complex double));
    for (int i = 0; i < n; i++)
    {
        p_complex[i] = p[i];
    }

    // Allocate memory for circulant matrix
    complex double *C = calloc((size_t)n, sizeof(complex double));
    C[0] = 2 + 0 * I;
    C[1] = 0.5 + 0 * I;
    C[n - 1] = 0.5 + 0 * I;

    // Allocate memory for solution vector
    complex double *c = calloc((size_t)n, sizeof(complex double));

    // Solve the circulant matrix system
    solve_circulant(C, p_complex, c, n);

    // Copy the solution back to the original array
    for (int i = 0; i < n; i++)
    {
        if (cimag(c[i]) > 1e-9)
        {
            fprintf(stderr, "Warning: Non-negligible imaginary part in solution.\n");
        }

        p[i] = creal(c[i]);
    }

    // Free allocated memory
    free(p_complex);
    free(C);
    free(c);
}

// Function to compute the FFT
void dft(complex double *x, int n, complex double *X)
{
    if (n <= 1)
        return;

    for (int k = 0; k < n; k++)
    {
        X[k] = 0;
        for (int t = 0; t < n; t++)
        {
            X[k] += x[t] * cexp(-2 * I * PI * k * t / n);
        }
    }
}

// Function to compute the Inverse FFT
void idft(complex double *X, int n, complex double *x)
{
    complex double *conj_X = malloc((size_t)n * sizeof(complex double));

    // Conjugate the input
    for (int i = 0; i < n; i++)
    {
        conj_X[i] = conj(X[i]);
    }

    // Compute the FFT
    dft(conj_X, n, x);

    // Conjugate the result and scale by 1/n
    for (int i = 0; i < n; i++)
    {
        x[i] = conj(x[i]) / n;
    }

    free(conj_X);
}

// Function to solve the circulant matrix system
void solve_circulant(complex double *c, complex double *b, complex double *x, int n)
{
    // Step 1: Compute the FFT of the first row of the circulant matrix
    complex double *lambda = malloc((size_t)n * sizeof(complex double));
    dft(c, n, lambda);

    // Step 2: Compute the FFT of the right-hand side vector
    complex double *b_hat = malloc((size_t)n * sizeof(complex double));
    for (int i = 0; i < n; i++)
    {
        b_hat[i] = b[i];
    }
    dft(b, n, b_hat);

    // Step 3: Solve for x_hat in the frequency domain
    complex double *x_hat = malloc((size_t)n * sizeof(complex double));
    for (int i = 0; i < n; i++)
    {
        if (cabs(lambda[i]) > 1e-9)
        { // Avoid division by zero
            x_hat[i] = b_hat[i] / lambda[i];
        }
        else
        {
            x_hat[i] = 0; // Handle singular case
        }
    }

    // Step 4: Compute the Inverse FFT to get the solution in the time domain
    for (int i = 0; i < n; i++)
    {
        x[i] = x_hat[i];
    }
    idft(x_hat, n, x);

    // Free allocated memory
    free(lambda);
    free(b_hat);
    free(x_hat);
}