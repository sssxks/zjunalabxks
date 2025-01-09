/**
12 23.64 17.39 12.77 16.62 10.67 14.85 12.68 26.90 28.30 15.59 37.99 23.18
---
9.20 5.58 3.24 7.00 1.99 6.36 2.25 10.01 11.52 0.50 17.65 4.88 
 */
#include <stdio.h>

#define Max_size 10000 /* max number of dishes */

void Price( int n, double p[] );

int main()
{
    int n, i;
    double p[Max_size];

    scanf("%d", &n);
    for (i=0; i<n; i++) 
        scanf("%lf", &p[i]);
    Price(n, p);
    for (i=0; i<n; i++)
        printf("%.2f ", p[i]);
    printf("\n");

    return 0;
}

/* Your function will be put here */

#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <stdint.h>
#include <string.h>
#include <stdbool.h>

#define PI 3.14159265358979323846

void solve_circulant(complex double *c, complex double *b, complex double *x, int n);
void dft(complex double *x, int n, complex double *X);
void idft(complex double *X, int n, complex double *x);


#ifndef nullptr
#define nullptr 0
#endif


bool is_pow_2(size_t n) {
    return (n & (n - 1)) == 0;
}

size_t next_pow_2(size_t x) {
    size_t h = 1;
    while (h < x)
        h *= 2;
    return h;
}

void swap_complex(float *a_re, float *a_im, float *b_re, float *b_im) {
    const float p_re = *a_re;
    const float p_im = *a_im;
    *a_re = *b_re;
    *a_im = *b_im;
    *b_re = p_re;
    *b_im = p_im;
}

static void bitrev(float *in_re, float *in_im, const size_t n) {
    for (size_t i = 0, j = 0; i < n - 1; i++) {
        if (i < j)
            swap_complex(&in_re[i], &in_im[i], &in_re[j], &in_im[j]);
        size_t m = n;
        do {
            m >>= 1;
            j ^= m;
        } while (!(j & m));
    }
}

static void fft_pow_2(float *in_re, float *in_im, const size_t n, int sign) {
    bitrev(in_re, in_im, n);
    float pi = (float) (sign * PI);
    for (size_t i = 1; i < n; i = 2 * i) {
        const float delta = pi / i;
        for (int j = 0; j < i; j++) {
            const float m = delta * j;
            const float w_re = cosf(m);
            const float w_im = sinf(m);
            for (int k = j; k < n; k += 2 * i) {
                const float v_re = in_re[k];
                const float v_im = in_im[k];
                const float x_re = in_re[k + i] * w_re - in_im[k + i] * w_im;
                const float x_im = in_re[k + i] * w_im + in_im[k + i] * w_re;
                in_re[k] = v_re + x_re;
                in_im[k] = v_im + x_im;
                in_re[k + i] = v_re - x_re;
                in_im[k + i] = v_im - x_im;
            }
        }
    }
}

void bluestein(float *in_re, float *in_im, const size_t n) {
    size_t n2 = 2 * n;
    const size_t nb = next_pow_2(n2);
    float *w_re = (float *) malloc(sizeof(float) * n);
    float *w_im = (float *) malloc(sizeof(float) * n);
    float *y_re = (float *) calloc(sizeof(float), nb);
    float *y_im = (float *) calloc(sizeof(float), nb);
    float *b_re = (float *) calloc(sizeof(float), nb);
    float *b_im = (float *) calloc(sizeof(float), nb);
    if (w_re == nullptr || w_im == nullptr ||
        y_re == nullptr || y_im == nullptr ||
        b_re == nullptr || b_im == nullptr) {
        if (w_re)free(w_re);
        if (w_im)free(w_im);
        if (y_re)free(y_re);
        if (y_im)free(y_im);
        if (b_re)free(b_re);
        if (b_im)free(b_im);
        return;
    }
    const float delta = (float) PI / n;
    w_re[0] = 1;
    w_im[0] = 0;
    y_re[0] = 1;
    y_im[0] = 0;
    for (int k = 1; k < n; k++) {
        const float m = delta * (int) ((k * k) % n2);
        w_re[k] = cosf(m);
        w_im[k] = sinf(m);
        y_re[k] = w_re[k];
        y_im[k] = w_im[k];
        y_re[nb - k] = w_re[k];
        y_im[nb - k] = w_im[k];
    }
    fft_pow_2(y_re, y_im, nb, -1);
    for (int i = 0; i < n; i++) {
        b_re[i] = w_re[i] * in_re[i] + w_im[i] * in_im[i];
        b_im[i] = w_re[i] * in_im[i] - w_im[i] * in_re[i];
    }
    fft_pow_2(b_re, b_im, nb, -1);
    for (int i = 0; i < nb; i++) {
        const float t = b_re[i];
        b_re[i] = t * y_re[i] - b_im[i] * y_im[i];
        b_im[i] = t * y_im[i] + b_im[i] * y_re[i];
    }
    fft_pow_2(b_re, b_im, nb, 1);
    for (int i = 0; i < n; i++) {
        in_re[i] = (w_re[i] * b_re[i] + w_im[i] * b_im[i]) / nb;
        in_im[i] = (w_re[i] * b_im[i] - w_im[i] * b_re[i]) / nb;
    }
    free(w_re);
    free(w_im);
    free(y_re);
    free(y_im);
    free(b_re);
    free(b_im);
}

static void fft(float *in_re, float *in_im, size_t n, int sign) {
    if (is_pow_2(n)) {
        fft_pow_2(in_re, in_im, n, sign);
    } else {
        if (sign == 1) {
            bluestein(in_im, in_re, n);
        } else {
            bluestein(in_re, in_im, n);
        }
    }
}

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
        if (cimag(c[i]) > 1e-3)
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
    float *in_re = (float *)malloc(n * sizeof(float));
    float *in_im = (float *)malloc(n * sizeof(float));
    float *out_re = (float *)malloc(n * sizeof(float));
    float *out_im = (float *)malloc(n * sizeof(float));

    for (int i = 0; i < n; i++) {
        in_re[i] = creal(x[i]);
        in_im[i] = cimag(x[i]);
    }

    fft(in_re, in_im, n, 1);

    for (int i = 0; i < n; i++) {
        X[i] = in_re[i] + I * in_im[i];
    }

    free(in_re);
    free(in_im);
    free(out_re);
    free(out_im);
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