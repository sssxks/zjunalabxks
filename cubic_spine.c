/**
2
0.0 1.0 2.0
0.0 1.0 2.0
1 1.0 1.0 0.0
0.0 3.0 2
---
0.00000000e+00 1.00000000e+00 0.00000000e+00 0.00000000e+00 
1.00000000e+00 1.00000000e+00 0.00000000e+00 0.00000000e+00 
f(0.00000000e+00) = 0.00000000e+00
f(1.50000000e+00) = 1.50000000e+00
f(3.00000000e+00) = 0.00000000e+00
 */
#include <stdio.h>

#define MAX_N 10

void Cubic_Spline(int n, double x[], double f[], int Type, double s0, double sn, double a[], double b[], double c[], double d[]);

double S( double t, double Fmax, int n, double x[], double a[], double b[], double c[], double d[] );

int main()
{
    int n, Type, m, i;
    double x[MAX_N], f[MAX_N], a[MAX_N], b[MAX_N], c[MAX_N], d[MAX_N];
    double s0, sn, Fmax, t0, tm, h, t;

    scanf("%d", &n);
    for (i=0; i<=n; i++) 
        scanf("%lf", &x[i]);
    for (i=0; i<=n; i++) 
        scanf("%lf", &f[i]);
    scanf("%d %lf %lf %lf", &Type, &s0, &sn, &Fmax);

    Cubic_Spline(n, x, f, Type, s0, sn, a, b, c, d);
    for (i=1; i<=n; i++)
        printf("%12.8e %12.8e %12.8e %12.8e \n", a[i], b[i], c[i], d[i]);

    scanf("%lf %lf %d", &t0, &tm, &m);
    h = (tm-t0)/(double)m;
    for (i=0; i<=m; i++) {
        t = t0+h*(double)i;
        printf("f(%12.8e) = %12.8e\n", t, S(t, Fmax, n, x, a, b, c, d));
    }

    return 0;
}

/**
 * @brief Computes the coefficients of a cubic spline interpolation.
 *
 * This function calculates the coefficients of a cubic spline interpolation
 * for a given set of data points. The spline can be either a natural spline
 * or a clamped spline, depending on the value of the Type parameter.
 *
 * @param n The number of data points minus one.
 * @param x An array of size n containing the x-coordinates of the data points.
 * @param f An array of size n containing the y-coordinates of the data points.
 * @param Type An integer indicating the type of spline:
 *             - 0 for a natural spline (second derivatives at endpoints are zero).
 *             - 1 for a clamped spline (first derivatives at endpoints are specified).
 * @param s0 The first derivative at the first data point (used only if Type is 1).
 * @param sn The first derivative at the last data point (used only if Type is 1).
 * @param a A 1-based array of size n-1 to store the coefficient a of each spline segment.
 * @param b A 1-based array of size n-1 to store the coefficient b of each spline segment.
 * @param c A 1-based array of size n-1 to store the coefficient c of each spline segment.
 * @param d A 1-based array of size n-1 to store the coefficient d of each spline segment.
 */
void Cubic_Spline(int n, double x[], double f[], int Type, double s0, double sn, double a[], double b[], double c[], double d[]) {
    double h[MAX_N], alpha[MAX_N], l[MAX_N], mu[MAX_N], z[MAX_N];

    // Step 1: Calculate h_i
    for (int i = 0; i < n; i++) {
        h[i] = x[i+1] - x[i];
    }

    // Step 2: Calculate alpha_i
    for (int i = 1; i < n; i++) {
        alpha[i] = (3.0 / h[i]) * (f[i+1] - f[i]) - (3.0 / h[i-1]) * (f[i] - f[i-1]);
    }

    // Step 3: Set up the tridiagonal system
    if (Type == 1) { // Clamped boundary condition
        alpha[0] = (3.0 / h[0]) * (f[1] - f[0]) - 3.0 * s0;
        alpha[n] = 3.0 * sn - (3.0 / h[n-1]) * (f[n] - f[n-1]);
        l[0] = 2.0 * h[0];
        mu[0] = 0.5;
        z[0] = alpha[0] / l[0];
    } else if (Type == 2) { // Natural boundary condition
        alpha[0] = 0.0;
        alpha[n] = 0.0;
        l[0] = 1.0;
        mu[0] = 0.0;
        z[0] = 0.0;
    }

    for (int i = 1; i < n; i++) {
        l[i] = 2.0 * (x[i+1] - x[i-1]) - h[i-1] * mu[i-1];
        mu[i] = h[i] / l[i];
        z[i] = (alpha[i] - h[i-1] * z[i-1]) / l[i];
    }

    if (Type == 1) {
        l[n] = h[n-1] * (2.0 - mu[n-1]);
        z[n] = (alpha[n] - h[n-1] * z[n-1]) / l[n];
        c[n] = z[n];
    } else if (Type == 2) {
        l[n] = 1.0;
        z[n] = 0.0;
        c[n] = 0.0;
    }

    // Step 4: Back substitution
    for (int i = n-1; i >= 0; i--) {
        c[i] = z[i] - mu[i] * c[i+1];
        b[i] = (f[i+1] - f[i]) / h[i] - h[i] * (c[i+1] + 2.0 * c[i]) / 3.0;
        d[i] = (c[i+1] - c[i]) / (3.0 * h[i]);
    }

    // Step 5: adjust for 1 based indexing
    for (int i = n; i > 0; i--)
    {
        a[i] = f[i-1];
        b[i] = b[i-1];
        c[i] = c[i-1];
        d[i] = d[i-1];
    }
}

/**
 * @brief Computes the value of a cubic spline function at a given point.
 *
 * This function evaluates the cubic spline function S(t) at the given point t.
 * The cubic spline is defined by the arrays x, a, b, c, and d, which contain
 * the spline coefficients and knot points. The function is typically used to
 * interpolate a smooth curve through a set of data points.
 *
 * @param t The point at which to evaluate the spline function.
 * @param Fmax Default value when t is out of range.
 * @param n The number of intervals (or segments) in the spline.
 * @param x An array of knot points (size n+1).
 * @param a An array of spline coefficients a (1-based, size n-1).
 * @param b An array of spline coefficients b (1-based, size n-1).
 * @param c An array of spline coefficients c (1-based, size n-1).
 * @param d An array of spline coefficients d (1-based, size n-1).
 * @return The value of the spline function at the point t.
 */
double S(double t, double Fmax, int n, double x[], double a[], double b[], double c[], double d[]) {
    if (t < x[0] || t > x[n]) {
        return Fmax;
    }

    for (int i = 0; i < n; i++) {
        if (t >= x[i] && t <= x[i+1]) {
            double dx = t - x[i];
            return a[i+1] + b[i+1] * dx + c[i+1] * dx * dx + d[i+1] * dx * dx * dx;
        }
    }

    return Fmax;
}