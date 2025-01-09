double int_term(int x);
double series_sum(double x, double phi[], int n, int m, int p);
void Series_Sum(double sum[]);

// theoretically, 128 is a closer upper bound for 10^(-10), but it doesn't AC in practice. 16x more is still acceptable
#define ITERS 2048
#define phi(x) phi[(x) * 10]
#define PI_SQUARED_OVER_6 1.64493406684823

void Series_Sum(double sum[])
{
    // "calculate" phi(0), compile time 😋
    sum[0] = PI_SQUARED_OVER_6; // 1/k^2 for i = 1 to oo

    // calculate phi(x) when x in Z 😙
    for (int i = 1; i <= 300; i++)
    {
        sum[i * 10] = int_term(i);
    }

    // evaluate at other point, using phi(x) x in Z 🤯
    for (int i = 0; i <= 2990; i++)
    {
        if (i % 10 == 0)
        {
            continue;
        }
        double x = i * 0.1;
        int n = i / 10;
        int m = i / 10 + 1;
        int p = i / 10 + 2;
        sum[i] = series_sum(x, sum, n, m, p);
    }
    // we can not use n, m, p as 👆
    // as p will become 301, array out of bound
    for (int i = 2991; i < 3000; i++)
    {
        double x = i * 0.1;
        int n = i / 10 - 1; // 298
        int m = i / 10;     // 299
        int p = i / 10 + 1; // 300
        sum[i] = series_sum(x, sum, n, m, p);
    }
}

double int_term(int x)
{
    double sum = 0;
    for (int i = x; i > 0; i--) 
    {
        sum += 1.0 / i;
    }

    return sum / x;
}

// calculate phi(x). assuming integer value already present
// n, m, p: reference value, order doesn't matter. should be close to x for better precision
double series_sum(double x, double phi[], int n, int m, int p)
{
    double main_part = 0.0;
    for (int k = ITERS; k > 0; k--) // hopefully compiler unroll this 🚀
    {
        main_part += 1 / (k * (k + x) * (k + n) * (k + m) * (k + p));
    }

    double ref_part = ((n - p) * (n - x) * (p - x) * phi(m) - (m - p) * (m - x) * (p - x) * phi(n) + (m - n) * (m - x) * (n - x) * phi(p)) / ((m - n) * (m - p) * (n - p));

    return (m - x) * (n - x) * (p - x) * main_part + ref_part;
}