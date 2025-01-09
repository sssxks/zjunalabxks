double f(int n, double c[], double x)
{
    double ret = 0;
    for (int i = n; i >= 0; i--)
    {
        ret = ret * x + c[i];
    }
    return ret;
}
double f_prime(int n, double c[], double x)
{
    double ret = 0;
    for (int i = n; i >= 1; i--)
    {
        ret = ret * x + c[i] * i;
    }
    return ret;
}
double f_double_prime(int n, double c[], double x)
{
    double ret = 0;
    for (int i = n; i >= 2; i--)
    {
        ret = ret * x + c[i] * i * (i - 1);
    }
    return ret;
}
double Polynomial_Root(int n, double c[], double a, double b, double EPS)
{
    if (a > b)
    {
        double temp = a;
        a = b;
        b = temp;
    }
    for (double j = 0; j < 2; j++)
    {
        double x = a + (b - a) * j / 2;
        for (int i = 1; i < 1000; i++)
        {
            double f_val = f(n, c, x);
            double f_prime_val = f_prime(n, c, x);
            double f_double_prime_val = f_double_prime(n, c, x);
            double xn = x - (f_prime_val * f_val) / (f_prime_val * f_prime_val - f_val * f_double_prime_val);

            if (xn > b || xn < a)
                break;
            
            if (fabs(x - xn) < EPS)
            {
                return xn;
            }
            x = xn;
        }
    }
}
