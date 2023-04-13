#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <vector>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <windows.h>
#include "log_duration.h"

std::vector<double> CreateRandomVector(size_t n)
{
    std::vector<double> result(n);
    for (int i = 0; i < n; i++)
    {
        result[i] = (double)(rand() % 10000 + 100) / 100;
    }
    return result;
}

///// <summary>
///// Прямой ход метода Прогонки
///// </summary>
///// <param name="matrix"></param>
///// <param name="b"></param>
//static void DirectCourse(TridiagonalMatrix matrix, double[] b)
//{
//    double a;
//    for (int i = 2; i < n; i++)
//    {
//        a = matrix[i, i - 1] / matrix[i - 1, i - 1];
//        matrix[i, i - 1] -= matrix[i - 1, i - 1] * a;
//        matrix[i, i] -= matrix[i - 1, i] * a;
//        b[i] -= b[i - 1] * a;
//    }
//}
///// <summary>
///// Обратный ход метода Прогонки
///// </summary>
///// <param name="matrix"></param>
///// <param name="b"></param>
///// <returns></returns>
//static double[] ReverseCourse(TridiagonalMatrix matrix, double[] b)
//{
//    double[] y = new double[b.Length];
//    y[n] = 0;
//    for (int i = n - 1; i > 0; i--)
//    {
//        y[i] = (b[i] - y[i + 1] * matrix[i, i + 1]) / matrix[i, i];
//    }
//    return y;
//}
///// <summary>
///// Метод Прогонки
///// </summary>
///// <param name="A"></param>
//static void Thomas(TridiagonalMatrix A)
//{
//    TridiagonalMatrix matrix = (TridiagonalMatrix)A.Clone();
//    double[] B = InitiateB();
//    DirectCourse(matrix, B);
//    double[] y = ReverseCourse(matrix, B);
//    string header = "Метод Прогонки";
//    reporter.Add(header, headers, DataToTable(y));
//}

std::vector<double> CyclicReduction(std::vector<double> a, std::vector<double> b, std::vector<double> c, std::vector<double> f, const int q, const int n)
{
    for (int k = 1; k < q; k++)
    {
        int twoInK = pow(2, k);
        int twoInKminusOne = pow(2, k - 1);
        omp_set_num_threads(6);
        double P;
        double Q;
        #pragma omp parallel for 
        for (int i = twoInK; i < n; i += twoInK)
        {
            P = a[i] / b[i - twoInKminusOne];
            Q = c[i] / b[i + twoInKminusOne];
            a[i] = P * a[i - twoInKminusOne];
            b[i] = b[i] - P * c[i - twoInKminusOne] - Q * a[i + twoInKminusOne];
            c[i] = Q * c[i + twoInKminusOne];
            f[i] = f[i] + P * f[i - twoInKminusOne] + Q * f[i + twoInKminusOne];
        }
    }
    std::vector<double> x(n + 1);
    x[0] = f[0];
    x[n] = f[n];
    for (int k = q; k > 0; k--)
    {
        int twoInK = pow(2, k);
        double twoInKminusOne = pow(2, k - 1);
        for (int i = twoInKminusOne; i <= n - twoInKminusOne; i += twoInK)
        {
            x[i] = (f[i] + a[i] * x[i - twoInKminusOne] + c[i] * x[i + twoInKminusOne]) / b[i];
        }
    }
    return x;
}

void TestReduction() 
{
    int q = 14; //n=2^q - количество уравнений
    int n = pow(2, q);
    std::vector<double> a = CreateRandomVector(n + 1);
    std::vector<double> b = CreateRandomVector(n + 1);
    std::vector<double> c = CreateRandomVector(n + 1);
    a[0] = 0;
    c[0] = 0;
    b[0] = 1;
    a[n] = 0;
    c[n] = 0;
    b[n] = 1;
    std::vector<double> x = CreateRandomVector(n + 1);
    std::vector<double> f(n + 1);
    f[0] = x[0];
    f[n] = x[n];
    for (int i = 1; i < n; i++)
    {
        f[i] = -a[i] * x[i - 1] + b[i] * x[i] - c[i] * x[i + 1];
    }
    std::vector<double> res = CyclicReduction(a, b, c, f, q, n);
    for (int i = 0; i < res.size(); i++)
    {
        //std::cout << x[i] << " " << res[i] << std::endl;
    }
}

const double a = -1;
const double b = 1;
const double epsilon = 0.05;

double q(double x)
{
    return 1.0 / epsilon;
}

double f(double x)
{
    return (1.0 / epsilon + M_PI * M_PI) * cos(M_PI * x);
}

double u(double x)
{
    return cos(M_PI_2 * x) + exp((x - 1) / sqrt(epsilon)) + exp(-(x + 1) / sqrt(epsilon));
}

//n=2^q - количество уравнений
void Task1(double degree)
{
    int n = pow(2, degree);
    double h = (b - a) / n;
    std::vector<double> c1(n + 1);
    std::vector<double> c2(n + 1);
    std::vector<double> c3(n + 1);
    std::vector<double> x(n + 1);
    std::vector<double> d(n + 1);
    for (int i = 1; i < n; i++)
    {
        x[i] = a + i * h;
        c1[i] = -1;
        c2[i] = 2 + h * h * q(x[i]);
        c3[i] = -1;
        d[i] = h * h * f(x[i]);
    }
    c1[0] = 0;
    c2[0] = 1;
    c3[0] = 0;
    c1[n] = 0;
    c2[n] = 1;
    c3[n] = 0;
    d[0] = u(a);
    d[n] = u(b);
    std::vector<double> res = CyclicReduction(c1, c2, c3, d, degree, n);
    
}

int main()
{
    {
        LOG_DURATION("Time");
        for (int i = 0; i < 10; i++)
            TestReduction();
    }
    //TestReduction();
}
