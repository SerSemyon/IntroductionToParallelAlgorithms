#define _USE_MATH_DEFINES
#include <omp.h>
#include <cmath>
#include <iostream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <windows.h>
#include <string>
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

int numThreads;

std::vector<double> CyclicReduction(const std::vector<double> A,const std::vector<double> B, const std::vector<double> C, const std::vector<double> F, const int q, const int n)
{
    double* P = new double[n];
    double* Q = new double[n]; 
    std::vector<double> a(A);
    std::vector<double> b(B);
    std::vector<double> c(C);
    std::vector<double> f(F);
    for (int k = 1; k < q; k++)
    {
        int twoInK = pow(2, k);
        int twoInKMinusOne = pow(2, k - 1);
#pragma omp parallel for num_threads(numThreads)
        for (int i = twoInK; i < n; i += twoInK)
        {
            P[i] = a[i] / b[i - twoInKMinusOne];
            Q[i] = c[i] / b[i + twoInKMinusOne];
            a[i] = P[i] * a[i - twoInKMinusOne];
            b[i] = b[i] - P[i] * c[i - twoInKMinusOne] - Q[i] * a[i + twoInKMinusOne];
            c[i] = Q[i] * c[i + twoInKMinusOne];
            f[i] = f[i] + P[i] * f[i - twoInKMinusOne] + Q[i] * f[i + twoInKMinusOne];
        }
    }
    std::vector<double> x(n + 1);
    x[0] = f[0];
    x[n] = f[n];
    for (int k = q; k > 0; k--)
    {
        int twoInK = pow(2, k);
        int twoInKMinusOne = pow(2, k - 1);
#pragma omp parallel for num_threads(numThreads)
        for (int i = twoInKMinusOne; i <= n - twoInKMinusOne; i += twoInK)
        {
            x[i] = (f[i] + a[i] * x[i - twoInKMinusOne] + c[i] * x[i + twoInKMinusOne]) / b[i];
        }
    }
    delete[] P;
    delete[] Q;
    return x;
}

std::vector<double> ThomasAlgorithm(std::vector<double> A, std::vector<double> B, std::vector<double> C, std::vector<double> F)
{
    std::vector<double> a(A);
    std::vector<double> b(B);
    std::vector<double> c(C);
    std::vector<double> f(F);
    int n = f.size() - 1;
    std::vector<double> P(n), Q(n), x(n+1);
    P[0] = c[0] / b[0]; Q[0] = f[0] / b[0];
    for (int i = 1; i < n; i++)
    {
        P[i] = c[i] / (b[i] - a[i] * P[i - 1]);
        Q[i] = (f[i] + a[i] * Q[i - 1]) / (b[i] - a[i] * P[i - 1]);
    }
    x[n] = (f[n] + a[n] * Q[n - 1]) / (b[n] - a[n] * P[n - 1]);
    for (int i = n - 1; i >= 0; i--)
    {
        x[i] = P[i] * x[i + 1] + Q[i];
    }
    return x;
}

void TestThomas()
{
    int q = 4; //n=2^q - количество уравнений
    int n = pow(2, q);
    std::vector<double> a = CreateRandomVector(n + 1);
    std::vector<double> b = CreateRandomVector(n + 1);
    std::vector<double> c = CreateRandomVector(n + 1);
    a[0] = 0;
    c[n] = 0;
    std::vector<double> x = CreateRandomVector(n + 1);
    std::vector<double> f(n + 1);
    f[0] = b[0] * x[0] - c[0] * x[1];
    for (int i = 1; i < n; i++)
    {
        f[i] = -a[i] * x[i - 1] + b[i] * x[i] - c[i] * x[i + 1];
    }
    f[n] = -a[n] * x[n - 1] + b[n] * x[n];
    //f[n] = 0;
    std::vector<double> res = ThomasAlgorithm(a, b, c, f);
    for (int i = 0; i < res.size(); i++)
    {
        //std::cout << x[i] << " " << res[i] << std::endl;
    }
}

void TestReduction() 
{
    int q = 15; //n=2^q - количество уравнений
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
    for (int i = 1; i < n; i++)
    {
        f[i] = -a[i] * x[i - 1] + b[i] * x[i] - c[i] * x[i + 1];
    }
    f[0] = x[0];
    f[n] = x[n];
    {
        LOG_DURATION("reduction");
        std::vector<double> res = CyclicReduction(a, b, c, f, q, n);
    }
    //std::vector<double> res = CyclicReduction(a, b, c, f, q, n);
    //for (int i = 0; i < res.size(); i++)
    //{
    //    //std::cout << x[i] << " " << res[i] << std::endl;
    //}
}

const double a = 0;
const double b = 1;
const double epsilon = 0.05;

double q(double x)
{
    return 1.0 / epsilon;
}

double f(double x)
{
    return 0;
}

double u(double x)
{
    return (exp(-x/sqrt(epsilon))-exp((x-2)/sqrt(epsilon)))/(1-exp(-2/sqrt(epsilon)));
}

double Inaccuracy(std::vector<double> x, std::vector<double> y)
{
    double inaccuracy = 0;
    for (int i = 0; i < y.size(); i++)
    {
        double diff = abs(y[i] - u(x[i]));
        if (diff > inaccuracy)
            inaccuracy = diff;
    }
    return inaccuracy;
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
        c1[i] = 1;
        c2[i] = 2 + h * h * q(x[i]);
        c3[i] = 1;
        d[i] = h * h * f(x[i]);
    }
    x[0] = a;
    x[n] = b;
    c1[0] = 0;
    c2[0] = 1;
    c3[0] = 0;
    c1[n] = 0;
    c2[n] = 1;
    c3[n] = 0;
    d[0] = u(x[0]);
    d[n] = u(x[n]);
    std::vector<double> res(n + 1);
    {
        LOG_DURATION("Thomas");
        res = ThomasAlgorithm(c1, c2, c3, d);
    }
    std::cout << "inaccuracy " << Inaccuracy(x, res) << std::endl;
    numThreads = 2;
    {
        //LOG_DURATION("Reduction "+std::to_string( numThreads));
        res = CyclicReduction(c1, c2, c3, d, degree, n);
    }
    {
        LOG_DURATION("Reduction " + std::to_string(numThreads));
        res = CyclicReduction(c1, c2, c3, d, degree, n);
    }
    std::cout << "inaccuracy " << Inaccuracy(x, res) << std::endl;
    numThreads = 4;
    {
        //LOG_DURATION("Reduction " + std::to_string(numThreads));
        res = CyclicReduction(c1, c2, c3, d, degree, n);
    }
    {
        LOG_DURATION("Reduction " + std::to_string(numThreads));
        res = CyclicReduction(c1, c2, c3, d, degree, n);
    }
    std::cout << "inaccuracy " << Inaccuracy(x, res) << std::endl;
    numThreads = 6;
    {
        //LOG_DURATION("Reduction " + std::to_string(numThreads));
        res = CyclicReduction(c1, c2, c3, d, degree, n);
    }
    {
        LOG_DURATION("Reduction " + std::to_string(numThreads));
        res = CyclicReduction(c1, c2, c3, d, degree, n);
    }
    std::cout << "inaccuracy " << Inaccuracy(x, res) << std::endl;
}

int main()
{
    /*{
        LOG_DURATION("Progonka");
        TestThomas();
    }
    {
        LOG_DURATION("TestReduction");
        for (int i = 0; i < 1; i++)
            TestReduction();
    }*/
    for (int i = 9; i < 20; i++)
    {
        std::cout << "number divisions " << pow(2,i) << std::endl;
        Task1(i);
        std::cout << std::endl;
    }
    /*for (int i = 0; i < 4; i++)
    {
        Task1(10);
    }*/
}
