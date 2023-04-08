﻿// IntroductionToParallelAlgorithms.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
#include <vector>
#include "TridiagonalMatrix.h"

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
        double twoInKminusOne = pow(2, k - 1);
        for (int i = twoInK; i < n; i += twoInK)
        {
            double P = a[i] / b[i - twoInKminusOne];
            double Q = c[i] / b[i + twoInKminusOne];
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

int main()
{
    int q = 4; //n=2^q - количество уравнений
    int n = pow(2, q);
    std::vector<double> a = CreateRandomVector(n+1);
    std::vector<double> b = CreateRandomVector(n+1);
    std::vector<double> c = CreateRandomVector(n+1);
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
    //std::cout << x[0] << '=' << f[0] << std::endl;
    for (int i = 1; i < n; i++)
    {
        f[i] = -a[i] * x[i - 1] + b[i] * x[i] - c[i]*x[i + 1];
        //std::cout << -a[i] << '*' << x[i - 1] << '+' << b[i] << '*' << x[i] << '-' << c[i] << '*' << x[i + 1] << '=' << f[i] << std::endl;
    }
    //std::cout << x[n] << '=' << f[n] << std::endl;

    std::vector<double> res = CyclicReduction(a, b, c, f, q, n);
    for (int i = 0; i < res.size(); i++)
    {
        std::cout << x[i] << " " << res[i] << std::endl;
    }
}
