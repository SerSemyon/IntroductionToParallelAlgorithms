#pragma once
/// <summary>
/// Класс квадратных трёхдиагональных матриц
/// </summary>
class TridiagonalMatrix
{
    double* diag;
    double* high_diag;
    double* low_diag;
    int size;
public:
    int GetSize()
    {
        return size;
    }
    /// <summary>
        /// Задание элемента верхней диагонали
        /// </summary>
        /// <param name="row"> Индекс строки </param>
        /// <param name="value"> Записываемое значение </param>
    void SetHighDiag(int row, double value)
    {
        high_diag[row] = value;
    }
    /// <summary>
    /// Задание элемента центральной диагонали
    /// </summary>
    /// <param name="row"> Индекс строки </param>
    /// <param name="value"> Записываемое значение </param>
    void SetDiag(int row, double value)
    {
        diag[row] = value;
    }
    /// <summary>
    /// Задание элемента нижней диагонали
    /// </summary>
    /// <param name="row"> Индекс строки </param>
    /// <param name="value"> Записываемое значение </param>
    void SetLowDiag(int row, double value)
    {
        low_diag[row - 1] = value;
    }
    TridiagonalMatrix(int n)
    {
        size = n;
        diag = new double[n];
        high_diag = new double[n - 1];
        low_diag = new double[n - 1];
    }
    TridiagonalMatrix(const TridiagonalMatrix& original) : TridiagonalMatrix(original.size)
    {
        for (int i = 0; i < size - 1; i++)
        {
            SetHighDiag(i, high_diag[i]);
        }
        for (int i = 0; i < size; i++)
        {
            SetDiag(i, diag[i]);
        }
        for (int i = 0; i < size - 1; i++)
        {
            SetLowDiag(i + 1, low_diag[i]);
        }
    }

    double& Element(int rows, int colums)
    {

        if (rows == colums)
        {
            return diag[colums];
        }
        else if (colums - rows == 1)
        {
            return high_diag[rows];
        }
        else if (rows - colums == 1)
        {
            return low_diag[colums];
        }
        double loc = 0;
        return loc;
    }

    void PrintMatrix()
    {
        for (int i = 0; i < size; i++)
        {
            for (int j = 0; j < size; j++)
            {
                std::cout << Element(i, j) << ' ';
            }
            std::cout << '\n';
        }

    }
    void SetMatrix(int row, int column, double value) {
        if (row == column)
        {
            diag[row] = value;
        }
        else if (column - row == 1)
        {
            high_diag[row] = value;
        }
        else if (row - column == 1)
        {
            low_diag[column] = value;
        }
    }
    double GetMatrix(int row, int column) {
        if (row == column)
        {
            return diag[row];
        }
        else if (column - row == 1)
        {
            return high_diag[row];
        }
        else if (row - column == 1)
        {
            return low_diag[column];
        }
        else return 0;

    }
};