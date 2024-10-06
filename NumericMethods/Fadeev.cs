using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace NumericMethods
{
    public class Fadeev
    {
        public static void CalculateEigenvaluesAndVectors(float[,] matrix, out float[] eigenvalues, out float[,] eigenvectors)
        {
            // Проверка размерности матрицы
            if (matrix.GetLength(0) != 4 || matrix.GetLength(1) != 4)
                throw new ArgumentException("Matrix must be 4x4");

            // Вычисление характеристического полинома
            float[] characteristic = CalculateCharacteristicPolynomial(matrix);

            // Вычисление собственных значений
            eigenvalues = FindEigenvalues(characteristic);

            // Вычисление собственных векторов
            eigenvectors = FindEigenvectors(matrix, eigenvalues);
        }

        private static float[] CalculateCharacteristicPolynomial(float[,] matrix)
        {
            float[] characteristic = new float[5];

            // Вычисление коэффициентов характеристического полинома
            characteristic[0] = 1;
            characteristic[1] = -matrix.GetLength(0);

            for (int i = 0; i < matrix.GetLength(0); i++)
            {
                float sum = 0;
                for (int j = 0; j < matrix.GetLength(1); j++)
                {
                    sum += matrix[i, j];
                }
                characteristic[1] += -sum;
            }

            for (int i = 0; i < matrix.GetLength(0); i++)
            {
                for (int j = 0; j < matrix.GetLength(1); j++)
                {
                    float product = 1;
                    for (int k = 0; k < matrix.GetLength(0); k++)
                    {
                        if (k != i && k != j)
                        {
                            product *= matrix[k, k];
                        }
                    }
                    characteristic[2] += product;
                }
            }

            for (int i = 0; i < matrix.GetLength(0); i++)
            {
                float product = 1;
                for (int j = 0; j < matrix.GetLength(1); j++)
                {
                    product *= matrix[j, j];
                }
                characteristic[3] += -product;
            }

            characteristic[4] = -matrix[0, 0] * matrix[1, 1] * matrix[2, 2] * matrix[3, 3];

            return characteristic;
        }

        private static float[] FindEigenvalues(float[] characteristic)
        {
            // Вычисление собственных значений с помощью усовершенствованного метода Фадеева
            float[] eigenvalues = new float[4];

            float p = characteristic[1];
            float q = characteristic[2];
            float r = characteristic[3];
            float s = characteristic[4];

            float discriminant = p * p - 3 * q;
            float a = (2 * p * p * p - 9 * p * q + 27 * s) / 54;
            float b = (discriminant) / 9;

            if (b > 0)
            {
                float theta = (float)Math.Acos(a / Math.Sqrt(b * b * b));
                eigenvalues[0] = -2 * (float)Math.Sqrt(b) * (float)Math.Cos(theta / 3) - p / 3;
                eigenvalues[1] = -2 * (float)Math.Sqrt(b) * (float)Math.Cos((theta + 2 * Math.PI) / 3) - p / 3;
                eigenvalues[2] = -2 * (float)Math.Sqrt(b) * (float)Math.Cos((theta - 2 * Math.PI) / 3) - p / 3;
                eigenvalues[3] = -p / 3;
            }
            else
            {
                float A = -a;
                float B = (float)Math.Sqrt(-b);
                eigenvalues[0] = 2 * B - p / 3;
                eigenvalues[1] = -B - p / 3;
                eigenvalues[2] = -B - p / 3;
                eigenvalues[3] = -p / 3;
            }

            return eigenvalues;
        }

        private static float[,] FindEigenvectors(float[,] matrix, float[] eigenvalues)
        {
            float[,] eigenvectors = new float[4, 4];

            for (int i = 0; i < 4; i++)
            {
                float[,] A = new float[4, 4];
                for (int j = 0; j < 4; j++)
                    for (int k = 0; k < 4; k++)
                        A[j, k] = matrix[j, k] - eigenvalues[i] * (j == k ? 1 : 0);

                float[] nullspace = FindNullspace(A);

                for (int j = 0; j < 4; j++)
                    eigenvectors[j, i] = nullspace[j];
            }

            return eigenvectors;
        }

        private static float[] FindNullspace(float[,] matrix)
        {
            float[] nullspace = new float[4];

            // Вычисление ядра матрицы с помощью метода Гаусса
            float[,] augmented = new float[4, 5];
            for (int i = 0; i < 4; i++)
                for (int j = 0; j < 4; j++)
                    augmented[i, j] = matrix[i, j];
            for (int i = 0; i < 4; i++)
                augmented[i, 4] = 0;

            int pivotCol = 0;
            for (int row = 0; row < 4; row++)
            {
                // Найти ведущий элемент
                float maxVal = 0;
                int maxRow = row;
                for (int i = row; i < 4; i++)
                {
                    if (Math.Abs(augmented[i, pivotCol]) > Math.Abs(maxVal))
                    {
                        maxVal = augmented[i, pivotCol];
                        maxRow = i;
                    }
                }

                // Поменять строки, если необходимо
                if (maxRow != row)
                {
                    for (int col = 0; col < 5; col++)
                    {
                        float temp = augmented[row, col];
                        augmented[row, col] = augmented[maxRow, col];
                        augmented[maxRow, col] = temp;
                    }
                }

                // Выполнить операцию вычитания
                for (int i = row + 1; i < 4; i++)
                {
                    float factor = augmented[i, pivotCol] / augmented[row, pivotCol];
                    for (int col = pivotCol; col < 5; col++)
                    {
                        augmented[i, col] -= factor * augmented[row, col];
                    }
                }

                pivotCol++;
            }

            // Найти вектор из ядра матрицы
            for (int i = 3; i >= 0; i--)
            {
                float sum = 0;
                for (int j = i + 1; j < 4; j++)
                {
                    sum += augmented[i, j] * nullspace[j];
                }
                nullspace[i] = -sum / augmented[i, i];
            }

            return nullspace;
        }
    }
}
