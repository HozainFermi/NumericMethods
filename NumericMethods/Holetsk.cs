using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace NumericMethods
{
    public static class Holetsk
    {

       public static void Holet(double[,] A1, ref double[] B1, int N)
        {
            int i, j, t;
            double[,] c = new double[N, N];
            double[,] L = new double[N, N + 1];
            double[] y = new double[N];
            double buf, summ;

            // Инициализация матриц
            for (i = 0; i < N; i++)
            {
                for (j = 0; j < N; j++)
                {
                    c[i, j] = 0;
                    L[i, j] = 0;
                    if (j < N) y[i] = 0;
                }
            }

            // Умножение матрицы на транспонированную
            for (i = 0; i < N; i++)
            {
                for (j = 0; j < N; j++)
                {
                    summ = 0.0;
                    for (t = 0; t < N; t++)
                    {
                        summ += A1[t, j] * A1[t, i];
                    }
                    c[i, j] = summ;
                }
            }

            //Умножение правой стороны на транспонированную матрицу 
            for (i = 0; i < N; i++)
            {
                for (j = 0; j < N; j++)
                {
                    y[i] += A1[j, i] * B1[j];
                }
            }

            for (i = 0; i < N; i++)
            {
                for (j = 0; j < N; j++)
                {
                    A1[i, j] = c[i, j];
                }
                B1[i] = y[i];
            }

            for (i = 0; i < N; i++)
            {
                for (j = 0; j <= i; j++)
                {
                    summ = 0;
                    for (t = 0; t < j; t++)
                    {
                        summ += L[i, t] * L[j, t];
                    }
                    if (i != j)
                    {
                        L[i, j] = (A1[i, j] - summ) / L[j, j];
                    }
                    else
                    {
                        L[i, i] = Math.Sqrt(A1[i, i] - summ);
                    }
                }
            }

            for (i = 0; i < N; i++)
            {
                L[i, N] = B1[i];
            }

            B1[0] = L[0, N] / L[0, 0];

            for (i = 1; i < N; i++)
            {
                for (j = 0; j < i; j++)
                {
                    L[i, N] -= L[i, j] * B1[j];
                }
                B1[i] = L[i, N] / L[i, i];
            }

            for (i = 0; i < N; i++)
            {
                for (j = i + 1; j < N; j++)
                {
                    L[i, j] = L[j, i];
                    L[j, i] = 0;
                }
                L[i, N] = B1[i];
            }

            B1[N - 1] = L[N - 1, N] / L[N - 1, N - 1];

            for (i = N - 2; i >= 0; i--)
            {
                for (j = i + 1; j < N; j++)
                {
                    L[i, N] -= L[i, j] * B1[j];
                }
                B1[i] = L[i, N] / L[i, i];
            }
        }

       public static void PrintMatrix(double[,] matrix, int N)
        {
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    Console.Write(matrix[i, j] + "\t");
                }
                Console.WriteLine();
            }
        }

       public static void PrintVector(double[] vector, int N)
        {
            for (int i = 0; i < N; i++)
            {
                Console.WriteLine(vector[i]);
            }
        }

    }
}
