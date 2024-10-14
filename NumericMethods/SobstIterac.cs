using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Text;
using System.Threading.Tasks;

namespace NumericMethods
{

    public class QrAlgorithm
    {
        public static void MatrixSolver()
        {

            Complex[,] matrix = {
            { new Complex(0.42016, 0), new Complex(-16.937, 0), new Complex(10.087, 0), new Complex(-2.8570, 0) },
            { new Complex(0.19439, 0), new Complex(-7.6571, 0), new Complex(4.5605, 0), new Complex(-1.3218, 0) },
            { new Complex(-6.1729, 0), new Complex(2.8952, 0), new Complex(-1.7253, 0), new Complex(0.41974, 0) },
            { new Complex(-0.20038, 0), new Complex(7.8932, 0), new Complex(-4.7011, 0), new Complex(1.3625, 0) }
        };

            int maxIterations = 100;
            double tolerance = 1e-10;
            FindEigenvalues(matrix, maxIterations, tolerance);
        }

        static void FindEigenvalues(Complex[,] matrix, int maxIterations, double tolerance)
        {
            int n = matrix.GetLength(0);
            Complex[,] currentMatrix = (Complex[,])matrix.Clone();

            for (int iter = 0; iter < maxIterations; iter++)
            {
                if (IsConverged(currentMatrix, tolerance))
                    break;

                Complex[,] q, r;
                QRDecomposition(currentMatrix, out q, out r);
                currentMatrix = MultiplyMatrices(r, q);
            }

            Console.WriteLine("Собственные значения:");
            for (int i = 0; i < n; i++)
            {
                Console.WriteLine(currentMatrix[i, i]);
            }
        }

        static bool IsConverged(Complex[,] matrix, double tolerance)
        {
            int n = matrix.GetLength(0);
            for (int i = 0; i < n - 1; i++)
            {
                if (Complex.Abs(matrix[i + 1, i]) > tolerance)
                    return false;
            }
            return true;
        }

        static void QRDecomposition(Complex[,] a, out Complex[,] q, out Complex[,] r)
        {
            int n = a.GetLength(0);
            q = new Complex[n, n];
            r = new Complex[n, n];

            Complex[,] u = new Complex[n, n];
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    u[j, i] = a[j, i];
                }

                for (int k = 0; k < i; k++)
                {
                    r[k, i] = 0;
                    for (int j = 0; j < n; j++)
                    {
                        r[k, i] += q[j, k] * a[j, i];
                    }

                    for (int j = 0; j < n; j++)
                    {
                        u[j, i] -= r[k, i] * q[j, k];
                    }
                }

                r[i, i] = 0;
                for (int j = 0; j < n; j++)
                {
                    r[i, i] += u[j, i] * Complex.Conjugate(u[j, i]);
                }
                r[i, i] = Complex.Sqrt(r[i, i]);

                for (int j = 0; j < n; j++)
                {
                    q[j, i] = u[j, i] / r[i, i];
                }
            }
        }

        static Complex[,] MultiplyMatrices(Complex[,] a, Complex[,] b)
        {
            int aRows = a.GetLength(0);
            int aCols = a.GetLength(1);
            int bRows = b.GetLength(0);
            int bCols = b.GetLength(1);

            if (aCols != bRows)
                throw new Exception("Матрицы не могут быть перемножены!");

            Complex[,] result = new Complex[aRows, bCols];
            for (int i = 0; i < aRows; i++)
            {
                for (int j = 0; j < bCols; j++)
                {
                    result[i, j] = 0;
                    for (int k = 0; k < aCols; k++)
                    {
                        result[i, j] += a[i, k] * b[k, j];
                    }
                }
            }
            return result;
        }

    }

}




