using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace NumericMethods
{
    public class Krilov
    {
        private float[,] matrix;

        public Krilov (float[,] matrix)
        {
            this.matrix = matrix;
        }

        public (float[], float[,]) CalculateEigenvaluesAndVectors()
        {
            int n = matrix.GetLength(0);
            float[] eigenvalues = new float[n];
            float[,] eigenvectors = new float[n, n];

            // Вычисление собственных значений
            float[] characteristic = CalculateCharacteristicPolynomial(matrix);
            eigenvalues = FindRoots(characteristic);

            // Вычисление собственных векторов
            for (int i = 0; i < n; i++)
            {
                float[] vector = CalculateEigenvector(matrix, eigenvalues[i]);
                for (int j = 0; j < n; j++)
                {
                    eigenvectors[j, i] = vector[j];
                }
            }

            return (eigenvalues, eigenvectors);
        }

        private float[] CalculateCharacteristicPolynomial(float[,] matrix)
        {
            int n = matrix.GetLength(0);
            float[] characteristic = new float[n + 1];

            characteristic[0] = 1;
            characteristic[1] = -matrix.GetLength(0);

            for (int i = 0; i < n; i++)
            {
                float sum = 0;
                for (int j = 0; j < n; j++)
                {
                    sum += matrix[i, j];
                }
                characteristic[1] += -sum;
            }

            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    float product = 1;
                    for (int k = 0; k < n; k++)
                    {
                        if (k != i && k != j)
                        {
                            product *= matrix[k, k];
                        }
                    }
                    characteristic[2] += product;
                }
            }

            for (int i = 0; i < n; i++)
            {
                float product = 1;
                for (int j = 0; j < n; j++)
                {
                    product *= matrix[j, j];
                }
                characteristic[3] += -product;
            }

            characteristic[4] = -matrix[0, 0] * matrix[1, 1] * matrix[2, 2] * matrix[3, 3];

            return characteristic;
        }

        private float[] FindRoots(float[] characteristic)
        {
            float[] roots = new float[characteristic.Length - 1];

            // Реализация метода Крылова для нахождения корней
            float p = characteristic[1] / characteristic[0];
            float q = characteristic[2] / characteristic[0];
            float r = characteristic[3] / characteristic[0];
            float s = characteristic[4] / characteristic[0];

            float discriminant = p * p - 4 * q;
            if (discriminant >= 0)
            {
                roots[0] = (-p + MathF.Sqrt(discriminant)) / 2;
                roots[1] = (-p - MathF.Sqrt(discriminant)) / 2;
            }
            else
            {
                float realPart = -p / 2;
                float imaginaryPart = MathF.Sqrt(-discriminant) / 2;
                roots[0] = realPart;
                roots[1] = realPart;
                roots[2] = realPart;
                roots[3] = realPart;
            }

            return roots;
        }

        private float[] CalculateEigenvector(float[,] matrix, float eigenvalue)
        {
            int n = matrix.GetLength(0);
            float[,] A = new float[n, n];

            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    A[i, j] = matrix[i, j] - eigenvalue * (i == j ? 1 : 0);
                }
            }

            float[] vector = new float[n];
            vector[0] = 1;

            for (int i = 1; i < n; i++)
            {
                float sum = 0;
                for (int j = 0; j < i; j++)
                {
                    sum += A[i, j] * vector[j];
                }
                vector[i] = -sum / A[i, i];
            }

            return vector;
        }
    }
}
