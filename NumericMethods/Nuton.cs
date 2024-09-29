using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace NumericMethods
{
    public class Nuton
    {

        const int n = 2;
        const double eps = 1e-11;

        static double[,] a = new double[n, n];
        static double[] x = new double[n];
        static double[] f = new double[n];
        static double[] dx = new double[n];

        public static void NutonMethod()
        {
            Console.Clear();
            x[0] = 0; // x[0] = x
            x[1] = 0; // x[1] = y
            int iter = 0;
            double max;

            do
            {
                VivodVectr(x, 3, 13);
                Console.WriteLine("nomer iterazii - " + iter);
                Console.WriteLine("=================");

                // Increment iteration count
                iter++;

                // Calculate Jacobian matrix elements
                for (int i = 0; i < n; i++)
                    for (int j = 0; j < n; j++)
                        a[i, j] = Jacobian(x, i + 1, j + 1);

                // Calculate the right-hand side of the system
                for (int i = 0; i < n; i++)
                    f[i] = -Func(x, i + 1);

                // Solve the system
                Simq(n, a, f);
                dx = f;
                max = Math.Abs(dx[0]);
                for (int i = 1; i < n; i++)
                    if (Math.Abs(dx[i]) > max)
                        max = Math.Abs(dx[i]);

                for (int i = 0; i < n; i++)
                    x[i] += dx[i];

            } while (max >= eps);

            Console.ReadLine();
        }

        // Function to calculate the value of the function
        static double Func(double[] x, int i)
        {
            switch (i)
            {
                case 1:
                    return Math.Cos(x[0]) + x[1] - 1.5; // First function
                case 2:
                    return 2*x[0]-Math.Sin(x[1]-0.5)-1; // Second function
                default:
                    throw new ArgumentOutOfRangeException();
            }
        }

        // Function to calculate the Jacobian matrix
        static double Jacobian(double[] x, int i, int j)
        {
            switch (i)
            {
                case 1:
                    switch (j)
                    {
                        case 1:
                            return -Math.Sin(x[0]); // Jacobian element (1,1)
                        case 2:
                            return 1; // Jacobian element (1,2)
                    }
                    break;
                case 2:
                    switch (j)
                    {
                        case 1:
                            return 2; // Jacobian element (2,1)
                        case 2:
                            return -Math.Cos(x[1]-0.5); // Jacobian element (2,2)
                    }
                    break;
            }
            throw new ArgumentOutOfRangeException();
        }

        // Function to display the vector
        static void VivodVectr(double[] vector, int N1, int N2)
        {
            for (int j = 0; j < n; j++)
                Console.WriteLine($"x{j + 1}= {vector[j]:N2}");
        }

        // Function to solve the system of linear equations
        static void Simq(int Nn, double[,] A, double[] Bb)
        {
            const double eps = 1e-21;
            double max, u, v;
            int k1;
            double[,] Aa = new double[Nn, Nn + 1];

            for (int i = 0; i < Nn; i++)
                Aa[i, Nn] = Bb[i];

            for (int i = 0; i < Nn; i++)
                for (int j = 0; j < Nn; j++)
                    Aa[i, j] = A[i, j];

            for (int i = 0; i < Nn; i++)
            {
                max = Math.Abs(Aa[i, i]);
                k1 = i;

                for (int l = i + 1; l < Nn; l++)
                    if (Math.Abs(Aa[l, i]) > max)
                    {
                        max = Math.Abs(Aa[l, i]);
                        k1 = l;
                    }

                if (max < eps)
                    return;

                if (k1 != i)
                    for (int j = 0; j < Nn + 1; j++)
                    {
                        u = Aa[i, j];
                        Aa[i, j] = Aa[k1, j];
                        Aa[k1, j] = u;
                    }

                v = Aa[i, i];
                for (int j = 0; j < Nn + 1; j++)
                    Aa[i, j] /= v;

                for (int l = i + 1; l < Nn; l++)
                {
                    v = Aa[l, i];
                    for (int j = i; j < Nn + 1; j++)
                        Aa[l, j] -= Aa[i, j] * v;
                }
            }

            Bb[Nn - 1] = Aa[Nn - 1, Nn];

            for (int i = Nn - 2; i >= 0; i--)
            {
                Bb[i] = Aa[i, Nn];
                for (int j = i + 1; j < Nn; j++)
                    Bb[i] -= Aa[i, j] * Bb[j];
            }

        }

    }
}
