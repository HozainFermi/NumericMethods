﻿

using System.Numerics;

namespace NumericMethods
{
    internal class Program
    {
        const int N = 3;

        
        class Matrix
        {
            public double[,] Data;

            public Matrix(int rows, int cols)
            {
                Data = new double[rows, cols];
            }

            public double this[int i, int j]
            {
                get { return Data[i - 1, j - 1]; }
                set { Data[i - 1, j - 1] = value; }
            }
        }

       
        class Vector
        {
            public double[] Data;

            public Vector(int size)
            {
                Data = new double[size];
            }

            public double this[int i]
            {
                get { return Data[i - 1]; }
                set { Data[i - 1] = value; }
            }
        }

        static Matrix A = new Matrix(N, N + 1)
        {
            Data = new double[,]
            {
                {3.88,3.78,3.45,0},
                {3.00,2.79,2.39,0},
                {2.67,2.37,1.96,0}
            }
        };

        static Vector B = new Vector(N)
        {
            Data = new double[] { 10.41, 8.36, 7.62 }
        };

        static void Main(string[] args)
        {
            //   double[] Bb = { 10.41, 8.36, 7.62 };
            //   double[,] Aa = new double[3, 4]
            //   {
            //       {3.88,3.78,3.45,0},
            //       {3.00,2.79,2.39,0},
            //       {2.67,2.37,1.96,0}
            //
            //   };
            //   int Ks = 0;
            //
            //   MatrixSolver.SIMQ(3, ref Aa, ref Bb, ref Ks);
            //
            //   Console.WriteLine("SIMQ");
            //   for (int i = 0; i < Bb.Length; i++)
            //   {
            //       Console.WriteLine($"x[{i}] = {Bb[i]:F3}");
            //       Console.WriteLine('\n');
            //   }
            //
            //   VivodMatr(A, N);
            //   Console.WriteLine("metod Holetckogo");
            //   Holet(A, B, N);
            //   VivodVectr(B, N);
            //   //
            //
            //   double[] By = { 9.7191, 10.5000, 10.9195, 10.9775 };
            //   double[,] Aaa = new double[4, 4]
            //   {
            //       {9.5000,0.0422, 0.0513, 0.0604},
            //       {0.0278,8.6000,0.0459,0.0550},
            //       {0.0224,0.0315,7.7000,0.0496},
            //       {0.0170,0.0261,0.0351,6.8000 }
            //
            //   };
            //   HoletskSELFMADE.HoletskSELF(Aa, Bb, 3);
            //
            //   int step = 0;
            //
            //
            //   Lab2.Yakoby(ref Aaa, ref By, 4, ref step);
            //   step = 0;
            //   double[] Bz = { 9.7191, 10.5000, 10.9195, 10.9775 };
            //   Lab2.Zeidel( Aaa, Bz, 4, ref step,0.01);
            //
            //   Console.WriteLine("Zeidel Ex");
            //   double w = 0.2;
            //   double[] Bzx = { 9.7191, 10.5000, 10.9195, 10.9775 };
            //   do
            //   {
            //       Lab2.ZeidelEx(Aaa, Bzx, 4, ref step, w);
            //       w = Math.Round(w, 1) + 0.2;
            //       w = Math.Round(w, 1);
            //       
            //   }
            //   while (w < 2.0);
            //   double[] bvozm = new double[2];
            //   
            //
            //   double[] b = { 0.502, 0.482 };
            //      double[,] a = new double[2, 2]
            //      {
            //          {1.03,0.998},
            //          {0.991,0.946},
            //      };
            //
            //   double[] X = new double[2];
            // 
            //   int p = 0;
            //
            //   Regul.RegulMethod(a, b, out X, out p);
            //   foreach (double x in X) { Console.WriteLine($"{x:F3}"); }
            //   Console.WriteLine('\n');
            //
            //   for (int i = 0; i < b.Length; i++)
            //   {
            //       b[i] += 0.005;
            //   }
            //
            //   Regul.RegulMethod(a, b, out X, out p);
            //   foreach (double x in X) { Console.WriteLine($"{x:F3}"); }
            //
            //   Console.WriteLine('\n');
            //   
            //
            //   double[] XN = new double[2];
            //   Gvines.Vrash(2, a, b, ref XN);
            //   foreach (double x in XN) Console.WriteLine ($"{x} "); 
            //   Console.WriteLine();
            //   
            //   for (int i = 0; i < bvozm.Length; i++)
            //   {
            //       bvozm[i] = b[i] + 0.005;
            //   }
            //
            //   Gvines.Vrash(2, a, bvozm, ref XN);
            //   foreach (double x in XN) Console.WriteLine($"{x} ");


            // Levere.KrilovMethod();
            //QrAlgorithm.MatrixSolver();

            Complex[,] matrix = {
            { 0.42016,-16.937,10.087,-2.8570  },
            { 0.19439,-7.6571,4.5605,-1.3218  },
            { -6.1729,2.8952,-1.7253, 0.41974 },
            {-0.20038,7.8932,-4.7011, 1.3625  }
        };

            Complex[] matrix2 = new Complex[2];

            SobstIterac2.Solve();

           

        }
       
        static void VivodMatr(Matrix A1, int N)
        {
            for (int i = 1; i <= N; i++)
            {
                for (int j = 1; j <= N; j++)
                {
                    Console.Write($"{A1[i, j],15:F6} ");
                }
                Console.WriteLine();
            }
            Console.WriteLine();
        }

        
        static void VivodVectr(Vector B1, int N)
        {
            for (int j = 1; j <= N; j++)
            {
                Console.WriteLine($"x{j}={B1[j]:F3}");
            }
            Console.WriteLine();
        }

        // Cholesky method
        static void Holet(Matrix A1, Vector B1, int N)
        {
            Matrix c = new Matrix(N, N + 1);
            Matrix L = new Matrix(N, N + 1);
            Vector y = new Vector(N);

            // Multiplying the matrix by its transpose
            for (int i = 1; i <= N; i++)
            {
                for (int j = 1; j <= N; j++)
                {
                    double summ = 0.0;
                    for (int t = 1; t <= N; t++)
                    {
                        summ += A1[t, j] * A1[t, i];
                    }
                    c[i, j] = summ;
                }
            }

            // Multiplying the right side by the transposed matrix
            for (int i = 1; i <= N; i++)
            {
                for (int j = 1; j <= N; j++)
                {
                    y[i] += A1[j, i] * B1[j];
                }
            }

            for (int i = 1; i <= N; i++)
            {
                for (int j = 1; j <= N; j++)
                {
                    A1[i, j] = c[i, j];
                }
                B1[i] = y[i];
            }

            for (int i = 1; i <= N; i++)
            {
                for (int j = 1; j <= i; j++)
                {
                    double summ = 0;
                    for (int t = 1; t <= j - 1; t++)
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

            for (int i = 1; i <= N; i++)
            {
                L[i, N + 1] = B1[i];
            }

            B1[1] = L[1, N + 1] / L[1, 1];
            for (int i = 2; i <= N; i++)
            {
                for (int j = 1; j <= i - 1; j++)
                {
                    L[i, N + 1] -= L[i, j] * B1[j];
                }
                B1[i] = L[i, N + 1] / L[i, i];
            }

            for (int i = 1; i <= N; i++)
            {
                for (int j = i + 1; j <= N; j++)
                {
                    L[i, j] = L[j, i];
                    L[j, i] = 0;
                }
                L[i, N + 1] = B1[i];
            }

            B1[N] = L[N, N + 1] / L[N, N];
            for (int i = N - 1; i >= 1; i--)
            {
                for (int j = i + 1; j <= N; j++)
                {
                    L[i, N + 1] -= L[i, j] * B1[j];
                }
                B1[i] = L[i, N + 1] / L[i, i];
            }
        }



    }
}
