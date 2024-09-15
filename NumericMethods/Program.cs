

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
            double[] Bb = { 10.41, 8.36, 7.62 };
            double[,] Aa = new double[3, 4]
            {
                {3.88,3.78,3.45,0},
                {3.00,2.79,2.39,0},
                {2.67,2.37,1.96,0}

            };
            int Ks = 0;

            MatrixSolver.SIMQ(3, ref Aa, ref Bb, ref Ks);

            Console.WriteLine("SIMQ");
            for (int i = 0; i < Bb.Length; i++)
            {
                Console.WriteLine($"x[{i}] = {Bb[i]:F3}");
                Console.WriteLine('\n');
            }

            VivodMatr(A, N);
            Console.WriteLine("metod Holetckogo");
            Holet(A, B, N);
            VivodVectr(B, N);
            //

            double[] By = { 9.7191, 10.5000, 10.9195, 10.9775 };
            double[,] Aaa = new double[4, 4]
            {
                {9.5000,0.0422, 0.0513, 0.0604},
                {0.0278,8.6000,0.0459,0.0550},
                {0.0224,0.0315,7.7000,0.0496},
                {0.0170,0.0261,0.0351,6.8000 }

            };
            int step = 0;


            Lab2.Yakoby(ref Aaa, ref By, 4, ref step);
            step = 0;
            double[] Bz = { 9.7191, 10.5000, 10.9195, 10.9775 };
            Lab2.Zeidel( Aaa, Bz, 4, ref step,0.01);

            Console.WriteLine("Zeidel Ex");
            double w = 0.2;
            double[] Bzx = { 9.7191, 10.5000, 10.9195, 10.9775 };
            do
            {
                Lab2.ZeidelEx(Aaa, Bzx, 4, ref step, w);
                w = Math.Round(w, 1) + 0.2;
                w = Math.Round(w, 1);
                
            }
            while (w < 2.0);

          






        }
        // Procedure for outputting a matrix
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

        // Procedure for outputting a vector
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
