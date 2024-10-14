using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Text;
using System.Threading.Tasks;

namespace NumericMethods
{
    using System;

   public class SobstIterac2
    {
        const int MaxSize = 4;
        public delegate void MatrixInput(int n, ref double[,] a);
        public delegate void VectorInput(int n, ref double[] b);

       static double[,] matrix = {
            { 0.42016,-16.937,10.087,-2.8570  },
            { 0.19439,-7.6571,4.5605,-1.3218  },
            { -6.1729,2.8952,-1.7253, 0.41974 },
            {-0.20038,7.8932,-4.7011, 1.3625  }
        };

        public  static void Solve()
        {
            double[,] a7 = new double[MaxSize, MaxSize];
            int n1 = 4;
            double[] y1 = new double[MaxSize];
            int k;

            while (true)
            {
                Console.Clear();
                Console.WriteLine("Wibirite:\n1- wwod matrici\n2- wwod wectora\n3- wichislit\n4- wihod");
                k = int.Parse(Console.ReadLine());

                if (k == 1) Vvod1(n1, ref a7);
                if (k == 2) Vvod2(n1, ref y1);
                if (k == 3) It(matrix, y1, 0.01, n1);
                if (k == 4) break;
            }
        }

        static void Vvod1(int n, ref double[,] a)
        {
            Console.WriteLine($"Wwedite {n * n} elementow");
            for (int i5 = 0; i5 < n; i5++)
            {
                for (int j5 = 0; j5 < n; j5++)
                {
                    a[i5, j5] = double.Parse(Console.ReadLine());
                }
            }
        }

        static void Vvod2(int n, ref double[] b)
        {
            Console.WriteLine($"Wwedite {n} elementow");
            for (int i5 = 0; i5 < n; i5++)
            {
                b[i5] = double.Parse(Console.ReadLine());
            }
        }

        static void Umnm(double[,] aa1, double[,] aa2, int n5, ref double[,] a)
        {
            for (int i1 = 0; i1 < n5; i1++)
            {
                for (int j1 = 0; j1 < n5; j1++)
                {
                    double r1 = 0;
                    for (int i = 0; i < n5; i++)
                    {
                        r1 += aa1[i1, i] * aa2[i, j1];
                    }
                    a[i1, j1] = r1;
                }
            }
        }

        static void Umnv(double[,] a5, double[] b5, int n5, ref double[] b)
        {
            for (int i1 = 0; i1 < n5; i1++)
            {
                double r1 = 0;
                for (int i = 0; i < n5; i++)
                {
                    r1 += a5[i1, i] * b5[i];
                }
                b[i1] = r1;
            }
        }

       public static void It(double[,] a, double[] y0, double eps, int n)
        {
            

            double[] dy = new double[MaxSize];
            double[] b1 = new double[MaxSize];
            double[] b2 = new double[MaxSize];
            double[,] a1 = new double[MaxSize, MaxSize];
            double[,] a2 = new double[MaxSize, MaxSize];
            double y, norm;
            bool b = false;
            int k = 0;

            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                    a1[i, j] = a[i, j];

            Array.Copy(y0, b2, n);

            while (!b && k < 14)
            {
                k++;
                b = true;

                Array.Copy(b2, b1, n);
                Umnv(a1, y0, n, ref b2);

                for (int i = 0; i < n; i++)
                {
                    dy[i] = b2[i] / b1[i];
                }

                y = 0;
                for (int i = 0; i < n; i++)
                {
                    Console.WriteLine($"{i + 1} : {dy[i]:F4}");
                    y += dy[i];
                }
                y /= n;
                Console.WriteLine($"na shage {k} chislo {y}");

                for (int i = 0; i < n; i++)
                    if (Math.Abs(y - dy[i]) > eps) b = false;

                Umnm(a1, a, n, ref a2);
                Array.Copy(a2, a1, n * n);
                for (int i = 0; i < n; i++)
                {
                    for (int j = 0; j < n; j++)
                    {
                        Console.Write($"{a1[i, j]:F4} ");
                    }
                    Console.WriteLine();
                }
                Console.ReadLine();
            }

            norm = 0;
            for (int i = 0; i < n; i++)
            {
                norm += Math.Pow(b2[i], 2);
            }
            norm = Math.Sqrt(norm);
            for (int i = 0; i < n; i++)
            {
                b2[i] /= norm;
                Console.WriteLine($"{b2[i]:F7}");
            }
            Console.ReadLine();
        }
    }

}
