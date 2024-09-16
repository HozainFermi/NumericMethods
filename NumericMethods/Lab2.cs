using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace NumericMethods
{
    public class Lab2
    {

        public static void Yakoby(ref double[,] a, ref double[] b, int n, ref int step)
        {
            double[] x0 = new double[n];
            double[] x = new double[n];
            double e;
            const double eps = 0.01;

            step = 0;
            for (int i = 0; i < n; i++)
            {
                x0[i] = b[i] / a[i, i];
            }

            do
            {
                for (int i = 0; i < n; i++)
                {
                    x[i] = b[i] / a[i, i];
                    for (int j = 0; j < i; j++)
                    {
                        x[i] -= a[i, j] * x0[j] / a[i, i];
                    }
                    for (int j = i + 1; j < n; j++)
                    {
                        x[i] -= a[i, j] * x0[j] / a[i, i];
                    }
                }

                e = 0;
                for (int i = 0; i < n; i++)
                {
                    if (Math.Abs(x[i] - x0[i]) > e)
                    {
                        e = Math.Abs(x[i] - x0[i]);
                    }
                    x0[i] = x[i];
                }
                step++;
            } while (e > eps);

            b = x0;
            Console.WriteLine("Yakoby");
            foreach(double i in b)
            {
                Console.WriteLine($"{i:F3}");
                
            }
            Console.WriteLine($"Количество итераций:{step}");
            Console.WriteLine();
        }

        public static void Zeidel(double[,] a, double[] b, int n, ref int step, double eps)
        {
            double[] x0 = new double[n];
            double[] x = new double[n];
            step = 0;

            // Initial approximation
            for (int i = 0; i < n; i++)
            {
                x0[i] = b[i] / a[i, i];
            }

            double e;
            do
            {
                for (int i = 0; i < n; i++)
                {
                    x[i] = b[i] / a[i, i];
                    for (int j = 0; j < i; j++)
                    {
                        x[i] -= a[i, j] * x[j] / a[i, i];
                    }
                    for (int j = i + 1; j < n; j++)
                    {
                        x[i] -= a[i, j] * x0[j] / a[i, i];
                    }
                }

                e = 0;
                for (int i = 0; i < n; i++)
                {
                    if (Math.Abs(x[i] - x0[i]) > e)
                    {
                        e = Math.Abs(x[i] - x0[i]);
                    }
                    x0[i] = x[i];
                }

                step++;
            } while (e > eps);

            // Update b with the solution
            for (int i = 0; i < n; i++)
            {
                b[i] = x0[i];
            }
            Console.WriteLine("Zeidel");
            foreach (double i in b)
            {
                Console.WriteLine($"{i:F3}");
               
            }
            Console.WriteLine($"Количество итераций:{step}");
            Console.WriteLine();

        }
        //
        public static void ZeidelEx(double[,] a, double[] b, int n, ref int step, double w)
        {
            double[] x0 = new double[n];
            double[] x = new double[n];
            double e;
            const double eps = 0.001;

            step = 0;

            for (int i = 0; i < n; i++)
            {
                x0[i] = b[i] / a[i, i];
            }

            do
            {
                for (int i = 0; i < n; i++)
                {
                    x[i] = w * b[i] / a[i, i] + (1 - w) * x0[i];

                    for (int j = 0; j < i; j++)
                    {
                        x[i] -= w * a[i, j] * x[j] / a[i, i];
                    }

                    for (int j = i + 1; j < n; j++)
                    {
                        x[i] -= w * a[i, j] * x0[j] / a[i, i];
                    }
                }

                e = 0;
                for (int i = 0; i < n; i++)
                {
                    if (Math.Abs(x[i] - x0[i]) > e)
                    {
                        e = Math.Abs(x[i] - x0[i]);
                    }
                    x0[i] = x[i];
                }

                step++;
            } while (e > eps);

            for (int i = 0; i < n; i++)
            {
                b[i] = x0[i];
            }
            
            for(int i=0;i<b.Length;i++)
            {
                Console.Write('\t'+$"x[{i}] "+$"{b[i]:F3}"+'\t'+"|");

            }
            Console.Write($"Количество итераций:{step}  |  W={w}");
            Console.WriteLine();
        }


    }

}

