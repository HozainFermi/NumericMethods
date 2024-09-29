using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace NumericMethods
{
    public class ProstIterac
    {
        const int n = 2;
        const double eps = 1e-4;
        static double[,] a = new double[n, n];
        static double[] x = new double[n];
        static double[] f = new double[n];
        static double[] x0 = new double[n];
        static int iter;

        public static void iterac()
        {

            x0[0] = 0;
            x0[1] = 0;
            iter = 0;
            double max;

            do
            {
                for (int i = 0; i < n; i++)
                {
                    for (int j = 0; j < n; j++)
                    {
                        a[i, j] = MatrJacobi(x0, i, j);
                    }
                }

                VivodVectr(x0, 3, 13);
                Console.WriteLine("норма = " + Norma(a, n));
                Console.WriteLine("номер итерации - " + iter);
                Console.WriteLine("=================");

                for (int i = 0; i < n; i++)
                {
                    x[i] = Func(x0, i);
                }

                max = Math.Abs(x[0] - x0[0]);
                for (int i = 1; i < n; i++)
                {
                    if (Math.Abs(x[i] - x0[i]) > max)
                    {
                        max = Math.Abs(x[i] - x0[i]);
                    }
                }

                Array.Copy(x, x0, n);
                iter++;
               // Console.ReadLine();
            } while (max >= eps && iter <= 20);
        }

        // вычисление нормы
        static double Norma(double[,] a, int n)
        {
            double res = 0;
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    res += a[i, j] * a[i, j];
                }
            }
            res = Math.Sqrt(res);
            return res;
        }

        // задание функции
        static double Func(double[] x, int i)
        {
            switch (i)
            {
                case 0:
                    return (1+Math.Sin(x[1]-0.5))/2; // вычисляем значение первой функции
                case 1:
                    return (1.5-Math.Cos(x[0])); // вычисляем значение второй функции
                default:
                    throw new ArgumentOutOfRangeException();
            }
        }

        // построение матрицы Якоби
        static double MatrJacobi(double[] x, int i, int j)
        {
            switch (i)
            {
                case 0:
                    switch (j)
                    {
                        case 0:
                            return Math.Sin(x[0]); // вычисляем значение элемента матрицы Якоби индексами 1,1
                        case 1:
                            return 0; // вычисляем значение элемента матрицы Якоби с индексами 1,2
                    }
                    break;
                case 1:
                    switch (j)
                    {
                        case 0:
                            return 0; // вычисляем значение элемента матрицы Якоби с индексами 2,1
                        case 1:
                            return 0.5 * Math.Cos(x[1]-0.5); // вычисляем значение элемента матрицы Якоби с индексами 2,2
                    }
                    break;
            }
            throw new ArgumentOutOfRangeException();
        }

        // вывод матрицы
        static void VivodMatr(double[,] mat, int N1, int N2)
        {
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    Console.Write(mat[i, j].ToString("F" + N2) + " ");
                }
                Console.WriteLine();
            }
        }

        // вывод вектора
        static void VivodVectr(double[] vector, int N1, int N2)
        {
            for (int j = 0; j < n; j++)
            {
                Console.WriteLine("x" + (j + 1) + "= " + vector[j].ToString("F" + N2));
            }
        }
    
    }

    }


