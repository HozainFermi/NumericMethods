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
        }

        public static void Zeidel(ref double[,] a, ref double[] b, int n, ref int step)
        {
            double[] x0 = new double[n + 1];
            double[] x = new double[n + 1];
            double e;
            step = 0;

            for (int i = 1; i <= n; i++)
            {
                x0[i] = b[i] / a[i, i];
            }

            do
            {
                for (int i = 1; i <= n; i++)
                {
                    x[i] = b[i] / a[i, i];
                    for (int j = 1; j < i; j++)
                    {
                        x[i] -= a[i, j] * x[j] / a[i, i];
                    }
                    for (int j = i + 1; j <= n; j++)
                    {
                        x[i] -= a[i, j] * x0[j] / a[i, i];
                    }
                }

                e = 0;
                for (int i = 1; i <= n; i++)
                {
                    if (Math.Abs(x[i] - x0[i]) > e)
                    {
                        e = Math.Abs(x[i] - x0[i]);
                    }
                    x0[i] = x[i];
                }

                step++;
            } while (e > 1e-10); // Assuming eps is 1e-10

            b = x0;
        }



    }

}

