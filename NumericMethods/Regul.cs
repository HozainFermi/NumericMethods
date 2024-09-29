using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace NumericMethods
{
    public static class Regul
    {
        const int Nn = 2;
        const int M = 3;
        const double Eps = 0.005;
        

        public static void RegulMethod(double[,] a, double[] b, out double[] X, out int p)
        {
            double[,] a1 = new double[Nn, Nn];
            double[,] a2 = new double[Nn, Nn];
            double[] b1 = new double[Nn];
            double[] b2 = new double[Nn];
            double[] x0 = new double[Nn];
            double alfa = 0, s, max;
            int i, j, k, l=0;

            X = new double[Nn];
            p = 0;

            // Initialize b2
            for (i = 0; i < Nn; i++)
            {
                b2[i] = b[i];
            }

            // Calculate a1 and b1
            for (i = 0; i < Nn; i++)
            {
                for (k = 0; k < Nn; k++)
                {
                    s = 0.0;
                    for (j = 0; j < Nn; j++)
                    {
                        s += a[j, i] * a[j, k];
                    }
                    a1[i, k] = s;
                }
            }

            for (i = 0; i < Nn; i++)
            {
                s = 0.0;
                for (j = 0; j < Nn; j++)
                {
                    s += a[j, i] * b[j];
                }
                b1[i] = s;
            }

            alfa = 0;
            k = 0;
           // Vozm(Nn, Eps, ref a2, ref b2);

            do
            {
                alfa += 1e-8;
                k++;
                Array.Copy(a1, a2, a1.Length);
                for (i = 0; i < Nn; i++)
                {
                    a2[i, i] = a1[i, i] + alfa;
                }
                for (i = 0; i < Nn; i++)
                {
                    b2[i] = b1[i] + alfa * x0[i];
                }
                Array.Copy(a2, a1, a2.Length);
                Array.Copy(b2, b1, b2.Length);
               MatrixSolver.SIMQ(Nn, a2, ref b2, ref l);
                Array.Copy(a1, a2, a1.Length);
                Array.Copy(b2, X, b2.Length);
                Array.Copy(X, x0, X.Length);
                Array.Copy(b1, b2, b1.Length);
                MatrixSolver.SIMQ(Nn, a2, ref b2, ref l);
                max = Math.Abs(b2[0] - X[0]);
                for (i = 1; i < Nn; i++)
                {
                    if (Math.Abs(b2[i] - X[i]) > max)
                    {
                        max = Math.Abs(b2[i] - X[i]);
                    }
                }
            } while (max >= Eps);
            p = k;
        }

        private static void Vozm(int nn, double eps, ref double[,] a, ref double[] b)
        {
            for (int i = 0; i < nn; i++)
            {
                b[i] += eps;
            }
        }

       


    }
}
