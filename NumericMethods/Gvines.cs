using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace NumericMethods
{
    public class Gvines
    {

        public static void Vrash(int Nn, double[,] Aa, double[] b, ref double[] x)
        {
            int I, J, K;
            double M, L, R;

            // Initialize the first column of Aa with vector b
            for (I = 0; I < Nn; I++)
            {
                Aa[I, 0] = b[I];
            }

            M = 0.0;

            for (I = 0; I < Nn - 1; I++)
            {
                for (K = I + 1; K < Nn; K++)
                {
                    if (Aa[I, I] != 0.0 || Aa[K, I] != 0.0)
                    {
                        M = Math.Sqrt(Aa[I, I] * Aa[I, I] + Aa[K, I] * Aa[K, I]);
                        L = -1.0 * Aa[K, I] / M;
                        M = Aa[I, I] / M;
                    }
                    else
                    {
                        M = 1.0;
                        L = 0.0;
                    }

                    for (J = 0; J < Nn; J++)
                    {
                        R = M * Aa[I, J] - L * Aa[K, J];
                        Aa[K, J] = L * Aa[I, J] + M * Aa[K, J];
                        Aa[I, J] = R;
                    }

                    R = M * Aa[I, 0] - L * Aa[K, 0];
                    Aa[K, 0] = L * Aa[I, 0] + M * Aa[K, 0];
                    Aa[I, 0] = R;
                }
            }

            for (I = Nn - 1; I >= 0; I--)
            {
                M = 0.0;
                for (K = 0; K < Nn - I - 1; K++)
                {
                    M += Aa[I, Nn - K - 1] * Aa[0, Nn - K - 1];
                }
                Aa[0, I] = (Aa[I, 0] - M) / Aa[I, I];
            }

            for (I = 0; I < Nn; I++)
            {
                x[I] = Aa[0, I];
            }
        }

    }
}
