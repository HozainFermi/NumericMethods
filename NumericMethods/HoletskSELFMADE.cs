using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace NumericMethods
{
    public  class HoletskSELFMADE
    {
        static void HoletskSELF(double[,] A, double[] B, int N )
        {
            int i, j, t;
            double[,] c = new double[N, N];
            double[,] L = new double[N, N + 1];
            double[] y = new double[N];
           // double buf, summ;

            // Initialize matrices
            for (i = 0; i < N; i++)
            {
                for (j = 0; j <= N; j++)
                {
                    c[i, j] = 0;
                    L[i, j] = 0;
                    if (j < N) y[i] = 0;
                }
            }

            L[0, 0] = Math.Sqrt(A[0,0]);
            for (i = 1; i < N; i++) 
            {
                L[i, 0] = A[i, 0] / L[0, 0];
            }
            L[1, 1] = Math.Sqrt(A[1, 1] - (L[1, 0] * L[1,0]));


            // 
            for (i = 1; i < N; i++)
            {
                for(j = i+1; j < N; j++)
                {
                    L[j, i] = (A[j, i] - L[j,(i-1)] * L[(j-1),i])/L[(j-1),i];
                    
                }

            }





        }
    }
}
