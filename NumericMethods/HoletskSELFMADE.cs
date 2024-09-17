using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace NumericMethods
{
    public class HoletskSELFMADE
    {
        static public void HoletskSELF(double[,] A, double[] B, int N)
        {
            double[,] lower = new double[N, N];
            double[,] upper = new double[N, N];

            // Decomposing a matrix
            // into Lower Triangular
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j <= i; j++)
                {
                    double sum = 0;

                    // summation for diagonals
                    if (j == i)
                    {
                        for (int k = 0; k < j; k++)
                            sum += (int)Math.Pow(lower[j, k],
                                                 2);
                        lower[j, j] = (int)Math.Sqrt(
                            A[j, j] - sum);
                    }

                    else
                    {

                        // Evaluating L(i, j)
                        // using L(j, j)
                        for (int k = 0; k < j; k++)
                            sum += (lower[i, k] * lower[j, k]);
                        lower[i, j] = (A[i, j] - sum)
                                      / lower[j, j];
                    }
                }
            }

            // Displaying Lower
            // Triangular and its Transpose
            Console.WriteLine(
                "  Lower Triangular\t   Transpose");
            for (int i = 0; i < N; i++)
            {

                // Lower Triangular
                for (int j = 0; j < N; j++)
                {
                    Console.Write($"{lower[i, j]:F3}" + "\t");
                }
                Console.Write("");

                // Transpose of
                // Lower Triangular
              //  for (int j = 0; j < N; j++)
              //  {
              //      upper[i, j] = lower[j, i];
              //      Console.Write($"{upper[i, j]:F3}" + "\t");
              //      
              //  }
                Console.WriteLine();

            }























                //    int i, j, t;
                //    double[,] c = new double[N, N];
                //    double[,] L = new double[N, N + 1];
                //    double[] y = new double[N];
                //   // double buf, summ;
                //
                //    // Initialize matrices
                //    for (i = 0; i < N; i++)
                //    {
                //        for (j = 0; j <= N; j++)
                //        {
                //            c[i, j] = 0;
                //            L[i, j] = 0;
                //            if (j < N) y[i] = 0;
                //        }
                //    }
                //
                //    L[0, 0] = Math.Sqrt(A[0,0]);
                //    for (i = 1; i < N; i++) 
                //    {
                //        L[i, 0] = A[i, 0] / L[0, 0];
                //    }
                //    L[1, 1] = Math.Sqrt(A[1, 1] - (L[1, 0] * L[1,0]));
                //
                //
                //    // 
                //    for (i = 1; i < N; i++)
                //    {
                //        for(j = i+1; j < N; j++)
                //        {
                //            L[j, i] = (A[j, i] - L[j,(i-1)] * L[(j-1),i])/L[(j-1),i];
                //            
                //        }
                //
                //    }
                //




            
        }
    }
}
