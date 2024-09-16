using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace NumericMethods;

public class MatrixSolver
{
    public static void SIMQ(int Nn, ref double[,] Aa, ref double[] Bb, ref int Ks)
    {
        const double Eps = 1e-21;
        double Max, U, V;
        int K1;

        // Копирование Bb в последний столбец Aa
        for (int I = 0; I < Nn; I++)
        {
            Aa[I, Nn] = Bb[I];
        }

        for (int I = 0; I < Nn; I++)
        {
            Max = Math.Abs(Aa[I, I]);
            K1 = I;

            // Нахождения максимального элемента в текущем столбе
            for (int L = I + 1; L < Nn; L++)
            {
                if (Math.Abs(Aa[L, I]) > Max)
                {
                    Max = Math.Abs(Aa[L, I]);
                    K1 = L;
                }
            }

            if (Max < Eps)
            {
                Ks = 1;
                goto M1; 
            }
            else
            {
                Ks = 0;
            }

            // Поменять строки если необходимо
            if (K1 != I)
            {
                for (int J = I; J <= Nn; J++)
                {
                    U = Aa[I, J];
                    Aa[I, J] = Aa[K1, J];
                    Aa[K1, J] = U;
                }
            }

            V = Aa[I, I];

            // Нормализовать текущую строку
            for (int J = I; J <= Nn; J++)
            {
                Aa[I, J] /= V;
            }

            // Eliminate the current column in the rows below
            for (int L = I + 1; L < Nn; L++)
            {
                V = Aa[L, I];
                for (int J = I + 1; J <= Nn; J++)
                {
                    Aa[L, J] -= Aa[I, J] * V;
                }
            }
        }

        // Обратная подстановка для решения для Bb
        Bb[Nn - 1] = Aa[Nn - 1, Nn];

        for (int I = Nn - 2; I >= 0; I--)
        {
            Bb[I] = Aa[I, Nn];
            for (int J = I + 1; J < Nn; J++)
            {
                Bb[I] -= Aa[I, J] * Bb[J];
            }
        }

    M1: // вых
        return;
    }
}
