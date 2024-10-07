using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Numerics;
using System.Text;
using System.Threading.Tasks;

namespace NumericMethods
{


    public class Levere
    {

       


        public static void LevereMethod()
        {
            float[,] A = new float[4, 4]
       {
            {0.42016f,-16.937f,10.087f,-2.8570f},
            {0.19439f,-7.6571f,4.5605f,-1.3218f},
            {-6.1729f,2.8952f,-1.7253f,0.41974f},
            {-0.20038f,7.8932f,-4.7011f,1.3625f}
       };

            float[,] Ak = new float[4, 4]
       {
            {0,0,0,0},
            {0,0,0,0},
            {0,0,0,0},
            {0,0,0,0}
       };



            float[] S = new float[4];
            float[] p = new float[4];
            float[] lambdas = new float[2];

            float summainp = 0;
            Ak = A;
            S[0] = TraceMatr(Ak);
            p[0] = -S[0];

            for (int i = 1; i < 4; i++)
            {
                Ak = Umnozenie(Ak, A);
                S[i] = TraceMatr(Ak);

                summainp += S[i];
                int index = i;
                for (int j = 0; j < i; j++)
                {
                    index = index - 1;
                    summainp += p[j] * S[index];

                }

                p[i] = (-summainp) / (i + 1);
                summainp = 0;
            }
            p[p.Length - 1] = 0.0000580709f;


            lambdas = IteracWithNewton(p);

           // double[] coeff = new double[5] { 5.80709E-05, -1.93418372, 61.2464867, 7.59974, 1 };
           // double[] coeff2 = new double[5] {1, 7.59974, 61.2464867, -1.93418372, 5.80709E-05 };

            Complex[] comproots = FindComplexRoots();
            CalculateEigenvectors(A, lambdas, comproots);

          //  float[,] AminLambdas = A;
          //  float[] B = new float[4] { 0, 0, 0, 0 };
          //  int Ks = 0;


            // for (int i = 0; i < lambdas.Length; i++)
            // {
            //     
            //
            //     for(int q = 0;q < 4; q++)
            //     {
            //         AminLambdas[q, q] = AminLambdas[q, q] - lambdas[i];
            //
            //     }
            //
            //     SIMQ(4,AminLambdas, ref B,ref Ks);
            //
            //     Console.WriteLine("Корни для лямбда"+i);
            //     foreach (float x in B)
            //     {
            //         Console.WriteLine(x);
            //     }
            //
            //     AminLambdas = A;
            //     B =[ 0, 0, 0, 0 ];
            // }
            //

        }

        public static void FadeevMethod()
        {
            float[,] A = new float[4, 4]
    {
            {0.42016f,-16.937f,10.087f,-2.8570f},
            {0.19439f,-7.6571f,4.5605f,-1.3218f},
            {-6.1729f,2.8952f,-1.7253f,0.41974f},
            {-0.20038f,7.8932f,-4.7011f,1.3625f}
    };

            float[,] A_new = new float[4, 4]
       {
            {0,0,0,0},
            {0,0,0,0},
            {0,0,0,0},
            {0,0,0,0}
       };

            float[,] B_new = new float[4, 4]
       {
            {0,0,0,0},
            {0,0,0,0},
            {0,0,0,0},
            {0,0,0,0}
       };
            float[] qs = new float[4];
            float[] ps = new float[4];

            float[] realroots = new float[2];
            Complex[] comproots = new Complex[2];

            A_new = (float[,])A.Clone();

            qs[0] = TraceMatr(A_new);
            for (int j = 0; j < 4; j++) {
            
                A_new[j,j] = A_new[j, j]-qs[0];
            }
            B_new = (float[,])A_new.Clone();
            ps[0] = -qs[0];

            for (int i = 1; i < 4; i++)
            {
                A_new = Umnozenie(A, B_new);
                qs[i]=TraceMatr(A_new);
                qs[i] = qs[i] /( i + 1);
                for (int j = 0; j < 4; j++)
                {

                    A_new[j, j] = A_new[j, j] - qs[i];
                }
                B_new = (float[,])A_new.Clone();
                ps[i]=-qs[i];
            }
            ps[ps.Length - 1] = 0.0000580709f;
            realroots = IteracWithNewton(ps);

            comproots = FindComplexRoots();
            CalculateEigenvectors(A, realroots, comproots);



        }

        public static void KrilovMethod()
        {
            float[,] A = new float[4, 4]
            {
            {0.42016f,-16.937f,10.087f,-2.8570f},
            {0.19439f,-7.6571f,4.5605f,-1.3218f},
            {-6.1729f,2.8952f,-1.7253f,0.41974f},
            {-0.20038f,7.8932f,-4.7011f,1.3625f}
            };

            float[] y0 = new float[4] { 1, 0, 1, 1 };
            float[] y_new = (float[])y0.Clone();
            float[] y_next = new float[4];
            float[] bvect = new float[4];
            
            float[,] ys = new float[4,4];
            float[] vecty;
            float temp = 0;



            for(int i = 0; i < 4; i++)
            {
                y_next = MultiplyMatrixByVector(A,y_new);
                y_new = (float[])y_next.Clone();
                
                for(int j = 0; j < 4; j++)
                {
                    ys[i, j] = y_new[j];
                }
            }

           
            float[,] YMatr = new float[4, 4]
                {
            {0,0,0,0},
            {0,0,0,0},
            {0,0,0,0},
            {0,0,0,0}
                 };
            
                for(int j = 0; j < 4; j++)
                {
                YMatr[j, 3] = y0[j];
                }
  
            for(int i = 0,yi=2;i < 3; i++,yi--)
            {
                for (int j = 0; j < 4; j++)
                {
                    YMatr[j, i] = ys[yi,j];
                    //Console.WriteLine(YMatr[j, i]+" ");
                }
                
            }
            int K = 0;
            for(int j = 0;j < 4; j++)
            {
                bvect[j] = -ys[3,j];
            }
            //float[] characticcoef = 
                SIMQ(4,YMatr, ref bvect,ref K);
            bvect[bvect.Length-1] = 0.0000580709f;

           float[] realroots = IteracWithNewton(bvect);

           Complex[] comproots = FindComplexRoots();
            CalculateEigenvectors(A, realroots, comproots);


        }


        public static float[,] Umnozenie(float[,] perv, float[,] vtor)
        {
            int k = perv.GetLength(0); //!!!
            float[,] result = new float[k, k];

            for (int i = 0; i < k; i++)
            {
                for (int j = 0; j < k; j++)
                {
                    float temp = 0;
                    for (int g = 0; g < k; g++)
                    {
                        temp += perv[i, g] * vtor[g, j];
                    }
                    result[i, j] = temp;
                }

            }
            return result;
        }
        public static float[] MultiplyMatrixByVector(float[,] matrix, float[] vector)
        {
            int rows = matrix.GetLength(0);
            int cols = matrix.GetLength(1);
            float[] result = new float[rows];

            for (int i = 0; i < rows; i++)
            {
                float sum = 0;
                for (int j = 0; j < cols; j++)
                {
                    sum += matrix[i, j] * vector[j];
                }
                result[i] = sum;
            }

            return result;
        }

        public static float TraceMatr(float[,] A)
        {

            float result = 0;

            for (int i = 0; i < A.GetLength(0); i++)
            {
                result += A[i, i];
            }
            return result;

        }

        public static float[] IteracWithNewton(float[] ps)
        {
            float[] result = new float[2];

            double eps = 1e-11;
            float y0 = 0;

            float y_new = 0;
            float y_pre = 0;
            Console.WriteLine("лямбды");
            for (int i = 0; i < 2; i++) {
                y_pre = y0;
                do
                {
                    y_new = y_pre - (float)((Math.Pow(y_pre, 4) + ps[0] * Math.Pow(y_pre, 3) + ps[1] * Math.Pow(y_pre, 2) + ps[2] * y_pre + ps[3]) / (4 * Math.Pow(y_pre, 3) + 3 * ps[0] * Math.Pow(y_pre, 2) + 2 * ps[1] * y_pre + ps[2]));
                    y_pre = y_new;

                } while (y_new - y_pre > eps);

                result[i] = y_new;
                Console.WriteLine($"y{i}=" + y_new);
                y0 += 0.05f;
            }

            return result;
        }


        public static void CalculateEigenvectors(float[,] matrix, float[] realEigenvalues, Complex[] complexEigenvalues)
        {
            int n = matrix.GetLength(0);
            Complex[,] eigenvectors = new Complex[n, n];

            // Вычисляем собственные векторы для действительных собственных значений
            for (int i = 0; i < realEigenvalues.Length; i++)
            {
                float eigenvalue = realEigenvalues[i];
                float[,] A = new float[n, n];

                // Заполняем матрицу A = A - λI
                for (int j = 0; j < n; j++)
                {
                    for (int k = 0; k < n; k++)
                    {
                        if (j == k)
                        {
                            A[j, k] = matrix[j, k] - eigenvalue;
                        }
                        else
                        {
                            A[j, k] = matrix[j, k];
                        }
                    }
                }

                // Находим ядро матрицы A
                float[,] nullspace = FindNullspace(A);

                // Нормализуем собственный вектор
                float norm = (float)Math.Sqrt(nullspace[0, 0] * nullspace[0, 0] + nullspace[1, 0] * nullspace[1, 0] + nullspace[2, 0] * nullspace[2, 0] + nullspace[3, 0] * nullspace[3, 0]);
                eigenvectors[0, i] = new Complex(nullspace[0, 0] / norm, 0);
                eigenvectors[1, i] = new Complex(nullspace[1, 0] / norm, 0);
                eigenvectors[2, i] = new Complex(nullspace[2, 0] / norm, 0);
                eigenvectors[3, i] = new Complex(nullspace[3, 0] / norm, 0);
            }

            // Вычисляем собственные векторы для комплексных собственных значений
            for (int i = 0; i < complexEigenvalues.Length; i++)
            {
                Complex eigenvalue = complexEigenvalues[i];
                float[,] A = new float[n, n];

                // Заполняем матрицу A = A - λI
                for (int j = 0; j < n; j++)
                {
                    for (int k = 0; k < n; k++)
                    {
                        if (j == k)
                        {
                            A[j, k] = (float)matrix[j, k] - (float)eigenvalue.Real;
                        }
                        else
                        {
                            A[j, k] = matrix[j, k];
                        }
                    }
                }

                // Находим ядро матрицы A
                float[,] nullspace = FindNullspace(A);

                // Нормализуем собственный вектор
                float norm = (float)Math.Sqrt(nullspace[0, 0] * nullspace[0, 0] + nullspace[1, 0] * nullspace[1, 0] + nullspace[2, 0] * nullspace[2, 0] + nullspace[3, 0] * nullspace[3, 0]);
                eigenvectors[0, i + realEigenvalues.Length] = new Complex(nullspace[0, 0] / norm, -nullspace[1, 0] / norm);
                eigenvectors[1, i + realEigenvalues.Length] = new Complex(nullspace[1, 0] / norm, nullspace[0, 0] / norm);
                eigenvectors[2, i + realEigenvalues.Length] = new Complex(nullspace[2, 0] / norm, -nullspace[3, 0] / norm);
                eigenvectors[3, i + realEigenvalues.Length] = new Complex(nullspace[3, 0] / norm, nullspace[2, 0] / norm);
            }

            // Выводим результаты
            Console.WriteLine("Собственные векторы:");
            for (int i = 0; i < n; i++)
            {
                if (i < realEigenvalues.Length)
                {
                    Console.WriteLine($"Для собственного значения {realEigenvalues[i]}:\n ({eigenvectors[0, i]},\n {eigenvectors[1, i]},\n{eigenvectors[2, i]},\n{eigenvectors[3, i]})");
                }
                else
                {
                    Console.WriteLine($"Для собственного значения {complexEigenvalues[i - realEigenvalues.Length]}:\n ({eigenvectors[0, i]},\n {eigenvectors[1, i]},\n {eigenvectors[2, i]},\n {eigenvectors[3, i]})");
                }
            }
        }

        private static float[,] FindNullspace(float[,] A)
        {
            int n = A.GetLength(0);
            float[,] nullspace = new float[n, 1];

            // Находим вектор в ядре матрицы A
            for (int i = 0; i < n; i++)
            {
                nullspace[i, 0] = 1;
                for (int j = 0; j < n; j++)
                {
                    nullspace[i, 0] -= A[i, j] * nullspace[j, 0];
                }
            }

            return nullspace;
        }








        public static Complex[] FindComplexRoots()
        {
            // Заданные коэффициенты
            double a4 = 1.0;
            double a3 = 7.59974;
            double a2 = 61.2464867;
            double a1 = -1.93418372;
            double a0 = 5.80709e-05;

            // Вычисляем коэффициенты приведенного квадратного уравнения
            double b2 = a4;
            double b1 = a3 / a4;
            double b0 = (a2 - (a1 * a1) / (4 * a4)) / a4;

            // Вычисляем дискриминант
            double discriminant = b1 * b1 - 4 * b2 * b0;

            // Находим комплексные корни
            Complex root1, root2;
            if (discriminant < 0)
            {
                double realPart = -b1 / (2 * b2);
                double imaginaryPart = Math.Sqrt(-discriminant) / (2 * b2);
                root1 = new Complex(realPart, imaginaryPart);
                root2 = new Complex(realPart, -imaginaryPart);
                Console.WriteLine($"Комплексные корни: {root1} и {root2}");
                return new Complex[2] { root1, root2 };
            }
           
            else
            {
                Console.WriteLine("Нет комплексных корней.");
                return new Complex[] {};
            }
        }

        public static float[] SolveLinearSystem(float[,] matrix, float[] vector)
        {
            int n = matrix.GetLength(0);
            float[] result = new float[n];

            // Прямой ход
            for (int i = 0; i < n; i++)
            {
                // Найти строку с максимальным элементом в i-ом столбце
                int maxRow = i;
                for (int j = i + 1; j < n; j++)
                {
                    if (Math.Abs(matrix[j, i]) > Math.Abs(matrix[maxRow, i]))
                    {
                        maxRow = j;
                    }
                }

                // Поменять местами строки, если необходимо
                if (maxRow != i)
                {
                    float temp = 0;
                    // Поменять местами строки в матрице
                    for (int j = 0; j < n; j++)
                    {
                        temp = matrix[i, j];
                        matrix[i, j] = matrix[maxRow, j];
                        matrix[maxRow, j] = temp;
                    }

                    // Поменять местами элементы в векторе свободных членов
                    temp = vector[i];
                    vector[i] = vector[maxRow];
                    vector[maxRow] = temp;
                }

                // Нормализовать i-ую строку
                float pivot = matrix[i, i];
                for (int j = i; j < n; j++)
                {
                    matrix[i, j] /= pivot;
                }
                vector[i] /= pivot;

                // Элиминировать элементы под главной диагональю
                for (int j = i + 1; j < n; j++)
                {
                    float factor = matrix[j, i];
                    for (int k = i; k < n; k++)
                    {
                        matrix[j, k] -= factor * matrix[i, k];
                    }
                    vector[j] -= factor * vector[i];
                }
            }

            // Обратный ход
            for (int i = n - 1; i >= 0; i--)
            {
                result[i] = vector[i];
                for (int j = i + 1; j < n; j++)
                {
                    result[i] -= matrix[i, j] * result[j];
                }
            }

            return result;
        }

        public static void SIMQ(int Nn, float[,] A, ref float[] Bb, ref int Ks)
        {
            const float Eps = 1e-21f;
            float Max, U, V;
            int K1;

            float[,] Aa = new float[Nn, Nn + 1];

            for (int i = 0; i < Nn; i++)
            {
                Aa[i, Nn] = 0;
            }

            for (int i = 0; i < Nn; i++)
            {
                for (int j = 0; j < Nn; j++)
                {
                    Aa[i, j] = A[i, j];
                    // Console.Write($"{Aa[i,j]:F3} ");

                }
                // Console.WriteLine();
            }



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




}
