using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Lab1_cm
{
    static class Approximator
    {
        private static string infinitCountOfSolutionsString = "Infinit number system's solutions";
        private static string noSolutionsString = "No system's solutions";
        private static string methodIsNotApplyableString = "Method is not applyable!";
        private static string incorrectMatrixString = "Matrix has incorrect lengths!";
        private static string matricesCannotBeMultipliedString = "Matricies cannot be multiplied!";

        public static double epsilan = Math.Pow(10, -4);
        public static double epsilanIter = Math.Pow(10, -3);
        public static int maxIterationCount = 10000;

        public static Random rnd = new Random();

        #region Accessory methods for equasion's methods

        #region Random for equasion's methods
        public static double RandomNumber(double a, double b)
        {
            return rnd.Next(Convert.ToInt32(a / epsilan), Convert.ToInt32(b / epsilan)) * epsilan;
        }
        public static double SuitableRandomArgumnent(Condition cond, double a, double b)
        {
            double arg;

            do
            {
                arg = RandomNumber(a, b);
            } while (!cond(arg));

            return arg;
        }
        #endregion

        #region Delegates for equasion's methods
        public delegate double Function(double x);
        public delegate double FunctionTwo(double x, double c);
        public delegate bool Condition(double func);
        #endregion

        private static bool IsGood(Condition cond, Function func, double a, double b)
        {
            try
            {
                for (double i = a; i <= b; i += epsilan)
                {
                    if (cond(func(i)))
                    {
                        return false;
                    }
                }

                return true;
            }
            catch (NotFiniteNumberException)
            {
                return false;
            }
            catch (DivideByZeroException)
            {
                return false;
            }
        }
        private static bool IsSuitable(Function func1, double a, double b)
        {
            return IsGood((double func) => { return double.IsInfinity(func) || double.IsNaN(func); },
                          func1, a, b) &&
                   (IsGood((double func) => { return func >= 0; },
                          func1, a, b) ||
                   IsGood((double func) => { return func <= 0; },
                          func1, a, b));
        }

        private static double MinValue(Function func, double a, double b)
        {
            double min = double.MaxValue;

            for (double x = a; x <= b; x += epsilan)
            {
                double f = Math.Abs(func(x));

                if (f < min) min = f;
            }

            return min;
        }
        private static double MaxValue(Function func, double a, double b)
        {
            double max = double.MinValue;

            for (double x = a; x <= b; x += epsilan)
            {
                double f = Math.Abs(func(x));

                if (f > max) max = f;
            }

            return max;
        }

        #region Funcs for equasion
        private static Function FuncSimpleIteration(Function func1, double a, double b)
        {
            double min = MinValue(func1, a, b),
                   max = MaxValue(func1, a, b);

            int k = IsGood((double func) => { return func <= 0; }, func1, a, b) ? -1 : 1;

            return (double x) => { return k * 2.0 / (min + max); };
        }
        private static Function FuncNuitonMethod(Function func1, double a, double b)
        {
            return (double x) => { return -1.0 / func1(x); };
        }
        private static Function FuncHordMethod(Function func, Function func1, Function func2, double a, double b, double c)
        {
            double funcc = func(c);
            return (double x) => { return -(x - c) / (func(x) - funcc); };
        }
        #endregion

        #region Initial X for equasion
        private static double InitialXSimpleIteration(double a, double b)
        {
            //return RandomNumber(a, b);
            return (a + b) / 2;
        }
        private static double InitialXNuitonMethod(Function func, Function func2, double a, double b)
        {
            return SuitableRandomArgumnent((double arg) => { return func(arg) * func2(arg) < 0; }, a, b);
        }
        private static double InitialXHordMethod(Function func, double a, double b, double c)
        {
            return SuitableRandomArgumnent((double x) => { return func(x) * func(c) < 0; }, a, b);
        }
        #endregion

        private static double GenericIteration(Function func, Function func1, Function ksi, double a, double b, double x0)
        {
            double xn = x0, counter = 0;

            do
            {
                xn = xn + ksi(xn) * func(xn);
                counter++;
            } while (counter < maxIterationCount && !(Math.Abs(func(xn)) / MinValue(func1, a, b) <= epsilanIter));

            if (counter == maxIterationCount) return xn;//  double.PositiveInfinity;

            return xn;
        }

        #endregion

        #region Main methods for equasion

        public static double SimpleIteration(Function func, Function func1, Function func2, double a, double b)
        {
            if (!IsSuitable(func1, a, b)) throw new Exception("Method can not be applied!");

            Function ksi = FuncSimpleIteration(func1, a, b);

            double x0 = InitialXSimpleIteration(a, b);

            return GenericIteration(func, func1, ksi, a, b, x0);
        }
        public static double NewtonMethod(Function func, Function func1, Function func2, double a, double b)
        {
            if (!IsSuitable(func1, a, b) ||
                !IsSuitable(func2, a, b))
                throw new Exception("Method can not be applied!");

            Function ksi = FuncNuitonMethod(func1, a, b);

            double x0 = InitialXNuitonMethod(func, func2, a, b);

            return GenericIteration(func, func1, ksi, a, b, x0);
        }
        public static double ChordsMethod(Function func, Function func1, Function func2, double a, double b)
        {
            if (!IsSuitable(func1, a, b) ||
                !IsSuitable(func2, a, b))
                throw new Exception("Method can not be applied!");

            double c = SuitableRandomArgumnent((double x) => { return func(x) * func2(x) > 0; }, a, b);

            Function ksi = FuncHordMethod(func, func1, func2, a, b, c);

            double x0 = InitialXHordMethod(func, a, b, c);

            return GenericIteration(func, func1, ksi, a, b, x0);
        }
        public static double HalfDivision(Function func, double a, double b)
        {
            if (func(a) * func(b) > 0)
                return double.PositiveInfinity;

            if (Math.Abs(func(a)) <= epsilanIter) return a;
            if (Math.Abs(func(b)) <= epsilanIter) return b;

            double xn = (a + b) / 2;

            while (Math.Abs(func(xn)) > epsilanIter)
            {
                if (func(xn) * func(a) < 0) b = xn;
                else if (func(xn) * func(b) < 0) a = xn;
                xn = (a + b) / 2;
            }

            return xn;
        }

        #endregion

        #region Gauss's method

        public static double[] GaussMethod(double[,] matrix)
        {
            int rowCount = matrix.GetLength(0),
                colCount = matrix.GetLength(1);
            int[] order = new int[colCount - 1];

            for (int i = 0; i < order.Length; i++)
                order[i] = i;

            for (int i = 0; i < rowCount; i++)
            {
                if (matrix[i, i] == 0)
                {
                    int index = FindNoZeroColumn(matrix, i);

                    if (index == -1)
                        throw new Exception(ResponseAboutExceptionalSolution(matrix[i, colCount - 1] == 0));

                    SwapColumns(ref matrix, i, index);

                    int tmp = order[i];
                    order[i] = order[index];
                    order[index] = i;
                }

                MakeCoeficientEqualToOne(ref matrix, i);
                MakeDownEqualToZero(ref matrix, i);
            }

            double[] solution = RetrieveSolution(matrix);
            for (int i = 0; i < colCount - 1; i++)
            {
                double tmp = solution[order[i]];
                solution[order[i]] = solution[i];
                solution[i] = tmp;
            }

            return solution;
        }

        private static int FindNoZeroColumn(double[,] matrix, int startIndex)
        {
            int colCount = matrix.GetLength(1);

            for (int i = startIndex; i < colCount - 1; i++)
            {
                if (matrix[startIndex, i] != 0)
                    return i;
            }

            return -1;
        }
        private static void SwapColumns(ref double[,] matrix, int firstIndex, int secondIndex)
        {
            int rowCount = matrix.GetLength(0),
                colCount = matrix.GetLength(1);

            if (firstIndex < 0 || firstIndex > colCount ||
                secondIndex < 0 || secondIndex > colCount)
                return;

            if (secondIndex == colCount - 1 ||
                firstIndex == colCount - 1)
                return;

            if (firstIndex == secondIndex)
                return;

            for (int i = 0; i < rowCount; i++)
            {
                Swap(ref matrix[i, firstIndex], ref matrix[i, secondIndex]);
            }
        }
        private static void MakeCoeficientEqualToOne(ref double[,] matrix, int row)
        {
            if (CheckIndex(matrix, row))
                return;

            if (matrix[row, row] == 1)
                return;

            double coef = 1 / matrix[row, row];
            int colCount = matrix.GetLength(1);

            for (int i = row; i < colCount; i++)
                matrix[row, i] *= coef;
        }
        private static void MakeDownEqualToZero(ref double[,] matrix, int row)
        {
            if (CheckIndex(matrix, row))
                return;

            int colCount = matrix.GetLength(1),
                rowCount = matrix.GetLength(0);

            for (int i = row + 1; i < rowCount; i++)
            {
                double coef = -matrix[i, row];

                for (int j = row; j < colCount; j++)
                    matrix[i, j] += matrix[row, j] * coef;
            }
        }

        private static bool CheckIndex(double[,] matrix, int index)
        {
            double colCount = matrix.GetLength(0),
                  rowCount = matrix.GetLength(1);

            if (index < 0 || index > colCount || index > rowCount)
                return false;

            return true;
        }

        private static string ResponseAboutExceptionalSolution(bool hasZeroRow)
        {
            return hasZeroRow ? infinitCountOfSolutionsString : noSolutionsString;
        }

        private static double[] RetrieveSolution(double[,] matrix)
        {
            int rowCount = matrix.GetLength(0),
                colCount = matrix.GetLength(1);

            double[] solution = new double[colCount - 1];
            bool[] isInitialize = new bool[colCount - 1];

            for (int i = 0; i < solution.Length; i++)
            {
                solution[i] = 1;
                isInitialize[i] = false;
            }

            for (int i = rowCount - 1; i >= 0; i--)
            {
                int index = FindNoZeroColumn(matrix, i);

                double newSolution = 1;

                for (int j = index + 1; j < colCount - 1; j++)
                {
                    newSolution -= matrix[i, j] * solution[j];
                }

                newSolution += matrix[i, colCount];
                newSolution /= matrix[i, index];

                if (newSolution != solution[index] && isInitialize[i])
                    throw new Exception(ResponseAboutExceptionalSolution(false));

                solution[index] = newSolution;
                isInitialize[i] = true;
            }

            return solution;
        }

        #endregion

        #region Holeckii's method

        public static double[] HoleckiiMethod(double[,] matrix)
        {
            if (!IsApplyable(matrix))
                throw new Exception(methodIsNotApplyableString);

            int colCount = matrix.GetLength(1),
                rowCount = matrix.GetLength(0);

            if (colCount - 1 != rowCount)
                throw new Exception(incorrectMatrixString);

            double[,] low, up;

            DetermineLowAndUpMatrix(matrix, out low, out up);

            double[] f = new double[rowCount];
            for (int i = 0; i < rowCount; i++)
                f[i] = matrix[i, colCount - 1];

            return GetXSolution(up, GetYSolution(low, f));
        }

        private static bool IsApplyable(double[,] matrix)
        {
            int colCount = matrix.GetLength(1),
                rowCount = matrix.GetLength(0);

            double determinant = 0;
            int[] indexes = new int[colCount - 1];

            for (int i = 0; i < indexes.Length; i++)
                indexes[i] = i;

            do
            {
                double newMember = Math.Pow(-1, InversionCount(indexes));
                for (int i = 0; i < indexes.Length; i++)
                    newMember *= matrix[i, indexes[i]];

                determinant += newMember;
            } while (NextSet(ref indexes));


            return determinant != 0;
        }
        private static bool NextSet(ref int[] set)
        {
            int n = set.Length;

            int j = n - 2;

            while (j != -1 && set[j] >= set[j + 1]) j--;

            if (j == -1)
                return false;

            int k = n - 1;

            while (set[j] >= set[k]) k--;

            Swap(ref set[j], ref set[k]);

            int l = j + 1, r = n - 1;

            while (l < r) Swap(ref set[l--], ref set[r--]);

            return true;
        }
        private static void Swap<T>(ref T val1, ref T val2)
        {
            T tmp = val1;
            val1 = val2;
            val2 = tmp;
        }
        private static int InversionCount(int[] indexes)
        {
            int length = indexes.Length;
            int count = 0;

            for (int i = 0; i < length; i++)
                for (int j = i + 1; j < length; j++)
                    if (indexes[i] > indexes[j])
                        count++;

            return count;
        }

        private static void DetermineLowAndUpMatrix(double[,] matrix, out double[,] low, out double[,] up)
        {
            int colCount = matrix.GetLength(1);

            low = new double[colCount - 1, colCount - 1];
            up = new double[colCount - 1, colCount - 1];

            for (int i = 0; i < colCount - 1; i++)
            {
                for (int j = 0; j < colCount - 1; j++)
                {
                    if (i >= j)
                    {
                        low[i, j] = matrix[i, j];
                        for (int k = 1; k < j - 1; k++)
                            low[i, j] -= low[i, k] * up[k, j];
                    }
                    else
                    {
                        up[i, j] = matrix[i, j];
                        for (int k = 1; k < i - 1; k++)
                            up[i, j] -= low[i, k] * up[k, j];

                        up[i, j] /= low[i, i];
                    }

                    if (i == j)
                        up[i, j] = 1;
                }
            }
        }
        private static double[] GetYSolution(double[,] low, double[] f)
        {
            int colCount = low.GetLength(1),
                rowCount = low.GetLength(0);

            double[] solution = new double[colCount];

            for (int i = 0; i < rowCount; i++)
            {
                solution[i] = f[i];
                for (int j = 0; j < i; j++)
                {
                    solution[i] -= low[i, j] * solution[j];
                }
                solution[i] /= low[i, i];
            }

            return solution;
        }
        private static double[] GetXSolution(double[,] up, double[] y)
        {
            int colCount = up.GetLength(1),
                rowCount = up.GetLength(0);

            double[] solution = new double[colCount];

            for (int i = rowCount - 1; i >= 0; i--)
            {
                solution[i] = y[i];
                for (int j = i + 1; j < rowCount; j++)
                {
                    solution[i] -= up[i, j] * solution[j];
                }
                solution[i] /= up[i, i];
            }

            return solution;
        }

        #endregion

       
    }
}
