using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Lab2_cm
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

        public static Matrix GaussMethod(Matrix matrix, Matrix f)
        {
            int[] order = new int[matrix.ColumnsCount];

            for (int i = 0; i < order.Length; i++)
                order[i] = i;

            for (int i = 0; i < matrix.RowsCount; i++)
            {
                if (matrix[i, i] == 0)
                {
                    int index = FindNoZeroColumn(matrix, i);

                    if (index == -1)
                        throw new Exception(ResponseAboutExceptionalSolution(matrix[i, matrix.ColumnsCount - 1] == 0));

                    SwapColumns(ref matrix, i, index);

                    int tmp = order[i];
                    order[i] = order[index];
                    order[index] = i;
                }

                MakeCoeficientEqualToOne(ref matrix, ref f, i);
                MakeDownEqualToZero(ref matrix, ref f, i);
            }

            double[] solution = RetrieveSolution(matrix);
            for (int i = 0; i < matrix.ColumnsCount; i++)
            {
                double tmp = solution[order[i]];
                solution[order[i]] = solution[i];
                solution[i] = tmp;
            }

            return solution;
        }

        private static int FindNoZeroColumn(Matrix matrix, int startIndex)
        {
            for (int i = startIndex; i < matrix.ColumnsCount; i++)
            {
                if (matrix[startIndex, i] != 0)
                    return i;
            }

            return -1;
        }
        private static void SwapColumns(ref Matrix matrix, int firstIndex, int secondIndex)
        {
            if (firstIndex < 0 || firstIndex > matrix.ColumnsCount - 1 ||
                secondIndex < 0 || secondIndex > matrix.ColumnsCount - 1)
                return;

            if (firstIndex == secondIndex)
                return;

            for (int i = 0; i < matrix.RowsCount; i++)
            {
                double val = matrix[i, firstIndex];

                matrix[i, firstIndex] = matrix[i, secondIndex];
                matrix[i, secondIndex] = val;
            }
        }
        private static void MakeCoeficientEqualToOne(ref Matrix matrix, ref Matrix f, int row)
        {
            if (CheckIndex(matrix, row))
                return;

            if (matrix[row, row] == 1)
                return;

            double coef = 1 / matrix[row, row];

            for (int i = row; i < matrix.ColumnsCount; i++)
                matrix[row, i] *= coef;

            f[row, 0] *= coef;
        }
        private static void MakeDownEqualToZero(ref Matrix matrix, ref Matrix f, int row)
        {
            if (CheckIndex(matrix, row))
                return;

            for (int i = row + 1; i < matrix.RowsCount; i++)
            {
                double coef = -matrix[i, row];

                for (int j = row; j < matrix.ColumnsCount; j++)
                    matrix[i, j] += matrix[row, j] * coef;

                f[i, 0] += f[row, 0] * coef;
            }
        }

        private static bool CheckIndex(Matrix matrix, int index)
        {
            if (index < 0 || index > matrix.ColumnsCount || index > matrix.RowsCount)
                return false;

            return true;
        }

        private static string ResponseAboutExceptionalSolution(bool hasZeroRow)
        {
            return hasZeroRow ? infinitCountOfSolutionsString : noSolutionsString;
        }

        private static Matrix RetrieveSolution(Matrix matrix)
        {
            double[] solution = new double[matrix.ColumnsCount];
            bool[] isInitialize = new bool[matrix.ColumnsCount];

            for (int i = 0; i < solution.Length; i++)
            {
                solution[i] = 1;
                isInitialize[i] = false;
            }

            for (int i = matrix.RowsCount - 1; i >= 0; i--)
            {
                int index = FindNoZeroColumn(matrix, i);

                double newSolution = 1;

                for (int j = index + 1; j < matrix.ColumnsCount; j++)
                {
                    newSolution -= matrix[i, j] * solution[j];
                }

                newSolution += matrix[i, matrix.ColumnsCount];
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

        public static Matrix HoleckiiMethod(Matrix matrix, Matrix f)
        {
            if (!matrix.IsNonDegenerate())
                throw new Exception(methodIsNotApplyableString);

            Matrix low, up;

            DetermineLowAndUpMatrix(matrix, out low, out up);

            return GetXSolution(up, GetYSolution(low, f));
        }

        private static void Swap<T>(ref T val1, ref T val2)
        {
            T tmp = val1;
            val1 = val2;
            val2 = tmp;
        }
        private static void DetermineLowAndUpMatrix(Matrix matrix, out Matrix low, out Matrix up)
        {
            low = new double[matrix.RowsCount, matrix.ColumnsCount];
            up = new double[matrix.RowsCount, matrix.ColumnsCount];

            for (int i = 0; i < matrix.RowsCount; i++)
            {
                for (int j = 0; j < matrix.ColumnsCount; j++)
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
        private static Matrix GetYSolution(Matrix low, Matrix f)
        {
            double[] solution = new double[low.ColumnsCount];

            for (int i = 0; i < low.RowsCount; i++)
            {
                solution[i] = f[i, 0];
                for (int j = 0; j < i; j++)
                {
                    solution[i] -= low[i, j] * solution[j];
                }
                solution[i] /= low[i, i];
            }

            return solution;
        }
        private static Matrix GetXSolution(Matrix up, Matrix y)
        {
            double[] solution = new double[up.ColumnsCount];

            for (int i = up.RowsCount - 1; i >= 0; i--)
            {
                solution[i] = y[i, 0];
                for (int j = i + 1; j < up.ColumnsCount; j++)
                {
                    solution[i] -= up[i, j] * solution[j];
                }
                solution[i] /= up[i, i];
            }

            return solution;
        }

        #endregion

        #region Simple Iteration method

        public static Matrix SimpleIteration(Matrix matrix, Matrix f)
        {
            if (!matrix.IsNonDegenerate())
                throw new Exception(methodIsNotApplyableString);

            double alpha = RandomNumber(0, 2.0 / matrix.Transposed().FirstNorm);

            Matrix bMatrix = Matrix.UnaryMatrix(matrix.ColumnsCount) - alpha * matrix.Transposed() * matrix,
                   gVector = alpha * matrix.Transposed() * f;

            Matrix xCurrent = new double[matrix.ColumnsCount],
                   xPrev;

            do
            {
                xPrev = xCurrent;
                xCurrent = bMatrix * xPrev + gVector;
            } while (bMatrix.FirstNorm > 0.5 && bMatrix.FirstNorm < 1 &&
                      bMatrix.FirstNorm / (1 - bMatrix.FirstNorm) * (xCurrent - xPrev).FirstNorm <= epsilan ||
                      (xCurrent - xPrev).FirstNorm <= epsilan);

            return xCurrent;
        }

        #endregion
    }
}
