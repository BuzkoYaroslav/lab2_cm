using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Lab2_cm
{
    class Matrix
    {
        private const string matrixWithoutDeterminantString = "Matrix does not have determinant!";
        private const string incorrectMatrixString = "Matrix has incorrect lengths!";
        private const string matriciesCannotBeMultipliedString = "Matricies cannot be multiplied!";
        private const string matriciesCannotBeAddedString = "Maticies cannot be added!";
        private const string matriciesCannotBeSubstrainedString = "Maticies cannot be substrained!";

        double[,] elements;

        public Matrix(double[,] matrix)
        {
            this.elements = matrix;
        }
        public Matrix(double[] vector)
        {
            elements = new double[vector.Length, 1];

            for (int i = 0; i < vector.Length; i++)
                elements[i, 0] = vector[i];
        }

        public double this[int i, int j]
        {
            get
            {
                if (i < 0 || i >= RowsCount ||
                    j < 0 || j >= ColumnsCount)
                    throw new IndexOutOfRangeException();
                return elements[i, j];
            }
            set
            {
                if (i < 0 || i >= RowsCount ||
                    j < 0 || j >= ColumnsCount)
                    throw new IndexOutOfRangeException();
                elements[i, j] = value;
            }
        }

        public int RowsCount
        {
            get
            {
                return elements.GetLength(0);
            }
        }
        public int ColumnsCount
        {
            get
            {
                return elements.GetLength(1);
            }
        }
        public double FirstNorm
        {
            get
            {
                double norm = -1;

                for (int i = 0; i < RowsCount; i++)
                {
                    double sum = 0;
                    for (int j = 0; j < ColumnsCount; j++)
                        sum += Math.Abs(this[i, j]);

                    if (norm < sum)
                        norm = sum;
                }

                return norm;
            }
        }

        public double Determinant
        {
            get
            {
                if (RowsCount != ColumnsCount)
                    throw new Exception(matrixWithoutDeterminantString);
                return CalculateDeterminant();
            }
        }

        private double CalculateDeterminant()
        {
            double determinant = 0;
            int[] indexes = new int[ColumnsCount];

            for (int i = 0; i < indexes.Length; i++)
                indexes[i] = i;

            do
            {
                double newMember = Math.Pow(-1, InversionCount(indexes));
                for (int i = 0; i < indexes.Length; i++)
                    newMember *= this[i, indexes[i]];

                determinant += newMember;
            } while (NextSet(ref indexes));

            return determinant;
        }
        private int InversionCount(int[] indexes)
        {
            int length = indexes.Length;
            int count = 0;

            for (int i = 0; i < length; i++)
                for (int j = i + 1; j < length; j++)
                    if (indexes[i] > indexes[j])
                        count++;

            return count;
        }
        private bool NextSet(ref int[] set)
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

            while (l < r) {
                Swap(ref set[l], ref set[r]);
                l++;
                r--;
            }

            return true;
        }
        private void Swap<T>(ref T val1, ref T val2)
        {
            T tmp = val1;
            val1 = val2;
            val2 = tmp;
        }

        public void SwapColumns(int firstIndex, int secondIndex)
        {
            if (firstIndex == secondIndex)
                return;

            for (int i = 0; i < RowsCount; i++)
            {
                double val = this[i, firstIndex];

                this[i, firstIndex] = this[i, secondIndex];
                this[i, secondIndex] = val;
            }
        }

        public bool IsNonDegenerate()
        {
            return Determinant != 0;
        }

        public Matrix Transposed()
        {
            double[,] transposed = new double[ColumnsCount, RowsCount];

            for (int i = 0; i < ColumnsCount; i++)
                for (int j = 0; j < RowsCount; j++)
                    transposed[i, j] = this[j, i];

            return new Matrix(transposed);
        }

        public static Matrix UnaryMatrix(int size)
        {
            double[,] elements = new double[size, size];

            for (int i = 0; i < size; i++)
                elements[i, i] = 1;

            return new Matrix(elements);
        }

        private static Matrix MultiplyMatricies(Matrix matrix1, Matrix matrix2)
        {
            double[,] result = new double[matrix1.RowsCount, matrix2.ColumnsCount];

            int rowCount = matrix1.RowsCount,
                colCount = matrix2.ColumnsCount;

            int length = matrix1.ColumnsCount;

            for (int i = 0; i < rowCount; i++)
            {
                for (int j = 0; j < colCount; j++)
                {
                    result[i, j] = 0;
                    for (int p = 0; p < length; p++)
                        result[i, j] += matrix1[i, p] * matrix2[p, j];
                }
            }

            return result;
        }
        private static Matrix AddMatricies(Matrix matrix1, Matrix matrix2)
        {
            double[,] matrix = new double[matrix1.RowsCount, matrix1.ColumnsCount];

            for (int i = 0; i < matrix1.RowsCount; i++)
                for (int j = 0; j < matrix1.ColumnsCount; j++)
                    matrix[i, j] = matrix1[i, j] + matrix2[i, j];

            return matrix;
        }

        public static Matrix operator*(double number, Matrix matrix)
        {
            Matrix copyMatrix = (double[,])matrix;

            for (int i = 0; i < copyMatrix.RowsCount; i++)
                for (int j = 0; j < copyMatrix.ColumnsCount; j++)
                    copyMatrix[i, j] *= number;

            return copyMatrix;
        }
        public static Matrix operator*(Matrix matrix, double number)
        {
            return number * matrix;
        }
        public static Matrix operator*(Matrix left, Matrix right)
        {
            if (left.ColumnsCount != right.RowsCount)
                throw new Exception(matriciesCannotBeMultipliedString);

            return MultiplyMatricies(left, right);
        }
        public static Matrix operator+(Matrix left, Matrix right)
        {
            if (left.RowsCount != right.RowsCount ||
                left.ColumnsCount != right.ColumnsCount)
                throw new Exception(matriciesCannotBeAddedString);

            return AddMatricies(left, right);
        }
        public static Matrix operator-(Matrix left, Matrix right)
        {
            if (left.RowsCount != right.RowsCount ||
                left.ColumnsCount != right.ColumnsCount)
                throw new Exception(matriciesCannotBeSubstrainedString);

            return AddMatricies(left, -right);
        }
        public static Matrix operator-(Matrix matrix)
        {
            Matrix copyMatrix = (double[,])matrix;

            for (int i = 0; i < copyMatrix.RowsCount; i++)
                for (int j = 0; j < copyMatrix.ColumnsCount; j++)
                    copyMatrix[i, j] = -copyMatrix[i, j];

            return copyMatrix;
        }

        public static explicit operator double[,](Matrix matrix)
        {
            return matrix.elements;
        }
        public static explicit operator double[](Matrix matrix)
        {
            if (matrix.RowsCount != 1 &&
                matrix.ColumnsCount != 1)
                throw new InvalidCastException();

            double[] result = new double[Math.Max(matrix.RowsCount, matrix.ColumnsCount)];
            for (int i = 0; i < result.Length; i++)
                if (matrix.RowsCount > 1)
                {
                    result[i] = matrix[i, 0];
                } else
                {
                    result[i] = matrix[0, i];
                }

            return result;
        }
        public static implicit operator Matrix(double[,] elements)
        {
            return new Matrix(elements);
        }
        public static implicit operator Matrix(double[] elements)
        {
            return new Matrix(elements);
        }
    }
}
