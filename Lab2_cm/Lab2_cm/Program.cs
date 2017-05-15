using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Lab2_cm
{
    class Program
    {
        static Matrix[] aMatricies = { new Matrix(new double[,] { { 13, 1, -1 }, 
                                                                 { -1, 10, 3 }, 
                                                                 { 2, 0, 8 } }),
                                       new Matrix(new double[,] { { 3, 1, -1 }, 
                                                                  { -5, 1, 3 }, 
                                                                  { 2, 0, 1 } }) };
        static Matrix[] fMatricies = { new Matrix(new double[] { 16, 16, -6}),
                                       new Matrix(new double[] { 6, -6, 1 }) };

        private const string gaussMethodKeySymbol = "G";
        private const string holeckiiMethodKeySymbol = "H";
        private const string simpleIterationKeySymbol = "S";
        private const string endProgramKeySymbol = "E";
              
        static void Main(string[] args)
        {
            while (true)
            {
                try
                {
                    Console.WriteLine("Input the system index: ");
                    int index = Convert.ToInt32(Console.ReadLine());

                    Console.WriteLine("Input operation symbol:\n" +
                                      "{0} - Gauss method,\n" +
                                      "{1} - Holeckii method,\n" +
                                      "{2} - Simple Iteration method,\n" +
                                      "{3} - Exit", gaussMethodKeySymbol,
                                                    holeckiiMethodKeySymbol,
                                                    simpleIterationKeySymbol,
                                                    endProgramKeySymbol);
                    switch (Console.ReadLine())
                    {
                        case endProgramKeySymbol:
                            return;
                        case gaussMethodKeySymbol:
                            Matrix xGauss = Approximator.GaussMethod(aMatricies[index], fMatricies[index]);

                            Console.WriteLine("Gauss method result: ");
                            WriteSolution(xGauss);
                            break;
                        case holeckiiMethodKeySymbol:
                            Matrix xHol = Approximator.HoleckiiMethod(aMatricies[index], fMatricies[index]);

                            Console.WriteLine("Holeckii method result: ");
                            WriteSolution(xHol);
                            break;
                        case simpleIterationKeySymbol:
                            Matrix xSI = Approximator.SimpleIteration(aMatricies[index], fMatricies[index]);

                            Console.WriteLine("Simple iteration method result: ");
                            WriteSolution(xSI);
                            break;
                    }
                }
                catch (Exception e)
                {
                    Console.WriteLine(e.Message);
                }
            }
        }

        static void WriteSolution(Matrix solution)
        {
            double[] sol = (double[])solution;
            for (int i = 0; i < solution.RowsCount; i++)
            {
                Console.WriteLine("x{0} = {1}", i + 1, sol[i]);
            }
        }
    }
}
