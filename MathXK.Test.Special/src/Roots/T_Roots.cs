//  (C) Copyright Rocco Romeo 2013.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using MathXK.Roots;


namespace MathXK.Test.Roots
{

    [TestClass]
    public class T_Roots
    {

        class TestData
        {
            public readonly string Description;
            public readonly Func<double, double> FuncValue;
            public readonly Func<double, double> FuncDer1;
            public readonly Func<double, double> FuncDer2;

            public readonly double Max; // these bracket the root;
            public readonly double Min; // these bracket the root;

            public readonly double Lmin; // limits for bracketing routine
            public readonly double Lmax; // limits for bracketing routine
            public readonly double Guess; // initial guess
            public readonly double Step; // step size


            public readonly double Expected;
            public readonly double RelTolerance; // tolerance from the real root
            public readonly double AbsTolerance;


            public TestData(string description,
                        Func<double, double> f, Func<double, double> fder, Func<double, double> fder2,
                        double min, double max,
                        double lmin, double lmax,
                        double guess, double step,
                        double expected, double tolerance, double absTolerance)
            {
                Description = description;
                FuncValue = f; FuncDer1 = fder; FuncDer2 = fder2;
                Min = min; Max = max;
                Lmin = lmin; Lmax = lmax;
                Guess = guess; Step = step;
                Expected = expected;
                RelTolerance = tolerance;
                AbsTolerance = absTolerance;
            }

            public (double, double) FuncVal_Der1(double x)
            {
                return (FuncValue(x), FuncDer1(x));
            }

            public ValueTuple<double, double, double> FuncVal_Der1_Der2(double x)
            {
                return (FuncValue(x), FuncDer1(x), FuncDer2(x));
            }

        };

        static readonly double _AbsTolerance = Math.Pow(2, -104); //Math.Pow(2,-52);

        // all these functions can be bracketed ( Sign(f(min))*Sign(f(max)) <= 0 )
        static readonly TestData[] _FunctionData = {
            new TestData("f(x) = x", 
                            x=>x, x=>1, x=>0, 
                            -1.0, 1.0,
                            double.NaN,double.NaN, 
                            1.0,0.1,
                            0.0,1e-14, _AbsTolerance),

            new TestData("f(x) = x*x-16", 
                            x=>x*x-16, x=>2*x, x=>2, 
                            0, 6,
                            double.NaN,double.NaN, 
                            1.0,0.1,
                            4.0,1e-14, 0),

            new TestData("f(x) = (x+2)*(x+1)*(x-1)", 
                            x=>x*x*x+2*x*x-x-2, x=>3*x*x+4*x-1, x=>6*x+4, 
                            0, 2,
                            double.NaN,double.NaN, 
                            0.5,0.1,
                            1.0,1e-14, 0),

            new TestData("f(x) = sin(x)", 
                            x=>Math.Sin(x), x=>Math.Cos(x), x=>-Math.Sin(x),
                            -Math.PI, Math.PI,
                            double.NaN,double.NaN, 
                            Math.PI/4,0.1,
                            0.0,1e-14, _AbsTolerance),

            new TestData("f(x) = cos(x)", 
                            x=>Math.Cos(x), x=>-Math.Sin(x), x=>-Math.Cos(x),
                            0, Math.PI,
                            double.NaN,double.NaN, 
                            Math.PI/4,0.1,
                            Math.PI/2,1e-14, 0)

                                        };

        // all these functions cannot be easily bracketed ( Sign(f(min))*Sign(f(max)) <= 0 )
        // they touch the axis only on a tangent
        static readonly TestData[] testUnbracketedData = {

            new TestData("f(x) = x*x", 
                x=>x*x, x=>2*x, x=>2, 
                -1, 1,
                double.NaN,double.NaN,
                1.0,0.1,
                0.0,1e-14, _AbsTolerance),

            new TestData("f(x) = sin(x) + 1", 
                x=>Math.Sin(x) + 1, x=>Math.Cos(x), x=>-Math.Sin(x),
                -Math.PI, -0.1,
                double.NaN,double.NaN, 
                -2, -1,
                -Math.PI/2,1e-14, 1e-8)

                                        };

        static bool AreNear(double x, double y, double relTolerance, double absTolerance)
        {
            return Math2.AreNearRel(x, y, relTolerance) || Math2.AreNearAbs(x, y, absTolerance);
        }


        [TestMethod]
        public void Bracket_Test()
        {

            foreach (TestData d in _FunctionData) {

                int iterations;
                RootBracketResult b = RootFinder.FindBracket(d.FuncValue, d.Guess, d.Step, FunctionShape.Unknown, d.Lmin, d.Lmax, RootFinder.DefaultIterations, out iterations);
                if (b == null) {
                    throw new AssertFailedException(" Bracket Solve: (" + d.Description + ") got null; Tolerance: " + d.RelTolerance);
                } if (!b.IsValid)
                    throw new AssertFailedException(" Bracket Solve: (" + d.Description + ") got " + b + " Tolerance: " + d.RelTolerance);
            }


        }


        [TestMethod]
        public void BisectionTest()
        {


            foreach (TestData d in _FunctionData) {

                double relTolerance = RootFinder.DefaultRelTolerance;
                double absTolerance = d.AbsTolerance;

                // test bisection
                RootResults results = RootFinder.Bisection(d.FuncValue, d.Min, d.Max, relTolerance, absTolerance, RootFinder.DefaultIterations);

                double solution = (results == null) ? double.NaN : results.SolutionX;
                if (!AreNear(solution, d.Expected, relTolerance, absTolerance))
                    throw new AssertFailedException(" Bisection Solve: (" + d.Description + ") = " + d.Expected + "; got " + solution + " Tolerance: " + d.RelTolerance);
            }



        }


        [TestMethod]
        public void BrentTest()
        {


            foreach (TestData d in _FunctionData) {

                double relTolerance = RootFinder.DefaultRelTolerance;
                double absTolerance = d.AbsTolerance;

                var tol = RootFinder.GetToleranceFunc(RootFinder.DefaultRelTolerance, d.AbsTolerance);

                RootResults results = RootFinder.Brent(d.FuncValue, d.Min, d.Max, relTolerance, absTolerance, RootFinder.DefaultIterations);
                double solution = (results == null) ? double.NaN : results.SolutionX;
                if (!AreNear(solution, d.Expected, relTolerance, absTolerance))
                    throw new AssertFailedException(" Brent Solve: (" + d.Description + ") = " + d.Expected + "; got " + solution + " Tolerance: " + d.RelTolerance);


                results = RootFinder.BrentBracket(d.FuncValue, d.Guess, d.Step, FunctionShape.Unknown, d.Lmin, d.Lmax, relTolerance, absTolerance, RootFinder.DefaultIterations);
                solution = (results == null) ? double.NaN : results.SolutionX;
                if (!AreNear(solution, d.Expected, relTolerance, absTolerance))
                    throw new AssertFailedException(" Brent Solve: (" + d.Description + ") = " + d.Expected + "; got " + solution + " Tolerance: " + d.RelTolerance);
            }



        }

        [TestMethod]
        public void Toms748Test()
        {


            foreach (TestData d in _FunctionData) {

                double relTolerance = RootFinder.DefaultRelTolerance;
                double absTolerance = d.AbsTolerance;

                RootResults results = RootFinder.Toms748(d.FuncValue, d.Min, d.Max, relTolerance, absTolerance, RootFinder.DefaultIterations);
                double solution = (results == null) ? double.NaN : results.SolutionX;
                if (!AreNear(solution, d.Expected, relTolerance, absTolerance))
                    throw new AssertFailedException(" Tom748 Solve: (" + d.Description + ") = " + d.Expected + "; got " + solution + " Tolerance: " + d.RelTolerance);


                results = RootFinder.Toms748Bracket(d.FuncValue, d.Guess, d.Step, FunctionShape.Unknown, d.Lmin, d.Lmax);
                solution = (results == null) ? double.NaN : results.SolutionX;
                if (!AreNear(solution, d.Expected, relTolerance, absTolerance))
                    throw new AssertFailedException(" Tom748 Solve: (" + d.Description + ") = " + d.Expected + "; got " + solution + " Tolerance: " + d.RelTolerance);
            }



        }

        [TestMethod]
        public void NewtonSafeTest()
        {


            foreach (TestData d in _FunctionData) {

                double relTolerance = RootFinder.DefaultRelTolerance;
                double absTolerance = d.AbsTolerance;


                RootResults results = RootFinder.NewtonSafe(d.FuncVal_Der1, d.Guess, d.Min, d.Max, relTolerance, absTolerance, RootFinder.DefaultIterations);
                double solution = (results == null) ? double.NaN : results.SolutionX;
                if (!AreNear(solution, d.Expected, relTolerance, absTolerance))
                    throw new AssertFailedException(" NewtonSafe Solve: (" + d.Description + ") = " + d.Expected + "; got " + solution + " Tolerance: " + d.RelTolerance);

            }


        }


        [TestMethod]
        public void NewtonTest()
        {


            foreach (TestData d in _FunctionData) {

                double relTolerance = RootFinder.DefaultRelTolerance;
                double absTolerance = d.AbsTolerance;

                RootResults results = RootFinder.Newton(d.FuncVal_Der1, d.Guess, d.Lmin, d.Lmax, relTolerance, absTolerance, RootFinder.DefaultIterations);
                double solution = (results == null) ? double.NaN : results.SolutionX;
                if (!AreNear(solution, d.Expected, relTolerance, absTolerance))
                    throw new AssertFailedException(" Newton Solve: (" + d.Description + ") = " + d.Expected + "; got " + solution + " Tolerance: " + d.RelTolerance);
            }


            foreach (TestData d in testUnbracketedData) {

                double relTolerance = RootFinder.DefaultRelTolerance;
                double absTolerance = d.AbsTolerance;

                RootResults results = RootFinder.Newton(d.FuncVal_Der1, d.Guess, d.Lmin, d.Lmax, relTolerance, absTolerance, RootFinder.DefaultIterations);
                double solution = (results == null) ? double.NaN : results.SolutionX;
                if (!AreNear(solution, d.Expected, relTolerance, absTolerance))
                    throw new AssertFailedException(" Newton Solve: (" + d.Description + ") = " + d.Expected + "; got " + solution + " Tolerance: " + d.RelTolerance);
            }

        }


        [TestMethod]
        public void HalleyTest()
        {


            foreach (TestData d in _FunctionData) {

                double relTolerance = RootFinder.DefaultRelTolerance;
                double absTolerance = d.AbsTolerance;

                RootResults results = RootFinder.Halley(d.FuncVal_Der1_Der2, d.Guess, d.Min, d.Max, relTolerance, absTolerance, RootFinder.DefaultIterations);
                double solution = (results == null) ? double.NaN : results.SolutionX;
                if (!AreNear(solution, d.Expected, relTolerance, absTolerance))
                    throw new AssertFailedException(" Halley Solve: (" + d.Description + ") = " + d.Expected + "; got " + solution + " Tolerance: " + d.RelTolerance);
            }


            foreach (TestData d in testUnbracketedData) {

                double relTolerance = RootFinder.DefaultRelTolerance;
                double absTolerance = d.AbsTolerance;

                RootResults results = RootFinder.Halley(d.FuncVal_Der1_Der2, d.Guess, d.Min, d.Max, relTolerance, absTolerance, RootFinder.DefaultIterations);
                double solution = (results == null) ? double.NaN : results.SolutionX;
                if (!AreNear(solution, d.Expected, relTolerance, absTolerance))
                    throw new AssertFailedException(" Halley Solve: (" + d.Description + ") = " + d.Expected + "; got " + solution + " Tolerance: " + d.RelTolerance);
            }

        }


    }
}
