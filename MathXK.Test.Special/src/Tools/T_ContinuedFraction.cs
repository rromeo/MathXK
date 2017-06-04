//  (C) Copyright Rocco Romeo 2013.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

using System;
using System.Collections.Generic;
using System.Diagnostics;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace MathXK.Test.Tools
{


    /// <summary>
    ///This is a test class for ContinuedFractionTest and is intended
    ///to contain all ContinuedFractionTest Unit Tests
    ///</summary>
    [TestClass()]
    public class T_ContinuedFraction
    {


        private TestContext testContextInstance;

        /// <summary>
        ///Gets or sets the test context which provides
        ///information about and functionality for the current test run.
        ///</summary>
        public TestContext TestContext
        {
            get
            {
                return testContextInstance;
            }
            set
            {
                testContextInstance = value;
            }
        }

        #region Additional test attributes
        // 
        //You can use the following additional attributes as you write your tests:
        //
        //Use ClassInitialize to run code before running the first test in the class
        //[ClassInitialize()]
        //public static void MyClassInitialize(TestContext testContext)
        //{
        //}
        //
        //Use ClassCleanup to run code after all tests in a class have run
        //[ClassCleanup()]
        //public static void MyClassCleanup()
        //{
        //}
        //
        //Use TestInitialize to run code before running each test
        //[TestInitialize()]
        //public void MyTestInitialize()
        //{
        //}
        //
        //Use TestCleanup to run code after each test has run
        //[TestCleanup()]
        //public void MyTestCleanup()
        //{
        //}
        //
        #endregion


        #region Test One Array


        struct TestFraction1
        {
            public readonly string FunctionName;
            public readonly Func<double, double[], double> EquivalentFunction;

            public TestFraction1(string name, Func<double, double[], double> equivalentFunction)
            {
                FunctionName = name;
                EquivalentFunction = equivalentFunction;
            }
        };

        static readonly TestFraction1[] _Fractions1 =  {
            new TestFraction1("b0",
                                (double b0,double[] b) => b0),
            new TestFraction1("b0 + 1.0/b[0]",
                                (double b0,double[] b) => b0 + 1.0/b[0]),
            new TestFraction1("b0 + 1.0/(b[0]+1.0/b[1])",
                                (double b0,double[] b) => b0 + 1.0/(b[0]+1.0/b[1])),
            new TestFraction1("b0 + 1.0/(b[0]+1.0/(b[1] + 1.0/b[2]))",
                                (double b0,double[] b) => b0 + 1.0/(b[0]+1.0/(b[1] + 1.0/b[2]))),
            new TestFraction1("b0 + 1.0/(b[0]+1.0/(b[1] + 1.0/(b[2]+1.0/b[3])))",
                                (double b0,double[] b) => b0 + 1.0/(b[0]+1.0/(b[1] + 1.0/(b[2]+1.0/b[3]))))
                                               };



        /// <summary>
        ///A test for Eval(double b0, double[] b)
        ///</summary>
        [TestMethod()]
        public void EvalTestOneArray()
        {

            double[] b0Values = { 0.0, 1.0 };
            double[][] bValues = {
                new double[] { },
                new double[] { 2.0 },
                new double[] { 2.0, 2.0 },
                new double[] { 2.0, 2.0, 2.0 },
                new double[] { 2.0, 2.0, 2.0, 2.0 }
                                 };

            foreach (double b0 in b0Values) {
                for (int i = 0; i < b0Values.Length; i++) {
                    double[] b = bValues[i];
                    TestFraction1 f = _Fractions1[i];
                    double result = ContinuedFraction.Eval(b0, 1.0, b);
                    double expected = f.EquivalentFunction(b0, b);

                    Assert2.AreNear(expected, result, 2);

                }
            }


        }

        #endregion

        #region Test Two Arrays

        struct TestFraction2
        {
            public readonly string FunctionName;
            public readonly Func<double, double[], double[], double> EquivalentFunction;

            public TestFraction2(string name, Func<double, double[], double[], double> equivalentFunction)
            {
                FunctionName = name;
                EquivalentFunction = equivalentFunction;
            }
        };

        static readonly TestFraction2[] _Fractions2 =  {
            new TestFraction2("b0",
                                (double b0, double[] a, double[] b) => b0),
            new TestFraction2("b0 + a[0]/b[0]",
                                (double b0,double[] a, double[] b) => b0 + a[0]/b[0]),
            new TestFraction2("b0 + a[0]/(b[0]+a[1]/b[1])",
                                (double b0,double[] a, double[] b) => b0 + a[0]/(b[0]+a[1]/b[1])),
            new TestFraction2("b0 + a[0]/(b[0]+a[1]/(b[1] + a[2]/b[2]))",
                                (double b0,double[] a, double[] b) => b0 + a[0]/(b[0]+a[1]/(b[1] + a[2]/b[2]))),
            new TestFraction2("b0 + a[0]/(b[0]+a[1]/(b[1] + a[2]/(b[2] + a[3]/b[3]))))",
                                (double b0,double[] a, double[] b) => b0 + a[0]/(b[0]+a[1]/(b[1] + a[2]/(b[2] + a[3]/b[3]))))
                                               };


        /// <summary>
        ///A test for Eval(double b0, double[] a, double[] b)
        ///</summary>
        [TestMethod()]
        public void EvalTestTwoArray()
        {
            double[] b0Values = { 0.0, 1.0 };

            Tuple<double[], double[]>[] abValues = {
                new Tuple<double[], double[]>(  new double[] { },
                                                new double[] { }),
                new Tuple<double[], double[]>(  new double[] { 1.0 },
                                                new double[] { 2.0 }),
                new Tuple<double[], double[]>(  new double[] { 1.0, 1.0 },
                                                new double[] { 2.0, 2.0 }),
                new Tuple<double[], double[]>(  new double[] { 1.0, 1.0, 1.0 },
                                                new double[] { 2.0, 2.0, 2.0 }),
                new Tuple<double[], double[]>(  new double[] { 1.0, 1.0, 1.0 },
                                                new double[] { 2.0, 2.0, 2.0 }),

                // infinite series would =1/(Sqrt(e)-1)
                // Source:http://mathworld.wolfram.com/ContinuedFraction.html
                new Tuple<double[], double[]>(  new double[] { 2.0, 4.0, 6.0, 8.0 },
                                                new double[] { 3.0, 5.0, 7.0, 9.0 }),

                // infinite series would =1/(e-1)
                // Source: http://mathworld.wolfram.com/ContinuedFraction.html
                new Tuple<double[], double[]>(  new double[] { 1.0, 2.0, 3.0, 4.0 },
                                                new double[] { 1.0, 2.0, 3.0, 4.0 })

                                                   };


            foreach (double b0 in b0Values) {
                foreach (var abValue in abValues) {
                    double[] a = abValue.Item1;
                    double[] b = abValue.Item2;
                    Debug.Assert(a.Length == b.Length);
                    TestFraction2 f = _Fractions2[a.Length];
                    double result = ContinuedFraction.Eval(b0, a, b);
                    double expected = f.EquivalentFunction(b0, a, b);

                    Assert2.AreNear(expected, result, 2);

                }
            }
        }

        #endregion


        #region Test Infinite Series Enumerable

        // infinite series 1/(Sqrt(e)-1)
        // Source:http://mathworld.wolfram.com/ContinuedFraction.html

        private static IEnumerable<(double, double)> Series1()
        {
            double an = 2.0;
            double bn = 3.0;
            for (; ; ) {
                yield return (an, bn);
                an += 2.0;
                bn += 2.0;
            }
        }

        // infinite series would =1/(e-1)
        // Source: http://mathworld.wolfram.com/ContinuedFraction.html
        private static IEnumerable<(double, double)> Series2()
        {
            double an = 1.0;
            double bn = 1.0;
            for (; ; ) {
                yield return (an, bn);
                an += 1.0;
                bn += 1.0;
            }
        }


        /// <summary>
        ///A test for Eval
        ///</summary>
        [TestMethod()]
        public void EvalTestInfinite()
        {
            double expected, actual;

            expected = 1.0 / (Math.Sqrt(Math.E) - 1.0);
            // actual = ContinuedFraction.Eval(1.0, Series1(),tolerance,100);
            actual = ContinuedFraction.Eval(1.0, Series1());
            Assert2.AreNear(expected, actual, 10);

            expected = 1.0 / (Math.E - 1.0);
            // actual = ContinuedFraction.Eval(0.0, Series2(),tolerance,100);
            actual = ContinuedFraction.Eval(0.0, Series2());
            Assert2.AreNear(expected, actual, 10);

        }

        #endregion

        #region Test Infinite Series Function

        // infinite series 1/(Sqrt(e)-1)
        // Source:http://mathworld.wolfram.com/ContinuedFraction.html

        private class Series1Function
        {
            double an = 0.0;
            double bn = 1.0;

            public (double, double) Next()
            {
                an += 2.0;
                bn += 2.0;

                return (an, bn);
            }
        }

        // infinite series would =1/(e-1)
        // Source: http://mathworld.wolfram.com/ContinuedFraction.html
        private class Series2Function
        {
            double an = 0.0;
            double bn = 0.0;

            public (double, double) Next()
            {
                an += 1.0;
                bn += 1.0;

                return (an, bn);
            }
        }


        /// <summary>
        ///A test for Eval
        ///</summary>
        [TestMethod()]
        public void EvalTestInfiniteF()
        {
            double expected, actual;

            Series1Function f1 = new Series1Function();
            expected = 1.0 / (Math.Sqrt(Math.E) - 1.0);
            // actual = ContinuedFraction.Eval(1.0, Series1(),tolerance,100);
            actual = ContinuedFraction.Eval(1.0, f1.Next);
            Assert2.AreNear(expected, actual, 10);

            Series2Function f2 = new Series2Function();
            expected = 1.0 / (Math.E - 1.0);
            // actual = ContinuedFraction.Eval(0.0, Series2(),tolerance,100);
            actual = ContinuedFraction.Eval(0.0, f2.Next);
            Assert2.AreNear(expected, actual, 10);

        }

        #endregion


    }
}
