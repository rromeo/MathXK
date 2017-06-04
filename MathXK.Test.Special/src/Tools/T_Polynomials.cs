//  (C) Copyright Rocco Romeo 2013.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

using Microsoft.VisualStudio.TestTools.UnitTesting;
using System;


namespace MathXK.Test.Tools
{

    /// <summary>
    ///This is a test class for Polynomials
    ///</summary>
    [TestClass]
    public class T_Polynomial
    {

        struct TestPolynomial
        {
            public string FunctionName;
            public double[] Polynomial;
            public Func<double, double> EquivalentFunction;

            public TestPolynomial(string name, double[] polynomial, Func<double, double> equivalentFunction)
            {
                FunctionName = name;
                Polynomial = polynomial;
                EquivalentFunction = equivalentFunction;
            }
        };

        static TestPolynomial[] _Polynomials =  {
            new TestPolynomial("1.0", new double[] { 1.0 }, x => 1.0),
            new TestPolynomial("x",new double[] { 0.0, 1.0 }, x => x),
            new TestPolynomial("x^2",new double[] { 0.0, 0.0, 1.0 },x => x*x),
            new TestPolynomial("x^3", new double[] { 0.0, 0.0, 0.0, 1.0 },x => x * x * x),
            new TestPolynomial("x^4",new double[] { 0.0, 0.0, 0.0, 0.0, 1.0 },x => x * x * x * x),
            new TestPolynomial("1.0+x+x^2+x^3+x^4", new double[] { 1.0, 1.0, 1.0, 1.0, 1.0 }, x => 1.0 + x + x*x + x*x*x + x * x * x * x)
                                               };

        static TestPolynomial[] _EvenPolynomials =  {
            new TestPolynomial("1.0", new double[] { 1.0 }, x => 1.0),
            new TestPolynomial("x^2",new double[] { 0.0, 1.0 },x => x*x),
            new TestPolynomial("x^4",new double[] { 0.0, 0.0, 1.0 },x => x * x * x * x),
            new TestPolynomial("1.0+x^2+x^4", new double[] { 1.0, 1.0, 1.0, }, x => 1.0 + x*x + x * x * x * x)
                                                    };

        static TestPolynomial[] _OddPolynomials =  {
            new TestPolynomial("x",new double[] { 1.0 }, x => x),
            new TestPolynomial("x^3", new double[] { 0.0, 1.0 },x => x * x * x),
            new TestPolynomial("x^5",new double[] { 0.0, 0.0, 1.0 },x => x * x * x * x * x),
            new TestPolynomial("x+x^3+x^5", new double[] { 1.0, 1.0, 1.0}, x => x + x*x*x + x * x * x * x * x)
                                                   };

        // match each entry of _Polynomials 
        static TestPolynomial[] _PolynomialDenominator =  {
            new TestPolynomial("2.0", new double[] { 2.0 }, x => 2.0),
            new TestPolynomial("2*x",new double[] { 0.0, 2.0 }, x => 2.0*x),
            new TestPolynomial("2*x^2",new double[] { 0.0, 0.0, 2.0 },x => 2.0*x*x),
            new TestPolynomial("2*x^3", new double[] { 0.0, 0.0, 0.0, 2.0 },x => 2.0*x * x * x),
            new TestPolynomial("2*x^4",new double[] { 0.0, 0.0, 0.0, 0.0, 2.0 },x => 2.0*x * x * x * x),
            new TestPolynomial("2.0+x+x^2+x^3+x^4", new double[] { 2.0, 1.0, 1.0, 1.0, 1.0 }, x => 2.0 + x + x*x + x*x*x + x * x * x * x)
                                               };


        /// <summary>
        ///A test for EvaluatePolynomial
        ///</summary>
        [TestMethod]
        public void EvaluatePolynomialTest()
        {
            double[] xValues = { -1.0, 0.0, 1.0, 2.0 };

            foreach (var item in _Polynomials) {
                foreach (double x in xValues) {
                    double result = Polynomial.Eval(item.Polynomial, x);
                    double expected = item.EquivalentFunction(x);
                    Assert.AreEqual(expected, result, $"f(x) = {item.FunctionName}; f({x}) = {result}; Expected: {expected}");
                }
            }

        }


        /// <summary>
        ///A test for EvaluateEvenPolynomial
        ///</summary>
        [TestMethod]
        public void EvaluateEvenPolynomialTest()
        {
            double[] xValues = { -1.0, 0.0, 1.0, 2.0 };

            foreach (var item in _EvenPolynomials) {
                foreach (double x in xValues) {
                    double result = Polynomial.EvalEven(item.Polynomial, x);
                    double expected = item.EquivalentFunction(x);
                    Assert.AreEqual(expected, result, $"f(x) = {item.FunctionName}; f({x}) = {result}; Expected: {expected}");
                }
            }

        }


        /// <summary>
        ///A test for EvaluateOddPolynomial
        ///</summary>
        [TestMethod]
        public void EvaluateOddPolynomialTest()
        {
            double[] xValues = { -1.0, 0.0, 1.0, 2.0 };

            foreach (var item in _OddPolynomials) {
                foreach (double x in xValues) {
                    double result = Polynomial.EvalOdd(item.Polynomial, x);
                    double expected = item.EquivalentFunction(x);
                    Assert.AreEqual(expected, result, $"f(x) = {item.FunctionName}; f({x}) = {result}; Expected: {expected}");
                }
            }

        }


        /// <summary>
        ///A test for EvaluateRational
        ///</summary>
        [TestMethod]
        public void EvaluateRationalTest()
        {
            double[] xValues = { -1.0, -0.5, 0.5, 1.0, 2.0 };

            for (int i = 0; i < _Polynomials.Length; i++) {
                TestPolynomial num = _Polynomials[i];
                TestPolynomial den = _PolynomialDenominator[i];
                foreach (double x in xValues) {
                    double result = Polynomial.EvalRational(num.Polynomial, den.Polynomial, x);
                    double expected = num.EquivalentFunction(x) / den.EquivalentFunction(x);
                    Assert.AreEqual(expected, result, $"f(x) = ({num.FunctionName})/({den.FunctionName}); f({x}) = {result}; Expected: {expected}");
                }
            }

        }

        /// <summary>
        ///A test for Evaluate Polynomial with two dimensions
        ///</summary>
        [TestMethod]
        public void EvaluatePolynomial2Test()
        {

            double[][] c = new double[][] {
                new double[] { 1.0 },
                new double[] { 0.0, 1.0 },
                new double[] { 0.0, 0.0, 1.0 },
                new double[] { 0.0, 0.0, 0.0, 1.0 },
                new double[] { 0.0, 0.0, 0.0, 0.0, 1.0 },
                new double[] { 1.0, 1.0, 1.0, 1.0, 1.0 },
            };

            Random rnd = new Random();
            for (int i = 0; i < 100; i++) {

                double x = rnd.NextDouble();
                double y = rnd.NextDouble();


                double a0 = 1;
                double a1 = x;
                double a2 = x*x;
                double a3 = x*x*x;
                double a4 = x*x*x*x;
                double a5 = 1.0 + x + x*x + x*x*x + x*x*x*x;

                double expected = a0 + y*(a1 + y*(a2 + y*(a3 + y*(a4 + y*a5))));
                double result = Polynomial.Eval(c, x, y);
                Assert2.AreNear(expected, result, 2, 
                    ()=> $"Polynomial.Eval(c, x, y) = {result}; Expected: {expected}");
            }

        }

    }
}
