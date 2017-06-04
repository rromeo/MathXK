//  Copyright (c) 2013 Rocco Romeo
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

//  Ported to .NET based on the original:
//      Copyright (c) John Maddock 2005-2006, Boost Software License, Version 1.0

using System;

namespace MathXK
{

    // For calculating fixed polynomials and a lack of an inline declaration, 
    // use these to copy and paste. :-( 
    // 0: c0;
    // 1: c0 + x * c1;
    // 2: c0 + x * (c1 + x * c2);
    // 3: c0 + x * (c1 + x * (c2 + x * c3));
    // 4: c0 + x * (c1 + x * (c2 + x * (c3 + x * c4)));
    // 5: c0 + x * (c1 + x * (c2 + x * (c3 + x * (c4 + x * c5))));
    // 6: c0 + x * (c1 + x * (c2 + x * (c3 + x * (c4 + x * (c5 + x * c6)))));
    // 7: c0 + x * (c1 + x * (c2 + x * (c3 + x * (c4 + x * (c5 + x * (c6 + x * c7))))));
    // 8: c0 + x * (c1 + x * (c2 + x * (c3 + x * (c4 + x * (c5 + x * (c6 + x * (c7 + x * c8)))))));
    // 9: c0 + x * (c1 + x * (c2 + x * (c3 + x * (c4 + x * (c5 + x * (c6 + x * (c7 + x * (c8 + x * c9))))))));
    //10: c0 + x * (c1 + x * (c2 + x * (c3 + x * (c4 + x * (c5 + x * (c6 + x * (c7 + x * (c8 + x * (c9 + x * c10)))))))));
    //11: c0 + x * (c1 + x * (c2 + x * (c3 + x * (c4 + x * (c5 + x * (c6 + x * (c7 + x * (c8 + x * (c9 + x * ( c10 + x * c11))))))))));
    //12: c0 + x * (c1 + x * (c2 + x * (c3 + x * (c4 + x * (c5 + x * (c6 + x * (c7 + x * (c8 + x * (c9 + x * ( c10 + x * (c11 + x * c12)))))))))));
    //13: c0 + x * (c1 + x * (c2 + x * (c3 + x * (c4 + x * (c5 + x * (c6 + x * (c7 + x * (c8 + x * (c9 + x * ( c10 + x * (c11 + x * (c12 + x * c13))))))))))));
    //14: c0 + x * (c1 + x * (c2 + x * (c3 + x * (c4 + x * (c5 + x * (c6 + x * (c7 + x * (c8 + x * (c9 + x * ( c10 + x * (c11 + x * (c12 + x * (c13 + x * c14)))))))))))));

    /// <summary>
    /// Routines for evaluating polynomials
    /// </summary>
    public static class Polynomial
    {

        /// <summary>
        /// Evaluates a polynomial with a given argument
        /// <para>Σ( c[k] * x^k ) k={0, c.Length-1}</para>
        /// </summary>
        /// <param name="c">The coefficient array</param>
        /// <param name="x">The polynomial argument</param>
        public static double Eval(double[] c, double x)
        {

            // Uses Horner's Method
            // see http://en.wikipedia.org/wiki/Horner_scheme

            if (c == null) {
                Policies.ReportDomainError("Requires Coefficients (c) != null");
                return double.NaN;
            }

            if (c.Length == 0)
                return 0;

            int n = c.Length - 1;
            double result = c[n];
            for (int i = n - 1; i >= 0; i--)
                result = result * x + c[i];
            return result;
        }

        /// <summary>
        /// Evaluates the first n terms of polynomial with a given argument
        /// <para>Σ( c[k] * x^k ) k={0, nTerms-1}</para>
        /// </summary>
        /// <param name="c">The coefficient array</param>
        /// <param name="x">The polynomial argument</param>
        /// <param name="nTerms">The maximum number of terms to evaluate</param>
        public static double Eval(double[] c, double x, int nTerms)
        {

            // Uses Horner's Method
            // see http://en.wikipedia.org/wiki/Horner_scheme

            if (c == null) {
                Policies.ReportDomainError("Requires Coefficients (c) != null");
                return double.NaN;
            }
            if (nTerms < 0 || nTerms > c.Length) {
                Policies.ReportDomainError("Requires nTerms in [0, {0}]; nTerms = {1}", c.Length - 1, nTerms);
                return double.NaN;
            }

            if (nTerms == 0)
                return 0;

            int n = nTerms - 1;
            double result = c[n];
            for (int i = n - 1; i >= 0; i--)
                result = result * x + c[i];
            return result;
        }


        /// <summary>
        /// Evaluates n terms of polynomial beginning at <paramref name="startIndex"/> with a specified argument
        /// <para>Σ( c[startIndex+k] * x^k ) k={0, nTerms-1}</para>
        /// </summary>
        /// <param name="c">The coefficient array</param>
        /// <param name="x">The polynomial argument</param>
        /// <param name="startIndex">The starting index</param>
        /// <param name="nTerms">The maximum number of terms to evaluate</param>
        public static double Eval(double[] c, double x, int startIndex, int nTerms)
        {

            // Uses Horner's Method
            // see http://en.wikipedia.org/wiki/Horner_scheme

            if (c == null) {
                Policies.ReportDomainError("Requires Coefficients (c) != null");
                return double.NaN;
            }
            if (nTerms < 0 || nTerms > c.Length) {
                Policies.ReportDomainError("Requires nTerms in [0, {0}]; nTerms = {1}", c.Length - 1, nTerms);
                return double.NaN;
            }
            if (startIndex < 0 || startIndex >= c.Length) {
                Policies.ReportDomainError("Requires start in [0, {0}]; startIndex = {1}", c.Length - 1, startIndex);
                return double.NaN;
            }
            if (startIndex + nTerms > c.Length) {
                Policies.ReportDomainError("Requires start + nTerms in [0, {0}]; startIndex = {1}; nTerms = {2}", c.Length - 1, startIndex, nTerms);
                return double.NaN;
            }

            if (nTerms == 0)
                return 0;

            int n = nTerms - 1;
            double result = c[n];
            for (int i = n - 1; i >= startIndex; i--)
                result = result * x + c[i];
            return result;
        }


        /// <summary>
        /// Evaluates an even polynomial with the given argument
        /// <para>Σ( c[k] * x^(2*k) ) k={0, c.Length-1}</para>
        /// </summary>
        /// <param name="c">The coefficient array</param>
        /// <param name="x">The polynomial argument</param>
        public static double EvalEven(double[] c, double x)
        {
            return Eval(c, x * x);
        }

        /// <summary>
        /// Evaluates an odd polynomial with the given argument
        /// <para>Σ( c[k] * x^(2*k+1) ) k={0, c.Length-1}</para>
        /// </summary>
        /// <param name="c">The coefficient array</param>
        /// <param name="x">The polynomial argument</param>
        public static double EvalOdd(double[] c, double x)
        {
            if (c == null) {
                Policies.ReportDomainError("Requires Coefficients (c) != null");
                return double.NaN;
            }
            if (c.Length == 0)
                return 0;
            // x*( c[0] + c[1]* x^2 + c[2]*x^4 ...)
            return x * Eval(c, x * x);
        }


        /// <summary>
        /// Evaluates a rational polynomial with the given argument
        /// <para>(num[0] + num[1]*x + num[2]*x^2 ...)/(denom[0] + denom[1]*x + denom[2]*x^2 ...)</para>
        /// </summary>
        /// <param name="num">The numerator coefficient array</param>
        /// <param name="denom">The denominator coefficient array</param>
        /// <param name="x">The argument</param>
        /// <remarks>
        /// When numerator and denominator are equal in order there
        /// are some tricks we can use to prevent overflow that might otherwise
        /// occur in direct polynomial evaluation, for example if <paramref name="x"/> is large
        /// </remarks>
        public static double EvalRational(double[] num, double[] denom, double x)
        {

            // Uses Horner's Method
            // see http://en.wikipedia.org/wiki/Horner_scheme

            if (num == null) {
                Policies.ReportDomainError("Requires num != null");
                return double.NaN;
            }
            if (denom == null) {
                Policies.ReportDomainError("Requires denom != null");
                return double.NaN;
            }
            if (denom.Length == 0) {
                Policies.ReportDomainError("Divide by zero. Requires denom != 0");
                return double.NaN;
            }
            if (num.Length == 0)
                return 0;

            if (num.Length != denom.Length) {
                return Polynomial.Eval(num, x) / Polynomial.Eval(denom, x);
            }


            if (Math.Abs(x) <= 1) {
                double z = x;

                int n = num.Length - 1;
                double numX = num[n];
                double denomX = denom[n];
                for (int i = n - 1; i >= 0; i--) {
                    numX = numX * z + num[i];
                    denomX = denomX * z + denom[i];
                }
                return numX / denomX;
            } else {

                // reverse the order
                // (1 + 1/x) / (2 + 2/x) = (x+1)/(2x+2)
                // can only do this when num.Length == denom.Length

                double z = 1 / x;
                double numX = num[0];
                double denomX = denom[0];
                for (int i = 1; i < num.Length; ++i) {
                    numX = numX * z + num[i];
                    denomX = denomX * z + denom[i];
                }
                return numX / denomX;
            }
        }

        /// <summary>
        /// Evaluates a two dimensional polynomial with the given arguments:
        /// <para>Σ( Eval(c[i],x) * y^i ) i={0, c[i].Length-1}</para>
        /// </summary>
        /// <param name="c">The array of coefficient arrays</param>
        /// <param name="x">The inner polynomial argument</param>
        /// <param name="y">The outer polynomial argument</param>
        public static double Eval(double[][] c, double x, double y)
        {
            if (c == null) {
                Policies.ReportDomainError("Requires coefficient array (c) != null");
                return double.NaN;
            }

            if (c.Length == 0)
                return 0;

            int n = c.Length - 1;
            double result = Polynomial.Eval(c[n],x);
            for (int i = n - 1; i >= 0; i--)
                result = result * y + Polynomial.Eval(c[i],x);
            return result;
        }
    }

} // namespace
