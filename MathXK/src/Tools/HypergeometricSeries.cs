//  Copyright (c) 2013 Rocco Romeo
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#define CMP_EQUAL // uses x == y
//#define CMP_FLOATING // uses |x| < |y| * _Tolerance 

using System;
using System.Diagnostics;

namespace MathXK
{


    // See: http://dlmf.nist.gov/16.2
    // R = (b1+b2+b3..bp) - (a1+a2+a3...aq)
    // 
    // Radius of convergence for the functions pFq: 
    // p &lt;= q, converges for finite z
    // p = q + 1,   on the circle |z| == 1, 
    //      if R &gt; 1, absolutely convergent at |z| &lt; 1,
    //      if -1 &lt; R &lt;= 0,  conditionally convergent if |z| != 1
    //      if R &lt;= -1, divergent



    /// <summary>
    /// These functions provide simple hypergeometric sums with minimal optimizations; 
    /// they are not Hypergeometric functions.
    /// </summary>
    internal static partial class HypergeometricSeries
    {

        // In the implementation of these routines, use
        // parameters a1 = a1 + n instead of a1++; b1 = b1+n instead of b1++, etc
        // so the routines have a better chance of working when any a_n or b_n > 1/ϵ



        /// <summary>
        /// Returns the sum of a hypergeometric series 0F0. 
        /// <para>Sum0F0 = addend + Σ( (z^n/n!)) n={0,∞} = e^z</para>
        /// </summary>
        /// <param name="z">The argument</param>
        /// <param name="addend">The addend to the series</param>
        /// <returns>
        /// The sum of the series, or NaN if the series does not converge
        /// within the policy maximum number of iterations. 
        /// </returns>
        public static double Sum0F0(double z, double addend = 0)
        {
            if (double.IsNaN(z)
            || double.IsNaN(addend)) {
                Policies.ReportDomainError("Sum0F0(z: {0}, addend: {1}): NaN not allowed", z, addend);
                return double.NaN;
            }

            double sum = 1.0 + addend;
            if (z == 0)
                return sum;

            double term = 1.0;
            int n = 0;
            for (; n < Policies.MaxSeriesIterations; n++) {

                double prevSum = sum;
                term *= (z / (n + 1));
                sum += term;

#if CMP_FLOAT
                if (Math.Abs(term) <= Math.Abs(prevSum) * Policies.SeriesTolerance)
                    return sum;
#else
                if (sum == prevSum)
                    return sum;
#endif

            }

            Policies.ReportConvergenceError("Sum0F0(z: {0}, addend: {1}): No convergence after {2} iterations", z, addend, n);
            return double.NaN;
        }



        /// <summary>
        /// Returns the sum of a hypergeometric series 0F1. 
        /// <para>Sum0F1 = addend + Σ( (z^n/n!)/b1_n ) n={0,∞}</para>
        /// </summary>
        /// <param name="b1">First denominator. Requires b1 &gt; 0</param>
        /// <param name="z">The argument</param>
        /// <param name="addend">The addend to the series</param>
        /// <returns>
        /// The sum of the series, or NaN if the series does not converge
        /// within the policy maximum number of iterations 
        /// </returns>
        public static double Sum0F1(double b1, double z, double addend = 0)
        {
            if ((double.IsNaN(z) || double.IsInfinity(z))
            || (double.IsNaN(b1) || double.IsInfinity(b1))
            || double.IsNaN(addend)) {
                Policies.ReportDomainError("Sum0F1(b1: {0}, z: {1}, addend: {2}): Requires finite arguments", b1, z, addend);
                return double.NaN;
            }
            if (b1 <= 0 && Math2.IsInteger(b1)) {
                Policies.ReportDomainError("Sum0F1(b1: {0}, z: {1}, addend: {2}): Requires b1 is not zero or a negative integer", b1, z, addend);
                return double.NaN;
            }

            double sum = 1.0 + addend;
            if (z == 0)
                return sum;

            double term = 1.0;
            int n = 0;
            for (; n < Policies.MaxSeriesIterations; n++) {

                double prevSum = sum;
                term *= (z / (n + 1)) / (b1 + n);
                sum += term;

#if CMP_FLOAT
                if (Math.Abs(term) <= Math.Abs(prevSum) * Policies.SeriesTolerance)
                    return sum;
#else
                if (sum == prevSum)
                    return sum;
#endif

            }

            Policies.ReportConvergenceError("Sum0F1(b1: {0}, z: {1}, addend: {2}): No convergence after {3} iterations", b1, z, n);
            return double.NaN;
        }

        /// <summary>
        /// Returns the sum of a hypergeometric series 1F1. 
        /// <para>Sum1F0 = addend + Σ( (a1)_n/n! * z^n) n={0,∞}</para>
        /// </summary>
        /// <param name="a1">First numerator</param>
        /// <param name="z">Argument</param>
        /// <param name="addend">The addend to the series</param>
        /// <returns>
        /// The sum of the series, or NaN if the series does not converge
        /// within the policy maximum number of iterations 
        /// </returns>
        public static double Sum1F0(double a1, double z, double addend = 0)
        {
            if ((double.IsNaN(a1) || double.IsInfinity(a1))
            || (double.IsNaN(z) || double.IsInfinity(z))
            || double.IsNaN(addend)) {
                Policies.ReportDomainError("Sum1F0(a1: {0}, z: {1}, addend: {2}): Requires finite arguments", a1, z, addend);
                return double.NaN;
            }

            double sum = 1.0 + addend;


            if (z == 0)
                return sum;
            if (a1 == 1) {
                if (Math.Abs(z) < 1)
                    return sum - z / (1 - z); // 1-1/(1-z)
                if (z >= 1) {
                    // Divergent Series
                    return double.PositiveInfinity;
                }

                Policies.ReportDomainError("Sum1F0(a1: {0}, z: {1}, addend: {2}): Divergent Series", a1, z, addend);
                return double.NaN;
            }


            double term = 1.0;
            int n = 0;
            for (; n < Policies.MaxSeriesIterations; n++) {

                double prevSum = sum;
                term *= ((a1 + n) / (n + 1)) * z;
                sum += term;

#if CMP_FLOAT
                if (Math.Abs(term) <= Math.Abs(prevSum) * Policies.SeriesTolerance)
                    return sum;
#else
                if (sum == prevSum)
                    return sum;
#endif

            }

            Policies.ReportConvergenceError("Sum1F0(a1: {0}, z: {1}, addend: {2}): No convergence after {3} iterations", a1, z, n);
            return double.NaN;
        }

        /// <summary>
        /// Optimization of Sum1F1 when a1 == 1. 
        /// <para>Sum1F1_A1 = addend + Σ( z^k/(b1)_k) k={0, ∞}</para>
        /// </summary>
        /// <param name="b1"></param>
        /// <param name="z"></param>
        /// <param name="addend">The addend to the series</param>
        /// <returns></returns>
        private static double Sum1F1_A1(double b1, double z, double addend)
        {

            double term = 1.0;
            double sum = 1.0 + addend;

            int n = 0;
            for (; n < Policies.MaxSeriesIterations; n++) {

                double prevSum = sum;
                term *= (z / (b1 + n));
                sum += term;

#if CMP_FLOAT
                if (Math.Abs(term) <= Math.Abs(prevSum) * Policies.SeriesTolerance)
                    return sum;
#else
                if (sum == prevSum)
                    return sum;
#endif
            }

            Policies.ReportConvergenceError("Sum1F1(a1: {0}, b1: {1}, z: {2}, addend: {3}): No convergence after {4} iterations", 1, b1, z, addend, n);
            return double.NaN;

        }

        /// <summary>
        /// Returns the sum of a hypergeometric series 1F1. 
        /// <para>Sum1F1 = addend + Σ( (a1)_n/(b1)_n) * (z^n/n!)) n={0,∞}</para>
        /// </summary>
        /// <param name="a1">First numerator</param>
        /// <param name="b1">First denominator. Requires b1 &gt; 0</param>
        /// <param name="z">Argument</param>
        /// <param name="addend">The addend to the series</param>
        /// <returns>
        /// The sum of the series, or NaN if the series does not converge
        /// within the policy maximum number of iterations 
        /// </returns>
        /// <remarks>
        /// The cases of a1 == b1 and a1 == 1 are optimized
        /// </remarks>
        public static double Sum1F1(double a1, double b1, double z, double addend = 0)
        {
            if ((double.IsNaN(a1) || double.IsInfinity(a1))
            || (double.IsNaN(z) || double.IsInfinity(z))
            || (double.IsNaN(b1) || double.IsInfinity(b1))
            || double.IsNaN(addend)) {
                Policies.ReportDomainError("Sum1F1(a1: {0}, b1: {1}, z: {2}, addend: {3}): Requires finite arguments", a1, b1, z, addend);
                return double.NaN;
            }
            if (b1 <= 0 && Math2.IsInteger(b1)) {
                Policies.ReportDomainError("Sum1F1(a1: {0}, b1: {1}, z: {2}, addend: {3}): Requires b1 is not zero or a negative integer", a1, b1, z, addend);
                return double.NaN;
            }

            double sum = 1.0 + addend;

            if (z == 0)
                return sum;
            if (a1 == b1)
                return Sum0F0(z, addend);
            if (a1 == 1)
                return Sum1F1_A1(b1, z, addend);

            double term = 1.0;
            int n = 0; 
            for (; n < Policies.MaxSeriesIterations; n++) {

                double prevSum = sum;
                term *= (z / (n + 1)) * ((a1 + n) / (b1 + n));
                sum += term;

#if CMP_FLOAT
                if (Math.Abs(term) <= Math.Abs(prevSum) * Policies.SeriesTolerance)
                    return sum;
#else
                if (sum == prevSum)
                    return sum;
#endif
            }

            Policies.ReportConvergenceError("Sum1F1(a1: {0}, b1: {1}, z: {2}, addend: {3}): No convergence after {4} iterations", a1, b1, z, addend, n);
            return double.NaN;
        }

        /// <summary>
        /// Optimization of Sum2F1 when a1 == 1. 
        /// <para>Sum2F1_A1 = addend + Σ( (z^k/(b1)_k)/(b2)_k) k={0, ∞}</para>
        /// </summary>
        /// <param name="b1"></param>
        /// <param name="b2"></param>
        /// <param name="z"></param>
        /// <param name="addend">The addend to the series</param>
        /// <returns></returns>
        private static double Sum1F2_A1(double b1, double b2, double z, double addend)
        {

            double term = 1.0;
            double sum = 1.0 + addend;

            int n = 0; 
            for (; n < Policies.MaxSeriesIterations; n++) {

                double prevSum = sum;
                term *= ((z / (b1 + n)) / (b2 + n));
                sum += term;

#if CMP_FLOAT
                if (Math.Abs(term) <= Math.Abs(prevSum) * Policies.SeriesTolerance)
                    return sum;
#else
                if (sum == prevSum)
                    return sum;
#endif

            }

            Policies.ReportConvergenceError("Sum1F2(a1: {0}, b1: {1}, b2: {2}, z: {3}, addend: {4}): No convergence after {0} iterations", 1, b1, b2, z, addend, n);
            return double.NaN;
        }


        /// <summary>
        /// Returns the sum of a hypergeometric series 1F2. 
        /// <para>Sum2F1 = addend + Σ( ((a1)_n)/((b1)_n * (b2)_n) * (z^n/n!)) n={0,∞}</para>
        /// </summary>
        /// <param name="a1">First numerator</param>
        /// <param name="b1">First denominator. Requires b1 &gt; 0</param>
        /// <param name="b2">Second denominator. Requires b2 &gt; 0</param>
        /// <param name="z">Argument</param>
        /// <param name="addend">The addend to the series</param>
        /// <returns>
        /// The sum of the series, or NaN if the series does not converge
        /// within the policy maximum number of iterations 
        /// </returns>
        public static double Sum1F2(double a1, double b1, double b2, double z, double addend = 0)
        {
            if ((double.IsNaN(a1) || double.IsInfinity(a1))
            || (double.IsNaN(z) || double.IsInfinity(z))
            || (double.IsNaN(b1) || double.IsInfinity(b1))
            || (double.IsNaN(b2) || double.IsInfinity(b2))
            || double.IsNaN(addend)) {
                Policies.ReportDomainError("Sum1F2(a1: {0}, b1: {1}, b2: {2}, z: {3}, addend: {4}): Requires finite arguments", a1, b1, b2, z, addend);
                return double.NaN;
            }
            if ((b1 <= 0 && Math2.IsInteger(b1))
                || (b2 <= 0 && Math2.IsInteger(b2))) {
                Policies.ReportDomainError("Sum1F2(a1: {0}, b1: {1}, b2: {2}, z: {3}, addend: {4}): Requires b1, b2 are not zero or negative integers", a1, b1, b2, z, addend);
                return double.NaN;
            }

            double sum = 1.0 + addend;

            if (z == 0)
                return sum;
            if (a1 == b1)
                return Sum0F1(b2, z, addend);
            if (a1 == b2)
                return Sum0F1(b1, z, addend);
            if (a1 == 1)
                return Sum1F2_A1(b1, b2, z, addend);

            double term = 1.0;
            int n = 0; 
            for (; n < Policies.MaxSeriesIterations; n++) {

                double prevSum = sum;
                term *= (z / (n + 1)) * (((a1 + n) / (b1 + n)) / (b2 + n));
                sum += term;

#if CMP_FLOAT
                if (Math.Abs(term) <= Math.Abs(prevSum) * Policies.SeriesTolerance)
                    return sum;
#else
                if (sum == prevSum)
                    return sum;
#endif

            }

            Policies.ReportConvergenceError("Sum1F2(a1: {0}, b1: {1}, b2: {2}, z: {3}, addend: {4}): No convergence after {0} iterations", a1, b1, b2, z, addend, n);
            return double.NaN;
        }


        /// <summary>
        /// Returns the sum of a hypergeometric series 2F0 when a1 == 1 or a2 == 2. 
        /// <para>Sum2F1_A1 = addend + Σ( (a1)_n * z^n ) n={0,∞}</para>
        /// </summary>
        /// <param name="a1">First numerator</param>
        /// <param name="z">Argument</param>
        /// <param name="addend">The addend to the series</param>
        /// <returns>
        /// The sum of the series, or NaN if the series does not converge
        /// within the policy maximum number of iterations 
        /// </returns>
        private static double Sum2F0_A1(double a1, double z, double addend)
        {
            double term = 1.0;
            double sum = 1.0 + addend;


            int n = 0; 
            for (; n < Policies.MaxSeriesIterations; n++) {

                double prevSum = sum;
                term *= z * (a1 + n);
                sum += term;

#if CMP_FLOAT
                if (Math.Abs(term) <= Math.Abs(prevSum) * Policies.SeriesTolerance)
                    return sum;
#else
                if (sum == prevSum)
                    return sum;
#endif

            }

            Policies.ReportConvergenceError("Sum2F0(a1: {0}, a2: {1}, z: {2}, addend: {3}): No converge after {4} iterations", 1, a1, z, addend, n);
            return double.NaN;
        }


        /// <summary>
        /// Returns the sum of a hypergeometric series 2F0. 
        /// <para>Sum2F0 = addend + Σ( (a1)_n/n! * (a2)_n * z^n) n={0,∞}</para>
        /// </summary>
        /// <param name="a1">First numerator</param>
        /// <param name="a2">Second numerator</param>
        /// <param name="z">Argument</param>
        /// <param name="addend">The addend to the series</param>
        /// <returns>
        /// The sum of the series, or NaN if the series does not converge
        /// within the policy maximum number of iterations 
        /// </returns>
        public static double Sum2F0(double a1, double a2, double z, double addend = 0)
        {
            if ((double.IsNaN(a1) || double.IsInfinity(a1))
            || (double.IsNaN(a2) || double.IsInfinity(a2))
            || (double.IsNaN(z) || double.IsInfinity(z))
            || double.IsNaN(addend)) {
                Policies.ReportDomainError("Sum2F0(a1: {0}, a2: {1}, z: {2}, addend: {3}): Requires arguments", a1, a2, z, addend);
                return double.NaN;
            }

            double sum = 1.0 + addend;

            if (z == 0)
                return sum;
            if (a1 == 1)
                return Sum2F0_A1(a2, z, addend);
            if (a2 == 1)
                return Sum2F0_A1(a1, z, addend);

            double term = 1.0;

            int n = 0; 
            for (; n < Policies.MaxSeriesIterations; n++) {

                double prevSum = sum;
                term *= z * ((a1 + n) / (n + 1)) * (a2 + n);
                sum += term;

#if CMP_FLOAT
                if (Math.Abs(term) <= Math.Abs(prevSum) * Policies.SeriesTolerance)
                    return sum;
#else
                if (sum == prevSum)
                    return sum;
#endif

            }

            Policies.ReportConvergenceError("Sum2F0(a1: {0}, a2: {1}, z: {2}, addend: {3}): No converge after {4} iterations", a1, a2, z, addend, n);
            return double.NaN;
        }


        /// <summary>
        /// Optimization of Sum2F1 when a1 == 1. 
        /// <para>Sum2F1_A1 = addend + Σ( z^k * ((a1)_k/(b1)_k)) k={0, ∞}</para>
        /// </summary>
        /// <param name="a1">First numerator</param>
        /// <param name="b1">Second numerator</param>
        /// <param name="z">Argument</param>
        /// <param name="addend">The addend to the series</param>
        /// <returns>
        /// The sum of the series, or NaN if the series does not converge
        /// within the policy maximum number of iterations 
        /// </returns>
        static double Sum2F1_A1(double a1, double b1, double z, double addend)
        {

            double term = 1.0;
            double sum = 1.0 + addend;

            int n = 0; 
            for (; n < Policies.MaxSeriesIterations; n++) {

                double prevSum = sum;
                term *= z * ((a1 + n) / (b1 + n));
                sum += term;

#if CMP_FLOAT
                if (Math.Abs(term) <= Math.Abs(prevSum) * Policies.SeriesTolerance)
                    return sum;
#else
                if (sum == prevSum)
                    return sum;
#endif
            }

            Policies.ReportConvergenceError("Sum2F1(a1: {0}, a2: {1}, b1: {2}, z: {3}, addend: {4}): No convergence after {5} iterations", 1, a1, b1, z, addend, n);
            return double.NaN;
        }



        /// <summary>
        /// Returns the sum of a hypergeometric series 2F1. 
        /// <para>Sum2F1 = addend + Σ(( (a1)_n * (a2)_n)/(b1)_n * (z^n/n!)) n={0,∞}</para>
        /// </summary>
        /// <param name="a1">First numerator</param>
        /// <param name="a2">Second numerator</param>
        /// <param name="b1">First denominator. Requires b1 > 0</param>
        /// <param name="z">Argument</param>
        /// <param name="addend">The addend to the series</param>
        /// <returns>
        /// The sum of the series, or NaN if the series does not converge
        /// within the policy maximum number of iterations 
        /// </returns>
        public static double Sum2F1(double a1, double a2, double b1, double z, double addend = 0)
        {
            if ((double.IsNaN(a1) || double.IsInfinity(a1))
            || (double.IsNaN(a2) || double.IsInfinity(a2))
            || (double.IsNaN(z) || double.IsInfinity(z))
            || (double.IsNaN(b1) || double.IsInfinity(b1))
            || double.IsNaN(addend)) {
                Policies.ReportDomainError("Sum2F1(a1: {0}, a2: {1}, b1: {2}, z: {3}, addend: {4}): Requires finite arguments", a1, a2, b1, z, addend);
                return double.NaN;
            }
            if (b1 <= 0 && Math2.IsInteger(b1)) {
                Policies.ReportDomainError("Sum2F1(a1: {0}, a2: {1}, b1: {2}, z: {3}, addend: {4}): Requires b1 is not zero or a negative integer", a1, a2, b1, z, addend);
                return double.NaN;
            }

            double sum = 1.0 + addend;

            if (z == 0)
                return sum;
            if (a1 == b1)
                return Sum1F0(a2, z, addend);
            if (a2 == b1)
                return Sum1F0(a1, z, addend);
            if (a1 == 1)
                return Sum2F1_A1(a2, b1, z, addend);
            if (a2 == 1)
                return Sum2F1_A1(a1, b1, z, addend);

            double term = 1.0;
            int n = 0; 
            for (; n < Policies.MaxSeriesIterations; n++) {

                double prevSum = sum;
                term *= (z / (n + 1)) * ((a1 + n) / (b1 + n)) * (a2 + n);
                sum += term;

#if CMP_FLOAT
                if (Math.Abs(term) <= Math.Abs(prevSum) * Policies.SeriesTolerance)
                    return sum;
#else
                if (sum == prevSum)
                    return sum;
#endif

            }

            Policies.ReportConvergenceError("Sum2F1(a1: {0}, a2: {1}, b1: {2}, z: {3}, addend: {4}): No convergence after {5} iterations", a1, a2, b1, z, addend, n);
            return double.NaN;
        }


        /// <summary>
        /// Optimization of Sum2F2 when a1 == 1 or a2 == 1. 
        /// <para>Sum2F2_A1 = addend + Σ( ((a1)_n)/((b1)_n * (b2)_n) * (z^n)) n={0,∞}</para>
        /// </summary>
        /// <param name="a1">First numerator</param>
        /// <param name="b1">First denominator. Requires b1 &gt; 0</param>
        /// <param name="b2">Second denominator. Requires b2 &gt; 0</param>
        /// <param name="z">Argument</param>
        /// <param name="addend">The addend to the series</param>
        /// <returns>
        /// The sum of the series, or NaN if the series does not converge
        /// within the policy maximum number of iterations 
        /// </returns>
        static double Sum2F2_A1(double a1, double b1, double b2, double z, double addend)
        {

            double sum = 1.0 + addend;

            double term = 1.0;
            int n = 0; 
            for (; n < Policies.MaxSeriesIterations; n++) {

                double prevSum = sum;
                term *= z * ((a1 + n)/(b1 + n))/(b2 + n);
                sum += term;

#if CMP_FLOAT
                if (Math.Abs(term) <= Math.Abs(prevSum) * Policies.SeriesTolerance)
                    return sum;
#else
                if (sum == prevSum)
                    return sum;
#endif
            }

            Policies.ReportConvergenceError("Sum2F2(a1: {0}, a2: {1}, b1: {2}, b2: {3}, z: {4}, addend: {5}): No converge after {6} iterations", 1, a1, b1, b2, z, addend, n);
            return double.NaN;
        }


        /// <summary>
        /// Returns the sum of a hypergeometric series 2F2. 
        /// <para>Sum2F2 = addend + Σ( ((a1)_n * (a2)_n)/((b1)_n * (b2)_n) * (z^n/n!)) n={0,∞}</para>
        /// </summary>
        /// <param name="a1">First numerator</param>
        /// <param name="a2">Numerator value 2</param>
        /// <param name="b1">First denominator. Requires b1 &gt; 0</param>
        /// <param name="b2">Second denominator. Requires b2 &gt; 0</param>
        /// <param name="z">Argument</param>
        /// <param name="addend">The addend to the series</param>
        /// <returns>
        /// The sum of the series, or NaN if the series does not converge
        /// within the policy maximum number of iterations 
        /// </returns>
        public static double Sum2F2(double a1, double a2, double b1, double b2, double z, double addend = 0)
        {
            if ((double.IsNaN(a1) || double.IsInfinity(a1))
            || (double.IsNaN(a2) || double.IsInfinity(a2))
            || (double.IsNaN(z) || double.IsInfinity(z))
            || (double.IsNaN(b1) || double.IsInfinity(b1))
            || (double.IsNaN(b2) || double.IsInfinity(b2))
            || double.IsNaN(addend)) {
                Policies.ReportDomainError("Sum2F2(a1: {0}, a2: {1}, b1: {2}, b2: {3}, z: {4}, addend: {5}): Requires finite arguments", a1, a2, b1, b2, z, addend);
                return double.NaN;
            }
            if ((b1 <= 0 && Math2.IsInteger(b1))
                || (b2 <= 0 && Math2.IsInteger(b2))) {
                Policies.ReportDomainError("Sum2F2(a1: {0}, a2: {1}, b1: {2}, b2: {3}, z: {4}, addend: {5}): Requires b1, b2 are not zero or negative integers", a1, a2, b1, b2, z, addend);
                return double.NaN;
            }


            double sum = 1.0 + addend;

            if (z == 0)
                return sum;
            if (a1 == b1)
                return Sum1F1(a2, b2, z, addend);
            if (a1 == b2)
                return Sum1F1(a2, b1, z, addend);
            if (a2 == b1)
                return Sum1F1(a1, b2, z, addend);
            if (a2 == b2)
                return Sum1F1(a1, b1, z, addend);
            if (a1 == 1)
                return Sum2F2_A1(a2, b1, b2, z, addend);
            if (a2 == 1)
                return Sum2F2_A1(a1, b1, b2, z, addend);


            double term = 1.0;
            int n = 0; 
            for (; n < Policies.MaxSeriesIterations; n++) {

                double prevSum = sum;
                term *= ((a1 + n) / (b1 + n)) * ((a2 + n) / (b2 + n)) * (z / (n + 1));
                sum += term;

#if CMP_FLOAT
                if (Math.Abs(term) <= Math.Abs(prevSum) * Policies.SeriesTolerance)
                    return sum;
#else
                if (sum == prevSum)
                    return sum;
#endif
            }

            Policies.ReportConvergenceError("Sum2F2(a1: {0}, a2: {1}, b1: {2}, b2: {3}, z: {4}, addend: {5}): No converge after {6} iterations", a1, a2, b1, b2, z, addend, n);
            return double.NaN;
        }


        /// <summary>
        /// Returns the sum of a hypergeometric series 3F0 with one numerator = 1. 
        /// <para>Sum3F0_A1 = addend + Σ(( (a1)_n * (a2)_n ) * (z^n)) n={0,∞}</para>
        /// </summary>
        /// <param name="a1">First numerator</param>
        /// <param name="a2">Second numerator</param>
        /// <param name="z">Argument</param>
        /// <param name="addend">The addend to the series</param>
        /// <returns>
        /// The sum of the series, or NaN if the series does not converge
        /// within the policy maximum number of iterations 
        /// </returns>
        private static double Sum3F0_A1(double a1, double a2, double z, double addend)
        {
            if ((double.IsNaN(a1) || double.IsInfinity(a1))
            || (double.IsNaN(a2) || double.IsInfinity(a2))
            || (double.IsNaN(z) || double.IsInfinity(z))
            || double.IsNaN(addend)) {
                Policies.ReportDomainError("Sum3F1(a1: {0}, a2: {1}, z: {2}, addend: {3}): Requires finite arguments", a1, a2, z, addend);
                return double.NaN;
            }

            double sum = 1.0 + addend;

            double term = 1.0;
            int n = 0; 
            for (; n < Policies.MaxSeriesIterations; n++) {

                double prevSum = sum;
                term *= z * (a1 + n) * (a2 + n);
                sum += term;

#if CMP_FLOAT
                if (Math.Abs(term) <= Math.Abs(prevSum) * Policies.SeriesTolerance)
                    return sum;
#else
                if (sum == prevSum)
                    return sum;
#endif
            }

            Policies.ReportConvergenceError("Sum3F1(a1: {0}, a2: {1}, z: {2}, addend: {3}): No convergence after {4} iterations", a1, a2, z, addend, n);
            return double.NaN;
        }



        /// <summary>
        /// Returns the sum of a hypergeometric series 3F0. 
        /// <para>Sum3F0 = addend + Σ(( (a1)_n * (a2)_n * (a3)_n) * (z^n/n!)) n={0,∞}</para>
        /// </summary>
        /// <param name="a1">First numerator</param>
        /// <param name="a2">Second numerator</param>
        /// <param name="a3">Third numerator</param>
        /// <param name="z">Argument</param>
        /// <param name="addend">The addend to the series</param>
        /// <returns>
        /// The sum of the series, or NaN if the series does not converge
        /// within the policy maximum number of iterations 
        /// </returns>
        public static double Sum3F0(double a1, double a2, double a3, double z, double addend = 0)
        {
            if ((double.IsNaN(a1) || double.IsInfinity(a1))
            || (double.IsNaN(a2) || double.IsInfinity(a2))
            || (double.IsNaN(a3) || double.IsInfinity(a3))
            || (double.IsNaN(z) || double.IsInfinity(z))
            || double.IsNaN(addend)) {
                Policies.ReportDomainError("Sum3F1(a1: {0}, a2: {1}, a3: {2}, z: {3}, addend: {4}): Requires finite arguments", a1, a2, a3, z, addend);
                return double.NaN;
            }

            double sum = 1.0 + addend;

            if (z == 0)
                return sum;
            if (a1 == 1)
                return Sum3F0_A1(a2, a3, z, addend);
            if (a2 == 1)
                return Sum3F0_A1(a1, a3, z, addend);
            if (a3 == 1)
                return Sum3F0_A1(a1, a2, z, addend);

            double term = 1.0;
            int n = 0; 
            for (; n < Policies.MaxSeriesIterations; n++) {

                double prevSum = sum;
                term *= (z / (n + 1)) * (a1 + n) * (a2 + n) * (a3 + n);
                sum += term;

#if CMP_FLOAT
                if (Math.Abs(term) <= Math.Abs(prevSum) * Policies.SeriesTolerance)
                    return sum;
#else
                if (sum == prevSum)
                    return sum;
#endif
            }

            Policies.ReportConvergenceError("Sum3F1(a1: {0}, a2: {1}, a3: {2}, z: {3}, addend: {4}): No convergence after {5} iterations", a1, a2, a3, z, addend, n);
            return double.NaN;
        }


        /// <summary>
        /// Returns the sum of a hypergeometric series 3F1. 
        /// <para>Sum3F1 = addend + Σ(( (a1)_n * (a2)_n * (a3)_n)/(b1)_n * (z^n/n!)) n={0,∞}</para>
        /// </summary>
        /// <param name="a1">First numerator</param>
        /// <param name="a2">Second numerator</param>
        /// <param name="a3">Third numerator</param>
        /// <param name="b1">First denominator. Requires b1 > 0</param>
        /// <param name="z">Argument</param>
        /// <param name="addend">The addend to the series</param>
        /// <returns>
        /// The sum of the series, or NaN if the series does not converge
        /// within the policy maximum number of iterations 
        /// </returns>
        public static double Sum3F1(double a1, double a2, double a3, double b1, double z, double addend = 0)
        {
            if ((double.IsNaN(a1) || double.IsInfinity(a1))
            || (double.IsNaN(a2) || double.IsInfinity(a2))
            || (double.IsNaN(a3) || double.IsInfinity(a3))
            || (double.IsNaN(z) || double.IsInfinity(z))
            || (double.IsNaN(b1) || double.IsInfinity(b1))
            || double.IsNaN(addend)) {
                Policies.ReportDomainError("Sum3F1(a1: {0}, a2: {1}, a3: {2}, b1: {3}, z: {4}, addend: {5}): Requires finite arguments", a1, a2, a3, b1, z, addend);
                return double.NaN;
            }
            if (b1 <= 0 && Math2.IsInteger(b1)) {
                Policies.ReportDomainError("Sum3F1(a1: {0}, a2: {1}, a3: {2}, b1: {3}, z: {4}, addend: {5}): Requires b1 is not zero or a negative integer", a1, a2, a3, b1, z, addend);
                return double.NaN;
            }

            double sum = 1.0 + addend;

            if (z == 0)
                return sum;
            if (a1 == b1)
                return Sum2F0(a2, a3, z, addend);
            if (a2 == b1)
                return Sum2F0(a1, a3, z, addend);
            if (a3 == b1)
                return Sum2F0(a1, a2, z, addend);
#if false
            // TODO:
            if (a1 == 1)
                return Sum3F1_A1(a2, a3, b1, z, addend);
            if (a2 == 1)
                return Sum3F1_A1(a1, a3, b1, z, addend);
            if (a3 == 1)
                return Sum3F1_A1(a1, a2, b1, z, addend);
#endif

            double term = 1.0;
            int n = 0; 
            for (; n < Policies.MaxSeriesIterations; n++) {

                double prevSum = sum;
                term *= (z / (n + 1)) * ((a1 + n) / (b1 + n)) * (a2 + n) * (a3 + n);
                sum += term;

#if CMP_FLOAT
                if (Math.Abs(term) <= Math.Abs(prevSum) * Policies.SeriesTolerance)
                    return sum;
#else
                if (sum == prevSum)
                    return sum;
#endif
            }

            Policies.ReportConvergenceError("Sum3F1(a1: {0}, a2: {1}, a3: {2}, b1: {3}, z: {4}, addend: {5}): No convergence after {6} iterations", a1, a2, a3, b1, z, addend, n);
            return double.NaN;
        }

        /// <summary>
        /// Returns the sum of a hypergeometric series 4F1. 
        /// <para>Sum4F1 = addend + Σ(( (a1)_n * (a2)_n * (a3)_n * (a4)_n)/(b1)_n * (z^n/n!)) n={0,∞}</para>
        /// </summary>
        /// <param name="a1">First numerator</param>
        /// <param name="a2">Second numerator</param>
        /// <param name="a3">Third numerator</param>
        /// <param name="a4">Fourth numerator</param>
        /// <param name="b1">First denominator. Requires b1 > 0</param>
        /// <param name="z">Argument</param>
        /// <param name="addend">The addend to the series</param>
        /// <returns>
        /// The sum of the series, or NaN if the series does not converge
        /// within the policy maximum number of iterations 
        /// </returns>
        public static double Sum4F1(double a1, double a2, double a3, double a4, double b1, double z, double addend = 0)
        {
            if ((double.IsNaN(a1) || double.IsInfinity(a1))
            || (double.IsNaN(a2) || double.IsInfinity(a2))
            || (double.IsNaN(a3) || double.IsInfinity(a3))
            || (double.IsNaN(a4) || double.IsInfinity(a4))
            || (double.IsNaN(z) || double.IsInfinity(z))
            || (double.IsNaN(b1) || double.IsInfinity(b1))
            || double.IsNaN(addend)) {
                Policies.ReportDomainError("Sum4F1(a1: {0}, a2: {1}, a3: {2}, a4: {3}, b1: {4}, z: {5}, addend: {6}): Requires finite arguments", a1, a2, a3, a4, b1, z, addend);
                return double.NaN;
            }
            if (b1 <= 0 && Math2.IsInteger(b1)) {
                Policies.ReportDomainError("Sum4F1(a1: {0}, a2: {1}, a3: {2}, a4: {3}, b1: {4}, z: {5}, addend: {6}): Requires b1 is not zero or a negative integer", a1, a2, a3, a4, b1, z, addend);
                return double.NaN;
            }

            double sum = 1.0 + addend;

            if (z == 0)
                return sum;

            if (a1 == b1)
                return Sum3F0(a2, a3, a4, z, addend);
            if (a2 == b1)
                return Sum3F0(a1, a3, a4, z, addend);
            if (a3 == b1)
                return Sum3F0(a1, a2, a4, z, addend);
            if (a4 == b1)
                return Sum3F0(a1, a2, a3, z, addend);
#if false
            if (a1 == 1)
                return Sum4F1_A1(a2, a3, a4, b1, z, addend);
            if (a2 == 1)
                return Sum4F1_A1(a1, a3, a4, b1, z, addend);
            if (a3 == 1)
                return Sum4F1_A1(a1, a2, a4, b1, z, addend);
            if (a4 == 1)
                return Sum4F1_A1(a1, a2, a3, b1, z, addend);
#endif

            double term = 1.0;
            int n = 0; 
            for (; n < Policies.MaxSeriesIterations; n++) {

                double prevSum = sum;
                term *= (z / (n + 1)) * ((a1 + n) / (b1 + n)) * (a2 + n) * (a3 + n) * (a4 + n);
                sum += term;

#if CMP_FLOAT
                if (Math.Abs(term) <= Math.Abs(prevSum) * Policies.SeriesTolerance)
                    return sum;
#else
                if (sum == prevSum)
                    return sum;
#endif

            }

            Policies.ReportConvergenceError("Sum4F1(a1: {0}, a2: {1}, a3: {2}, a4: {3}, b1: {4}, z: {5}, addend: {6}): Not convergence after {7} iterations", a1, a2, a3, a4, b1, z, addend, n);
            return double.NaN;
        }


        /// <summary>
        /// Returns the sum of a hypergeometric series PFQ, which has a variable number of numerators and denominators. 
        /// <para> PFQ = addend + Σ(( Prod(a[j]_n)/Prod(b[j])_n * (z^n/n!)) n={0,∞}</para>
        /// </summary>
        /// <param name="a">The array of numerators</param>
        /// <param name="b">The array of denominators. (b != 0 or b != negative integer)</param>
        /// <param name="z">Argument</param>
        /// <param name="addend">The addend to the series</param>
        /// <returns>
        /// The sum of the series, or NaN if the series does not converge
        /// within the policy maximum number of iterations 
        /// </returns>
        public static double SumPFQ(double[] a, double[] b, double z, double addend = 0)
        {
            if ( (a == null || b == null)
                || (double.IsNaN(z) || double.IsInfinity(z))
                || double.IsNaN(addend)) {
                string aStr = (a == null) ? "null" : a.ToString();
                string bStr = (b == null) ? "null" : b.ToString();
                Policies.ReportDomainError("SumPFQ(a: {0}; b: {1}, z: {2}, addend: {3}): Requires finite non-null arguments", aStr, bStr, z, addend);
                return double.NaN;
            }

            if (a.Length == 0 && b.Length == 0)
                return Sum0F0(z, addend);

            // check for infinities or non-positive integer b

            bool hasBadParameter = false;
            for (int i = 0; i < a.Length; i++) {
                if (double.IsInfinity(a[i])) {
                    hasBadParameter = true;
                    break;
                }
            }

            if (!hasBadParameter) {
                for (int i = 0; i < b.Length; i++) {
                    if (double.IsInfinity(b[i]) || (b[i] <= 0 && Math2.IsInteger(b[i]))) {
                        hasBadParameter = true;
                        break;
                    }


                }
            }

            if (hasBadParameter) {
                Policies.ReportDomainError("SumPFQ(a: {0}; b: {1}, z: {2}, addend: {3}): Requires finite arguments and b not zero or a negative integer", a, b, z, addend);
                return double.NaN;
            }

            double sum = 1.0 + addend;

            if (z == 0)
                return sum;

            double term = 1.0;
            int n = 0; 
            for (; n < Policies.MaxSeriesIterations; n++) {


                double term_part = (z / (n + 1));

                // use a/b where we can for max(a.Length, b.Length)

                int j = 0;
                for (; j < a.Length && j < b.Length; j++) {
                    term_part *= (a[j] + n) / (b[j] + n);
                }

                // otherwise multiply or divide the rest of the factors 
                // either a or b

                for (; j < a.Length; j++) {
                    term_part *= (a[j] + n);
                }

                for (; j < b.Length; j++) {
                    term_part /= (b[j] + n);
                }

                term *= term_part;

                double prevSum = sum;
                sum += term;

#if CMP_FLOAT
                if (Math.Abs(term) <= Math.Abs(prevSum) * Policies.SeriesTolerance)
                    return sum;
#else
                if (sum == prevSum)
                    return sum;
#endif

            }

            Policies.ReportConvergenceError("No convergence after {0} iterations", n);
            return double.NaN;

        }

    }


} // namespace