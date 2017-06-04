//  Copyright (c) 2013 Rocco Romeo
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

//  Based on the original:
//      Copyright (c) John Maddock 2005-2006, Boost Software License, Version 1.0
// History:
//      Source originally by JM. 
//      RR ported code to .NET and added code to fit the framework better.


using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;

namespace MathXK
{

    /// <summary>
    /// Routines for evaluating Continued Fractions
    /// </summary>
    public static class ContinuedFraction
    {

        private const double _DefaultTolerance = 2 * DoubleLimits.MachineEpsilon;


        //
        // continued_fraction
        // Evaluates:
        //
        //            a0
        // c + ---------------
        //      b0 +     a1
        //           ----------
        //           b1 +   a2
        //                -----
        //                b2 + ...
        //
        //

        // For more information see: http://en.wikipedia.org/wiki/Generalized_continued_fraction



        /// <summary>
        /// Evaluates a finite continued fraction
        /// <para>Return = initialValue + a/(b[0]+(a/b[1]+...(b[n-1] + a/b[n])))</para>
        /// </summary>
        /// <param name="initialValue">The initial value that is added to the finite fraction</param>
        /// <param name="a">The numerator</param>
        /// <param name="b">The array of denominators</param>
        /// <returns></returns>
        public static double Eval(double initialValue, double a, double[] b)
        {
            if (double.IsNaN(initialValue)) {
                Policies.ReportDomainError("Requires finite initialValue: initialValue = {0}", initialValue);
                return double.NaN;
            }
            if (b == null) {
                Policies.ReportDomainError("Requires b != null");
                return double.NaN;
            }

            IEnumerable<(double, double)> series = Enumerable.Repeat(a, b.Length).Zip(b, (x, y) => (x, y));

            return Eval(initialValue, series, _DefaultTolerance, b.Length);
        }

        /// <summary>
        /// Evaluates a finite continued fraction
        /// <para>Return = initialValue + a[0]/(b[0]+(a[1]/b[1]+...(b[n-1]+a[n]/b[n])))</para>
        /// </summary>
        /// <param name="initialValue">The initial value is added to the finite fraction</param>
        /// <param name="a">The array of numerators</param>
        /// <param name="b">The array of denominators</param>
        /// <returns></returns>
        public static double Eval(double initialValue, double[] a, double[] b)
        {
            if (double.IsNaN(initialValue)) {
                Policies.ReportDomainError("Requires finite initialValue: initialValue = {0}", initialValue);
                return double.NaN;
            }
            if (a == null) {
                Policies.ReportDomainError("Requires a != null");
                return double.NaN;
            }
            if (b == null) {
                Policies.ReportDomainError("Requires b != null");
                return double.NaN;
            }
            if (a.Length != b.Length) {
                Policies.ReportDomainError("Requires equal sized arrays: a.Length = {0}; b.Length = {1}", a.Length, b.Length);
                return double.NaN;
            }

            IEnumerable<(double, double)> series = a.Zip(b, (x, y) => (x, y)); 

            return Eval(initialValue, series, _DefaultTolerance, a.Length);
        }

        /// <summary>
        /// Evaluates the continued fraction: initialValue +  K(i=0,inf, a[i]/b[i] )
        /// </summary>
        /// <param name="initialValue">The initial value that is added to the fraction</param>
        /// <param name="series">Each pair of (a[i],b[i])</param>
        /// <param name="tolerance">Stop when changes go below |tolerance|</param>
        /// <param name="maxIterations">The maximum number of iterations before stopping</param>
        /// <returns></returns>
        public static double Eval(double initialValue, IEnumerable<(double, double)> series, double tolerance, int maxIterations)
        {
            if (double.IsNaN(initialValue)) {
                Policies.ReportDomainError("Requires finite initialValue: initialValue = {0}", initialValue);
                return double.NaN;
            }
            if (series == null) {
                Policies.ReportDomainError("Requires series != null");
                return double.NaN;
            }
            if (maxIterations < 0) {
                Policies.ReportDomainError("Requires maxIterations >=0: maxIterations = {0}", maxIterations);
                return double.NaN;
            }


            // Uses the modified Lentz's Algroritm to compute the continued fraction.
            // Developed using psuedocode from "Numerical Recipies: the Art of Scientific Computing", By William H. Press
            // Chapter 5

            if (maxIterations == 0)
                return initialValue;

            // do _not_ use double.Epsilon here -- the routine will fail
            // DoubleLimits.MinNormalizedValue is the minimum safe floating point denominator
            const double tiny = DoubleLimits.MinNormalValue;
            double f = initialValue;
            if (f == 0)
                f = tiny;
            double C = f;
            double D = 0;

            int nIterations = 0;
            foreach (var v in series) {
                double a = v.Item1;
                double b = v.Item2;

                D = b + a * D;
                if (D == 0)
                    D = tiny;

                C = b + a / C;
                if (C == 0)
                    C = tiny;

                D = 1.0 / D;
                double delta = C * D;
                double prevf = f;
                f = f * delta;

                if (Math.Abs(f - prevf) <= Math.Abs(prevf) * tolerance)
                    return f;

                //if ( Math2.AreNearUlps(f,prevf,tolerance) )
                //  return f;

                if (nIterations++ >= maxIterations) {
                    Policies.ReportConvergenceError("No convergence : {0} iterations exceeded", maxIterations);
                    return double.NaN;
                }
            }

            return f;

        }

        /// <summary>
        /// Evaluates the continued fraction: initialValue +  K(i=0,inf, a[i]/b[i] )
        /// </summary>
        /// <param name="initialValue">The initial value that is added to the finite fraction</param>
        /// <param name="series">each pair of (a[i],b[i])</param>
        /// <param name="ulpsTolerance">Stop when changes go below |tolerance|</param>
        /// <returns></returns>
        public static double Eval(double initialValue, IEnumerable<(double, double)> series, uint ulpsTolerance)
        {
            return Eval(initialValue, series, ulpsTolerance, Policies.MaxFractionIterations);
        }

        /// <summary>
        /// Evaluates the continued fraction: initialValue +  K(i=0,inf, a[i]/b[i] )
        /// </summary>
        /// <param name="initialValue">The initial value that is added to the fraction</param>
        /// <param name="series">each pair of (a[i],b[i])</param>
        /// <returns></returns>
        public static double Eval(double initialValue, IEnumerable<(double, double)> series)
        {
            return Eval(initialValue, series, _DefaultTolerance, Policies.MaxFractionIterations);
        }

        /// <summary>
        /// Evaluates the continued fraction: initialValue +  K(i=0,inf, a[i]/b[i] )
        /// </summary>
        /// <param name="initialValue">The initial value that is added to the fraction</param>
        /// <param name="seriesF">The function that returns each (a[i],b[i])</param>
        /// <param name="tolerance">Stop when changes go below |tolerance|</param>
        /// <param name="maxIterations">The maximum number of iterations before stopping</param>
        /// <returns></returns>
        public static double Eval(double initialValue, Func<(double, double)> seriesF, double tolerance, int maxIterations)
        {
            if (double.IsNaN(initialValue)) {
                Policies.ReportDomainError("Requires finite initialValue: initialValue = {0}", initialValue);
                return double.NaN;
            }
            if (seriesF == null) {
                Policies.ReportDomainError("Requires seriesF != null");
                return double.NaN;
            }
            if (maxIterations < 0) {
                Policies.ReportDomainError("Requires maxIterations >=0: maxIterations = {0}", maxIterations);
                return double.NaN;
            }

            // Uses the modified Lentz's Algroritm to compute the continued fraction.
            // Developed using psuedocode from "Numerical Recipies: the Art of Scientific Computing", By William H. Press
            // Chapter 5

            if (maxIterations == 0)
                return initialValue;

            // do _not_ use double.Epsilon here -- the routine will fail
            // DoubleLimits.MinNormalizedValue is the minimum safe floating point denominator
            const double tiny = DoubleLimits.MinNormalValue;
            double f = initialValue;
            if (f == 0)
                f = tiny;
            double C = f;
            double D = 0;

            for (int nIterations = 0; nIterations < maxIterations; nIterations++) {
                (double, double) v = seriesF();
                double a = v.Item1;
                double b = v.Item2;

                D = b + a * D;
                if (D == 0)
                    D = tiny;

                C = b + a / C;
                if (C == 0)
                    C = tiny;

                D = 1.0 / D;
                double delta = C * D;
                double prevf = f;
                f = f * delta;

                if (Math.Abs(f - prevf) <= Math.Abs(prevf) * tolerance)
                    return f;

                //if ( Math2.AreNearUlps(f,prevf,ulpsTolerance) )
                //  return f;

            }

            Policies.ReportConvergenceError("No convergence : {0} iterations exceeded", maxIterations);
            return double.NaN;
        }

        /// <summary>
        /// Evaluates the continued fraction: initialValue +  K(i=0,inf, a[i]/b[i] )
        /// </summary>
        /// <param name="initialValue">The initial value that is added to the fraction</param>
        /// <param name="seriesF">Each pair of (a[i],b[i])</param>
        /// <returns></returns>
        public static double Eval(double initialValue, Func<(double, double)> seriesF)
        {
            return Eval(initialValue, seriesF, _DefaultTolerance, Policies.MaxFractionIterations);
        }

    }


} // namespace




