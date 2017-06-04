//  Copyright (c) 2013 Rocco Romeo
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

//  Ported to .NET based on the original:
//      Copyright (c) John Maddock 2006, Boost Software License, Version 1.0


using System;

namespace MathXK
{

    public static partial class Math2
    {

        /// <summary>
        /// Returns the value of the Laguerre Polynomial of order <paramref name="n"/> at point <paramref name="x"/>
        /// </summary>
        /// <param name="n">Polynomial degree. Requires <paramref name="n"/> ≥ 0. Limited to <paramref name="n"/> ≤ 127</param>
        /// <param name="x">Polynomial argument</param>
        public static double LaguerreL(int n, double x)
        {
            if (n < 0) {
                Policies.ReportDomainError("LaguerreL(n: {0}, x: {1}): Requires n >= 0", n, x);
                return double.NaN;
            }

            if (n >= 128) {
                Policies.ReportNotImplementedError("LaguerreL(n: {0}, x: {1}): Not implemented for |n| >= 128", n, x);
                return double.NaN;
            }

            if (n == 0)
                return 1;

            if (double.IsNaN(x)) {
                Policies.ReportDomainError("LaguerreL(n: {0}, x: {1}): NaN not allowed", n, x);
                return double.NaN;
            }
            if (double.IsInfinity(x)) 
                return (x > 0 && IsOdd(n)) ? double.NegativeInfinity : double.PositiveInfinity;


            // use recurrence
            // p0 - current
            // p1 - next

            double p0 = 1;
            double p1 = 1 - x;

            for (uint k = 1; k < (uint)n; k++) {
                double next = ((2 * k + 1 - x) * p1 - k * p0) / (k + 1);
                p0 = p1;
                p1 = next;

            }
            return p1;
        }


        /// <summary>
        /// Returns the value of the Associated Laguerre polynomial of degree <paramref name="n"/> and order <paramref name="m"/> at point <paramref name="x"/>
        /// </summary>
        /// <param name="n">Polynomial degree. Requires <paramref name="n"/> ≥ 0. Limited to <paramref name="n"/> ≤ 127</param>
        /// <param name="m">Polynomial order. Requires <paramref name="m"/> ≥ 0</param>
        /// <param name="x">Polynomial argument</param>
        public static double LaguerreL(int n, int m, double x)
        {
            if ((n < 0) || (m < 0)) {
                Policies.ReportDomainError("LaguerreL(n: {0}, m: {1}, x: {2}): Requires n,m >= 0", n, m, x);
                return double.NaN;
            }
            if ( n >= 128) {
                Policies.ReportNotImplementedError("LaguerreL(n: {0}, m: {1}, x: {2}): Not implemented for n >= 128", n, m, x);
                return double.NaN;
            }

            // Special cases:
            if (m == 0)
                return Math2.LaguerreL(n, x);
            if (n == 0)
                return 1;

            if (double.IsNaN(x)) {
                Policies.ReportDomainError("LaguerreL(n: {0}, m: {1}, x: {2}): NaN not allowed", n, m, x);
                return double.NaN;
            }
            if (double.IsInfinity(x))
                return (x > 0 && IsOdd(n)) ? double.NegativeInfinity : double.PositiveInfinity;

            // use recurrence
            // p0 - current
            // p1 - next

            double p0 = 1;
            double p1 = m + 1 - x;

            for (uint k = 1; k < (uint)n; k++) {
                double next = ((2 * k + (uint)m + 1 - x) * p1 - (k + (uint)m) * p0) / (k + 1);
                p0 = p1;
                p1 = next;
            }
            return p1;
        }
    }


} // namespaces





