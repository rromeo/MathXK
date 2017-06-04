//  Copyright (c) 2013 Rocco Romeo
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

//  Ported to .NET based on the original:
//      Copyright (c) John Maddock 2006, Boost Software License, Version 1.0


namespace MathXK
{

    public static partial class Math2
    {

        /// <summary>
        /// Returns the value of the Hermite Polynomial of order <paramref name="n"/> at point <paramref name="x"/>
        /// </summary>
        /// <param name="n">The polynomial order</param>
        /// <param name="x">The polynomial argument</param>
        public static double HermiteH(int n, double x)
        {
            if (n < 0) {
                Policies.ReportDomainError("Hermite(n: {0}, x: {1}): Requires n >= 0", n, x);
                return double.NaN;
            }

            // H_{0}(x) = 1, 
            if (n == 0)
                return 1;

            if (double.IsNaN(x)) {
                Policies.ReportDomainError("Hermite(n: {0}, x: {1}): NaN not allowed", n, x);
                return double.NaN;
            }

            if (double.IsInfinity(x)) 
                return (x < 0 && IsOdd(n)) ? double.NegativeInfinity : double.PositiveInfinity;
            
            // use recurrence
            // p0 - current
            // p1 - next

            double p0 = 1;
            double p1 = 2 * x;

            for (int k = 1; k < n; k++) {
                double next = 2 * (x * p1 - k * p0);
                p0 = p1;
                p1 = next;
            }
            return p1;
        }
    }

} // namespaces





