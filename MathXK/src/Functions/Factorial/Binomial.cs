//  Copyright (c) 2013 Rocco Romeo
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

//  Ported to .NET based on the original:
//      Copyright (C) John Maddock 2006, Boost Software License, Version 1.0

using System;

namespace MathXK
{

    public static partial class Math2
    {

        /// <summary>
        /// Returns the binomial coefficient
        /// <para>BinomialCoefficient(n,k) = n! / (k!(n-k)!)</para>
        /// </summary>
        /// <param name="n">Requires n ≥ 0</param>
        /// <param name="k">Requires k in [0,n]</param>
        public static double BinomialCoefficient(int n, int k)
        {
            if ((n < 0)
            || (k < 0 || k > n)) {
                Policies.ReportDomainError("BinomialCoefficient(n: {0},k: {1}): Requires n >= 0; k in [0,n]", n, k);
                return double.NaN;
            }

            if (k == 0 || k == n)
                return 1;
            if (k == 1 || k == n - 1)
                return n;

            double result;
            if (n < Math2.FactorialTable.Length) {
                // Use fast table lookup:
                result = (FactorialTable[n] / FactorialTable[n - k]) / FactorialTable[k];
            } else {
                // Use the beta function:
                if (k < n - k)
                    result = k * Math2.Beta(k, n - k + 1);
                else
                    result = (n - k) * Math2.Beta(k + 1, n - k);

                if (result == 0) 
                    return double.PositiveInfinity;
                
                result = 1 / result;
            }
            // convert to nearest integer:
            return Math.Ceiling(result - 0.5);
        }
    }


} // namespace 






