//  Copyright (c) 2013 Rocco Romeo
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

//  Ported to .NET based on the original:
//      Copyright (c) John Maddock 2005-2006, Boost Software License, Version 1.0

using System;

namespace MathXK
{

    public static partial class Math2
    {

        /// <summary>
        /// Returns Log(1 + x) with improved accuracy for |x| ≤ 0.5
        /// </summary>
        /// <param name="x">The function argument. Requires x ≥ -1</param>
        public static double Log1p(double x)
        {
            if (!(x >= -1)) {
                Policies.ReportDomainError("Log1p(x: {0}): Requires x >= -1", x);
                return double.NaN;
            }
            if (x == -1)
                return double.NegativeInfinity;
            if (double.IsInfinity(x))
                return double.PositiveInfinity;


            double u = 1 + x;
            if (u == 1.0)
                return x;

            return Math.Log(u) * (x / (u - 1.0));
        }
    }


} // namespaces





