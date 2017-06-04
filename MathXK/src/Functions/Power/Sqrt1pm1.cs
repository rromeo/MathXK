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
        /// Returns Sqrt(x+1)-1 with improved accuracy for |x| ≤ 0.75
        /// </summary>
        /// <param name="x">The function argument</param>
        public static double Sqrt1pm1(double x)
        {
            if (!(x >= -1)) {
                Policies.ReportDomainError("Sqrt1pm1(x: {0}): Requires x >= -1", x);
                return double.NaN;
            }
            if (double.IsInfinity(x))
                return x;

            if (Math.Abs(x) <= 0.75)
                return Math2.Expm1(Math2.Log1p(x) / 2);

            return Math.Sqrt(1 + x) - 1;
        }


    }

}