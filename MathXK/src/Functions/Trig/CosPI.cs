//  Copyright (c) 2013 Rocco Romeo
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

//  Ported to .NET based on the original:
//      Copyright (c) 2007 John Maddock, Boost Software License, Version 1.0

using System;
using System.Diagnostics;

namespace MathXK
{

    public static partial class Math2
    {

        /// <summary>
        /// Returns Cos(π * x)
        /// </summary>
        /// <param name="x">Argument</param>
        public static double CosPI(double x)
        {
            if (double.IsNaN(x) || double.IsInfinity(x)) {
                Policies.ReportDomainError("CosPI(x: {0}): Requires finite x", x);
                return double.NaN;
            }

            // cos(-x) == cos(x)
            if (x < 0)
                x = -x;

            // Reduce angles over 2*π
            if (x >= 2)
                x -= 2 * Math.Floor(0.5 * x);

            Debug.Assert(x >= 0 && x < 2);


            // reflect so that x in [0,1]*π
            // cos(x) = cos(2π-x)
            if (x >= 1.0)
                x = 2.0 - x;

            // -cos(x) = cos(π-x)
            bool neg = false;
            if (x > 0.5) {
                x = 1 - x;
                neg = true;
            }

            if (x == 0.5)
                return 0;

            // using the identity: cos(x) = sin(π/2-x)
            // when x > 0.25 significantly improves the precision of this routine in .NET

            double result;
            if (x > 0.25)
                result = Math.Sin((0.5 - x) * Math.PI);
            else
                result = Math.Cos(Math.PI * x);

            return (neg) ? -result : result;
        }



    }


}
