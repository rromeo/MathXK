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
        /// Returns Sin(π * x)
        /// </summary>
        /// <param name="x">Argument</param>
        public static double SinPI(double x)
        {
            if (double.IsNaN(x) || double.IsInfinity(x)) {
                Policies.ReportDomainError("SinPI(x: {0}): Requires finite x", x);
                return double.NaN;
            }

            // sin(-x) == -sin(x)
            bool neg = false;

            if (x < 0) {
                x = -x;
                neg = !neg;
            }

            // Reduce angles over 2*PI
            if (x >= 2)
                x -= 2 * Math.Floor(x / 2);

            Debug.Assert(x >= 0 && x < 2);

            // reflect in the x-axis
            // and negate so that x in [0,1]*π
            if (x > 1) {
                x = 2.0 - x;
                neg = !neg;
            }

            // reflect in the y-axis so that x in [0,0.5]*π
            if (x > 0.5)
                x = 1.0 - x;

            double result;
            if (x == 0.5) {
                result = 1.0;
            } else if (x > 0.25) {
                result = Math.Cos((0.5 - x) * Math.PI);
            } else {
                result = Math.Sin(Math.PI * x);
            }

            return (neg) ? -result : result;

        }



    }


}