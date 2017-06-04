//  Copyright (c) 2013 Rocco Romeo
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)


using System;
using System.Diagnostics;

namespace MathXK
{

    public static partial class Math2
    {

        /// <summary>
        /// Returns Tan(π * x)
        /// </summary>
        /// <param name="x">Argument</param>
        public static double TanPI(double x)
        {
            if (double.IsNaN(x) || double.IsInfinity(x)) {
                Policies.ReportDomainError("TanPI(x: {0}): Requires finite x", x);
                return double.NaN;
            }

            // save x for error messages
            double originalX = x;

            // tan(-x) == -tan(x)
            bool neg = false;
            if (x < 0) {
                x = -x;
                neg = !neg;
            }

            // the period for tan(x) is π
            if (x >= 1)
                x -= Math.Floor(x);

            Debug.Assert(x >= 0 && x < 1);

            // shift x to [-0.5,0.5]*π
            if (x > 0.5)
                x -= 1.0;


            if (x == 0.5) {
                Policies.ReportPoleError("TanPI(x: {0}): Requires x not a half integer", originalX);
                return double.NaN;
            }

            double result = Math.Tan(Math.PI * x);

            return (neg) ? -result : result;
        }




    }



}
