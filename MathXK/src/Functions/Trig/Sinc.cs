//  Copyright (c) 2013 Rocco Romeo
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

//  Ported to .NET based on the original:
//      Copyright (c) Hubert Holin 2001, Boost Software License, Version 1.0

using System;

namespace MathXK
{

    public static partial class Math2
    {
        /// <summary>
        /// Returns Sin(x)/x
        /// </summary>
        /// <param name="x">The Sinc argument</param>
        /// <returns></returns>
        public static double Sinc(double x)
        {
            if (double.IsNaN(x)) {
                Policies.ReportDomainError("Sinc(x: {0}): Requires x not NaN", x);
                return double.NaN;
            }
            if (double.IsInfinity(x))
                return 0;

            // for |x| < eps^(1/6) use a taylor series
            // http://functions.wolfram.com/ElementaryFunctions/Sin/06/01/02/01/
            // Sinc(x) = 1 - x^2/3! + x^4/5! - x^6/7!
            //
            // otherwise use
            // Sinc(x) = Sin(x)/x

            double absX = Math.Abs(x);
            if (absX < DoubleLimits.RootMachineEpsilon._2)
                return 1;

            if (absX < DoubleLimits.RootMachineEpsilon._6) {
                double x2 = x * x;

                const double c0 = 1.0;
                const double c1 = -1.0 / 6;
                const double c2 = 1.0 / 120;

                return c0 + x2 * (c1 + x2 * c2);
            }

            return (Math.Sin(x) / x);
        }

    }

} //namespace



