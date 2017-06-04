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
        /// Returns Sinh(x)/x
        /// </summary>
        /// <param name="x">Argument</param>
        /// <returns></returns>
        public static double Sinhc(double x)
        {
            if (double.IsNaN(x)) {
                Policies.ReportDomainError("Sinhc(x: {0}): NaNs not allowed", x);
                return double.NaN;
            }
            if (double.IsInfinity(x))
                return double.PositiveInfinity;

            // for |x| < (5040*eps)^(1/6) use a taylor series
            // http://functions.wolfram.com/ElementaryFunctions/Sinh/06/01/02/01/
            // Sinhc(x) = 1 + x^2/3! + x^4/5! + x^6/7!
            //
            // otherwise use
            // Sinhc(x) = Sinh(x)/x

            double absX = Math.Abs(x);
            if (absX < DoubleLimits.RootMachineEpsilon._2)
                return 1;
            if (absX < 4*DoubleLimits.RootMachineEpsilon._6) {
                double x2 = x * x;

                const double c0 = 1.0;
                const double c1 = 1.0 / 6;
                const double c2 = 1.0 / 120;

                return c0 + x2 * (c1 + x2 * c2);
            }

            return (Math.Sinh(x) / x);
        }
    }

} // namespace



