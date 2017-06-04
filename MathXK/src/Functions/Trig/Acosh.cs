//  Copyright (c) 2013 Rocco Romeo
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

//  Ported to .NET based on the original:
//      Copyright (c) John Maddock 2008, Boost Software License, Version 1.0
//      Copyright (c) Eric Ford & Hubert Holin 2001

using System;

namespace MathXK
{

    public static partial class Math2
    {

        /// <summary>
        /// Returns the Inverse Hyperbolic Cosine
        /// <para>Acosh(x) = ln(x + Sqrt(x^2 - 1))</para>
        /// </summary>
        /// <param name="x">Acosh argument. Requires x ≥ 1</param>
        public static double Acosh(double x)
        {
            if (!(x >= 1)) {
                Policies.ReportDomainError("Acosh(x:{0}): Requires x >= 1", x);
                return double.NaN;
            }
            if (double.IsInfinity(x))
                return double.PositiveInfinity;

            // for 1 < x < 1 + sqrt(eps), use a taylor series
            // Acosh(x) = Sqrt(2 * y) * (1 - y / 12 + 3 * y^2 / 160 ...) where y = x-1
            // http://functions.wolfram.com/ElementaryFunctions/ArcCosh/06/01/04/01/
            // 
            // otherwise use variations of the standard form
            // Acosh(x) = ln(x + sqrt(x^2 - 1))
            // http://functions.wolfram.com/ElementaryFunctions/ArcCosh/02/

            if (x < 1.5) {
                double y = x - 1;

                if (y < DoubleLimits.RootMachineEpsilon._2) {
                    const double c0 = 1.0;
                    const double c1 = -1 / 12.0;
                    const double c2 = 3 / 160.0;
                     
                    return Math.Sqrt(2 * y) * (c0 + y*(c1 + y*c2));
                }

                // This is just a rearrangement of the standard form 
                // devised to minimize loss of precision when x ~ 1:
                return Math2.Log1p(y + Math.Sqrt(y * ( y + 2 )));
            }

            // for large x, result = Math.Log(2*x), but watch for overflows
            if (x > 1 / DoubleLimits.RootMachineEpsilon._2)
                return Constants.Ln2 + Math.Log(x);

            return (Math.Log(x + Math.Sqrt(x * x - 1)));

        }
    }

}



