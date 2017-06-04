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
        /// Returns the Inverse Hyperbolic Sine
        /// <para>Asinh(x) = ln(x + Sqrt(x^2 + 1))</para>
        /// </summary>
        /// <param name="x">Asinh argument</param>
        public static double Asinh(double x)
        {
            if (double.IsNaN(x)) {
                Policies.ReportDomainError("Asinh(x: {0}): NaNs not allowed", x);
                return double.NaN;
            }
            if (double.IsInfinity(x))
                return x;

            // odd function
            if (x < 0)
                return -Math2.Asinh(-x);

            // For |x| < eps^(1/4), use a taylor series
            // Asinh(x) = z - z^3/6 + 3 * z^5/40 ...
            // http://functions.wolfram.com/ElementaryFunctions/ArcSinh/06/01/03/01/      
            //
            // otherwise use varying forms of
            // Ln(x + Sqrt(x * x + 1))
            // http://functions.wolfram.com/ElementaryFunctions/ArcSinh/02/


            if (x < 0.75) {
                if (x < DoubleLimits.RootMachineEpsilon._4)
                    return x * (1.0 - x * x / 6);

                // rearrangment of log(x + sqrt(x^2 + 1)) to preserve digits:
                return Math2.Log1p(x + Math2.Sqrt1pm1(x * x));
            }

            // for large x, result = Math.Log(2*x), but watch for overflows
            if (x > 1 / DoubleLimits.RootMachineEpsilon._2)
                return Constants.Ln2 + Math.Log(x);


            return (Math.Log(x + Math.Sqrt(x * x + 1)));

        }
    }



} //namespaces

