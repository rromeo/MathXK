//  Copyright (c) 2013 Rocco Romeo
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

//  Ported to .NET based on the original:
//      Copyright (c) John Maddock 2008, Boost Software License, Version 1.0
//      Copyright (c) Hubert Holin 2001

using System;

namespace MathXK
{

    public static partial class Math2
    {

        /// <summary>
        /// Returns the Inverse Hyperbolic Tangent
        /// <para>Atanh(x) = 1/2 * ln((1+x)/(1-x))</para>
        /// </summary>
        /// <param name="x">Atanh argument. Requires |x| &lt;= 1</param>
        public static double Atanh(double x)
        {
            if (!(x >= -1 && x <= 1)) {
                Policies.ReportDomainError("Atanh(x: {0}): Requires |x| <= 1", x);
                return double.NaN;
            }
            if (x == -1) 
                return double.NegativeInfinity;
            if (x == 1) 
                return double.PositiveInfinity;

            // For |x| < eps^(1/4), use a taylor series
            // Atanh(x) = z + z^3/3 + z^5/5 ...
            // http://functions.wolfram.com/ElementaryFunctions/ArcTanh/06/01/03/01/
            //
            // for  eps^(1/4) <= |x| < 0.5
            // Atanh(x) = 1/2 * ( ln(1+x) - ln(1-x) )
            // http://functions.wolfram.com/ElementaryFunctions/ArcTanh/02/
            //
            // otherwise use 
            // Atanh(x) = 1/2 * ln((1+x)/(1-x))

            double absX = Math.Abs(x);
            if (absX < 0.5) {
                if (absX < DoubleLimits.RootMachineEpsilon._4)
                    return x * (1.0 + x * x / 3);
                return (Math2.Log1p(x) - Math2.Log1p(-x)) / 2;
            }

            return (Math.Log((1 + x) / (1 - x)) / 2);

        }
    }



} //namespaces




