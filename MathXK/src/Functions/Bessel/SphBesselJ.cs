//  Copyright (c) 2013 Rocco Romeo
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

//  Ported to .NET based on the original:
//      Copyright (c) 2007 John Maddock, Boost Software License v1.0

using System;

namespace MathXK
{

    public partial class Math2
    {
        /// <summary>
        /// Returns the Spherical Bessel function of the first kind: j<sub>n</sub>(x)
        /// </summary>
        /// <param name="n">Order. Requires n &gt;= 0</param>
        /// <param name="x">Argument. Requires x &gt;= 0</param>
        /// <returns></returns>
        public static double SphBesselJ(int n, double x)
        {
            if (!(n >= 0) || !(x >= 0)) {
                Policies.ReportDomainError("SphBesselJ(n: {0}, x: {1}): Requires n >= 0, x >= 0", n, x);
                return double.NaN;
            }
            if (double.IsInfinity(x))
                return 0;

            //
            // Special case, n == 0 resolves down to the sinus cardinal of x:
            //
            if (n == 0)
                return Math2.Sinc(x);
            if ( x < 1) {
                // the straightforward evaluation below can underflow too quickly when x is small
                double result = (0.5 * Constants.SqrtPI);
                result *= (Math.Pow(x / 2, n) / Math2.Tgamma(n + 1.5));
                if (result != 0)
                    result *= HypergeometricSeries.Sum0F1(n + 1.5, -x * x / 4);
                return result;
            }

            // Default case is just a straightforward evaluation of the definition:
            return Constants.SqrtHalfPI * BesselJ(n + 0.5, x) / Math.Sqrt(x);
        }
    }

} // namespaces
