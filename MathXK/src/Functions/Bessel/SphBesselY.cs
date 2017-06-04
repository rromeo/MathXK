//  Copyright (c) 2013 Rocco Romeo
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

//  Ported to .NET based on the original:
//      Copyright (c) 2007 John Maddock, Boost Software License v1.0

using System;

namespace MathXK {

    public partial class Math2
    {
        /// <summary>
        /// Returns the Spherical Bessel function of the second kind: y<sub>n</sub>(x) 
        /// </summary>
        /// <param name="v">Order. Requires n &gt;= 0</param>
        /// <param name="x">Argument. Requires x &gt;= 0</param>
        /// <returns></returns>
        public static double SphBesselY(int v, double x)
        {
            const double LowerLimit = 2 * DoubleLimits.MinNormalValue;

            if (!(v >= 0) || !(x >= 0)) {
                Policies.ReportDomainError("SphBesselY(v: {0}, x: {1}): Requires v >= 0, x >= 0", v, x);
                return double.NaN;
            }
            if (x < LowerLimit) 
                return double.NegativeInfinity;
            if (double.IsInfinity(x))
                return 0;

            // The code below is the following
            // double result = (Constants.SqrtHalfPI/Math.Sqrt(x)) * Math2.CylBesselY(v+0.5, x);
            // but for x << v, spherical Bessel y can overflow too soon.
            var (J, Y) = _Bessel.JY(v + 0.5, x, false, true);
            double prefix = (Constants.SqrtHalfPI / Math.Sqrt(x));
            double result = (double)(prefix * Y);

            return result;
        }
    }

}