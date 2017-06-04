//  Copyright (c) 2013 Rocco Romeo
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

//  Ported to .NET based on the original:
//      Copyright (c) 2007 John Maddock, Boost Software License v1.0
//      Most of the implementation details copyright Xiaogang Zhang.
//

using System;
using System.Diagnostics;
using MathXK.Numerics;

namespace MathXK
{


    static partial class _Bessel
    {

        /// <summary>
        /// Series evaluation for Jv(x), for small x
        /// </summary>
        /// <param name="v"></param>
        /// <param name="x"></param>
        /// <returns></returns>
        public static double J_SmallArg(double v, double x)
        {
            Debug.Assert(v >= 0);

            // See http://functions.wolfram.com/Bessel-TypeFunctions/BesselJ/06/01/04/01/01/0003/
            // Converges rapidly for all x << v.
            // Most of the error occurs calculating the prefix

            double prefix;
            if (v <= DoubleLimits.MaxGamma - 1)
                prefix = Math.Pow(x / 2, v) / Math2.Tgamma(v + 1);
            else
                prefix = StirlingGamma.PowDivGamma(x / 2, v)/v;

            if (prefix == 0)
                return prefix;

            double series = HypergeometricSeries.Sum0F1(v + 1, -x * x / 4);

            return prefix * series;
        }

    }

    public partial class Math2
    {


        /// <summary>
        /// Returns the Cylindrical Bessel J function: J<sub>v</sub>(x)
        /// </summary>
        /// <param name="v">Order</param>
        /// <param name="x">Argument</param>
        /// <returns></returns>
        public static double BesselJ(double v, double x)
        {
            if (double.IsNaN(x) || double.IsNaN(v)) {
                Policies.ReportDomainError("BesselJ(v: {0}, x: {1}): Requires finite arguments", v, x);
                return double.NaN;
            }
            if (x < 0 && !IsInteger(v)) {
                Policies.ReportDomainError("BesselJ(v: {0}, x: {1}): Requires integer v when x < 0", v, x);
                return double.NaN;
            }
            if (x == 0 && v < 0) {
                Policies.ReportDomainError("BesselJ(v: {0}, x: {1}): Complex results not supported. Requires v >= 0 when x = 0", v, x);
                return double.NaN;
            }

            if (double.IsInfinity(x))
                return 0;

            if (x == 0) {
                Debug.Assert(v >= 0);

                if (v == 0)
                    return 1;
                return 0;
            }

            if (IsInteger(v))
                return _Bessel.JN(v, x);

            Debug.Assert(x > 0);

            return _Bessel.JY(v, x, true, false).J;
        }
    };

} //namespaces