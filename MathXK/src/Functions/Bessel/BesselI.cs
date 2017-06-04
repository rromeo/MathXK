//  Copyright (c) 2013 Rocco Romeo
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

//  Ported to .NET based on the original:
//      Copyright (c) 2007 John Maddock, Boost Software License v1.0
//      Most of the implementation details copyright Xiaogang Zhang.
//

using System;

namespace MathXK
{

    internal partial class _Bessel
    {

        /// <summary>
        /// Returns Iv(x) for small <paramref name="x"/>
        /// </summary>
        /// <param name="v"></param>
        /// <param name="x"></param>
        /// <returns></returns>
        public static double I_SmallArg(double v, double x)
        {
            double prefix;
            if (v <= DoubleLimits.MaxGamma - 1) {
                prefix = Math.Pow(x / 2, v) / Math2.Tgamma(v + 1);
            } else {
                prefix = StirlingGamma.PowDivGamma(x / 2, v) / v;
            }
            if (prefix == 0)
                return prefix;

            double series = HypergeometricSeries.Sum0F1(v + 1, x * x / 4);
            return prefix * series;
        }


        /// <summary>
        /// Returns I{0.5}(x) = Sqrt(2/PI)*Sinh(x)/Sqrt(x) 
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        public static double I0_5(double x)
        {
            // common special case, note try to avoid overflow in Math.Exp(x):
            if (x >= DoubleLimits.MaxLogValue) {
                double e = Math.Exp(x / 2);
                return e * (e * (Constants.RecipSqrt2PI / Math.Sqrt(x)));
            }

            // Avoid overflows near 0
            return Constants.RecipSqrtHalfPI * Math.Sqrt(x) * Math2.Sinhc(x);
        }


    }

    public partial class Math2
    {
        /// <summary>
        /// Returns the Cylindrical Bessel I function: I<sub>v</sub>(x)
        /// </summary>
        /// <param name="v">Order</param>
        /// <param name="x">Argument</param>
        /// <returns></returns>
        public static double BesselI(double v, double x)
        {
            //
            // This handles all the bessel I functions, note that we don't optimise
            // for integer v, other than the v = 0 or 1 special cases, as Millers
            // algorithm is at least as inefficient as the general case (the general
            // case has better error handling too).
            //

            if (double.IsNaN(v) || double.IsInfinity(v) || double.IsNaN(x)) {
                Policies.ReportDomainError("BesselI(v: {0}, x: {1}): NaNs not allowed; Requres v != Infinity", v, x);
                return double.NaN;
            }

            if (x < 0) {
                // better have integer v:
                if (!IsInteger(v)) {
                    Policies.ReportDomainError("BesselI(v: {0}, x: {1}): Requires integer v when x < 0", v, x);
                    return double.NaN;
                }

                double r = Math2.BesselI(v, -x);
                if (IsOdd(v))
                    r = -r;
                return r;
            }

            if (x == 0)
                return (v == 0) ? 1 : 0;

            // at this point x > 0
            if (double.IsInfinity(x)) 
                return double.PositiveInfinity;

            if (v == 0)
                return _Bessel.I0(x);

            if (v == 1)
                return _Bessel.I1(x);

            if (v == 0.5)
                return _Bessel.I0_5(x);

            if (v > 0 && (x < 2 || 3*(v + 1) > x * x))
                return _Bessel.I_SmallArg(v, x);

            double Iv = _Bessel.IK(v, x, true, false).I;

            return Iv;
        }
    }



} // namespaces