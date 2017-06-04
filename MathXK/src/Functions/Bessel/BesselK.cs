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
        /// Returns the Cylindrical Bessel K function: K<sub>v</sub>(x)
        /// </summary>
        /// <param name="v">Order</param>
        /// <param name="x">Argument</param>
        /// <returns></returns>
        public static double BesselK(double v, double x)
        {
            if (!(x >= 0)) {
                Policies.ReportDomainError("BesselK(v: {0}, x: {1}): Requires x >= 0", v, x);
                return double.NaN;
            }
            if (x == 0) {
                if (v == 0) 
                    return double.PositiveInfinity;

                Policies.ReportDomainError("BesselK(v: {0}, x: {1}): Requires x > 0 when v != 0", v, x);
                return double.NaN;
            }


            // K{-v} = K{v}
            // negative orders are handled in each of the following routines
            if (IsInteger(v) && (v > int.MinValue && v <= int.MaxValue) )
                return _Bessel.KN((int)v, x);

            double kv = _Bessel.IK(v, x, false, true).K;
           
            return kv;
        }
    }

}