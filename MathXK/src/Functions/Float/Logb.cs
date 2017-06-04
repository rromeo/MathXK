//  (C) Copyright Rocco Romeo 2017.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

using System;
using System.Diagnostics;

namespace MathXK {

    public static partial class Math2 {


        /// <summary>
        /// Returns Floor(Log(2,|<paramref name="x"/>|)).
        /// Note: this is the exponent field of a double precision number.
        /// <para>If x is NaN then x is returned</para>
        /// <para>If x == ±∞ then +∞ is returned</para>
        /// <para>If x == ±0 then -∞ is returned</para>
        /// </summary>
        /// <param name="x">Argument</param>
        public static double Logb(double x)
        {
            IEEEDouble rep = new IEEEDouble(x);

            // Check for +/- Inf or NaN
            if (!rep.IsFinite) {
                if ( rep.IsInfinity )
                    return double.PositiveInfinity;
                Policies.ReportDomainError("Logb(x: {0}): NaNs not allowed", x);
                return x;
            }

            if ( x == 0 ) {
                Policies.ReportPoleError("Logb(x: {0}): Logb(0) == -∞", x);
                return double.NegativeInfinity;

            }

            int exponent = 0;
            if (rep.IsSubnormal) {
                // Multiply by 2^53 to normalize the number
                const double normalMult = (1L << IEEEDouble.MantissaBits);
                rep = new IEEEDouble(x * normalMult);
                exponent = -IEEEDouble.MantissaBits;

                Debug.Assert(!rep.IsSubnormal);
            }

            exponent += rep.ExponentValue; 
            return exponent;
        }


    }
} // namespace