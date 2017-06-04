//  (C) Copyright Rocco Romeo 2013.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

using System;
using System.Diagnostics;

namespace MathXK
{

    public static partial class Math2
    {


        /// <summary>
        /// Returns the fraction where <paramref name="x" /> == fraction * 2 ^ exponent.
        /// <para>Note: |fraction| is in [0.5,1.0) or is 0</para>
        /// <para>if x is NaN or ±∞ then x is returned and the exponent is unspecified</para>
        /// </summary>
        /// <param name="x">Argument</param>
        /// <param name="exponent">Output the exponent value</param>
        /// <returns>System.Double.</returns>
        public static double Frexp(double x, out int exponent)
        {
            exponent = 0;

            // Check for +/- Inf or NaN
            if (double.IsNaN(x)) {
                Policies.ReportDomainError("Frexp(x: {0},...): NaNs not allowed", x);
                return x;
            }

            if (double.IsInfinity(x) || x == 0 ) 
                return x; // +/-0, +/-Inf


            IEEEDouble bits = new IEEEDouble(x);
            if (bits.IsSubnormal) {
                // Multiply by 2^52 to normalize the number
                const double normalMult = (1L << IEEEDouble.MantissaBits);
                bits = new IEEEDouble(x * normalMult);
                exponent = -IEEEDouble.MantissaBits;

                Debug.Assert(!bits.IsSubnormal);
            }

            // x is normal from here on

            // must return x in [0.5,1.0) = [2^-1,2^0) 
            // so: x_exponent must = -1 while keeping the mantissa and sign the same
            exponent += bits.ExponentValue + 1; // 2^0 = 1 = 0.5 * 2^1; so increment exponent

            // replace the exponent with 2^-1 so that x is in [0.5,1.0)
            const Int64 ExpHalf = 0x3fe0000000000000;
            return BitConverter.Int64BitsToDouble( (bits & ~IEEEDouble.ExponentMask)| ExpHalf );
        }


    }
} // namespace