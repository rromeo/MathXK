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
        /// Returns <paramref name="x"/> * 2^n
        /// </summary>
        /// <param name="x">Argument</param>
        /// <param name="n">The binary exponent (2^n)</param>
        /// <returns>
        /// Returns <paramref name="x"/> * 2^n
        /// <para>If <paramref name="x"/> is NaN, returns NaN.</para> 
        /// <para>If <paramref name="x"/> is ±0 or ±∞, returns x.</para> 
        /// <para>If n is 0, returns x.</para>
        /// </returns> 
        public static double Ldexp(double x, int n)
        {
            const double toNormal = (1L << IEEEDouble.MantissaBits);
            const double fromNormal = 1.0/(1L << IEEEDouble.MantissaBits);

            // handle special cases
            if (double.IsNaN(x)) {
                Policies.ReportDomainError("Ldexp(x: {0}, n: {1}): Requires x not NaN", x, n);
                return x;
            }
            if (double.IsInfinity(x))
                return x; // +/-Inf

            if (x == 0 || n == 0)
                return x; // +/-0, x

            IEEEDouble bits = new IEEEDouble(x);
            int xExp = 0;
            if ( bits.IsSubnormal ) {
                bits = new IEEEDouble(x * toNormal);
                xExp -= IEEEDouble.MantissaBits;
            }
            xExp += bits.ExponentValue + n;

            // there will be an overflow
            if (xExp > IEEEDouble.MaxBinaryExponent) 
                return (x > 0) ? double.PositiveInfinity : double.NegativeInfinity;
            
            // cannot represent with a denorm -- underflow
            if (xExp < IEEEDouble.MinBinaryDenormExponent) 
                return (x > 0) ? 0.0 : -0.0;
            
            if (xExp >= IEEEDouble.MinBinaryExponent) {
                // in the normal range - just set the exponent to the sum
                bits.ExponentValue = xExp;
                return bits;
            }

            // We have a denormalized number

#if true
            // Use C99 definition x*2^exp
            xExp += IEEEDouble.MantissaBits;
            Debug.Assert(xExp >= IEEEDouble.MinBinaryExponent);

            bits.ExponentValue = xExp;
            return ((double)bits) * fromNormal;
#else
            // Make it compatible with MS C library - no rounding for denorms
            // kept for future reference

            int shift = IEEEDouble.MinBinaryExponent - xExp;
            Debug.Assert(shift > 0);

            // Recover the hidden mantissa bit and shift
            const Int64 impliedbit = 1L << IEEEDouble.MantissaBits; // 1.MANTISSA = 1 << 52
            Int64 mantissa = (bits.Significand | impliedbit) >> shift;
            return BitConverter.Int64BitsToDouble(bits.Sign | mantissa);
#endif
        }

    }




} // namespace