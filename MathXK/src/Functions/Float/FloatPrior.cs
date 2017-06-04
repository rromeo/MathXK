//  Copyright (c) 2013 Rocco Romeo
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

//  Ported to .NET based on the original:
//      Copyright (C) John Maddock 2008, Boost Software License, Version 1.0

using System;
using System.Diagnostics;

namespace MathXK
{

    public static partial class Math2
    {

        /// <summary>
        /// Returns the previous representable value which is less than x.
        /// </summary>
        /// <param name="x">FloatPrior argument</param>
        /// <returns>The previous floating point value
        /// <para>If x == NaN || x == Infinity, returns NaN</para>
        /// <para>If x == MinValue, returns NegativeInfinity</para>
        /// </returns>
        public static double FloatPrior(double x)
        {

            if (double.IsNaN(x) || double.IsInfinity(x)) {
                Policies.ReportDomainError("FloatPrior(x: {0}): Requires finite x", x);
                return double.NaN;
            }

            if (x <= double.MinValue)
                return double.NegativeInfinity;
            if (x == 0)
                return -double.Epsilon;


            Int64 bits = BitConverter.DoubleToInt64Bits(x);

            // check if we should increase or decrease the mantissa
            if (x < 0)
                bits++;
            else
                bits--;

            return BitConverter.Int64BitsToDouble(bits);

        }
    }

}