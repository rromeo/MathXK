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
        /// Returns a floating point number that is the given <paramref name="distance"/> from <paramref name="x"/>. 
        /// </summary>
        /// <param name="x">The value</param>
        /// <param name="distance">The number of distinct values to increment/decrement x</param>
        public static double FloatAdvance(double x, int distance)
        {
            if (double.IsNaN(x) || double.IsInfinity(x)) {
                Policies.ReportDomainError("FloatAdvance(x: {0}, distance: {1}): Requires finite x", x, distance);
                return double.NaN;
            }


            if (distance == 0)
                return x;

            Int64 bits = BitConverter.DoubleToInt64Bits(x);
            Int64 expon = bits & IEEEDouble.ExponentMask;
            Int64 significand = bits & IEEEDouble.SignificandMask;


            if (expon == 0) {
                // x is either 0.0 or subnormal
                // there could be wraparound here


                // treat -0.0 as 0.0
                // Note: that Int64.MinValue has the same bit pattern as -0.0.

                if (bits == Int64.MinValue)
                    bits = 0;

                if (bits >= 0) {
                    // wraparound?
                    if (distance >= 0 || bits >= -distance)
                        bits += distance;
                    else
                        bits = Int64.MinValue - (bits + distance);


                } else {

                    // wraparound?
                    if (distance < 0 || significand >= distance)
                        bits -= distance;
                    else
                        bits = distance - significand;

                }


                // check for 0.0 or -0.0
                if (bits == Int64.MinValue)
                    return 0.0;


                return BitConverter.Int64BitsToDouble(bits);

            }

            if (bits >= 0)
                bits += distance;
            else
                bits -= distance;

            // if the exponent was 
            //      Max Exponent and the distance was positive
            // or
            //      Min Exponent and the distance was negative
            // there could be an overflow to bitwise NaN -- return infinity instead
            if ((bits & IEEEDouble.InfinityExponentBits) == IEEEDouble.InfinityExponentBits) // infinity or NaN
                bits &= ~IEEEDouble.SignificandMask; // set the significand to 0 

            return BitConverter.Int64BitsToDouble(bits);


        }

    }


} // namespaces



