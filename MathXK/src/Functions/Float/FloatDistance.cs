//  Copyright (c) 2013 Rocco Romeo
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

using System;
using System.Diagnostics;

namespace MathXK
{

    public static partial class Math2
    {
        // this routine adapted from 
        // David LeBlanc, Integer Handling with the C++ SafeInt Class
        // http://msdn.microsoft.com/en-us/library/ms972705


        /// <summary>
        /// Return the difference as a double
        /// Note: if the difference is large, there will be a loss of precision
        /// </summary>
        /// <param name="lhs">The left hand side argument</param>
        /// <param name="rhs">The right hand side argument</param>
        /// <returns>System.Double.</returns>
        private static double SafeSubtractD(Int64 lhs, Int64 rhs)
        {

            //test for +/- combo
            if ((rhs ^ lhs) < 0) {

                //mixed positive and negative
                //two cases - +X - -Y => X + Y - check for overflow against MaxInt
                //            -X - +Y - check for overflow against MinInt


                if (lhs >= 0) {

                    // first case: X >= 0, Y < 0
                    // test is     X - Y > MaxInt, 
                    // equivalent to X > MaxInt + Y

                    // Y == MinInt creates special case because -MinInt overflows
                    // However, the equation still works due to the fact that
                    // MaxInt - MinInt == -1, and lhs is non-negative
                    if (lhs > Int64.MaxValue + rhs) {
                        //remember that rhs is negative
                        // TODO: Use int64 math for better precision
                        return (double)(lhs) - (double)(rhs);
                    }

                } else {

                    // second case: X < 0, Y >= 0
                    // test is X - Y < MinInt
                    // or      X < MinInt + Y

                    // Y== MinInt does not cause any problems here because 
                    // abs(MinInt) > MaxInt
                    if (lhs < Int64.MinValue + rhs) {
                        // TODO: Use int64 math for better precision
                        return (double)(lhs) - (double)(rhs);
                    }
                }
            }

            //both negative, or both positive
            //no possible overflow

            return (double)(lhs - rhs);
        }


        /// <summary>
        /// Returns the number of floating point representations between <paramref name="x"/> and <paramref name="y"/> 
        /// </summary>
        /// <param name="x">The first value</param>
        /// <param name="y">The second value</param>
        /// <returns>
        /// The result is always a signed integer value (stored in floating-point) representing the number of distinct representations between x and y.
        /// </returns>
        public static double FloatDistance(double x, double y)
        {

            if ((double.IsNaN(x) || double.IsInfinity(x))
            || (double.IsNaN(y) || double.IsInfinity(y))) {
                Policies.ReportDomainError("FloatDistance(x: {0}, y: {1}): Requires finite x, y", x, y);
                return double.NaN;
            }

            // Uses Bruce Dawson's technique
            // See http://www.cygnus-software.com/papers/comparingfloats/comparingfloats.htm

            // Make xInt lexicographically ordered as a twos-complement int
            Int64 xInt = BitConverter.DoubleToInt64Bits(x);
            if (xInt < 0)
                xInt = Int64.MinValue - xInt;

            // Make yInt lexicographically ordered as a twos-complement int
            Int64 yInt = BitConverter.DoubleToInt64Bits(y);
            if (yInt < 0)
                yInt = Int64.MinValue - yInt;

            // Now we can compare xInt and yInt to find out how far apart X and Y are.
            // use SafeSubtract, could be an overflow
            double result = SafeSubtractD(xInt, yInt);
            result = Math.Abs(result);

            if (y > x)
                result = -result;
            return result;


        }

    }

}