//  (C) Copyright Rocco Romeo 2013.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)


using System;

namespace MathXK
{
    /// <summary>
    /// Base class for special functions
    /// </summary>
    public static partial class Math2
    {

        // ============================================
        // IsOdd
        // ============================================

        /// <summary>Returns true if <paramref name="x"/> is odd</summary>
        /// <param name="x">The argument</param>
        internal static bool IsOdd(char x)
        {
            return (x & 1) != 0;
        }

        /// <summary>Returns true if <paramref name="x"/> is odd</summary>
        /// <param name="x">The argument</param>
        internal static bool IsOdd(byte x)
        {
            return (x & 1) != 0;
        }

        /// <summary>Returns true if <paramref name="x"/> is odd</summary>
        /// <param name="x">The argument</param>
        internal static bool IsOdd(short x)
        {
            return (x & 1) != 0;
        }

        /// <summary>Returns true if <paramref name="x"/> is odd</summary>
        /// <param name="x">The argument</param>
        internal static bool IsOdd(int x)
        {
            return (x & 1) != 0;
        }

        /// <summary>Returns true if <paramref name="x"/> is odd</summary>
        /// <param name="x">The argument</param>
        internal static bool IsOdd(long x)
        {
            return (x & 1) != 0;
        }

        /// <summary>Returns true if <paramref name="x"/> is odd</summary>
        /// <param name="x">The argument</param>
        internal static bool IsOdd(ushort x)
        {
            return (x & 1) != 0;
        }

        /// <summary>Returns true if <paramref name="x"/> is odd</summary>
        /// <param name="x">The argument</param>
        internal static bool IsOdd(uint x)
        {
            return (x & 1) != 0;
        }

        /// <summary>Returns true if <paramref name="x"/> is odd</summary>
        /// <param name="x">The argument</param>
        internal static bool IsOdd(ulong x)
        {
            return (x & 1) != 0;
        }

        /// <summary>Returns true if <paramref name="x"/> is odd</summary>
        /// <param name="x">The argument</param>
        internal static bool IsOdd(float x)
        {
            return (x - 2 * Math.Floor(x / 2)) != 0;
        }

        /// <summary>Returns true if <paramref name="x"/> is odd</summary>
        /// <param name="x">The argument</param>
        internal static bool IsOdd(double x)
        {
            return (x - 2 * Math.Floor(x / 2)) != 0;
        }

        // ============================================
        // IsInteger
        // ============================================

        /// <summary>Returns true if <paramref name="x"/> is an integer</summary>
        /// <param name="x">The argument</param>
        internal static bool IsInteger(float x)
        {
            return (x == Math.Floor(x));
        }

        /// <summary>Returns true if <paramref name="x"/> is an integer</summary>
        /// <param name="x">The argument</param>
        internal static bool IsInteger(double x)
        {
            return (x == Math.Floor(x));
        }


        // ============================================
        // Powers
        // ============================================

        /// <summary>Returns x^2</summary>
        /// <param name="x">The argument</param>
        internal static double Squared(double x) { return x * x; }

        /// <summary>Returns x^3</summary>
        /// <param name="x">The argument</param>
        internal static double Cubed(double x) { return x * x * x; }


    }


}
