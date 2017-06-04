//  (C) Copyright Rocco Romeo 2013.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

using System;

namespace MathXK
{

    internal static class Utility
    {

        /// <summary>
        /// Swaps two variables
        /// </summary>
        /// <typeparam name="T"></typeparam>
        /// <param name="ref1"></param>
        /// <param name="ref2"></param>
        public static void Swap<T>(ref T ref1, ref T ref2)
        {
            T temp = ref1;
            ref1 = ref2;
            ref2 = temp;
        }


        /// <summary>
        /// Returns the maximum value of a, b, c
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="c"></param>
        /// <returns></returns>
        public static double Max(double a, double b, double c)
        {
            if (a > b) {
                return (a > c) ? a : c;
            }
            return (b > c) ? b : c;
        }

        /// <summary>
        /// Returns the maximum value of a, b, c, d
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="c"></param>
        /// <param name="d"></param>
        /// <returns></returns>
        public static double Max(double a, double b, double c, double d)
        {
            return Math.Max(Math.Max(a, b), Math.Max(c, d));
        }

        /// <summary>
        /// Order parameters x, y, z such that x &lt;= y &lt;= z
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="z"></param>
        public static void OrderAscending(ref double x, ref double y, ref double z)
        {
            if (x > y)
                Utility.Swap(ref x, ref y);
            if (y > z)
                Utility.Swap(ref y, ref z);
            if (x > y)
                Utility.Swap(ref x, ref y);
        }

    }


} // namespace