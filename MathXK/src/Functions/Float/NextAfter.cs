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
        /// Returns the next representable value after x in the direction of y. 
        /// </summary>
        /// <param name="x">the value to get the next of</param>
        /// <param name="y">direction</param>
        /// <returns>
        /// The next floating point value in the direction of y
        /// <para>If x == NaN || x == Infinity, returns NaN</para>
        /// <para>If y == NaN, returns NaN</para>
        /// <para>If y &gt; x, returns FloatNext(x)</para>
        /// <para>If y = x, returns x</para>
        /// <para>If y &lt; x, returns FloatPrior(x)</para>
        /// </returns>
        public static double Nextafter(double x, double y)
        {
            if ((double.IsNaN(x) || double.IsInfinity(x))
            || (double.IsNaN(y))) {
                Policies.ReportDomainError("Nextafter(x: {0}, y: {1}): Requires finite x, y not NaN", x, y);
                return double.NaN;
            }

            if (y > x)
                return Math2.FloatNext(x);
            if (y == x)
                return x;
            return Math2.FloatPrior(x);
        }
    }

}