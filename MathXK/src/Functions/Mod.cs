//  Copyright (c) 2013 Rocco Romeo
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)


using System;

namespace MathXK
{

    public partial class Math2
    {

        /// <summary>
        /// Returns <paramref name="x"/> mod <paramref name="y"/>
        /// <para>Mod(x,y) = Sign(x) * (|x| - Floor(|x| / |y|) * |y|)</para>
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y">Requires y != 0 </param>
        public static double Mod(double x, double y)
        {
            if ((double.IsNaN(x))
            || (double.IsNaN(y) || y == 0)) {
                Policies.ReportDomainError("Mod(x: {0}, y: {1}): Requires finite x, y; y != 0", x, y);
                return double.NaN;
            }

            double absX = Math.Abs(x);
            double absY = Math.Abs(y);

            return Math.Sign(x) * (absX - Math.Floor(absX / absY) * absY);
        }
    }

} //namespace