//  Copyright (c) 2013 Rocco Romeo
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

//  Ported to .NET based on the original:
//      Copyright (c) John Maddock 2005-2006, Boost Software License, Version 1.0

using System;

namespace MathXK
{

    public partial class Math2
    {
        /// <summary>
        /// Returns Sqrt(x*x + y*y) while trying to avoid overflows
        /// </summary>
        /// <param name="x">Hypot function argument 1</param>
        /// <param name="y">Hypot function argument 2</param>
        /// <returns></returns>
        public static double Hypot(double x, double y)
        {
            if (double.IsNaN(x) || double.IsNaN(y)) {
                Policies.ReportDomainError("Hypot(x: {0}, y: {1}) NaN not allowed", x, y);
                return double.NaN;
            }
            if (double.IsInfinity(x) || double.IsInfinity(y)) 
                return double.PositiveInfinity;

            double absX = Math.Abs(x);
            double absY = Math.Abs(y);

            if (absX == 0)
                return absY;

            if (absY == 0)
                return absX;

            if (absX > absY) {
                double h = absY / absX; // h < 1
                return absX * Math.Sqrt(1.0 + h * h);
            } else {
                double h = absX / absY; // h < 1
                return absY * Math.Sqrt(1.0 + h * h);
            }

        }
    }

} //namespace





