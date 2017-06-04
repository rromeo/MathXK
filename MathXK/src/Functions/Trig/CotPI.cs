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


        /// <summary>
        /// Returns Cot(π * x)
        /// </summary>
        /// <param name="x">Argument</param>
        public static double CotPI(double x)
        {
            if (double.IsNaN(x) || double.IsInfinity(x)) {
                Policies.ReportDomainError("CotPI(x: {0}): Requires finite x", x);
                return double.NaN;
            }

            double absX = Math.Abs(x);

            // the period for cot(x) is π
            if (absX >= 1)
                absX -= Math.Floor(absX);

            if (absX == 0) {
                Policies.ReportPoleError("CotPI(x: {0}): Requires non-integer value for x", x);
                return double.NaN;
            }

            Debug.Assert(absX > 0 && absX < 1);

            // cot(x) == tan(π/2 - x)
            double result;
            if (absX >= 0.5)
                result = -Math.Tan((absX - 0.5) * Math.PI);
            else
                result = 1.0 / Math.Tan(absX * Math.PI);

            // cot(-x) == -cot(x)
            return (x < 0) ? -result : result;
        }




    }



}
