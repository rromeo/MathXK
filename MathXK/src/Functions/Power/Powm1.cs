//  Copyright (c) 2013 Rocco Romeo
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

//  Ported to .NET based on the original:
//      Copyright (c) John Maddock 2006, Boost Software License, Version 1.0

using System;
using System.Diagnostics;

namespace MathXK
{

    public static partial class Math2
    {

        /// <summary>
        /// Returns (x^y)-1 with improved accuracy when x is close to 1, or when y is very small
        /// </summary>
        /// <param name="x">Powm1 function value</param>
        /// <param name="y">Powm1 function exponent</param>
        /// <returns></returns>
        /// <remarks>
        /// Note: The .NET treatment of Math.Pow for NaN parameters is different from that of ISO C.
        /// Math.Pow: for any parameter that is NaN, return NaN. 
        /// ISO C: 1^y == 1, x^0 = 1, even if x,y are NaN
        /// see: http://msdn.microsoft.com/en-us/library/system.math.pow.aspx
        /// see also: http://pubs.opengroup.org/onlinepubs/009695399/functions/pow.html
        /// for the conditions on Math.Pow.
        /// This implementation treats NaN similarly to .NET - treat every NaN as an error condition. 
        /// </remarks>
        public static double Powm1(double x, double y)
        {
            if (double.IsNaN(x) || double.IsNaN(y)) {
                Policies.ReportDomainError("Powm1(x: {0}, y: {1}): Requires x, y not NaN", x, y);
                return double.NaN;
            }

            if (y == 0 || x == 1.0)
                return 0; // x^0 = 1, 1^y = 1
            if (x == 0) {
                return (y > 0) ? -1 : double.PositiveInfinity;
            }

            double xOrig = x;
            double yOrig = y;

            if (x < 0) {
                if (!IsInteger(y)) {
                    Policies.ReportDomainError("Powm1(x: {0}, y: {1}): Requires integer y when x < 0", x, y);
                    return double.NaN;
                }

                if (IsOdd(y))
                    return Math.Pow(x, y) - 1;

                x = -x;
            }

            Debug.Assert(x > 0, "Expecting x > 0");

            // log(x) = log((z+1)/(z-1)) = 2z * (1+z^2/3 .. ), where z = (x-1)/(x+1)
            // see https://en.wikipedia.org/wiki/Natural_logarithm
            // Series converges fastest when z is closest to 1

            if ((x - 1) * y <= 0.25 * (x + 1)) {
                double p = Math.Log(x) * y;
                if (p <= 1)
                    return Math2.Expm1(p);
                // otherwise fall though:
            }
            return Math.Pow(x, y) - 1;
        }
    }

}
