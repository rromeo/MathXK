//  Copyright (c) 2013 Rocco Romeo
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

//  Ported to .NET based on the original:
//        Copyright (c) 2006 Xiaogang Zhang, Boost Software License v1.0
//        Copyright (c) 2006, 2015 John Maddock, Boost Software License v1.0
//  History:
//      * XZ wrote the original of this file as part of the Google Summer of Code 2006.  
//      * JM modified it to fit into the Boost.Math conceptual framework better, and to correctly
//      handle the y < 0 case.
//      * JM updated to use Carlson's latest methods (2015)

using System;
using System.Diagnostics;

namespace MathXK
{

    public partial class Math2
    {

        // Carlson's degenerate elliptic integral
        // Carlson, Numerische Mathematik, vol 33, 1 (1979)


        /// <summary>
        /// Carlson's symmetric form of elliptical integrals 
        /// <para>R<sub>C</sub>(x,y) = R<sub>F</sub>(x,y,y) = (1/2) * ∫ dt/((t+y)*Sqrt(t+x)), t={0,∞}</para>
        /// </summary>
        public static double EllintRC(double x, double y)
        {
            if ((!(x >= 0) || double.IsInfinity(x)) 
                || (y == 0 || double.IsInfinity(y) || double.IsNaN(y))) {
                Policies.ReportDomainError("EllintRC(x: {0}, y: {1}): Requires finite x >= 0; finite y != 0", x, y);
                return double.NaN;
            }

            // Taylor series: Max Value
            const double tolerance = 0.003100392679625389599124425076703727071077135406054400409781;
            Debug.Assert(Math2.AreNearUlps(tolerance, Math.Pow(4 * DoubleLimits.MachineEpsilon, 1.0 / 6), 5), "Incorrect constant");

            // for y < 0, the integral is singular, return Cauchy principal value
            double prefix = 1;
            if (y < 0) {
                prefix = Math.Sqrt(x / (x - y));
                x -= y;
                y = -y;
            }

            if (x == 0) 
                return prefix * ((Math.PI / 2) / Math.Sqrt(y));
            if (x == y) 
                return prefix / Math.Sqrt(x);

            int k = 1;
            for(; k < Policies.MaxSeriesIterations; k++ ) {
                double u = (x + y + y) / 3;
                double S = y / u - 1; // 1 - x / u = 2 * S

                if (2 * Math.Abs(S) < tolerance) {
                    // Taylor series expansion to the 5th order
                    double value = (1 + S * S * (3.0 / 10 + S * (1.0 / 7 + S * (3.0 / 8 + S * (9.0 / 22))))) / Math.Sqrt(u);
                    return value * prefix;
                }

                double sx = Math.Sqrt(x);
                double sy = Math.Sqrt(y);
                double lambda = 2 * sx * sy + y;
                x = (x + lambda) / 4;
                y = (y + lambda) / 4;

            } 

            Policies.ReportDomainError("EllintRC(x: {0}, y: {1}): No convergence after {2} iterations", x, y, k);
            return double.NaN;
        }
    }

} // namespaces



