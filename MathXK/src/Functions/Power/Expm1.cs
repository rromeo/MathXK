//  Copyright (c) 2013 Rocco Romeo
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

//  Ported to .NET based on the original:
//      Copyright (c) John Maddock 2006, Boost Software License, Version 1.0

using System;

namespace MathXK
{

    public partial class Math2
    {


        /// <summary>
        /// Returns Exp(x) - 1 with improved accuracy for |x| &lt; 0.5
        /// </summary>
        /// <param name="x">Expm1 function argument</param>
        public static double Expm1(double x)
        {
            if (double.IsNaN(x)) {
                Policies.ReportDomainError("Expm1(x: {0}): Requires x not NaN", x);
                return double.NaN;
            }

            double absX = Math.Abs(x);
            if (absX > 0.5)
                return Math.Exp(x) - 1.0;

            if (absX < DoubleLimits.MachineEpsilon)
                return x;

            const double Y = 0.10281276702880859e1;

            const double p0 = -0.28127670288085937e-1;
            const double p1 = 0.51278186299064534e0;
            const double p2 = -0.6310029069350198e-1;
            const double p3 = 0.11638457975729296e-1;
            const double p4 = -0.52143390687521003e-3;
            const double p5 = 0.21491399776965688e-4;

            const double q0 = 1;
            const double q1 = -0.45442309511354755e0;
            const double q2 = 0.90850389570911714e-1;
            const double q3 = -0.10088963629815502e-1;
            const double q4 = 0.63003407478692265e-3;
            const double q5 = -0.17976570003654402e-4;

            double p = p0 + x * (p1 + x * (p2 + x * (p3 + x * (p4 + x * p5))));
            double q = q0 + x * (q1 + x * (q2 + x * (q3 + x * (q4 + x * q5))));

            return x * (Y + p/q);
        }

    }


} // namespaces






