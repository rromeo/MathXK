//  Copyright (c) 2015 Rocco Romeo
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

//  Ported to .NEdouble based on the original:
//      Copyright (c) John Maddock 2015, Boost Software License v1.0


using System;
using System.Diagnostics;

namespace MathXK
{

    public partial class Math2
    {
        /// <summary>
        /// Returns the Jacobi Zeta function
        /// Z(k, φ) = E(k, φ) - E(k)/K(k)*F(k, φ)
        /// </summary>
        /// <param name="k">The modulus. Requires |k| ≤ 1</param>
        /// <param name="phi">The amplitude</param>
        /// <returns></returns>
        public static double JacobiZeta(double k, double phi)
        {
            if ((!(k >= -1 && k <= 1)) || double.IsNaN(phi) || double.IsInfinity(phi) ) {
                Policies.ReportDomainError("JacobiZeta(k: {0}, phi: {1}): Requires |k| <= 1; finite phi", k, phi);
                return double.NaN;
            }

            if (phi == 0)
                return 0;

            double sign = 1;
            if(phi < 0) {
               phi = Math.Abs(phi);
               sign = -1;
            }

            double k2 = k * k;
            if (k2 < 8 * DoubleLimits.MachineEpsilon) {
                // See: http://functions.wolfram.com/EllipticIntegrals/JacobiZeta/06/01/08/0001/
                if (k2 == 0)
                    return 0;
                return sign * Math.Sin(2 * phi) * k2/4;
            }

            double sinp = Math.Sin(phi);
            double cosp = Math.Cos(phi);

            if (k == 1) {
                // We get here by simplifying JacobiZeta[w, 1] in Mathematica, and the fact that 0 <= phi.
                return sign * sinp * Math.Sign(cosp);  
            }

            double x = 0;
            double y = 1 - k2;
            double z = 1;

            double p = 1 - k2 * sinp * sinp;
            if ( p < 0.125 ) // second form is more accurate for small p
                p = (1 - k2) + k2 * cosp * cosp;

            double result = sign * k2 * sinp * cosp * Math.Sqrt(p) * EllintRJ(x, y, z, p)/ (3 * EllintK(k));
            return result;
        }

    }


}