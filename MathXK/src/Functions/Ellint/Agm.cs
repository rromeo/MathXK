//  Copyright (c) 2015 Rocco Romeo
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)


using System;

namespace MathXK
{

    public partial class Math2
    {
        /// <summary>
        /// Returns the Arithmetic-Geometric Mean of a and b
        /// </summary>
        /// <param name="a">First Agm argument</param>
        /// <param name="b">Second Agm argument</param>
        /// <returns></returns>
        public static double Agm(double a, double b) 
        {
            if (double.IsNaN(a) || double.IsNaN(b) ) {
                Policies.ReportDomainError("Agm(a: {0}, b: {1}): NaNs not allowed", a, b);
                return double.NaN;
            }

            // Agm(x, 0) = Agm(0, x) = 0
            if (a == 0 || b == 0)
                return 0;

            // We don't handle complex results
            if (Math.Sign(a) != Math.Sign(b)) {
                Policies.ReportDomainError("Agm(a: {0}, b: {1}): Requires a, b have the same sign", a, b);
                return double.NaN;
            }

            // Agm(± a, ± ∞) = ± ∞
            if ( double.IsInfinity(a) )
                return a;
            if ( double.IsInfinity(b) )
                return b;

            if (a == b)
                return a; //or b

            // save a, b for error messages
            double aOriginal = a;
            double bOriginal = b;

            // Agm(-a, -b) = -Agm(a, b)
            double mult = 1;
            if (a < 0) {
                a = -a; 
                b = -b;
                mult = -1;
            } 

            // Set a > b: Agm(a, b) == Agm(b, a)
            if (a < b)
                Utility.Swap(ref a, ref b);


            // If b/a < eps^(1/4), we can transform Agm(a, b) to a * Agm(1, b/a)
            // and use the Taylor series Agm(1, y) at y = 0
            // See: http://functions.wolfram.com/EllipticFunctions/ArithmeticGeometricMean/06/01/01/

            double y = b / a;
            if (y <= DoubleLimits.RootMachineEpsilon._4) {

                // Compute log((b/a)/4) being careful not to underflow on b/a
                double logyDiv4;
                if (y < 4 * DoubleLimits.MinNormalValue) {
                    logyDiv4 = (y < DoubleLimits.MinNormalValue) ? Math.Log(b) - Math.Log(a) : Math.Log(y);
                    logyDiv4 -= 2 * Constants.Ln2;
                } else {
                    logyDiv4 = Math.Log(y/4);
                }

                double r = (Math.PI / 8) * (-4.0 + (1 + 1 / logyDiv4) * y * y) / logyDiv4;
                double result = mult * a * r;
                return result;
            }

            // scale by a to guard against overflows/underflows in a*b
            mult *= a;
            b = y;
            a = 1;

            // Use a series of arithmetic and geometric means.
            int n = 0;
            for (; n < Policies.MaxSeriesIterations; n++) {

#if !USE_GENERAL
                // If |b/a-1| < eps^(1/9), we can transform Agm(a, b) to a*Agm(1, b/a)
                // and use the Taylor series Agm(1, z) at z = 1:
                // Agm(1, z) = 1 + (z-1)/2 - (z-1)^2/16 + ....
                // http://functions.wolfram.com/EllipticFunctions/ArithmeticGeometricMean/06/01/02/0001/

                double z = (b-a);
                if (Math.Abs(z) < a * DoubleLimits.RootMachineEpsilon._9) {
                    const double c0 = 1;
                    const double c1 = 1.0 / 2;
                    const double c2 = -1.0 / 16;
                    const double c3 = 1.0 / 32;
                    const double c4 = -21.0 / 1024;
                    const double c5 = 31.0 / 2048;
                    const double c6 = -195.0 / 16384;
                    const double c7 = 319.0 / 32768;
                    const double c8 = -34325.0 / 4194304;

                    z /= a; // z=(b-a)/a = b/a-1
                    double r = c0 + z * (c1 + z * (c2 + z * (c3 + z * (c4 + z * (c5 + z * (c6 + z * (c7 + z * c8)))))));
                    double result = mult * a * r;
                    return result;
                }

                // take a step

                double nextA = 0.5 * (a + b);
                b = Math.Sqrt(a*b);
                a = nextA;
#else
                // the double specialized series above is a little faster, 
                // so this more general approach is disabled

                // tol = Sqrt(8 * eps)
                const double tol = 2 * Constants.Sqrt2 * DoubleLimits.RootMachineEpsilon._2;
                double nextA = 0.5 * (a + b);
                if (Math.Abs(b - a) < a * tol)
                    return mult * nextA;

                // take a step

                b = Math.Sqrt(a*b);
                a = nextA;

#endif

            }

            Policies.ReportConvergenceError("Agm(a: {0}, b: {1}): Did not converge after {2} iterations", aOriginal, bOriginal, n);
            return double.NaN; 
        }

    }
}