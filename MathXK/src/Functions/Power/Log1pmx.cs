//  Copyright (c) 2013, 2015 Rocco Romeo
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

//  Ported to .NET based on the original:
//      Copyright (c) John Maddock 2005-2006, Boost Software License, Version 1.0

using System;
using System.Diagnostics;

namespace MathXK {

    public static partial class Math2 {


        /// <summary>
        /// Returns Log(1 + x) - x with improved accuracy for |x| ≤ 0.5
        /// </summary>
        /// <param name="x">The function argument. Requires x ≥ -1</param>
        public static double Log1pmx(double x)
        {
            if (!(x >= -1)) {
                Policies.ReportDomainError("Log1pmx(x: {0}): Requires x >= -1", x);
                return double.NaN;
            }
            if (x == -1)
                return double.NegativeInfinity;
            if (double.IsInfinity(x))
                return double.PositiveInfinity;

            // if x is not in [-0.5, 0.625] use the function
            if (x < -0.5 || x > 0.625)
                return Math.Log(1 + x) - x;

            // Log1p works fine
            if (x < -7.0 / 16)
                return Math2.Log1p(x) - x;

            double absX = Math.Abs(x);
            if (absX < (3.0 / 4) * DoubleLimits.MachineEpsilon)
                return -x * x / 2;

            // Based on the techniques used in DiDonato & Morris  
            // Significant Digit Computation of the Incomplete Beta function
            // Function: Rlog1 

            // The following taylor series converges very slowly:
            // log(1+x)-x = -x^2 * (1/2 - x/3 + x^2/4 - x^3/5 + x^4/6 - x^5/7 ... )
            //
            // So, use a combination of preset points and argument reduction:
            // Step #1: Create preset intervals
            // let a = preset constant near x
            // log((1+x)/(1+a)) = log(1+(x-a)/(1+a)) = log(1+u) where u = (x-a)/(1+a); x=a+u+au
            // log(1+x) = log(1+u) + log(1+a)
            // log(1+x)-x = log(1+u)  + log(1+a) - (a+u+a*u)
            // log(1+x)-x = (log(1+u)-u)  + (log(1+a)-a) - a*u
            //
            // Step #2: Argument reduction
            // Using argument reduction: http://dlmf.nist.gov/4.6#E6
            // log(1+u) = 2 * (y + y^3/3 + y^5/5 + y^7/7 ...), where y = u/(2+u)
            // Thus, u maps to u = 2*y/(1-y) 
            //
            // log(1 + u) - u 
            // = 2 * y * (1 + y^2/3 + y^4/5 + y^6/7 ...) - 2y/(1-y)
            // = 2 * y * (1 - 1/(1-y) + y^2/3 + y^4/5 + y^6/7 ...)
            // = 2 * y * y * (-1/(1-y) + y/3 + y^3/5 + y^5/7 ...)
            // = -2 * y * y * (1 + u/2 - y*(1/3 + y^2/5 + y^4/7 ...))
            //
            // Step#3: Together
            // log(1+x)-x = (log(1+a)-a) - a*u - 2 * y * y * (1 + u/2 - y*(1/3 + y^2/5 + y^4/7 ...))
            // where 
            // a = an arbitrary point near x
            // u = (x-a)/(1+a)
            // y = u/(2+u) 



#if !USE_GENERAL
            // Approach for double precision:
            // With some preset points, 
            // Set u, y for a = 0;

            double result = 0;
            double u = x;
            double y = x / (x + 2);

            // Set maximum limit for y = +/-0.125 = 1/8 => x = [-2/9, 2/7] = [ -0.2222..., 0.2857... ]
            const double ymax = 0.125;
            if (y < -ymax) {
                // Try to set the magnitude of "a" as low as possible to avoid cancellation errors
                // but ensure that |(x-a)/(1+a)| < ymax for all x

                // a, Log[1 + a] - a
                const double a = -9.0 / 32;
                const double log1pma = -0.0489916868705768562794077754807;
                const double minX = (a - ymax * (2 + a)) / (1 + ymax); // ~= -0.44

                Debug.Assert(x > minX, "Expecting x > minX, a: " + a + " ymax: " + ymax + " minX: " + minX + " x: " + x);

                // adjust u, y for a;
                u = (x - a) / (1 + a);
                y = u / (u + 2);
                result = log1pma - a * u;

            } else if (y > ymax) {
                // Try to set the magnitude of "a" as low as possible to avoid cancellation errors
                // but ensure that |(x-a)/(1+a)| < ymax for all x

                // a, Log[1 + a] - a
                const double a = 9.0 / 32;
                const double log1pma = -0.0334138360954187432193972342535;
                const double maxX = (a + ymax * (2 + a)) / (1 - ymax); // ~= 0.645

                Debug.Assert(x < maxX, "Expecting x < maxX, a: " + a + " ymax: " + ymax + " maxX: " + maxX + " x: " + x);

                // adjust u, y for a;
                u = (x - a) / (1 + a);
                y = u / (u + 2);
                result = log1pma - a * u;
            }


            Debug.Assert(Math.Abs(y) <= ymax, "Y too large for the series. x: " + x + " y: " + y);


            double y2 = y * y;

            // Set upper limit to y = +/- 0.125
            const double c0 = 1.0 / 3;
            const double c1 = 1.0 / 5;
            const double c2 = 1.0 / 7;
            const double c3 = 1.0 / 9;
            const double c4 = 1.0 / 11;
            const double c5 = 1.0 / 13;
            const double c6 = 1.0 / 15;
            const double c7 = 1.0 / 17;

            double series = (c0 + y2 * (c1 + y2 * (c2 + y2 * (c3 + y2 * (c4 + y2 * (c5 + y2 * (c6 + y2 * c7)))))));
            result += -2 * y2 * (1 + u / 2 - y * series);
            return result;

#else
            // The generalized approach

            double y = x / (x + 2);
            double y2 = y * y;

            double sum = 1.0 / 3;
            double mult = y2;
            double k = 5;
            for (int n = 1; n < Policies.MaxSeriesIterations; n++) {

                double prevSum = sum;
                sum += mult / k;

                if (prevSum == sum)
                    return -2 * y2 * (1 + x / 2 - y * sum);

                mult *= y2;
                k += 2;
            }

            Policies.ReportConvergenceError("Series did not converge after {0} iterations", Policies.MaxSeriesIterations);
            return double.NaN;


#endif

        }
    }


} // namespaces





