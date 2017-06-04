//  Copyright (c) 2013 Rocco Romeo
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

//  Ported to .NET based on the original:
//      Copyright (c) 2006 Xiaogang Zhang, Boost Software License v1.0
//      Copyright (c) 2006, 2015 John Maddock, Boost Software License v1.0
//  History:
//      * XZ wrote the original of this file as part of the Google Summer of Code 2006.  
//      * JM modified it to fit into the Boost.Math conceptual framework better, and to ensure
//      that the code continues to work no matter how many digits type double has.
//      * JM updated to use Carlson's latest methods (2015)

using System;
using System.Diagnostics;

namespace MathXK
{

    public partial class Math2
    {
        // Carlson's elliptic integral of the second kind
        // Carlson, Numerische Mathematik, vol 33, 1 (1979)


        /// <summary>
        /// Carlson's symmetric form of elliptical integrals 
        /// <para>R<sub>D</sub>(x,y,z) = R<sub>J</sub>(x,y,z,z) = (3/2) * ∫ dt/((t+z)*Sqrt((t+x)(t+y)(t+z))), t={0,∞}</para>
        /// </summary>
        /// <param name="x">Argument. Requires x >= 0</param>
        /// <param name="y">Argument. Requires y >= 0</param>
        /// <param name="z">Argument. Requires z > 0</param>
        public static double EllintRD(double x, double y, double z)
        {
            if ((!(x >= 0) || double.IsInfinity(x))
            || (!(y >= 0) || double.IsInfinity(y))
            || (!(z > 0) || double.IsInfinity(z))
            || (x + y == 0)) {
                Policies.ReportDomainError("EllintRD(x: {0}, y: {1}, z: {2}) requires: finite x,y >= 0; finite z > 0; at most one argument == 0", x, y, z);
                return double.NaN;
            }

            // Special cases from http://dlmf.nist.gov/19.20#iv

            if ( x > y )
                Utility.Swap(ref x, ref y);
            if (x == 0) {
                if (y == 0)
                    return double.PositiveInfinity; // RD(0, 0, z)
                if (y == z)
                    return 3 * (Math.PI / 4) / (y * Math.Sqrt(y)); // RD(0, y, y)
                //
                // Special handling for common case, from
                // Numerical Computation of Real or Complex Elliptic Integrals, eq.47
                //
                double xn = Math.Sqrt(y);
                double yn = Math.Sqrt(z);
                double x0 = xn;
                double y0 = yn;
                double sum = 0;
                double sum_pow = 0.25;

                double c = xn - yn;
                while (Math.Abs(c) >= 2.7 * DoubleLimits.RootMachineEpsilon._2 * xn) {
                    double t = Math.Sqrt(xn * yn);
                    xn = (xn + yn) / 2;
                    yn = t;

                    c = xn - yn;
                    sum_pow *= 2;
                    sum += sum_pow * c * c;
                }
                double RF = Math.PI / (xn + yn);
                //
                // This following calculation suffers from serious cancellation when y ~ z
                // unless we combine terms.  We have:
                //
                // ( ((x0 + y0)/2)^2 - z ) / (z(y-z))
                //
                // Substituting y = x0^2 and z = y0^2 and simplifying we get the following:
                //
                double pt = (x0 + 3 * y0) / (4 * z * (x0 + y0));
                //
                // Since we've moved the demoninator from eq.47 inside the expression, we
                // need to also scale "sum" by the same value:
                //
                pt -= sum / (z * (y - z));
                return pt * RF * 3;
            }

            if ( x == y && x == z )
                return 1 / (x * Math.Sqrt(x));


            // Taylor series upper limit
            const double tolerance = 0.002049052858245406612358440409548690638254492642397575328974;
            Debug.Assert(Math2.AreNearUlps(Math.Pow(DoubleLimits.MachineEpsilon / 3, 1.0 / 6.0), tolerance, 5), "Incorrect constant");

            // duplication
            double sigma = 0;
            double factor = 1;
            int k = 1;
            for (; k < Policies.MaxSeriesIterations; k++) {

                double u = (x + y + z + z + z) / 5;
                double X = (u - x);
                double Y = (u - y);
                double Z = (u - z);

                // Termination condition: 
                double utol = u * tolerance;
                if (Math.Abs(X) < utol && Math.Abs(Y) < utol && Math.Abs(Z) < utol) {

                    X /= u;
                    Y /= u;
                    Z /= u;

                    // Taylor series expansion to the 5th order
                    double EA = X * Y;
                    double EB = Z * Z;
                    double EC = EA - EB;
                    double ED = EA - 6 * EB;
                    double EE = ED + EC + EC;
                    double S1 = ED * (ED * (9.0 / 88) - Z * EE * (9.0 / 52) - 3.0 / 14);
                    double S2 = Z * (EE / 6 + Z * (-EC * (9.0 / 22) + Z * EA * (3.0 / 26)));
                    double value = 3 * sigma + factor * (1 + S1 + S2) / (u * Math.Sqrt(u));

                    return value;
                }

                double sx = Math.Sqrt(x);
                double sy = Math.Sqrt(y);
                double sz = Math.Sqrt(z);
                double lambda = sy * (sx + sz) + sz * sx;
                sigma += factor / (sz * (z + lambda));
                factor /= 4;
                x = (x + lambda) / 4;
                y = (y + lambda) / 4;
                z = (z + lambda) / 4;
            }

            Policies.ReportDomainError("EllintRD(x: {0}, y: {1}, z: {2}) No convergence after {3} iterations", x, y, z, k);
            return double.NaN;
        }
    }

} // namespaces



