//  Copyright (c) 2013 Rocco Romeo
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

//  Ported to .NET based on the original:
//      Copyright (c) 2006 Xiaogang Zhang, Boost Software License v1.0
//      Copyright (c) 2006 John Maddock, Boost Software License v1.0
//  History:
//      * XZ wrote the original of this file as part of the Google Summer of Code 2006.  
//      * JM modified it to fit into the Boost.Math conceptual framework better, and to correctly
//      handle the p < 0 case.

using System;
using System.Diagnostics;

namespace MathXK
{

    public static partial class Math2
    {


        /// <summary>
        /// Carlson's symmetric form of elliptical integrals 
        /// <para>R<sub>J</sub>(x,y,z,p) = (3/2) * ∫ dt/((t+p)*Sqrt((t+x)(t+y)(t+z))), t={0,∞}</para>
        /// </summary>
        /// <param name="x">Argument. Requires finite x >= 0</param>
        /// <param name="y">Argument. Requires finite y >= 0</param>
        /// <param name="z">Argument. Requires finite z >= 0</param>
        /// <param name="p">Argument. Requires finite p != 0</param>
        public static double EllintRJ(double x, double y, double z, double p)
        {
            if ((!(x >= 0) || double.IsInfinity(x))
                || (!(y >= 0) || double.IsInfinity(y))
                || (!(z >= 0) || double.IsInfinity(z))
                || (p == 0 || double.IsInfinity(p) || double.IsNaN(p))
                || (x + y == 0 || y + z == 0 || z + x == 0)) {
                Policies.ReportDomainError("EllintRJ(x: {0}, y: {1}, z: {2}, p: {3}) requires: finite x,y,z >= 0; finite p != 0; at most one of x,y,z == 0", x, y, z, p);
                return double.NaN;
            }

            // See Carlson, Numerische Mathematik, vol 33, 1 (1979)

            double value;

            // error scales as the 6th power of tolerance
            const double tolerance = 0.002049052858245406612358440409548690638254492642397575328974;
            Debug.Assert(Math2.AreNearUlps(Math.Pow(DoubleLimits.MachineEpsilon / 3, 1.0 / 6.0), tolerance, 5), "Incorrect constant");

            // reorder x,y,z such that x <= y <= z 
            // so that 0 shows up in x first
            // and so that the p<0 case below doesn't suffer cancellation errors
            if (x > y)
                Utility.Swap(ref x, ref y);
            if (y > z)
                Utility.Swap(ref y, ref z);
            if (x > y)
                Utility.Swap(ref x, ref y);

            // for p < 0, the integral is singular, return Cauchy principal value
            if (p < 0) {
                // We must ensure that (z - y) * (y - x) is positive.
                Debug.Assert(x <= y && y <= z, "x, y, z must be in ascending order");

                double q = -p;
                double pmy = (z - y) * (y - x) / (y + q); // p - y

                Debug.Assert(pmy >= 0);

                double pAdj = pmy + y;
                value = pmy * Math2.EllintRJ(x, y, z, pAdj);
                value -= 3 * Math2.EllintRF(x, y, z);
                value += 3 * Math.Sqrt((x * y * z) / (x * z + pAdj * q)) * Math2.EllintRC(x * z + pAdj * q, pAdj * q);
                value /= (y + q);
                return value;
            }

            if (x == 0) {
                if (y == 0) //RJ(0, 0, z, p)
                    return double.PositiveInfinity;
                if (y == z) //RJ(0, y, y, p)
                    return 3 * (Math.PI / 2) / (y * Math.Sqrt(p) + p * Math.Sqrt(y));
            }

            if (p == z) {
                if (z == y) { 
                    if ( x == y )
                        return 1 / (x * Math.Sqrt(x)); // RJ(x, x, x, x)
                    return EllintRD(x, y, y); // RJ(x, y, y, y)
                }
                return EllintRD(x, y, z); // RJ(x, y, z, z)
            }

            if (x == y && y == z ) 
                return EllintRD(p, p, x); // RJ(x, x, x, p)

            // duplication
            double sigma = 0;
            double factor = 1;
            int k = 1;
            for (; k < Policies.MaxSeriesIterations; k++) {
                double u = (x + y + z + p + p) / 5;
                double X = (u - x);
                double Y = (u - y);
                double Z = (u - z);
                double P = (u - p);

                // Termination condition: 
                double utol = u * tolerance;
                if (Math.Abs(X) < utol && Math.Abs(Y) < utol && Math.Abs(Z) < utol && Math.Abs(P) < utol) {

                    X /= u;
                    Y /= u;
                    Z /= u;
                    P /= u;

                    // Taylor series expansion to the 5th order
                    double EA = X * Y + Y * Z + Z * X;
                    double EB = X * Y * Z;
                    double EC = P * P;
                    double E2 = EA - 3 * EC;
                    double E3 = EB + 2 * P * (EA - EC);
                    double S1 = 1 + E2 * (E2 * (9.0 / 88) - E3 * (9.0 / 52) - 3.0 / 14);
                    double S2 = EB * (1.0 / 6 + P * (-6.0 / 22 + P * (3.0 / 26)));
                    double S3 = P * ((EA - EC) / 3 - P * EA * (3.0 / 22));
                    value = 3 * sigma + factor * (S1 + S2 + S3) / (u * Math.Sqrt(u));

                    return value;
                }

                double sx = Math.Sqrt(x);
                double sy = Math.Sqrt(y);
                double sz = Math.Sqrt(z);

                double lambda = sy * (sx + sz) + sz * sx;
                double alpha = p * (sx + sy + sz) + sx * sy * sz;
                alpha *= alpha;
                double beta = p * (p + lambda) * (p + lambda);
                sigma += factor * Math2.EllintRC(alpha, beta);
                factor /= 4;
                x = (x + lambda) / 4;
                y = (y + lambda) / 4;
                z = (z + lambda) / 4;
                p = (p + lambda) / 4;
            }

            Policies.ReportDomainError("EllintRJ(x: {0}, y: {1}, z: {2}, p: {3}) No convergence after {4} iterations", x, y, z, p, k);
            return double.NaN;
        }
    }

} // namespaces


