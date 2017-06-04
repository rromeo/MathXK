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
        // Carlson's elliptic integral of the first kind
        // Carlson, Numerische Mathematik, vol 33, 1 (1979)

        /// <summary>
        /// Carlson's symmetric form of elliptical integrals 
        /// <para>R<sub>F</sub>(x,y,z) = (1/2) * ∫ dt/((t+x)(t+y)(t+z)), t={0,∞}</para>
        /// </summary>
        /// <param name="x">Argument. Requires x >= 0</param>
        /// <param name="y">Argument. Requires y >= 0</param>
        /// <param name="z">Argument. Requires z >= 0</param>
        public static double EllintRF(double x, double y, double z)
        {
            if ((!(x >= 0) || double.IsInfinity(x))
            || (!(y >= 0) || double.IsInfinity(y))
            || (!(z >= 0) || double.IsInfinity(z))
            || (x + y == 0 || y + z == 0 || z + x == 0)) {
                Policies.ReportDomainError("EllintRF(x: {0}, y: {1}, z: {2}) requires: finite x, y, z >= 0; at most one of x,y,z == 0", x, y, z);
                return double.NaN;
            }

            double xOrig = x;
            double yOrig = y;
            double zOrig = z;

            // reorder x,y,z such that x <= y <= z 
            if (x > y)
                Utility.Swap(ref x, ref y);
            if (y > z)
                Utility.Swap(ref y, ref z);
            if (x > y)
                Utility.Swap(ref x, ref y);

            // Special cases from http://dlmf.nist.gov/19.20#i
            if (x == 0) {
                //Debug.Assert(y != 0 || z != 0, "Expecting only one zero");
                if (y == z)
                    return (Math.PI / 2) / Math.Sqrt(y);

                // http://dlmf.nist.gov/19.22#E8
                return (Math.PI / 2) / Agm(Math.Sqrt(y), Math.Sqrt(z));
            }

            if (x == y) {
                if (x == z)
                    return 1 / Math.Sqrt(x); // EllintRF(x, x, x);
                return EllintRC(z, x);  // EllintRF(x, x, z);
            }

            Debug.Assert(x != z, "Parameters are not ordered properly");

            if (y == z)
                return EllintRC(x, y); // EllintRF(x, y, y)

            // Carlson scales error as the 6th power of tolerance,
            // but this seems not to work for types larger than
            // 80-bit reals, this heuristic seems to work OK:
            const double tolerance = 0.003100392679625389599124425076703727071077135406054400409781;
            Debug.Assert(Math2.AreNearUlps(tolerance, Math.Pow(4 * DoubleLimits.MachineEpsilon, 1.0 / 6), 5), "Incorrect constant");

            // duplication
            int k = 1;
            for (; k < Policies.MaxSeriesIterations; k++) {
                double u = (x + y + z) / 3;
                double X = (u - x);
                double Y = (u - y);
                double Z = (u - z);

                // Termination condition: 
                double utol = u*tolerance;
                if ( Math.Abs(X) < utol && Math.Abs(Y) < utol && Math.Abs(Z) < utol) {

                    // Taylor series expansion to the 5th order
                    X /= u;
                    Y /= u;
                    Z /= u;

                    double E2 = X * Y - Z * Z;
                    double E3 = X * Y * Z;
                    double value = (1 + E2 * (E2 / 24 - E3 * (3.0 / 44) - 0.1) + E3 / 14) / Math.Sqrt(u);

                    return value;
                }

                double sx = Math.Sqrt(x);
                double sy = Math.Sqrt(y);
                double sz = Math.Sqrt(z);
                double lambda = sy * (sx + sz) + sz * sx;
                x = (x + lambda) / 4;
                y = (y + lambda) / 4;
                z = (z + lambda) / 4;
            }


            Policies.ReportDomainError("EllintRF(x: {0}, y: {1}, z: {2}) No convergence after {3} iterations", xOrig, yOrig, zOrig, k);
            return double.NaN;

        }

    }


} // namespaces


