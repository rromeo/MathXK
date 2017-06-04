﻿//  Copyright (c) 2013 Rocco Romeo
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

//  Ported to .NET based on the original:
//      Copyright (c) 2006 Xiaogang Zhang, Boost Software License v1.0
//      Copyright (c) 2006 John Maddock, Boost Software License v1.0
//  History:
//      * XZ wrote the original of this file as part of the Google Summer of Code 2006.  
//      * JM modified it to fit into the Boost.Math conceptual framework better, and to ensure
//      that the code continues to work no matter how many digits type double has.

using System;

namespace MathXK
{

    public partial class Math2
    {
        // Elliptic integrals (complete and incomplete) of the second kind
        // Carlson, Numerische Mathematik, vol 33, 1 (1979)


        /// <summary>
        /// Returns the incomplete elliptic integral of the second kind E(φ, k)
        /// <para>E(φ, k) =  ∫ sqrt(1-k^2*sin^2(θ)) dθ, θ={0,φ}</para>
        /// </summary>
        /// <param name="k">The modulus. Requires |k| ≤ 1</param>
        /// <param name="phi">The amplitude</param>
        public static double EllintE(double k, double phi)
        {
            if ((!(k >= -1 && k <= 1)) || (double.IsNaN(phi))) {
                Policies.ReportDomainError("EllintE(k: {0}, phi: {1}): Requires |k| <= 1; phi not NaN", k, phi);
                return double.NaN;
            }
            if (double.IsInfinity(phi))
                return phi;

            if (Math.Abs(phi) > Trig.PiReductionLimit) {
                Policies.ReportNotImplementedError("EllintE(k: {0}, phi: {1}): |phi| > {2} not implemented", k, phi, Trig.PiReductionLimit);
                return double.NaN;
            }

            // special values
            if (k == 0)
                return phi;
            if (phi == 0)
                return 0;
            if (phi == Math.PI / 2)
                return Math2.EllintE(k);

            // Carlson's algorithm works only for |phi| <= π/2,
            // use the integrand's periodicity to normalize phi
            // E(k, phi + π*mult) = E(k, phi) + 2*mult*E(k)

            double result = 0;
            double rphi = Math.Abs(phi);
            if (rphi > Math.PI / 2) {
                // Normalize periodicity so that |rphi| <= π/2
                var (angleMultiple, angleRemainder) = Trig.RangeReducePI(rphi);
                double mult = 2 * angleMultiple;
                rphi = angleRemainder;
                if (mult != 0)
                    result += mult * EllintE(k);
            }

            if (k == 1) {

                result += Math.Sin(rphi);

            } else {

                double k2 = k * k;
                double sinp = Math.Sin(rphi);
                double cosp = Math.Cos(rphi);
                double x = cosp * cosp;
                double t = k2 * sinp * sinp;
                double y = (t < 0.875) ? 1 - t : (1 - k2) + k2 * x;
                double z = 1;

                result += sinp * (Math2.EllintRF(x, y, z) - t * Math2.EllintRD(x, y, z) / 3);
            }

            return (phi < 0) ? -result : result;
        }

        /// <summary>
        /// Returns the complete elliptic integral of the second kind E(k) = E(π/2,k) 
        /// <para>E(k) = ∫ sqrt(1-k^2*sin^2(θ)) dθ, θ={0,π/2}</para>
        /// </summary>
        /// <param name="k">The modulus. Requires |k| ≤ 1</param>
        public static double EllintE(double k)
        {
            if (!(k >= -1 && k <= 1)) {
                Policies.ReportDomainError("EllintE(k: {0}): Requires |k| <= 1", k);
                return double.NaN;
            }

            // Small series at k == 0
            // Note that z = k^2 in the following link
            // http://functions.wolfram.com/EllipticIntegrals/EllipticE/06/01/03/01/0004/

            double k2 = k * k;
            if (k2 < DoubleLimits.RootMachineEpsilon._13) {
                if (k2 < 4*DoubleLimits.MachineEpsilon)
                    return Math.PI / 2; // E(0) = Pi/2
                return (Math.PI / 2) * HypergeometricSeries.Sum2F1(-0.5, 0.5, 1, k2);
            }

            if (Math.Abs(k) == 1)
                return 1;

            double x = 0;
            double y = 1 - k2;
            double z = 1;

            double value = 2 * Math2.EllintRG(x, y, z);
            return value;
        }
    }

} // namespaces



