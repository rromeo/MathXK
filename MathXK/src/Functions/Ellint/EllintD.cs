//  Copyright (c) 2013 Rocco Romeo
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
        /// Returns the incomplete elliptic integral of the second kind D(φ, k)
        /// <para>D(φ, k) =  ∫ sin^2(θ)/sqrt(1-k^2*sin^2(θ)) dθ, θ={0,φ}</para>
        /// </summary>
        /// <param name="k">The modulus. Requires |k| ≤ 1</param>
        /// <param name="phi">The amplitude</param>
        public static double EllintD(double k, double phi)
        {
            if ((!(k >= -1 && k <= 1)) || (double.IsNaN(phi))) {
                Policies.ReportDomainError("EllintD(k: {0}, phi: {1}): Requires |k| <= 1; phi not NaN", k, phi);
                return double.NaN;
            }
            if (double.IsInfinity(phi))
                return phi; // Note: result != phi, only +/- infinity

            if (Math.Abs(phi) > Trig.PiReductionLimit) {
                Policies.ReportNotImplementedError("EllintD(k: {0}, phi: {1}): |phi| > {2} not implemented", k, phi, Trig.PiReductionLimit);
                return double.NaN;
            }

            // special values
            //if (k == 0) // too much cancellation error near phi=0, general case works fine
            //    return (phi - Cos(phi)*Sin(phi))/2
            if (phi == 0)
                return 0;
            if (phi == Math.PI / 2)
                return Math2.EllintD(k);

            // Carlson's algorithm works only for |phi| <= π/2,
            // use the integrand's periodicity to normalize phi
            // D(k, phi + π*mult) = D(k, phi) + 2*mult*D(k)

            double result = 0;
            double rphi = Math.Abs(phi);
            if (rphi > Math.PI / 2) {
                // Normalize periodicity so that |rphi| <= π/2
                var (angleMultiple, angleRemainder) = Trig.RangeReducePI(rphi);
                double mult = 2 * angleMultiple;
                rphi = angleRemainder;
                if (mult != 0)
                    result += mult * EllintD(k);
            }

            double k2 = k * k;
            double sinp = Math.Sin(rphi);
            double cosp = Math.Cos(rphi);
            double x = cosp * cosp;
            double t = k2 * sinp * sinp;
            double y = (t < 0.875) ? 1 - t : (1 - k2) + k2 * x;
            double z = 1;

            // http://dlmf.nist.gov/19.25#E13
            // and RD(lambda*x, lambda*y, lambda*z) = lambda^(-3/2) * RD(x, y, z) 
            result += EllintRD(x, y, z) * sinp * sinp * sinp / 3;

            return (phi < 0) ? -result : result;
        }

        /// <summary>
        /// Returns the complete elliptic integral of the second kind D(k) = D(π/2,k) 
        /// <para>D(k) = ∫ sin^2(θ)/sqrt(1-k^2*sin^2(θ)) dθ, θ={0,π/2}</para>
        /// </summary>
        /// <param name="k">The modulus. Requires |k| ≤ 1</param>
        public static double EllintD(double k)
        {
            if (!(k >= -1 && k <= 1)) {
                Policies.ReportDomainError("EllintD(k: {0}): Requires |k| <= 1", k);
                return double.NaN;
            }

            // special values
            if (k == 0)
                return Math.PI / 4;

            if (Math.Abs(k) == 1)
                return Double.PositiveInfinity;

            // http://dlmf.nist.gov/19.5#E3
            if (k * k < DoubleLimits.RootMachineEpsilon._13)
                return (Math.PI / 4) * HypergeometricSeries.Sum2F1(1.5, 0.5, 2, k * k);

            double x = 0;
            double y = 1 - k * k;
            double z = 1;

            double value = Math2.EllintRD(x, y, z) / 3;
            return value;
        }
    }

} // namespaces



