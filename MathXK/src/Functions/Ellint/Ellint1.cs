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
//      * RR added small series, better performing agm, and better angle reduction 

#if DEBUG
//#define EXTRA_DEBUG
#endif

using System;
using System.Diagnostics;

namespace MathXK
{

    public partial class Math2
    {

        // Elliptic integrals (complete and incomplete) of the first kind
        // Carlson, Numerische Mathematik, vol 33, 1 (1979)


        /// <summary>
        /// Returns the incomplete elliptic integral of the first kind 
        /// <para>F(φ, k) = F(φ | k^2) = F(sin(φ); k)</para>
        /// <para>F(φ, k) = ∫ dθ/sqrt(1-k^2*sin^2(θ)), θ={0,φ}</para>
        /// <para>F(x; k) = ∫ dt/sqrt((1-t^2)*(1-k^2*t^2)), t={0,x}</para>
        /// </summary>
        /// <param name="k">The modulus. Requires |k| ≤ 1</param>
        /// <param name="phi">The amplitude</param>
        public static double EllintF(double k, double phi)
        {
            if ((!(k >= -1 && k <= 1)) || (double.IsNaN(phi))) {
                Policies.ReportDomainError("EllintF(k: {0}, phi: {1}): Requires |k| <= 1; phi not NaN", k, phi);
                return double.NaN;
            }
            if (double.IsInfinity(phi))
                return phi;

            if (Math.Abs(phi) > Trig.PiReductionLimit) {
                Policies.ReportNotImplementedError("EllintF(k: {0}, phi: {1}): |phi| > {2} not implemented", k, phi, Trig.PiReductionLimit);
                return double.NaN;
            }

#if EXTRA_DEBUG 
            Debug.WriteLine("EllintF(phi = {0}; k = {1})", phi, k);
#endif

            // special values
            if (k == 0)
                return phi;
            if (phi == 0)
                return 0;
            if (phi == Math.PI / 2)
                return Math2.EllintK(k);


            // Carlson's algorithm works only for |phi| <= π/2,
            // use the integrand's periodicity to normalize phi
            // F(k, phi + π*mult) = F(k, phi) + 2*mult*K(k)

            double result = 0;
            double rphi = Math.Abs(phi);
            if (rphi > Math.PI / 2) {
                // Normalize periodicity so that |rphi| <= π/2
                var (angleMultiple, angleRemainder) = Trig.RangeReducePI(rphi);
                double mult = 2 * angleMultiple;
                rphi = angleRemainder;
                if (mult != 0)
                    result += mult * EllintK(k);
            }


            double sinp = Math.Sin(rphi);
            double cosp = Math.Cos(rphi);

            double k2 = k * k;
            double x = cosp * cosp;
            double y = 1 - k2 * sinp * sinp;
            if ( y > 0.125 ) // use a more accurate form with less cancellation
                y = (1 - k2) + k2 * x;
            double z = 1;

            result += sinp * EllintRF(x, y, z);

            return (phi < 0) ? -result : result;

        }

        /// <summary>
        /// Returns the complete elliptic integral of the first kind K(k) = F(π/2,k) 
        /// <para>K(k) = ∫ dθ/sqrt(1-k^2*sin^2(θ)), θ={0,π/2}</para>
        /// <para>K(k) = ∫ dt/sqrt((1-t^2)*(1-k^2*t^2)), t={0,1}</para>
        /// </summary>
        /// <param name="k">The modulus. Requires |k| ≤ 1</param>
        public static double EllintK(double k)
        {

            if (!(k >= -1 && k <= 1)) {
                Policies.ReportDomainError("EllintK(k: {0}): Requires |k| <= 1", k);
                return double.NaN;
            }

            // See: http://functions.wolfram.com/EllipticIntegrals/EllipticK/06/01/03/01/0004/
            // K(k) = 2F1(1/2, 1/2; 1; k^2)
            // Be careful of notational differences
            double k2 = k * k;
            if ( k2 < DoubleLimits.RootMachineEpsilon._13) {
                if (k2 < 4 *DoubleLimits.MachineEpsilon )
                    return Math.PI / 2;
                return (Math.PI / 2) * HypergeometricSeries.Sum2F1(0.5, 0.5, 1, k2);
            }

            if (Math.Abs(k) == 1)
                return double.PositiveInfinity;

#if false
            // Save this for future reference
            double x = 0;
            double y = (1 - k)*(1 + k); // 1-k^2
            double z = 1;
            double value = EllintRF(x, y, z);
            return value;
#else
            return (Math.PI / 2) / Agm(1 - k, 1 + k);
#endif
        }
    }

} // namespaces


