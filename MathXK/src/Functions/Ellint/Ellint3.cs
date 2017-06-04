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
//      * RR added angle reduction and asymptotic routine for large n.

#if DEBUG
//#define EXTRA_DEBUG
#endif

using System;
using System.Diagnostics;

namespace MathXK
{

    internal static class _EllintPi
    {

        // Elliptic integrals (complete and incomplete) of the third kind
        // Carlson, Numerische Mathematik, vol 33, 1 (1979)

        public static double Imp(double k, double n, double nc)
        {
            // Note arg nc = 1-n, possibly without cancellation errors

            Debug.Assert(k >= -1 && k <= 1, "Requires |k| <= 1: k = " + k);

            // special values 
            // See: http://dlmf.nist.gov/19.6.E3

            if (n == 0) {
                if (k == 0)
                    return Math.PI / 2;
                return Math2.EllintK(k);
            }

            // special values
            // See: http://functions.wolfram.com/EllipticIntegrals/EllipticPi/03/01/01/

            if (k == 0)
                return (Math.PI / 2) / Math.Sqrt(nc);

            if (k == 1)
                return double.NegativeInfinity / Math.Sign(n - 1);

            if (n < 0) {
                // when the magnitude of n gets very large, N below becomes ~= 1, so
                // there can be some error, so use the series at n=Infinity
                // http://functions.wolfram.com/EllipticIntegrals/EllipticPi/06/01/05/
                if (n < -1 / DoubleLimits.MachineEpsilon)
                    return (Math.PI / 2) / Math.Sqrt(-n) + (Math2.EllintE(k) - Math2.EllintK(k)) / n;

                // A&S 17.7.17:
                // shift to k^2 < N < 1

                double k2 = k * k;
                double N = (k2 - n) / nc;
                double Nc = (1 - k2) / nc; // Nc = 1-N

                double result = k2 / (k2 - n) * Math2.EllintK(k);
                result -= (n / nc) * ((1 - k2) / (k2 - n)) * _EllintPi.Imp(k, N, Nc);
                return result;

            }

            double x = 0;
            double y = 1 - k * k;
            double z = 1;
            double p = nc;
            
            double value = Math2.EllintK(k) + n * Math2.EllintRJ(x, y, z, p) / 3;

            return value;
        }

        public static double Imp(double k, double n, double nc, double phi)
        {
            // Note arg nc = 1-n, possibly without cancellation errors

            Debug.Assert(k >= -1 && k <= 1, "Requires |k| <= 1: k = " + k);

            // See: http://dlmf.nist.gov/19.6#iv
            if (phi == 0)
                return 0;
            if (phi == Math.PI / 2)
                return Math2.EllintPi(k, n);

            if (n == 0) {
                if (k == 0)
                    return phi;
                return Math2.EllintF(k, phi);
            }


            // Carlson's algorithm works only for |phi| <= π/2,
            // use the integrand's periodicity to normalize phi
            // Π(k, n, phi + π*mult) = Π(k, n, phi) + 2*mult*Π(k,n)


            double result = 0;
            double rphi = Math.Abs(phi);
            if (rphi > Math.PI / 2) {
                // Normalize periodicity so that |rphi| <= π/2
                var (angleMultiple, angleRemainder) = Trig.RangeReducePI(rphi);
                double mult = 2 * angleMultiple;
                rphi = angleRemainder;
                if ((mult > 0) && (nc > 0))
                    result = mult * _EllintPi.Imp(k, n, nc);
            }


            if (k == 0) {
                double ncr;

                // A&S 17.7.20:
                if (n < 1) {
                    if (nc == 1)
                        return phi;
                    ncr = Math.Sqrt(nc);
                    result += Math.Atan(ncr * Math.Tan(rphi)) / ncr;
                } else if (n == 1) {
                    result += Math.Tan(rphi);
                } else {
                    // n > 1:
                    ncr = Math.Sqrt(-nc);
                    result += Math2.Atanh(ncr * Math.Tan(rphi)) / ncr;
                }

                return (phi < 0) ? -result : result;

            }


            double sphi = Math.Sin(Math.Abs(phi));

            if (n > 1 / (sphi * sphi)) {
                // Complex result is a domain error:
                Policies.ReportDomainError("EllintPi(k: {0}, nu: {1}, phi: {2}): Complex results not supported. Requires n > 1 / sin^2(phi)", k, n, phi);
                return double.NaN;
            }


            // Special cases first:

            if (n < 0) {
                //
                // If we don't shift to 0 <= n <= 1 we get
                // cancellation errors later on.  Use
                // A&S 17.7.15/16 to shift to n > 0:
                //
                double k2 = k * k;
                double N = (k2 - n) / nc;
                double Nc = (1 - k2) / nc; // Nc = 1-N = (1-k^2)/nc

                // check to see if a rounding error occurred
                //  k^2 <= N <= 1  
                if (N < k2) {
#if EXTRA_DEBUG
                    Debug.WriteLine("Rounding error: EllintPi(k: {0} , nu: {1} , phi: {2})",k, n, phi);
#endif
                    N = k2;
                    Nc = 1 - N;
                }

                double p2 = Math.Sqrt(-n * N);
                double delta = Math.Sqrt(1 - k2 * sphi * sphi);


                // Reduce A&S 17.7.15/16
                // Mathematica eqns below
                // N is protected in Mathematica, so use V 
                // V = (k2 - n)/(1 - n)
                // Assuming[(k2 >= 0 && k2 <= 1) && n < 0,  FullSimplify[Sqrt[(1 - V)*(1 - k2/V)]/Sqrt[((1 - n)*(1 - k2/n))]]]
                // Result: ((-1 + k2) n)/((-1 + n) (-k2 + n))

                // Assuming[(k2 >= 0 && k2 <= 1) && n < 0,  FullSimplify[k2/(Sqrt[-n*(k2 - n)/(1 - n)]*Sqrt[(1 - n)*(1 - k2/n)])]]
                // Result: k2/(k2 - n)

                // Assuming[(k2 >= 0 && k2 <= 1) && n < 0,  FullSimplify[Sqrt[1/((1 - n)*(1 - k2/n))]]]
                // Result: Sqrt[n/((k2 - n) (-1 + n))]

                double nResult = -(n / nc) * ((1 - k2) / (k2 - n)) * _EllintPi.Imp(k, N, Nc, rphi);
                nResult += k2 / (k2 - n) * Math2.EllintF(k, rphi);
                nResult += Math.Sqrt((n / nc) / (n - k2)) * Math.Atan((p2 / (2 * delta)) * Math.Sin(2 * rphi));

                result += nResult;

                return (phi < 0) ? -result : result;
            }


            double sinp = Math.Sin(rphi);
            double cosp = Math.Cos(rphi);
            double x = cosp * cosp;
            double t = sinp * sinp;
            double y = 1 - k * k * t;
            double z = 1;
            double p = (n * t < 0.5) ? 1 - n * t : x + nc * t;

            result += sinp * (Math2.EllintRF(x, y, z) + n * t * Math2.EllintRJ(x, y, z, p) / 3);

            return (phi < 0) ? -result : result;
        }

    };


    public partial class Math2
    {
        /// <summary>
        /// Returns the incomplete elliptic integral of the third kind Π(n, φ, k)
        /// <para>Π(n, φ, k) = ∫ dθ/((1-n^2*sin^2(θ)) * sqrt(1-k^2*sin^2(θ))), θ={0,φ}</para>
        /// </summary>
        /// <param name="k">The modulus. Requires |k| ≤ 1</param>
        /// <param name="nu">The characteristic. Requires nu &lt; 1/sin^2(φ)</param>
        /// <param name="phi">The amplitude</param>
        /// <remarks>
        /// Perhaps more than any other special functions the elliptic integrals are expressed in a variety of different ways. 
        /// In particular, the final parameter k (the modulus) may be expressed using a modular angle α, or a parameter m. 
        /// These are related by: 
        ///    <para>k = sin α</para> 
        ///    <para>m = k^2 = (sin α)^2</para> 
        /// So that the incomplete integral of the third kind may be expressed as either: 
        /// <para>Π(n, φ, k)</para> 
        /// <para>Π(n, φ \ α)</para> 
        /// <para>Π(n, φ| m)</para> 
        /// To further complicate matters, some texts refer to the complement of the parameter m, or 1 - m, where: 
        /// 1 - m = 1 - k^2 = (cos α)^2 
        /// </remarks>
        public static double EllintPi(double k, double nu, double phi)
        {
            if ((!(k >= -1 && k <= 1)) || (double.IsNaN(nu) || double.IsNaN(phi))) {
                Policies.ReportDomainError("EllintPi(k: {0}, nu: {1}, phi: {2}): Requires |k| <= 1; No NaNs", k, nu, phi);
                return double.NaN;
            }
            if (double.IsInfinity(phi))
                return phi;

            if (Math.Abs(phi) > Trig.PiReductionLimit) {
                Policies.ReportNotImplementedError("EllintPi(k: {0}, nu: {1}, phi: {2}): |phi| > {3} not implemented", k, nu, phi, Trig.PiReductionLimit);
                return double.NaN;
            }

            return _EllintPi.Imp(k, nu, 1 - nu, phi);
        }

        /// <summary>
        /// Returns the complete elliptic integral of the first kind Π(n, k) = Π(n, π/2, k)
        /// <para>Π(n, k) = ∫ dθ/((1-n^2*sin^2(θ)) * sqrt(1-k^2*sin^2(θ))), θ={0,π/2}</para>
        /// </summary>
        /// <param name="k">The modulus. Requires |k| ≤ 1</param>
        /// <param name="nu">The characteristic. Requires nu ≤ 1</param>
        /// <remarks>
        /// Perhaps more than any other special functions the elliptic integrals are expressed in a variety of different ways. 
        /// In particular, the final parameter k (the modulus) may be expressed using a modular angle α, or a parameter m. 
        /// These are related by: 
        ///    <para>k = sin α</para> 
        ///    <para>m = k^2 = (sin α)^2</para> 
        /// So that the complete integral of the third kind may be expressed as either: 
        /// <para>Π(n, k)</para> 
        /// <para>Π(n \ α)</para> 
        /// <para>Π(n | m)</para> 
        /// To further complicate matters, some texts refer to the complement of the parameter m, or 1 - m, where: 
        /// 1 - m = 1 - k^2 = (cos α)^2 
        /// </remarks>
        public static double EllintPi(double k, double nu)
        {
            if ((!(k >= -1 && k <= 1)) || (!(nu <= 1))) {
                Policies.ReportDomainError("EllintPi(k: {0}, nu: {1}): Requires |k| <= 1; nu <= 1", k, nu);
                return double.NaN;
            }

            if (nu == 1)
                return double.PositiveInfinity;
            if (double.IsInfinity(nu))
                return 0;

            return _EllintPi.Imp(k, nu, 1 - nu);
        }
    }
} // namespaces



