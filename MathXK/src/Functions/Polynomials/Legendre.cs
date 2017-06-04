//  Copyright (c) 2013 Rocco Romeo
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

//  Ported to .NET based on the original:
//      Copyright (c) John Maddock 2006, Boost Software License, Version 1.0

using System;
using System.Diagnostics;

namespace MathXK
{

    public static partial class Math2
    {

        /// <summary>
        /// Returns the value of the Legendre Polynomial of the first kind
        /// </summary>
        /// <param name="n">Polynomial degree. Limited to |n| ≤ 127</param>
        /// <param name="x">Requires |x| ≤ 1</param>
        public static double LegendreP(int n, double x)
        {
            if (!(x >= -1 && x <= 1)) {
                Policies.ReportDomainError("LegendreP(n: {0}, x: {1}): Requires |x| <= 1", n, x);
                return double.NaN;
            }
            if (n <= -128 || n >= 128) {
                Policies.ReportNotImplementedError("LegendreP(n: {0}, x: {1}): Not implemented for |n| >= 128", n, x);
                return double.NaN;
            }

            // LegendreP(n,x) == LegendreP(-n-1,x)
            if (n < 0)
                n = -(n + 1);


            // Special values
            if (n == 0)
                return 1;
            if (x == 1)
                return 1;
            if (x == -1)
                return (IsOdd(n)) ? -1 : 1;

            // Implement Legendre P polynomials via recurrance
            // p0 - current
            // p1 - next

            double p0 = 1;
            double p1 = x;

            for (uint k = 1; k < (uint)n; k++) {
                double next = ((2 * k + 1) * x * p1 - k * p0) / (k + 1);
                p0 = p1;
                p1 = next;
            }
            return p1;
        }


        /// <summary>
        /// Returns the value of the Legendre polynomial that is the second solution to the Legendre differential equation
        /// </summary>
        /// <param name="n">Polynomial degree. Requires n ≥ 0. Limited to n ≤ 127</param>
        /// <param name="x">Requires |x| ≤ 1</param>
        public static double LegendreQ(int n, double x)
        {
            if ((n < 0)
            || (!(x >= -1 && x <= 1))) {
                Policies.ReportDomainError("LegendreQ(n: {0}, x: {1}): Requires n >= 0; |x| <= 1", n, x);
                return double.NaN;
            }

            if (x == 1) 
                return double.PositiveInfinity;

            if (x == -1)
                return IsOdd(n) ? double.PositiveInfinity : double.NegativeInfinity;

            if (n >= 128) {
                Policies.ReportNotImplementedError("LegendreQ(n: {0}, x: {1}): Not implemented for n >= 128", n, x);
                return double.NaN;
            }

            // Implement Legendre Q polynomials via recurrance
            // p0 - current
            // p1 - next

            double p0 = Math2.Atanh(x);
            double p1 = x * p0 - 1;
            if (n == 0)
                return p0;

            for (uint k = 1; k < (uint)n; k++) {
                double next = ((2 * k + 1) * x * p1 - k * p0) / (k + 1);
                p0 = p1;
                p1 = next;
            }
            return p1;
        }


        /// <summary>
        /// Returns the Associated Legendre Polynomial of the first kind
        /// </summary>
        /// <param name="n">Polynomial degree. Limited to |<paramref name="n"/>| ≤ 127</param>
        /// <param name="m">Polynomial order</param>
        /// <param name="x">Requires |<paramref name="x"/>| ≤ 1</param>
        public static double LegendreP(int n, int m, double x)
        {
            if (!(x >= -1 && x <= 1)) {
                Policies.ReportDomainError("LegendreP(n: {0}, m: {1}, x: {2}): Requires |x| <= 1", n, m, x);
                return double.NaN;
            }
            if (n <= -128 || n >= 128) {
                Policies.ReportNotImplementedError("LegendreP(n: {0}, m: {1}, x: {2}): Not implemented for |n| >= 128", n, m, x);
                return double.NaN;
            }

            // P{l,m}(x) == P{-l-1,m}(x); watch overflow
            if (n < 0)
                n = -(n + 1);

            if (m == 0)
                return Math2.LegendreP(n, x);

            // by definition P{n,m}(x) == 0 when |m| > l
            if (m < -n || m > n)
                return 0;

            // P{n,-m} = (-1)^m * (n-m)!/(n+m)! * P{n,m}
            if (m < 0) {
                double legP = LegendreP(n, -m, x);
                if (legP == 0)
                    return legP;
                if (IsOdd(m))
                    legP = -legP;

                int MaxFactorial = Math2.FactorialTable.Length - 1;
                if (n - m <= MaxFactorial) {
                    double ratio = Math2.FactorialTable[n + m] / Math2.FactorialTable[n - m];
                    return ratio * legP;
                }

                if (n - m <= 300) {
                    // break up the Tgamma ratio into two separate factors
                    // to prevent an underflow: a!/b! = (a!/170!) * (170!/b!)
                    // 170!/300! ~= 2.37e-308 > double.MinNormalizedValue

                    return (legP * (Math2.FactorialTable[n + m] / Math2.FactorialTable[MaxFactorial])) * Math2.TgammaRatio(MaxFactorial + 1, n + 1 - m);

                }

#if false
                // since we are setting Nmax = 128, so Mmax = 128, and Nmax+Mmax = 256
                // so this should never be called, but left here for completeness
                // 170!/300! ~= 2.37e-308 > double.MinNormalizedValue

                double ratioL = Math2.TgammaRatio(n+m+1, n+1-m);
                if ( ratioL >= DoubleLimits.MinNormalizedValue ) 
                    return  ratioL * legP;

                double value = Math.Exp(Math2.Lgamma(n+m+1) - Math2.Lgamma(n+1-m) + Math.Log(Math.Abs(legP)));
                return Math.Sign(legP)*value;
#else
                Policies.ReportNotImplementedError("LegendreP(n: {0}, m: {1}, x: {2}): Not implemented for |n| >= 128", n, m, x);
                return double.NaN;
#endif

            }

            // at this point n >= 0 && 0 <= m <= n
            Debug.Assert(n >= 0 && (m >= 0 && m <= n));

            // use the following recurrence relations:
            // P{n, n} = (-1)^n * (2*n - 1)!! * (1 - x^2)^(n/2) 
            // P{n+1, n} = (2*n + 1) * x * P{n, n}
            // (n-m+1)*P{n+1,m} = (2*n + 1)* x * P{n,m} - (n+m)*P{n-1, m}
            // see: http://en.wikipedia.org/wiki/Associated_Legendre_polynomials


            // the following is: 
            // pmm = Math2.Factorial2(2 * m - 1) * Math.Pow(1.0 - x*x, m/2.0);
            // being careful not to underflow too soon

            double pmm = 0;
            if (Math.Abs(x) > Constants.RecipSqrt2) {
                double inner = (1 - x) * (1 + x);
                int exp = 0;
                double mant = Math2.Frexp(inner, out exp);
                // careful: use integer division 
                double pmm_scaled = Math.Pow(mant, m / 2) * Math2.Factorial2(2 * m - 1);
                pmm = Math2.Ldexp(pmm_scaled, (m / 2) * exp);
                if (IsOdd(m))
                    pmm *= Math.Sqrt(inner);
            } else {
                double p_sin_theta = Math.Exp(m / 2.0 * Math2.Log1p(-x * x)); //Math.Pow(1.0 - x*x, m/2.0);
                pmm = Math2.Factorial2(2 * m - 1) * p_sin_theta;
            }

            if (IsOdd(m))
                pmm = -pmm;

            if (n == m)
                return pmm;

            // p0 - current
            // p1 - next

            double p0 = pmm;
            double p1 = (2.0 * m + 1.0) * x * p0; //P{m+1,m}(x)

            for (uint k = (uint)m + 1; k < (uint)n; k++) {
                double next = ((2 * k + 1) * x * p1 - (k + (uint)m) * p0) / (k - (uint)m + 1);
                p0 = p1;
                p1 = next;
            }


            return p1;
        }

        //
        // Spherical
        //


        /// <summary>
        /// Returns the associated legendre polynomial with angle parameterization P{l,m}( cos(theta) ) 
        /// </summary>
        /// <param name="n"></param>
        /// <param name="m"></param>
        /// <param name="theta"></param>
        /// <returns></returns>
        static double LegendrePAngle(int n, int m, double theta)
        {
            if (n < 0)
                n = -(n + 1); // P{n,m}(x) = P{-n-1,m}(x)

            if (m == 0)
                return Math2.LegendreP(n, Math.Cos(theta));

            // if ( abs(m) < n ) return 0 by definition; handling int.MinValue too
            if (m < -n || m > n)
                return 0;

            // P{n,-m} = (-1)^m * (n-m)!/(n+m)! * P{n,m}
            if (m < 0) {
                double legP = LegendrePAngle(n, -m, theta);
                if (legP == 0)
                    return legP;
                if (IsOdd(m))
                    legP = -legP;

                int MaxFactorial = Math2.FactorialTable.Length - 1;
                if (n - m <= MaxFactorial) {
                    double ratio = Math2.FactorialTable[n + m] / Math2.FactorialTable[n - m];
                    return ratio * legP;
                }

                if (n - m <= 300) {
                    // break up the Tgamma ratio into two separate factors
                    // to prevent an underflow: a!/b! = (a!/170!) * (170!/b!)
                    // 170!/300! ~= 2.37e-308 > double.MinNormalizedValue

                    return (legP * (Math2.FactorialTable[n + m] / Math2.FactorialTable[MaxFactorial])) * Math2.TgammaRatio(MaxFactorial + 1, n + 1 - m);

                }


#if false
                // since we are setting Lmax = 128, so Mmax = 128, and Lmax+Mmax = 256
                // so this should never be called, but left here for completeness
                // 170!/300! ~= 2.37e-308 > double.MinNormalizedValue

                double ratio = Math2.TgammaRatio(l+m+1, l+1-m);
                if ( ratio >= DoubleLimits.MinNormalizedValue ) 
                    return  ratio * legP;

                double value = Math.Exp(Math2.Lgamma(l+m+1) - Math2.Lgamma(l+1-m) + Math.Log(Math.Abs(legP)));
                return Math.Sign(legP)*value;
#else
                Policies.ReportNotImplementedError("LegendrePAngle(l: {0}, m: {1}, theta: {2}): Not implemented for |l| >= 128", n, m, theta);
                return double.NaN;
#endif

            }

            // at this point n >= 0 && 0 <= m <= n
            Debug.Assert(n >= 0 && (m >= 0 && m <= n));

            // use the recurrence relations
            // P{n, n} = (-1)^n * (2*n - 1)!! * (1 - x^2)^(n/2) 
            // P{n+1, n} = (2*n + 1) * x * P{n, n}
            // (n-m+1)*P{n+1,m} = (2*n + 1)* x * P{n,m} - (n+m)*P{n-1, m}
            // see: http://en.wikipedia.org/wiki/Associated_Legendre_polynomials

            double x = Math.Cos(theta);
            double sin_theta = Math.Sin(theta);

            // the following computes (2m-1)!! * sin_theta^m
            // taking special care not to underflow too soon
            // double pmm = Math.Pow(sin_theta, m) * Math2.Factorial2(2 * m - 1)

            int exp = 0;
            double mant = Math2.Frexp(sin_theta, out exp);
            double pmmScaled = Math.Pow(mant, m) * Math2.Factorial2(2 * m - 1);
            double pmm = Math2.Ldexp(pmmScaled, m * exp);

            if (IsOdd(m))
                pmm = -pmm;

            if (n == m)
                return pmm;

            // p0 - current
            // p1 - next

            double p0 = pmm;
            double p1 = (2.0 * m + 1.0) * x * p0; //P{m+1,m}(x)

            for (uint k = (uint)m + 1; k < (uint)n; k++) {

                double next = ((2 * k + 1) * x * p1 - (k + (uint)m) * p0) / (k - (uint)m + 1);

                p0 = p1;
                p1 = next;
            }

            return p1;
        }


        /// <summary>
        /// Returns the Spherical Harmonic Y{n,m}(theta,phi = 0)
        /// <para>Y{n,m}(theta,0) = Sqrt((2*n+1)/(4*π) * (n-m)!/(n+m)!)* P{n,m}(cos(theta))</para>
        /// </summary>
        /// <param name="n">Polynomial degree. Requires <paramref name="n"/> ≥ 0. Limited to <paramref name="n"/> ≤ 127</param>
        /// <param name="m">Polynomial order</param>
        /// <param name="theta">Colatitude, or polar angle, [0,π] where 0=North Pole, π=South Pole, π/2=Equator</param>
        public static double SphLegendre(int n, int m, double theta)
        {
            if ((n < 0)
            || (double.IsNaN(theta) || double.IsInfinity(theta))) {
                Policies.ReportDomainError("SphLegendre(n: {0}, m: {1}, theta: {2}): Requires n >= 0; finite theta", n, m, theta);
                return double.NaN;
            }

            // by definition: if ( abs(m) < n ) return 0
            if (m < -n || m > n)
                return 0;

            if (n >= 128) {
                Policies.ReportNotImplementedError("SphLegendre(n: {0}, m: {1}, theta: {2}): Not implemented for n >= 128", n, m, theta);
                return double.NaN;
            }

            bool sign = false;
            if (m < 0) {
                sign = IsOdd(m);
                m = Math.Abs(m);
            }

            // at this point n >= 0 && 0 <= m <= n
            Debug.Assert(n >= 0 && (m >= 0 && m <= n));
            Debug.Assert(n < 128);

            double result;
            double legP = LegendrePAngle(n, m, theta);

            int MaxFactorial = Math2.FactorialTable.Length - 1;
            if (n + m <= MaxFactorial) {
                double ratio = Math2.FactorialTable[n - m] / Math2.FactorialTable[n + m];
                result = 0.5 * Constants.RecipSqrtPI * Math.Sqrt((2 * n + 1) * ratio) * legP;
            } else if (n + m <= 300) {
                // break up the Tgamma ratio into two separate factors
                // to prevent an underflow: a!/b! = (a!/170!) * (170!/b!)
                // 170!/300! ~= 2.37e-308 > double.MinNormalizedValue

                double ratio1 = Math2.FactorialTable[n - m] / Math2.FactorialTable[MaxFactorial];
                double ratio2 = Math2.TgammaRatio(MaxFactorial + 1, n + m + 1);
                result = 0.5 * Constants.RecipSqrtPI * Math.Sqrt(ratio1) * Math.Sqrt((2 * n + 1) * ratio2) * legP;
            } else {
                // use logs 
                // note: with a limit of n=128, we should never reach here, but it is left for completeness
                double ratio = Math2.TgammaRatio(n - m + 1, n + m + 1);
                if (ratio < DoubleLimits.MinNormalValue) {
                    double mult = (2 * n + 1) / (4 * Math.PI);
                    result = Math.Exp((Math2.Lgamma(n - m + 1) - Math2.Lgamma(n + m + 1) + Math.Log(mult)) / 2 + Math.Log(Math.Abs(legP)));
                    result *= Math.Sign(legP);
                } else {
                    result = 0.5 * Constants.RecipSqrtPI * Math.Sqrt((2 * n + 1) * ratio) * legP;
                }
            }

            return sign ? -result : result;
        }

    }

} // namespaces





