//  Copyright (c) 2013 Rocco Romeo
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

//  Ported to .NET based on the original:
//      Copyright (c) 2006 Xiaogang Zhang, Boost Software License v1.0
// History:
//      Source originally by XZ. 
//      RR restructured code, added scaling, uniform asymptotics, and improved hankel asymptotics.

using System;
using System.Diagnostics;
using MathXK.Numerics;

namespace MathXK
{

    internal static partial class _Bessel
    {


        // Modified Bessel functions of the first and second kind of fractional order


        /// <summary>
        /// Returns (K{v}(x), K{v+1}(x)) 
        /// <para>|x| ≤ 2, the Temme series converges rapidly</para>
        /// <para>|x| &gt; 2, convergence slows as |x| increases</para>
        /// </summary>
        /// <param name="v"></param>
        /// <param name="x"></param>
        /// <returns></returns>
        static (double Kv, double Kvp1) K_Temme(double v, double x)
        {
            // see: Temme, Journal of Computational Physics, vol 21, 343 (1976)

            const double tolerance = DoubleLimits.MachineEpsilon;

            Debug.Assert(Math.Abs(x) <= 2);
            Debug.Assert(Math.Abs(v) <= 0.5);

            double gp = Math2.Tgamma1pm1(v);
            double gm = Math2.Tgamma1pm1(-v);

            double a = Math.Log(x / 2);
            double b = Math.Exp(v * a);
            double sigma = -a * v;
            double c = Math.Abs(v) < DoubleLimits.MachineEpsilon ? 1 : (Math2.SinPI(v) / (v * Math.PI));
            double d = Math2.Sinhc(sigma);
            double gamma1 = Math.Abs(v) < DoubleLimits.MachineEpsilon ? -Constants.EulerMascheroni : ((0.5 / v) * (gp - gm) * c);
            double gamma2 = (2 + gp + gm) * c / 2;

            // initial values
            double p = (gp + 1) / (2 * b);
            double q = (1 + gm) * b / 2;
            double f = (Math.Cosh(sigma) * gamma1 + d * (-a) * gamma2) / c;
            double h = p;
            double coef = 1;
            double sum = coef * f;
            double sum1 = coef * h;

            // series summation
            int k;
            for (k = 1; k < Policies.MaxSeriesIterations; k++) {
                f = (k * f + p + q) / (k * k - v * v);
                p /= k - v;
                q /= k + v;
                h = p - k * f;
                coef *= x * x / (4 * k);
                sum += coef * f;
                sum1 += coef * h;
                if (Math.Abs(coef * f) < Math.Abs(sum) * tolerance) {
                    break;
                }
            }


            if (k >= Policies.MaxSeriesIterations) {
                Policies.ReportConvergenceError("K_Temme(v:{0}, x:{1}) did not converge after {2} iterations", v, x, Policies.MaxSeriesIterations);
                return (double.NaN, double.NaN);
            }

            return (sum, 2 * sum1 / x);
        }


        /// <summary>
        /// Evaluate I{v+1}(x) / I{v}(x) using continued fractions
        /// <para>|x| ≤ |v|, CF1_ik converges rapidly</para>
        /// <para>|x| &gt; |v|, CF1_ik needs O(|x|) iterations to converge</para>
        /// </summary>
        /// <param name="v"></param>
        /// <param name="x"></param>
        /// <returns></returns>
        static double I_CF1(double v, double x)
        {
            const double tolerance = 2 * DoubleLimits.MachineEpsilon;

            double f;
            int k;

            // See Abramowitz and Stegun, Handbook of Mathematical Functions, 1972, 9.1.73

            // modified Lentz's method, see
            // Lentz, Applied Optics, vol 15, 668 (1976)
            double tiny = Math.Sqrt(DoubleLimits.MinNormalValue);

            double C = f = tiny;
            double D = 0;
            for (k = 1; k < Policies.MaxSeriesIterations; k++) {
                double a = 1;
                double b = 2 * (v + k) / x;
                C = b + a / C;
                D = b + a * D;
                if (C == 0) { C = tiny; }
                if (D == 0) { D = tiny; }
                D = 1 / D;
                double delta = C * D;
                f *= delta;
                if (Math.Abs(delta - 1) <= tolerance) {
                    break;
                }
            }

            if (k >= Policies.MaxSeriesIterations) {
                Policies.ReportConvergenceError("CF1(v:{0}, x:{1}) did not converge after {2} iterations", v, x, Policies.MaxSeriesIterations);
                return double.NaN;
            }

            return f;
        }


        /// <summary>
        /// Calculate (K{v}(x), K{v+1}(x)) by evaluating continued fraction
        /// z1 / z0 = U(v+1.5, 2v+1, 2x) / U(v+0.5, 2v+1, 2x)
        /// <para>|x| &gt;= |v|, CF2_ik converges rapidly</para>
        /// <para>|x| -&gt; 0, CF2_ik fails to converge</para>
        /// </summary>
        /// <param name="v"></param>
        /// <param name="x"></param>
        /// <param name="expScale">if true, exponentially scale the results</param>
        /// <returns>
        /// If <paramref name="expScale"/> is false K{v}(x) and K{v+1}(x); 
        /// otherwise e^x * K{v}(x) and e^x * K{v+1}(x)
        /// </returns>
        static (double Kv, double Kvp1) K_CF2(double v, double x, bool expScale)
        {
            // See Thompson and Barnett, Computer Physics Communications, vol 47, 245 (1987)

            const double tolerance = DoubleLimits.MachineEpsilon;
            double C, f, q, delta;
            int k;

            Debug.Assert(Math.Abs(x) > 1);

            // Steed's algorithm, see Thompson and Barnett,
            // Journal of Computational Physics, vol 64, 490 (1986)
            double a = v * v - 0.25;
            double b = 2 * (x + 1);
            double D = 1 / b;
            f = delta = D;                                // f1 = delta1 = D1, coincidence
            double prev = 0;
            double current = 1;
            double Q = C = -a;
            double S = 1 + Q * delta;

            // starting from 2
            for (k = 2; k < Policies.MaxSeriesIterations; k++) {
                // continued fraction f = z1 / z0
                a -= 2 * (k - 1);
                b += 2;
                D = 1 / (b + a * D);
                delta *= b * D - 1;
                f += delta;

                q = (prev - (b - 2) * current) / a;
                prev = current;
                current = q;                        // forward recurrence for q
                C *= -a / k;
                Q += C * q;
                S += Q * delta;

                // S converges slower than f

                if (Math.Abs(Q * delta) < Math.Abs(S) * tolerance) {
                    break;
                }
            }

            if (k >= Policies.MaxSeriesIterations) {
                Policies.ReportConvergenceError("K_CF2(v:{0}, x:{1}) did not converge after {2} iterations", v, x, Policies.MaxSeriesIterations);
                return (double.NaN, double.NaN);
            }

            double Kv, Kv1;
            if (expScale)
                Kv = (Constants.SqrtHalfPI / Math.Sqrt(x)) / S;
            else
                Kv = (Constants.SqrtHalfPI / Math.Sqrt(x)) * Math.Exp(-x) / S;

            Kv1 = Kv * (0.5 + v + x + (v * v - 0.25) * f) / x;

            return (Kv, Kv1);
        }


        /// <summary>
        /// Compute K{v}(x) and K{v+1}(x). Approximately O(v).
        /// <para>Multiply K{v}(x) and K{v+1}(x) by 2^binaryScale to get the true result</para>
        /// </summary>
        /// <param name="v"></param>
        /// <param name="x"></param>
        /// <returns></returns>
        static (double Kv, double Kvp1, int BScale) K_CF(double v, double x)
        {
            int binaryScale = 0;

            // Ku = K{u}, Ku1 = K{u+1}
            double Ku, Kup1;

            int n = (int)Math.Floor(v + 0.5);
            double u = v - n;                              // -1/2 <= u < 1/2

            // for x in (0, 2], use Temme series
            // otherwise use scaled continued fraction K_CF2
            // to prevent Kv from underflowing too quickly for large x
            bool expScale = false;
            if (x <= 2) {
                (Ku, Kup1) = K_Temme(u, x);
            } else {
                expScale = true;
                (Ku, Kup1) = K_CF2(u, x, expScale);
            }

            var (Kvp1, Kv, bScale) = Recurrence.ForwardK_B(u + 1, x, n, Kup1, Ku);


            if (expScale) {
                var sf = DoubleX.Ldexp(DoubleX.Exp(-x), bScale);
                Kv = Kv * sf.Mantissa;
                Kvp1 = Kvp1 * sf.Mantissa;
                binaryScale = sf.Exponent;
            } else {
                binaryScale = bScale;
            }


            return (Kv, Kvp1, binaryScale);
        }


        /// <summary>
        /// Simultaneously compute I{v}(x) and K{v}(x) 
        /// </summary>
        /// <param name="v"></param>
        /// <param name="x"></param>
        /// <param name="needI"></param>
        /// <param name="needK"></param>
        /// <returns>A tuple of (I, K). If needI == false or needK == false, I = NaN or K = NaN  respectively</returns>
        public static (double I, double K) IK(double v, double x, bool needI, bool needK)
        {
            Debug.Assert(needI || needK, "NeedI and NeedK cannot both be false");

            // set initial values for parameters
            double I = double.NaN; 
            double K = double.NaN;

            // check v first so that there are no integer throws later
            if (Math.Abs(v) > int.MaxValue) {
                Policies.ReportNotImplementedError("BesselIK(v: {0}, x: {1}): Requires |v| <= {2}", v, x, int.MaxValue);
                return (I, K);
            }


            if (v < 0) {
                v = -v;

                // for integer orders only, we can use the following identities:
                //  I{-n}(x) = I{n}(v) 
                //  K{-n}(x) = K{n}(v) 
                // for any v, use reflection rule:
                //  I{-v}(x) = I{v}(x) + (2/PI) * sin(pi*v)*K{v}(x)
                //  K{-v}(x) = K{v}(x)

                if (needI && !Math2.IsInteger(v)) {
                    var (IPos, KPos) = IK(v, x, true, true);

                    I = IPos + (2 / Math.PI) * Math2.SinPI(v) * KPos;
                    if (needK)
                        K = KPos;
                    return (I, K);
                }
            }

            if (x < 0) {
                Policies.ReportDomainError("BesselIK(v: {0}, x: {1}): Complex result not supported. Requires x >= 0", v, x);
                return (I, K);
            }

            // both x and v are non-negative from here
            Debug.Assert(x >= 0 && v >= 0);

            if (x == 0) {
                if (needI)
                    I = (v == 0) ? 1.0 : 0.0;
                if (needK) {
                    K = double.PositiveInfinity;
                }
                return (I, K);
            }

            if (needI && (x < 2 || 3*(v + 1) > x * x )) {
                I = I_SmallArg(v, x);
                needI = false;
                if ( !needK )
                    return (I, K);
            }

            // Hankel is fast, and reasonably accurate.
            // at x == 32, it will converge in about 15 iterations
            if (x >= HankelAsym.IKMinX(v)) { 
                if ( needK )
                    K = HankelAsym.K(v, x);
                if ( needI )
                    I = HankelAsym.I(v, x);
                return (I, K);
            }

            // the uniform expansion is here as a last resort 
            // to limit the number of recurrences, but it is less accurate.
            if (UniformAsym.IsIKAvailable(v, x)) {
                if (needK)
                    K = UniformAsym.K(v, x);
                if (needI)
                    I = UniformAsym.I(v, x);
                return (I, K);
            }


            // K{v}(x) and K{v+1}(x), binary scaled
            var (Kv, Kv1, binaryScale) = K_CF(v, x);
            if ( needK )
                K = Math2.Ldexp(Kv, binaryScale);
            if (needI) {
                // use the Wronskian relationship

                // Since CF1 is O(x) for x > v, try to weed out obvious overflows
                // Note: I{v+1}(x)/I{v}(x) is in [0, 1].
                // I{0}(713. ...) == MaxDouble
                const double I0MaxX = 713;
                if (x > I0MaxX &&  x > v) {
                    I = Math2.Ldexp(1.0 / (x * (Kv + Kv1)), -binaryScale);
                    if (double.IsInfinity(I))
                        return (I, K);
                }

                double W = 1 / x;                                 
                double fv = I_CF1(v, x);
                I = Math2.Ldexp(W / (Kv * fv + Kv1), -binaryScale);
            }


            return (I, K);
        }

    }

} // namespaces


