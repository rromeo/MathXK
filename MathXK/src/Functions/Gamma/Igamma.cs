//  Copyright (c) 2013 Rocco Romeo
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

//  Ported to .NET based on the original:
//      Copyright (c) John Maddock 2006-7, Boost Software License v1.0
//      Copyright (c) Paul A. Bristow 2007, Boost Software License v1.0

//  History:
//      JM & PB created original C++ code
//      RR ported to C# and restructured. Also added asymptotic series.

#if DEBUG
//#define EXTRA_DEBUG
#endif

using System;
using System.Diagnostics;
using MathXK.Numerics;

namespace MathXK
{

    static partial class _Igamma
    {
        enum Routine
        {
            Upper = 0,
            Lower = 1,
            Q = 2,
            P = 3
        };

        /// <summary>
        /// Upper incomplete gamma fraction
        /// <para>Γ(a,z) = (z^a * e^-z) * UpperFraction(a,z)</para>
        /// <para>Q(a,z) = (z^a * e^-z)/Γ(a) * UpperFraction(a,z)</para>
        /// </summary>
        /// <param name="a"></param>
        /// <param name="z"></param>
        /// <returns></returns>
        public static double UpperFraction(double a, double z)
        {
            double z0 = z - a + 1;
            int k = 0;
            Func<(double, double)> f = () => {
                k++;
                return (k * (a - k), z0 + 2 * k);
            };
            return 1.0 / ContinuedFraction.Eval(z0, f);

        }

        /// <summary>
        /// Lower gamma series sum: Σ( z^k/(a+1)_k) k={0, inf}
        /// <para>γ(a,z) = ((z^a) * (e^-z) / a) * LowerSeriesSum(a,z)</para>
        /// </summary>
        /// <param name="a"></param>
        /// <param name="z"></param>
        /// <returns></returns>
        /// <see href="http://dlmf.nist.gov/8.5#E1"/>
        public static double LowerSeries(double a, double z)
        {
            // Sum1F1 is optimized for a1 == 1
            return HypergeometricSeries.Sum1F1(1, 1 + a, z);
        }

        /// <summary>
        /// Returns the minimum z value to use the small a large z asymptotic series.
        /// <para>Currently set to Max(50, a*10)</para>
        /// </summary>
        /// <param name="a"></param>
        /// <returns></returns>
        static double Asym_MinLargeZ(double a)
        {
            return Math.Max(50, a * 10);
        }

        /// <summary>
        /// Asymptotic series for Γ(a, z) with z &gt; a
        /// <para>Currently set to z &gt;= Max(50, a * 10)</para>
        /// <para>Σ( (-1)^k * (1-a)_{k} * z^-k, k={0, inf})</para>
        /// <para>Γ(a, z) = e^-z * z^(a-1) * TgammaAsymSeries(a,z)</para>
        /// </summary>
        /// <param name="a"></param>
        /// <param name="z"></param>
        /// <returns></returns>
        /// <seealso href="http://functions.wolfram.com/GammaBetaErf/Gamma2/06/02/02/"/>
        /// <remarks>
        /// Targeting n=20: As a-> 0, n!/z^n &lt; eps, so z = 50
        /// Targeting n=16: As a-> inf, (a/z)^n &lt; eps, z ~= a*2^(52/16)
        /// Result ~= 1
        /// </remarks>
        public static double Asym_SeriesLargeZ(double a, double z)
        {
            //Sum2F0 is optimized for a0=1
            return HypergeometricSeries.Sum2F0(1, 1 - a, -1 / z);
        }


        /// <summary>
        /// Series representation for γ(a,x)-1/a for small a
        /// <para>LowerSmallASeriesSum = Σ((-x)^k/((a+k)*k!)) k={1, inf}</para>
        /// <para>γ(a,x) = z^a * LowerSmallASeriesSum(a, x, 1/a)</para>
        /// </summary>
        /// <param name="a"></param>
        /// <param name="z"></param>
        /// <param name="initValue"></param>
        /// <returns></returns>
        /// <seealso href="http://functions.wolfram.com/GammaBetaErf/Gamma2/06/01/04/01/01/0003/"/>
        public static double LowerSmallASeries(double a, double z, double initValue)
        {
            double term = 1.0;
            double sum = initValue;

            for (int n = 1; n < Policies.MaxSeriesIterations; n++) {

                double prevSum = sum;
                term *= (-z / n);
                sum += term / (a + n);

                if (Math.Abs(term) <= Math.Abs(prevSum) * Policies.SeriesTolerance)
                    return sum;
            }

            Policies.ReportConvergenceError("Series did not converge after {0} iterations", Policies.MaxSeriesIterations);
            return double.NaN;
        }


        /// <summary>
        /// Calculate mult * (z^a)(e^-z). Used in the non-regularized incomplete gammas
        /// </summary>
        /// <param name="a"></param>
        /// <param name="z"></param>
        /// <param name="mult"></param>
        /// <returns></returns>
        public static double Prefix(double a, double z, double mult = 1)
        {
            Debug.Assert(a >= 0);
            Debug.Assert(z >= 0);

            if (z == 0)
                return 0;
            if (double.IsInfinity(z))
                return 0;

            double prefix;
            if ((a <= 1 && z < -DoubleLimits.MinLogValue) || (z <= 1 && mult <= 1)) {
                // can't overflow
                prefix = Math.Pow(z, a) * Math.Exp(-z) * mult;

            } else {

                double alz = a * Math.Log(z);

                if ((alz > DoubleLimits.MinLogValue && alz < DoubleLimits.MaxLogValue) && (z < -DoubleLimits.MinLogValue)) {
                    prefix = Math.Pow(z, a) * Math.Exp(-z) * mult;
                } else if ((alz > 2 * DoubleLimits.MinLogValue && alz < 2 * DoubleLimits.MaxLogValue) && (z < -2 * DoubleLimits.MinLogValue)) {
                    double rp = Math.Pow(z, a / 2) * Math.Exp(-z / 2);
                    prefix = (mult * rp) * rp;
                } else if ((alz > 4 * DoubleLimits.MinLogValue && alz < 4 * DoubleLimits.MaxLogValue) && (z < -4 * DoubleLimits.MinLogValue)) {
                    double rp = Math.Pow(z, a / 4) * Math.Exp(-z / 4);
                    prefix = (mult * rp) * rp * rp * rp;
                } else {
                    prefix = Math.Exp(alz - z + Math.Log(mult));
                }
            }

            return prefix;

        }


        /// <summary>
        /// Compute (z^a)(e^-z)/Γ(a) used in the regularized incomplete gammas
        /// </summary>
        /// <param name="a"></param>
        /// <param name="z"></param>
        /// <returns></returns>
        /// <remarks>Note: most of the error occurs in this function</remarks>
        public static double PrefixRegularized(double a, double z)
        {
            Debug.Assert(a > 0);
            Debug.Assert(z >= 0);

            // if we can't compute Tgamma directly, use Stirling
            if (a > DoubleLimits.MaxGamma || (a >= StirlingGamma.LowerLimit && z > -DoubleLimits.MinLogValue ))
                return StirlingGamma.PowExpDivGamma(a, z);

            // avoid overflows from dividing small Gamma(a)
            // for z > 1, e^(eps*ln(z)-z) ~= e^-z, so underflow is ok
            if (a < DoubleLimits.MachineEpsilon/Constants.EulerMascheroni)
                return a * Math.Pow(z, a) * Math.Exp(-z);

            double alz = a * Math.Log(z);
            double g = Math2.Tgamma(a);

            double prefix;
            if ((alz > DoubleLimits.MinLogValue && alz < DoubleLimits.MaxLogValue) && (z < -DoubleLimits.MinLogValue)) {
                prefix = Math.Pow(z, a) * Math.Exp(-z) / g;
            } else if ((alz > 2 * DoubleLimits.MinLogValue && alz < 2 * DoubleLimits.MaxLogValue) && (z < -2 * DoubleLimits.MinLogValue)) {
                double rp = Math.Pow(z, a / 2) * Math.Exp(-z / 2);
                prefix = (rp/g) * rp;
            } else if ((alz > 4 * DoubleLimits.MinLogValue && alz < 4 * DoubleLimits.MaxLogValue) && (z < -4 * DoubleLimits.MinLogValue)) {
                double rp = Math.Pow(z, a / 4) * Math.Exp(-z / 4);
                prefix = (rp/g) * rp * rp * rp;
            } else {
                prefix = Math.Exp(alz - z - Math.Log(g));
            }

            return prefix;
        }

        /// <summary>
        /// Computes the specified Gamma function when a is small. Uses:
        /// <para>LowerSmallASeriesSum = Σ((-x)^k/((a+k)*k!)) k={1, inf}</para>
        /// <para>γ(a,x) = z^a * LowerSmallASeriesSum(a, x, 1/a)</para>
        /// </summary>
        /// <param name="r"></param>
        /// <param name="a"></param>
        /// <param name="x"></param>
        /// <returns></returns>
        static double SmallA(Routine r, double a, double x)
        {
            if (r == Routine.P || r == Routine.Lower) {
                double l = Math.Pow(x, a) * LowerSmallASeries(a, x, 1.0 / a);
                return (r == Routine.Lower) ? l : l / Math2.Tgamma(a);
            } else {
                double p = Math2.Powm1(x, a);
                double g = Math2.Tgamma1pm1(a);
                double initValue = (g - p) / a;

                double u = initValue - (p + 1) * LowerSmallASeries(a, x, 0.0);
                return (r == Routine.Upper) ? u : u * a / (g + 1);
            }
        }


        //
        // Upper gamma fraction for integer a:
        //

        /// <summary>
        /// Returns initValue + Σ(x^k/k!) k={0, n} 
        /// </summary>
        /// <param name="n"></param>
        /// <param name="x"></param>
        /// <param name="initValue"></param>
        /// <returns></returns>
        static double ExponentialSeries(int n, double x, double initValue = 0)
        {
            Debug.Assert(n >= 0);

            double sum = 1.0 + initValue;
            double term = 1.0;
            for (int k = 1; k <= n; k++) {
                term *= x / k;
                double lastSum = sum;
                sum += term;
                if (sum == lastSum)
                    break;
            }

            return sum;
        }


        /// <summary>
        /// Calculates the specified Gamma function when a is an integer
        /// <para>Q(a,x) = e^-x * (Σ(x^k/k!) k={0,a-1})</para>
        /// </summary>
        /// <param name="r"></param>
        /// <param name="a"></param>
        /// <param name="x"></param>
        /// <returns></returns>
        static double IntegerA(Routine r, double a, double x)
        {
            if (r == Routine.Q || r == Routine.Upper) {
                double q = Math.Exp(-x) * ExponentialSeries((int)(a - 1), x);
                return (r == Routine.Q) ? q : q * Math2.Tgamma(a);
            } else {
                // avoid this routine for small x, because there's too much cancellation
                // unsuccessfuly tried:
                //p = -Math2.Expm1(-x)  - Math.Exp(-x)*ExponentialSeries(a - 1, x, -1);
                double p = 1.0 - Math.Exp(-x) * ExponentialSeries((int)(a - 1), x);
                return (r == Routine.P) ? p : p * Math2.Tgamma(a);
            }

        }


        /// <summary>
        /// Calculates the specified Gamma function when <paramref name="a"/> is a half integer
        /// </summary>
        /// <param name="r"></param>
        /// <param name="a"></param>
        /// <param name="x"></param>
        /// <returns></returns>
        static double HalfIntegerA(Routine r, double a, double x)
        {
            Debug.Assert(a - Math.Floor(a) == 0.5);

            double sum = 0;
            if (a > 1) {
                sum = 1.0;
                double term = 1.0;
                for (int n = 2; n < a; ++n) {
                    term *= x / (n - 0.5);
                    sum += term;
                }
            }

            double result;
            double mult = 2.0 * Constants.RecipSqrtPI * Math.Sqrt(x);
            if (r == Routine.Q || r == Routine.Upper) {
                result = Math2.Erfc(Math.Sqrt(x)) + sum * mult * Math.Exp(-x);
                return (r == Routine.Q) ? result : result * Math2.Tgamma(a);
            } else {
                result = Math2.Erf(Math.Sqrt(x)) - sum * mult * Math.Exp(-x);
                return (r == Routine.P) ? result : result * Math2.Tgamma(a);
            }
        }


        /// <summary>
        /// Compute Γ(a,x) = Γ(a) - γ(a,x). Use this function when x is small
        /// </summary>
        /// <param name="a"></param>
        /// <param name="x"></param>
        /// <returns></returns>
        static double TgammaMinusLowerSeries(double a, double x)
        {
            double lowerSeriesSum = LowerSeries(a, x);
            if (a <= DoubleLimits.MaxGamma )
                return Math2.Tgamma(a) - Prefix(a, x) / a * lowerSeriesSum;

            // Since we could overflow too soon, 
            // or get  (inf - inf) = NaN, instead of inf
            // use logs: Log(Γ(a)*(1-P(a,x)))

            return Math.Exp(Math2.Lgamma(a) + Math2.Log1p(-PrefixRegularized(a, x) / a * lowerSeriesSum));
        }

        /// <summary>
        /// Compute γ(a,x) = Γ(a) - Γ(a,x). Use this function when x is large
        /// </summary>
        /// <param name="a"></param>
        /// <param name="x"></param>
        /// <returns></returns>
        static double TgammaMinusUpperFraction(double a, double x)
        {
            double upperFraction = UpperFraction(a, x);
            if (a <= DoubleLimits.MaxGamma )
                return Math2.Tgamma(a) - Prefix(a, x) * upperFraction;

            // Since we could overflow too soon, 
            // or get  (inf - inf) = NaN, instead of inf
            // use logs: Log(Γ(a)*(1-Q(a,x)))

            return Math.Exp(Math2.Lgamma(a) + Math2.Log1p(-PrefixRegularized(a, x) * upperFraction));
        }



        public static double TgammaImp(double a, double x)
        {
            const Routine routine = Routine.Upper;

            // special values
            if (x == 0)
                return Math2.Tgamma(a);
            if (double.IsInfinity(x))
                return 0.0;

            // Γ(0,x) = E_1(x)
            if (a == 0)
                return Math2.Expint(1, x);
            if (a == 1.0)
                return Math.Exp(-x);

            // Process small a, large x with asymptotic series, which is faster than CF
            if (x >= Asym_MinLargeZ(a))
                return Asym_SeriesLargeZ(a, x) * Prefix(a, x)/x;

            if (a == 0.5)
                return Constants.SqrtPI * Math.Exp(-x) * Math2.Erfcx(Math.Sqrt(x));


            // is a small
            if ((a < 30) && (x >= a + 1) && (x < DoubleLimits.MaxLogValue)) {
                double frac = a - Math.Floor(a);
                if (frac == 0) {
                    return IntegerA(routine, a, x);
                } else if (frac == 0.5) {
                    return HalfIntegerA(routine, a, x);
                }
            }

            // Use the series 
            if (x < (a + 1)) {
                if (a < 1)
                    return SmallA(routine, a, x);
                return TgammaMinusLowerSeries(a, x);
            }

            // Use CF for x > a+1 
            return Prefix(a, x, UpperFraction(a, x));

        }

        public static double TgammaLowerImp(double a, double x)
        {
            const Routine routine = Routine.Lower;

            // special values
            if (x == 0.0)
                return 0.0;
            if (double.IsInfinity(x))
                return Math2.Tgamma(a);

            if (a == 1.0)
                return -Math2.Expm1(-x); // 1.0 - Math.Exp(-x);
            if (a == 0.5)
                return Constants.SqrtPI * Math2.Erf(Math.Sqrt(x));

            // Process small a, large x with asymptotic series, which is faster than CF
            if (x >= Asym_MinLargeZ(a))
                return Math2.Tgamma(a) - Asym_SeriesLargeZ(a, x) * Prefix(a, x) / x;

            // is a small
            if ((a < 30) && (x >= a + 1) && (x < DoubleLimits.MaxLogValue)) {
                double frac = a - Math.Floor(a);
                if (frac == 0) {
                    return IntegerA(routine, a, x);
                } else if (frac == 0.5) {
                    return HalfIntegerA(routine, a, x);
                }
            }

            if (x < (a + 1)) {
                if (a < 1)
                    return SmallA(routine, a, x);
                return Prefix(a, x, 1/a) * LowerSeries(a, x);
            }

            // Use CF for x > a+1 
            return TgammaMinusUpperFraction(a, x);
        }

        public static double GammaQImp(double a, double x)
        {
            const Routine routine = Routine.Q;

            // special values
            if (x == 0)
                return 1.0;
            if (double.IsInfinity(x))
                return 0.0;

            if (a == 0)
                return 0;
            if (a == 1.0)
                return Math.Exp(-x);

            // Process small a, large x with asymptotic series, which is faster than CF
            if (x >= Asym_MinLargeZ(a))
                return Asym_SeriesLargeZ(a, x) * PrefixRegularized(a, x) / x;

            if (a == 0.5)
                return Math.Exp(-x) * Math2.Erfcx(Math.Sqrt(x));

            // is a small
            if ((a < 30) && (x >= a + 1) && (x < DoubleLimits.MaxLogValue)) {
                double frac = a - Math.Floor(a);
                if (frac == 0) {
                    return IntegerA(routine, a, x);
                } else if (frac == 0.5)
                    return HalfIntegerA(routine, a, x);
            }



            //
            // Begin by testing whether we're in the "bad" zone
            // where the result will be near 0.5 and the usual
            // series and continued fractions are slow to converge:
            //
            if (a > 20) {

                // The second limit below is chosen so that we use Temme's expansion
                // only if the result would be larger than about 10^-6.
                // Below that the regular series and continued fractions
                // converge OK, and if we use Temme's method we get increasing
                // errors from the dominant erfc term as it's (inexact) argument
                // increases in magnitude.

                double sigma = Math.Abs((x - a) / a);
                if (sigma < 0.4 || (a > 200 && 20 / a > sigma * sigma)) {
                    double t = TemmeSymmetricAsym.GammaP(a, x);
                    return (x >= a) ? t : 1 - t;
                }

            }

            if (x < (a + 1)) {
                if (a < 1)
                    return SmallA(routine, a, x);
                return 1.0 - PrefixRegularized(a, x) * LowerSeries(a, x) / a;
            }

            // Use CF for x > a+1 
            return PrefixRegularized(a, x) * UpperFraction(a, x);

        }

        public static double GammaPImp(double a, double x)
        {
            const Routine routine = Routine.P;

            // special values
            if (x == 0.0)
                return 0.0;
            if (double.IsInfinity(x))
                return 1.0;

            if (a == 0)
                return 1;
            if (a == 1.0)
                return -Math2.Expm1(-x); // 1.0 - Math.Exp(-x)
            if (a == 0.5)
                return Math2.Erf(Math.Sqrt(x));

            // Process small a, large x with asymptotic series, which is faster than CF
            if (x >= Asym_MinLargeZ(a))
                return 1 - Asym_SeriesLargeZ(a, x) * PrefixRegularized(a, x) / x;

            // is a small
            if ((a < 30) && (x >= a + 1) && (x < DoubleLimits.MaxLogValue)) {
                double frac = a - Math.Floor(a);
                if (frac == 0) {
                    return IntegerA(routine, a, x);
                } else if (frac == 0.5) {
                    return HalfIntegerA(routine, a, x);
                }
            }

            //
            // Begin by testing whether we're in the "bad" zone
            // where the result will be near 0.5 and the usual
            // series and continued fractions are slow to converge:
            //
            if (a > 20) {

                // This second limit below is chosen so that we use Temme's expansion
                // only if the result would be larger than about 10^-6.
                // Below that the regular series and continued fractions
                // converge OK, and if we use Temme's method we get increasing
                // errors from the dominant erfc term as it's (inexact) argument
                // increases in magnitude.

                double sigma = Math.Abs((x - a) / a);
                if (sigma < 0.4 || (a > 200 && 20 / a > sigma * sigma)) {
                    double t = TemmeSymmetricAsym.GammaP(a, x);
                    return (x >= a) ? 1 - t : t;
                }

            }

            if (x < (a + 1)) {
                if (a < 1)
                    return SmallA(routine, a, x);
                return PrefixRegularized(a, x) * LowerSeries(a, x) / a;
            }

            // Use CF for x > a+1 
            return 1.0 - PrefixRegularized(a, x) * UpperFraction(a, x);

        }



    }


    public static partial class Math2
    {


        /// <summary>
        /// Returns ∂P(a,x)/∂x = e^(-x) * x^(a-1) / Γ(a)
        /// <para>Note: ∂Q(a,x)/∂x = -∂P(a,x)/∂x</para>
        /// </summary>
        /// <param name="a">Requires finite a &gt; 0</param>
        /// <param name="x">Requires x ≥ 0</param>
        public static double GammaPDerivative(double a, double x)
        {
            if ((!(a > 0) || double.IsInfinity(a))
            || (!(x >= 0))) {
                Policies.ReportDomainError("GammaPDerivative(a: {0}, x: {1}): Requires finite a > 0; x >= 0", a, x);
                return double.NaN;
            }

            //
            // Now special cases:
            //
            if (x == 0) {
                if (a > 1)
                    return 0;
                if (a == 1)
                    return 1;
                return double.PositiveInfinity;
            }

            // lim {x->+inf} e^(-x + (a-1)*ln(x)) / Γ(a) = 0/Γ(a) = 0 
            if (double.IsInfinity(x))
                return 0;

            // when x is small, we can underflow too quickly
            // For x < eps, e^-x ~= 1, so use 
            //      x^(a-1)/Γ(a)
            // and for eps < x < 1, use:
            //      IgammaPrefixRegularized(a, x) = (x/(a-1)) * IgammaPrefixRegularized(a-1, x)

            if (x <= DoubleLimits.MachineEpsilon) {
                if (a < DoubleLimits.MachineEpsilon/Constants.EulerMascheroni)
                    return (a / x) * Math.Pow(x, a);
                if (a < 1)
                    return Math.Pow(x, a) / (x * Math2.Tgamma(a));
                return Math.Pow(x, a - 1) / Math2.Tgamma(a);
            }

            if (x < 1 && a >= 2)
                return _Igamma.PrefixRegularized(a - 1, x) / (a - 1);

            // Normal case
            double result = _Igamma.PrefixRegularized(a, x) / x;

            return result;
        }

        /// <summary>
        /// Returns the full (non-normalised) upper incomplete gamma function 
        /// <para>Γ(a,x) = ∫ t^(a-1) * e^-t dt, t={x, ∞}</para>
        /// </summary>
        /// <param name="a">Requires finite a ≥ 0</param>
        /// <param name="x">Requires x ≥ 0</param>
        public static double Tgamma(double a, double x)
        {
            if ((!(a >= 0) || double.IsInfinity(a))
            || (!(x >= 0))) {
                Policies.ReportDomainError("Tgamma(a: {0}, x: {1}): Requires finite a >= 0; x >= 0", a, x);
                return double.NaN;
            }
            if (a == 0 && x == 0) {
                Policies.ReportDomainError("Tgamma(a: {0}, x: {1}): Requires that both a,x are not zero", a, x);
                return double.NaN;
            }

            return _Igamma.TgammaImp(a, x);

        }


        /// <summary>
        /// Returns the full (non-normalised) lower incomplete gamma function
        /// <para> γ(a,x) = ∫ t^(a-1) * e^-t dt, t={0, x}</para>
        /// </summary>
        /// <param name="a">Requires finite a &gt; 0</param>
        /// <param name="x">Requires x ≥ 0</param>
        public static double TgammaLower(double a, double x)
        {
            if ((!(a > 0) || double.IsInfinity(a))
            || (!(x >= 0))) {
                Policies.ReportDomainError("TgammaLower(a: {0}, x: {1}): Requires finite a > 0; x >= 0", a, x);
                return double.NaN;
            }

            return _Igamma.TgammaLowerImp(a, x);
        }

        /// <summary>
        /// Returns the regularized upper incomplete gamma function
        /// <para>Q(a,x) = Γ(a,x) / Γ(a)</para>
        /// </summary>
        /// <param name="a">Requires finite a &gt; 0</param>
        /// <param name="x">Requires x ≥ 0</param>
        public static double GammaQ(double a, double x)
        {
            if ((!(a > 0) || double.IsInfinity(a))
            || (!(x >= 0))) {
                Policies.ReportDomainError("GammaQ(a: {0}, x: {1}): Requires finite a > 0; x >= 0", a, x);
                return double.NaN;
            }

            return _Igamma.GammaQImp(a, x);

        }

        /// <summary>
        /// Returns the regularized lower incomplete gamma function 
        /// <para>P(a,x) = γ(a,x) / Γ(a)</para>
        /// </summary>
        /// <param name="a">Requires finite a &gt; 0</param>
        /// <param name="x">Requires x ≥ 0</param>
        public static double GammaP(double a, double x)
        {
            if ((!(a > 0) || double.IsInfinity(a))
            || (!(x >= 0))) {
                Policies.ReportDomainError("GammaP(a: {0}, x: {1}): Requires finite a > 0; x >= 0", a, x);
                return double.NaN;
            }

            return _Igamma.GammaPImp(a, x);
        }
    }



} // namespaces





