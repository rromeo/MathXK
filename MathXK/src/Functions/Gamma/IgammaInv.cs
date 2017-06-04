//  Copyright (c) 2013 Rocco Romeo
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

//  Ported to .NET based on the original:
//      Copyright (c) John Maddock 2006, Boost Software License v1.0

#if DEBUG
//#define EXTRA_DEBUG
#endif

using System;
using System.Runtime.InteropServices;
using MathXK;
using MathXK.Roots;
using System.Diagnostics;

namespace MathXK {

    static class _IgammaInv {
        /// <summary>
        /// Use the root finder to solve P(x,a) == p, for x
        /// </summary>
        /// <param name="a"></param>
        /// <param name="p"></param>
        /// <param name="guess"></param>
        /// <param name="min"></param>
        /// <param name="max"></param>
        /// <returns></returns>
        public static double SolveGivenAP(double a, double p, double guess, double min, double max)
        {
            Func<double, ValueTuple<double, double, double>> halleyF = (double x) => {
                // Calculate P(x) - p and the first two derivatives

                double f = Math2.GammaP(a, x) - p;
                double f1 = Math2.GammaPDerivative(a, x);
                double f2 = -f1 * (1.0 + (1.0 - a) / x);
                return (f, f1, f2);
            };

            const double relTolerance = RootFinder.DefaultRelTolerance;
            const double absTolerance = 0;
            RootResults rr = RootFinder.Halley(halleyF, guess, min, max, relTolerance, absTolerance, Policies.MaxRootIterations);
            if (rr == null) {
                Policies.ReportRootNotFoundError("Invalid parameter in root solver");
                return double.NaN;
            }

#if EXTRA_DEBUG
            Debug.WriteLine("Halley iterations: {0}",rr.Iterations);
#endif

            if (!rr.Success) {
                Policies.ReportRootNotFoundError("Root not found after {0} iterations", rr.Iterations);
                return double.NaN;
            }

            return rr.SolutionX;
        }

        /// <summary>
        /// Use the root finder to solve Q(a,x) == q, for x
        /// </summary>
        /// <param name="a"></param>
        /// <param name="q"></param>
        /// <param name="guess"></param>
        /// <param name="min"></param>
        /// <param name="max"></param>
        /// <returns></returns>
        public static double SolveGivenAQ(double a, double q, double guess, double min, double max)
        {
            Func<double, ValueTuple<double, double, double>> halleyF = (double x) => {
                // Calculate Q(x) - q and the first two derivatives

                double f = Math2.GammaQ(a, x) - q;
                double f1 = -Math2.GammaPDerivative(a, x);
                double f2 = -f1 * (1.0 + (1.0 - a) / x);
                return (f, f1, f2);
            };


            const double relTolerance = RootFinder.DefaultRelTolerance;
            const double absTolerance = 0;
            RootResults rr = RootFinder.Halley(halleyF, guess, min, max, relTolerance, absTolerance, Policies.MaxRootIterations);
            if (rr == null) {
                Policies.ReportRootNotFoundError("Invalid parameter in root solver");
                return double.NaN;
            }
#if EXTRA_DEBUG
            Debug.WriteLine("Halley iterations: {0}",rr.Iterations);
#endif
            if (!rr.Success) {
                Policies.ReportRootNotFoundError("Root not found after {0} iterations", rr.Iterations);
                return double.NaN;
            }

            return rr.SolutionX;
        }

        static double FindInverseS(double p, double q)
        {
            //
            // Computation of the Incomplete Gamma Function Ratios and their Inverse
            // ARMIDO R. DIDONATO and ALFRED H. MORRIS, JR.
            // ACM Transactions on Mathematical Software, Vol. 12, No. 4,
            // December 1986, Pages 377-393.
            //
            // See equation 32.
            //

            double t;
            if (p < 0.5) {
                t = Math.Sqrt(-2 * Math.Log(p));
            } else {
                t = Math.Sqrt(-2 * Math.Log(q));
            }

            const double p0 = 3.31125922108741;
            const double p1 = 11.6616720288968;
            const double p2 = 4.28342155967104;
            const double p3 = 0.213623493715853;

            const double q0 = 1;
            const double q1 = 6.61053765625462;
            const double q2 = 6.40691597760039;
            const double q3 = 1.27364489782223;
            const double q4 = 0.3611708101884203e-1;

            double P = p0 + t * (p1 + t * (p2 + t * p3));
            double Q = q0 + t * (q1 + t * (q2 + t * (q3 + t * q4)));

            double s = t - P / Q;
            if (p < 0.5)
                s = -s;
            return s;
        }

        static double DidonatoSN(double a, double x, int N, double tolerance = 0)
        {
            Debug.Assert(N >= 0);

            //
            // Computation of the Incomplete Gamma Function Ratios and their Inverse
            // ARMIDO R. DIDONATO and ALFRED H. MORRIS, JR.
            // ACM Transactions on Mathematical Software, Vol. 12, No. 4,
            // December 1986, Pages 377-393.
            //
            // See equation 34.
            //
            double sum = 1;
            if (N >= 1) {
                double partial = x / (a + 1);
                sum += partial;
                for (int i = 2; i <= N; ++i) {
                    partial *= x / (a + i);
                    sum += partial;
                    if (partial < tolerance)
                        break;
                }
            }
            return sum;
        }

#if false
        // unreferenced, but kept for future reference

        static double DidonatoFN(double p, double a, double x, uint N, double tolerance)
        {
           //
           // Computation of the Incomplete Gamma Function Ratios and their Inverse
           // ARMIDO R. DIDONATO and ALFRED H. MORRIS, JR.
           // ACM Transactions on Mathematical Software, Vol. 12, No. 4,
           // December 1986, Pages 377-393.
           //
           // See equation 34.
           //
   
           double u = Math.Log(p) + Math2.Lgamma(a + 1);
           return Math.Exp((u + x - Math.Log(didonato_SN(a, x, N, tolerance))) / a);
        }
#endif

        public static double GammaPInvGuess(double a, double p, double q, out bool p_has_10_digits)
        {
            //
            // In order to understand what's going on here, you will
            // need to refer to:
            //
            // Computation of the Incomplete Gamma Function Ratios and their Inverse
            // ARMIDO R. DIDONATO and ALFRED H. MORRIS, JR.
            // ACM Transactions on Mathematical Software, Vol. 12, No. 4,
            // December 1986, Pages 377-393.
            //


            double result;
            p_has_10_digits = false;

            if (a == 1) {
                result = -Math.Log(q);
            } else if (a < 1) {
                double g = Math2.Tgamma(a);
                double b = q * g;
                if ((b > 0.6) || ((b >= 0.45) && (a >= 0.3))) {
                    // DiDonato & Morris Eq 21: 
                    //
                    // There is a slight variation from DiDonato and Morris here:
                    // the first form given here is unstable when p is close to 1,
                    // making it impossible to compute the inverse of Q(a,x) for small
                    // q.  Fortunately the second form works perfectly well in this case.
                    //
                    double u;
                    if ((b * q > 1e-8) && (q > 1e-5)) {
                        u = Math.Pow(p * g * a, 1 / a);
                    } else {
                        u = Math.Exp((-q / a) - Constants.EulerMascheroni);
                    }
                    result = u / (1 - (u / (a + 1)));
                } else if ((a < 0.3) && (b >= 0.35)) {
                    // DiDonato & Morris Eq 22:
                    double t = Math.Exp(-Constants.EulerMascheroni - b);
                    double u = t * Math.Exp(t);
                    result = t * Math.Exp(u);
                } else if ((b > 0.15) || (a >= 0.3)) {
                    // DiDonato & Morris Eq 23:
                    double y = -Math.Log(b);
                    double u = y - (1 - a) * Math.Log(y);
                    result = y - (1 - a) * Math.Log(u) - Math.Log(1 + (1 - a) / (1 + u));
                } else if (b > 0.1) {
                    // DiDonato & Morris Eq 24:
                    double y = -Math.Log(b);
                    double u = y - (1 - a) * Math.Log(y);
                    result = y - (1 - a) * Math.Log(u) - Math.Log((u * u + 2 * (3 - a) * u + (2 - a) * (3 - a)) / (u * u + (5 - a) * u + 2));
                } else {
                    // DiDonato & Morris Eq 25:
                    double y = -Math.Log(b);
                    double c1 = (a - 1) * Math.Log(y);
                    double c1_2 = c1 * c1;
                    double c1_3 = c1_2 * c1;
                    double c1_4 = c1_2 * c1_2;
                    double a_2 = a * a;
                    double a_3 = a_2 * a;

                    double c2 = (a - 1) * (1 + c1);
                    double c3 = (a - 1) * (-(c1_2 / 2) + (a - 2) * c1 + (3 * a - 5) / 2);
                    double c4 = (a - 1) * ((c1_3 / 3) - (3 * a - 5) * c1_2 / 2 + (a_2 - 6 * a + 7) * c1 + (11 * a_2 - 46 * a + 47) / 6);
                    double c5 = (a - 1) * (-(c1_4 / 4)
                                      + (11 * a - 17) * c1_3 / 6
                                      + (-3 * a_2 + 13 * a - 13) * c1_2
                                      + (2 * a_3 - 25 * a_2 + 72 * a - 61) * c1 / 2
                                      + (25 * a_3 - 195 * a_2 + 477 * a - 379) / 12);

                    double y_2 = y * y;
                    double y_3 = y_2 * y;
                    double y_4 = y_2 * y_2;
                    result = y + c1 + (c2 / y) + (c3 / y_2) + (c4 / y_3) + (c5 / y_4);
                    if (b < 1e-28)
                        p_has_10_digits = true;
                }
            } else {
                // DiDonato and Morris Eq 31:
                double s = FindInverseS(p, q);

                double s_2 = s * s;
                double s_3 = s_2 * s;
                double s_4 = s_2 * s_2;
                double s_5 = s_4 * s;
                double ra = Math.Sqrt(a);


                double w = a + s * ra + (s * s - 1) / 3;
                w += (s_3 - 7 * s) / (36 * ra);
                w -= (3 * s_4 + 7 * s_2 - 16) / (810 * a);
                w += (9 * s_5 + 256 * s_3 - 433 * s) / (38880 * a * ra);


                if ((a >= 500) && (Math.Abs(1 - w / a) < 1e-6)) {
                    result = w;
                    p_has_10_digits = true;
                } else if (p > 0.5) {
                    if (w < 3 * a) {
                        result = w;
                    } else {
                        double D = Math.Max(2.0, (a * (a - 1)));
                        double lg = Math2.Lgamma(a);
                        double lb = Math.Log(q) + lg;
                        if (lb < -D * 2.3) {
                            // DiDonato and Morris Eq 25:
                            double y = -lb;
                            double c1 = (a - 1) * Math.Log(y);
                            double c1_2 = c1 * c1;
                            double c1_3 = c1_2 * c1;
                            double c1_4 = c1_2 * c1_2;
                            double a_2 = a * a;
                            double a_3 = a_2 * a;

                            double c2 = (a - 1) * (1 + c1);
                            double c3 = (a - 1) * (-(c1_2 / 2) + (a - 2) * c1 + (3 * a - 5) / 2);
                            double c4 = (a - 1) * ((c1_3 / 3) - (3 * a - 5) * c1_2 / 2 + (a_2 - 6 * a + 7) * c1 + (11 * a_2 - 46 * a + 47) / 6);
                            double c5 = (a - 1) * (-(c1_4 / 4)
                                              + (11 * a - 17) * c1_3 / 6
                                              + (-3 * a_2 + 13 * a - 13) * c1_2
                                              + (2 * a_3 - 25 * a_2 + 72 * a - 61) * c1 / 2
                                              + (25 * a_3 - 195 * a_2 + 477 * a - 379) / 12);

                            double y_2 = y * y;
                            double y_3 = y_2 * y;
                            double y_4 = y_2 * y_2;
                            result = y + c1 + (c2 / y) + (c3 / y_2) + (c4 / y_3) + (c5 / y_4);
                        } else {
                            // DiDonato and Morris Eq 33:
                            double u = -lb + (a - 1) * Math.Log(w) - Math.Log(1 + (1 - a) / (1 + w));
                            result = -lb + (a - 1) * Math.Log(u) - Math.Log(1 + (1 - a) / (1 + u));
                        }
                    }
                } else {
                    double z = w;
                    double ap1 = a + 1;
                    double ap2 = a + 2;
                    if (w < 0.15 * ap1) {
                        // DiDonato and Morris Eq 35:
                        double v = Math.Log(p) + Math2.Lgamma(ap1);
                        double s_ = 1;
                        z = Math.Exp((v + w) / a);
                        s_ = Math2.Log1p(z / ap1 * (1 + z / ap2));
                        z = Math.Exp((v + z - s_) / a);
                        s_ = Math2.Log1p(z / ap1 * (1 + z / ap2));
                        z = Math.Exp((v + z - s_) / a);
                        s_ = Math2.Log1p(z / ap1 * (1 + z / ap2 * (1 + z / (a + 3))));
                        z = Math.Exp((v + z - s_) / a);
                    }

                    if ((z <= 0.01 * ap1) || (z > 0.7 * ap1)) {
                        result = z;
                        if (z <= 0.002 * ap1)
                            p_has_10_digits = true;
                    } else {
                        // DiDonato and Morris Eq 36:
                        double ls = Math.Log(DidonatoSN(a, z, 100, 1e-4));
                        double v = Math.Log(p) + Math2.Lgamma(ap1);
                        z = Math.Exp((v + z - ls) / a);
                        result = z * (1 - (a * Math.Log(z) - z - v + ls) / (a - z));

                    }
                }
            }
            return result;
        }


        // get the upper and lower limits on x for the GammaInv functions
        static (double upperLimit, double lowerLimit) GammaInvLimits(double a, double p, double q, bool getQ)
        {
            Debug.Assert(a > 0 && a != 1);
            Debug.Assert(p >= 0 && p <= 1);
            Debug.Assert(q >= 0 && q <= 1);

            // set the upper and lower limits
            double lower = 0;
            double upper = double.MaxValue;

            // Find the upper or lower limits using http://dlmf.nist.gov/8.10.E11
            // For Q, the limits are:
            //      limit1 = -ln(1-(1-q)^(1/a))
            //      limit2 = limit1 * (Gamma(1+a)^(1/a))
            // For P, the limits are
            //      limit1 = -ln(1-p^(1/a))
            //      limit2 = limit1 * (Gamma(1+a)^(1/a))

            double limit1;
            if (getQ) {

                // if (1-q)^(1/a) > 1/2 
                double y = Math2.Log1p(-q) / a;
                if (y > -Constants.Ln2)
                    limit1 = -Math.Log(-Math2.Expm1(y));
                else
                    limit1 = -Math2.Log1p(-Math.Exp(y));
            } else {
                //      double y = Math.Exp(Math.Log(p)/a);
                double y = Math.Pow(p, 1 / a);
                if (y > 0.5)
                    limit1 = -Math.Log(-Math2.Expm1(Math.Log(p) / a));
                else
                    limit1 = -Math2.Log1p(-y);
            }

            if (a < 1) {
                upper = limit1;
                lower = limit1 * Math.Exp(Math2.Log1p(Math2.Tgamma1pm1(a)) / a);

            } else if (a > 1) {
                lower = limit1;
                upper = limit1 * Math.Exp(Math2.Lgamma(1 + a) / a);
            }

            Debug.Assert(upper >= lower, "Upper < lower");

            // Our computations for upper and lower can be slightly off in the last few digits
            // which can be problematic when the solution is close to max or min
            // So, give ourselves some wiggle room
            const double factor = 1.125;

            return (lower / factor, upper * factor);

        }

        public static (double upperLimit, double lowerLimit) GammaPInvLimits(double a, double p)
        {
            return GammaInvLimits(a, p, 1 - p, false);
        }

        public static (double upperLimit, double lowerLimit) GammaQInvLimits(double a, double q)
        {
            return GammaInvLimits(a, 1 - q, q, true);
        }
    }

    public static partial class Math2 {

        /// <summary>
        /// Returns the value "x" such that: p == GammaP(a, x);
        /// </summary>
        /// <param name="a">Requires a > 0</param>
        /// <param name="p">Requires 0 ≤ p ≤ 1</param>
        public static double GammaPInv(double a, double p)
        {
            if ((!(a > 0) || double.IsInfinity(a))
            || (!(p >= 0 && p <= 1))) {
                Policies.ReportDomainError("GammaPInv(a: {0}, p: {1}): Requires finite a > 0; p in [0,1]", a, p);
                return double.NaN;
            }

            if (p == 1)
                return double.PositiveInfinity;
            if (p == 0)
                return 0;

            // Since: P(1,x) = p = 1.0 - Math.Exp(-x), therefore: x = ln(1-p)
            if (a == 1.0)
                return -Math2.Log1p(-p);
            // Since P(1/2, x) = p = erf(sqrt(x)), therefore x = ErfInv(p)^2
            if (a == 0.5)
                return Squared(Math2.ErfInv(p));

            const double cutoff = 0.8125; // chosen through experimentation
            if (p > cutoff)
                return Math2.GammaQInv(a, 1.0 - p);

            double guess = _IgammaInv.GammaPInvGuess(a, p, 1 - p, out bool has10Digits);

            // set the upper and lower limits
            var (lower, upper) = _IgammaInv.GammaPInvLimits(a, p);

            // there can be some numerical uncertainties in the limits,
            // particularly for denormalized values, so... adjust the limits
            if (upper == 0 && lower == 0)
                return 0.0;

            if (guess <= DoubleLimits.MinNormalValue)
                return guess;

            if (guess <= lower)
                lower = DoubleLimits.MinNormalValue;

            if (guess > upper)
                upper = guess;

            return _IgammaInv.SolveGivenAP(a, p, guess, lower, upper);
        }


        /// <summary>
        /// Returns the value "x" such that q == GammaQ(a, x);
        /// </summary>
        /// <param name="a">Requires a > 0</param>
        /// <param name="q">Requires 0 ≤ q ≤ 1</param>
        public static double GammaQInv(double a, double q)
        {
            if ((!(a > 0) || double.IsInfinity(a))
            || (!(q >= 0 && q <= 1))) {
                Policies.ReportDomainError("GammaQInv(a: {0}, q: {1}): Requires finite a > 0; q in [0,1]", a, q);
                return double.NaN;
            }

            if (q == 0)
                return double.PositiveInfinity;
            if (q == 1)
                return 0;

            // Q(1,x) = e^-x
            if (a == 1)
                return -Math.Log(q);
            // Q(1/2, x) = erfc(sqrt(x))
            if (a == 0.5)
                return Squared(Math2.ErfcInv(q));


            double guess = _IgammaInv.GammaPInvGuess(a, 1 - q, q, out bool has10Digits);
            var (lower, upper) = _IgammaInv.GammaQInvLimits(a, q);

            // there can be some numerical uncertainties in the limits,
            // particularly for denormalized values, so adjust the limits
            if (upper == 0 && lower == 0)
                return 0.0;
            if (guess <= DoubleLimits.MinNormalValue)
                return guess;
            if (guess <= lower)
                lower = DoubleLimits.MinNormalValue;
            if (guess > upper)
                upper = guess;


            return _IgammaInv.SolveGivenAQ(a, q, guess, lower, upper);
        }

    }


} // namespace




