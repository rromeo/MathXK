//  Copyright (c) 2013 Rocco Romeo
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

//  Ported to .NET based on the original:
//      Copyright (c) 2006 Xiaogang Zhang, Boost Software License v1.0
//  History:
//      Source originally by XZ. 
//      RR restructured code and made extensive changes to improve accuracy.

#if DEBUG
// uncomment to get more detailed debug messages
//#define EXTRA_DEBUG

#endif

using System;
using System.Diagnostics;
using System.Runtime.InteropServices;
using MathXK.Numerics;

namespace MathXK
{

    static partial class _Bessel
    {


        /// <summary>
        /// Computes J{v}(x), Y{v}(x) using a combination of asymptotic approximation and recurrence
        /// </summary>
        /// <param name="v"></param>
        /// <param name="x"></param>
        /// <param name="needJ"></param>
        /// <param name="needY"></param>
        /// <returns></returns>
        static (double J, DoubleX Y) JY_AsymRecurrence(double v, double x, bool needJ, bool needY)
        {
            Debug.Assert(v >= 0 && x >= 1, "Requires positive values for v,x");
            Debug.Assert(v < int.MaxValue, "v too large: v = " + v);

            var J = double.NaN;
            var Y = DoubleX.NaN;

            int nPos = (int)Math.Floor(v);
            double uPos = v - nPos;

            // Using Hankel, find:
            // J{v-Floor(v)}(x) and J{v-Floor(v)+1}(x) 
            // Y{v-Floor(v)}(x) and Y{v-Floor(v)+1}(x) 
            // then use recurrence to find J{v}(x), Y{v}(x)

            double u0, u1;;       
            int n;

            if (x >= 9) {

                // set the start of the recurrence near sqrt(x)
                double maxV = Math.Floor(Math.Sqrt(x));
                u1 = (maxV - 1) + uPos;
                u0 = u1 - 1;
                n = (int)Math.Floor(v - u1 + 0.5);
                Debug.Assert(n >= 0);
            } else {
                u0 = uPos;
                u1 = uPos + 1;
                n = nPos - 1;
            }

            Debug.Assert(x >= HankelAsym.JYMinX(u1), "x is too small for HankelAsym");

            var (Ju, Yu) = HankelAsym.JY(u0, x);
            var (Jup1, Yup1) = HankelAsym.JY(u1, x);

            if (needJ) {
                if (v < x) {
                    J = Recurrence.ForwardJY(u1, x, n, Jup1, Ju).JYvpn;
                } else {
                    // Use fv = J{v+1}(x) / J{v}(x)
                    // and backward recurrence to find (J{v-n+1}/J{v-n})
                    var (fv, s) = J_CF1(v, x);
                    var (Jvmn, Jvmnp1, scale) = Recurrence.BackwardJY_B(v, x, n, s, fv * s);

                    var Jv = Math2.Ldexp(Jup1 / Jvmn, -scale);

                    J = s * Jv;      // normalization
                }
            }

            if (needY) {
                var (Yv, Yvm1, YScale) = Recurrence.ForwardJY_B(u1, x, n, Yup1, Yu);
                Y = DoubleX.Ldexp(Yv, YScale);
            }

            return (J, Y);
        }


        /// <summary>
        /// Returns (Y{v}(x), Y{v+1}(x)) by Temme's method
        /// </summary>
        /// <param name="v"></param>
        /// <param name="x"></param>
        /// <returns></returns>
        static (double Yv, double Yvp1) Y_Temme(double v, double x)
        {
            //  see Temme, Journal of Computational Physics, vol 21, 343 (1976)

            Debug.Assert(Math.Abs(v) <= 0.5);  // precondition for using this routine

            double gp = Math2.Tgamma1pm1(v);
            double gm = Math2.Tgamma1pm1(-v);
            double spv = Math2.SinPI(v);
            double spv2 = Math2.SinPI(v / 2);
            double xp = Math.Pow(x / 2, v);

            double a = Math.Log(x / 2);
            double sigma = -a * v;
            double d = Math2.Sinhc(sigma);
            double e = Math.Abs(v) < DoubleLimits.MachineEpsilon ? ((Math.PI * Math.PI / 2) * v) : (2 * spv2 * spv2 / v);

            double g1 = (v == 0) ? -Constants.EulerMascheroni : ((gp - gm) / ((1 + gp) * (1 + gm) * 2 * v));
            double g2 = (2 + gp + gm) / ((1 + gp) * (1 + gm) * 2);
            double vspv = (Math.Abs(v) < DoubleLimits.MachineEpsilon) ? 1 / Math.PI : v / spv;
            double f = (g1 * Math.Cosh(sigma) - g2 * a * d) * 2 * vspv;

            double p = vspv / (xp * (1 + gm));
            double q = vspv * xp / (1 + gp);

            double g = f + e * q;
            double h = p;
            double coef = 1;
            double sum = coef * g;
            double sum1 = coef * h;

            double v2 = v * v;
            double coef_mult = -x * x / 4;


            int k = 1;
            for (; k < Policies.MaxSeriesIterations; k++) {
                f = (k * f + p + q) / (k * k - v2);
                p /= k - v;
                q /= k + v;
                g = f + e * q;
                h = p - k * g;
                coef *= coef_mult / k;
                sum += coef * g;
                sum1 += coef * h;
                if (Math.Abs(coef * g) < Math.Abs(sum) * Policies.SeriesTolerance)
                    return (-sum, -2 * sum1 / x);
            }

            Policies.ReportConvergenceError("Series did not converge after {0} iterations", Policies.MaxSeriesIterations);
            return (double.NaN, double.NaN);
        }

        /// <summary>
        /// Evaluate fv = J{v+1}(x)/J{v}(x) using continued fractions
        /// <para>|x| &lt;= |v|, function converges rapidly</para>
        /// <para>|x| &gt; |v|, function needs O(|x|) iterations to converge</para>
        /// </summary>
        /// <param name="v"></param>
        /// <param name="x"></param>
        /// <returns></returns>
        static (double JRatio, int sign) J_CF1(double v, double x)
        {
            //  See Abramowitz and Stegun, Handbook of Mathematical Functions, 1972, 9.1.73

            // note that the interations required to converge are approximately order(x)
            const int MaxIterations = 100000;
            const double tolerance = 2 * DoubleLimits.MachineEpsilon;
            double tiny = Math.Sqrt(DoubleLimits.MinNormalValue);

            // modified Lentz's method, see
            // Lentz, Applied Optics, vol 15, 668 (1976)

            int s = 1; // sign of denominator
            double f;
            double C = f = tiny;
            double D = 0;

            int k;
            for (k = 1; k < MaxIterations; k++) {

                double a = -1;
                double b = 2 * (v + k) / x;

                C = b + a / C;
                D = b + a * D;
                if (C == 0) { C = tiny; }
                if (D == 0) { D = tiny; }
                D = 1 / D;
                double delta = C * D;
                f *= delta;
                if (D < 0) { s = -s; }
                if (Math.Abs(delta - 1) < tolerance)
                    return (-f, s);
            }

            Policies.ReportConvergenceError("Series did not converge after {0} iterations", Policies.MaxSeriesIterations);
            return (double.NaN, 0);

        }

        /// <summary>
        /// p+iq = (J'{v}+iY'{v})/(J{v}+iY{v})
        /// <para>|x| >= |v|, function converges rapidly</para>
        /// <para>|x| -> 0, function fails to converge</para>
        /// </summary>
        /// <param name="v"></param>
        /// <param name="x"></param>
        /// <returns></returns>
        static (double P, double Q) JY_CF2(double v, double x)
        {
            //
            // This algorithm was originally written by Xiaogang Zhang
            // using complex numbers to perform the complex arithmetic.
            // However, that turns out to 10x or more slower than using
            // all real-valued arithmetic, so it's been rewritten using
            // real values only.
            //

            Debug.Assert(Math.Abs(x) > 1);

            // modified Lentz's method, complex numbers involved, see
            // Lentz, Applied Optics, vol 15, 668 (1976)
            const double tolerance = 2 * DoubleLimits.MachineEpsilon;
            const double tiny = DoubleLimits.MinNormalValue;
            double Cr = -0.5 / x, fr = Cr;
            double Ci = 1, fi = Ci;
            
            double v2 = v * v;
            double a = (0.25 - v2) / x;
            double br = 2 * x;
            double bi = 2;
            double temp = Cr * Cr + 1;
            Ci = bi + a * Cr / temp;
            Cr = br + a / temp;
            double Dr = br;
            double Di = bi;
            if (Math.Abs(Cr) + Math.Abs(Ci) < tiny) { Cr = tiny; }
            if (Math.Abs(Dr) + Math.Abs(Di) < tiny) { Dr = tiny; }
            temp = Dr * Dr + Di * Di;
            Dr = Dr / temp;
            Di = -Di / temp;
            double delta_r = Cr * Dr - Ci * Di;
            double delta_i = Ci * Dr + Cr * Di;
            temp = fr;
            fr = temp * delta_r - fi * delta_i;
            fi = temp * delta_i + fi * delta_r;
            int k;
            for (k = 2; k < Policies.MaxSeriesIterations; k++) {
                a = k - 0.5;
                a *= a;
                a -= v2;
                bi += 2;
                temp = Cr * Cr + Ci * Ci;
                Cr = br + a * Cr / temp;
                Ci = bi - a * Ci / temp;
                Dr = br + a * Dr;
                Di = bi + a * Di;
                if (Math.Abs(Cr) + Math.Abs(Ci) < tiny) { Cr = tiny; }
                if (Math.Abs(Dr) + Math.Abs(Di) < tiny) { Dr = tiny; }
                temp = Dr * Dr + Di * Di;
                Dr = Dr / temp;
                Di = -Di / temp;
                delta_r = Cr * Dr - Ci * Di;
                delta_i = Ci * Dr + Cr * Di;
                temp = fr;
                fr = temp * delta_r - fi * delta_i;
                fi = temp * delta_i + fi * delta_r;
                if (Math.Abs(delta_r - 1) + Math.Abs(delta_i) < tolerance)
                    break;

            }

            if (k >= Policies.MaxSeriesIterations) {
                Policies.ReportConvergenceError("JY_CF2(v:{0}, x:{1}) did not converge after {2} iterations", v, x, Policies.MaxSeriesIterations);
                return (double.NaN, double.NaN);
            }

            return (fr, fi);
        }

        /// <summary>
        /// Compute (J, Y) using Steeds method
        /// </summary>
        /// <param name="v"></param>
        /// <param name="x"></param>
        /// <param name="needJ"></param>
        /// <param name="needY"></param>
        /// <returns></returns>
        static (double J, DoubleX Y) JY_Steed(double v, double x, bool needJ, bool needY)
        {
            Debug.Assert(x > 0 && v >= 0);
            Debug.Assert(needJ || needY);

            double J = double.NaN;
            DoubleX Y = DoubleX.NaN;

            int n = (int)Math.Floor(v + 0.5);
            double u = v - n;                              
            Debug.Assert(u >= -0.5 && u < 0.5); // Ensure u in [-1/2, 1/2)

            // compute J{v+1}(x) / J{v}(x)
            var (fv, s) = J_CF1(v, x);
            var (Jvmn, Jvmnp1, scale) = Recurrence.BackwardJY_B(v, x, n, s, fv * s);

            double ratio = s / Jvmn;       // normalization
            double fu = Jvmnp1 / Jvmn;

            // prev/current, can also call CF1_jy() to get fu, not much difference in precision
            //double fuCF;
            //int sfu;
            //Bessel.CF1_jy(u,x,fuCF,sfu);
            //fu = fuCF;

            var (p, q) = JY_CF2(u, x);          // continued fraction JY_CF2
            double t = u / x - fu;              // t = J'/J
            double gamma = (p - t) / q;

            //
            // We can't allow gamma to cancel out to zero competely as it messes up
            // the subsequent logic.  So pretend that one bit didn't cancel out
            // and set to a suitably small value.  The only test case we've been able to
            // find for this, is when v = 8.5 and x = 4*PI.
            //
            if (gamma == 0)
                gamma = u * DoubleLimits.MachineEpsilon / x;

            double W = (2 / Math.PI) / x;               // Wronskian
            double Ju = Math.Sign(Jvmn) * Math.Sqrt(W / (q + gamma * (p - t)));

            double Jv = Ju * ratio;                    // normalization
            J = Math2.Ldexp(Jv, -scale);

            if (needY) {
                double Yu = gamma * Ju;
                double Yu1 = Yu * (u / x - p - q / gamma);

                var (JYvpn, JYvpnm1, YScale) = Recurrence.ForwardJY_B(u + 1, x, n, Yu1, Yu);
                Y = DoubleX.Ldexp(JYvpnm1, YScale);
            }

            return (J, Y);


        }
        
        /// <summary>
        /// Try to use asymptotics to generate the result
        /// </summary>
        /// <param name="v"></param>
        /// <param name="x"></param>
        /// <param name="J"></param>
        /// <param name="Y"></param>
        /// <param name="needJ"></param>
        /// <param name="needY"></param>
        /// <returns></returns>
        public static bool JY_TryAsymptotics(double v, double x, out double J, out DoubleX Y, bool needJ, bool needY)
        {
            J = double.NaN;
            Y = DoubleX.NaN;

            // Try asymptotics directly:

            if (x > v) {

                // x > v*v
                if ( x >= HankelAsym.JYMinX(v) ) {
                    var result = HankelAsym.JY(v, x);
                    if (needJ)
                        J = result.J;
                    if (needY)
                        Y = result.Y;
                    return true;
                }

                // Try Asymptotic Phase for x > 47v
                if (x >= MagnitudePhase.MinX(v)) {
                    var result = MagnitudePhase.BesselJY(v, x);
                    if (needJ)
                        J = result.J;
                    if (needY)
                        Y = result.Y;
                    return true;
                }
            }

            if (UniformAsym.IsJYPrecise(v, x)) {
                var (Jv, Yv) = UniformAsym.JY(v, x);
                if (needJ)
                    J = Jv;
                if (needY)
                    Y = Yv;
                return true;
            }

            // Try asymptotics with recurrence
            if (x > v && x >= HankelAsym.JYMinX(v - Math.Floor(v) + 1)) {
                var (Jv, Yv) = JY_AsymRecurrence(v, x, needJ, needY);
                if (needJ)
                    J = Jv;
                if (needY)
                    Y = Yv;
                return true;
            }

            return false;

        }


        /// <summary>
        /// Compute J{v}(x) and Y{v}(x)
        /// </summary>
        /// <param name="v"></param>
        /// <param name="x"></param>
        /// <param name="needJ"></param>
        /// <param name="needY"></param>
        /// <returns></returns>
        public static (double J, DoubleX Y) JY(double v, double x, bool needJ, bool needY)
        {
            Debug.Assert(needJ || needY);

            // uses  Steed's method
            // see Barnett et al, Computer Physics Communications, vol 8, 377 (1974)

            // set [out] parameters
            double J = double.NaN;
            DoubleX Y = DoubleX.NaN;

            // check v first so that there are no integer throws later
            if (Math.Abs(v) > int.MaxValue) {
                Policies.ReportNotImplementedError("BesselJY(v: {0}): Large |v| > int.MaxValue not yet implemented", v);
                return (J, Y);
            }


            if (v < 0) {

                v = -v;

                if (Math2.IsInteger(v)) {

                    // for integer orders only, we can use the following identities:
                    //      J{-n}(x) = (-1)^n * J{n}(x) 
                    //      Y{-n}(x) = (-1)^n * Y{n}(x) 

                    if (Math2.IsOdd(v)) {

                        var (JPos, YPos) = JY(v, x, needJ, needY);

                        if (needJ)
                            J = -JPos;
                        if (needY)
                            Y = -YPos;

                        return (J, Y);
                    }


                } if (v - Math.Floor(v) == 0.5) {

                    Debug.Assert(v >= 0);

                    // use reflection rule:
                    // for integer m >= 0 
                    // J{-(m+1/2)}(x) = (-1)^(m+1) * Y{m+1/2}(x)
                    // Y{-(m+1/2)}(x) = (-1)^m * J{m+1/2}(x)

                    // call the general bessel functions with needJ and needY reversed
                    var (JPos, YPos) = JY(v, x, needY, needJ);

                    double m = v - 0.5;
                    bool isOdd = Math2.IsOdd(m);

                    if (needJ) {
                        double y = (double)YPos;
                        J = isOdd ? y : -y;
                    }
                    if (needY)
                        Y = isOdd ? -JPos : JPos;

                    return (J, Y);

                } else {

                    // use reflection rule:
                    // J{-v}(x) = cos(pi*v)*J{v}(x) - sin(pi*v)*Y{v}(x)
                    // Y{-v}(x) = sin(pi*v)*J{v}(x) + cos(pi*v)*Y{v}(x)

                    var (JPos, YPos) = JY(v, x, true, true);

                    double cp = Math2.CosPI(v);
                    double sp = Math2.SinPI(v);

                    J = cp * JPos - (double)(sp * YPos);
                    Y = sp * JPos + cp * YPos;

                    return (J, Y);
                }
            }

            // both x and v are positive from here
            Debug.Assert(x >= 0 && v >= 0);

            if (x == 0) {
                // For v > 0 
                if (needJ)
                    J = 0;
                if (needY)
                    Y = DoubleX.NegativeInfinity;
                return (J, Y);
            }

            int n = (int)Math.Floor(v + 0.5);
            double u = v - n;                              // -1/2 <= u < 1/2

            // is it an integer?
            if (u == 0) {
                if (v == 0) {
                    if (needJ)
                        J = J0(x);
                    if (needY)
                        Y = Y0(x);
                    return (J, Y);
                }

                if (v == 1) {
                    if (needJ)
                        J = J1(x);
                    if (needY)
                        Y = Y1(x);
                    return (J, Y);
                }

                // for integer order only
                if (needY && x < DoubleLimits.RootMachineEpsilon._2) {
                    Y = YN_SmallArg(n, x);
                    if (!needJ)
                        return (J, Y);
                    needY = !needY;
                }
            }


            if (needJ && ((x < 5) || (v > x * x / 4))) {
                // always use the J series if we can
                J = J_SmallArg(v, x);
                if (!needY)
                    return (J, Y);
                needJ = !needJ;
            }

            if (needY && x <= 2) {
                // J should have already been solved above
                Debug.Assert(!needJ);

                // Evaluate using series representations.
                // Much quicker than Y_Temme below.
                // This is particularly important for x << v as in this
                // area Y_Temme may be slow to converge, if it converges at all.

                // for non-integer order only
                if ( u != 0 ) {
                    if ((x < 1) && (Math.Log(DoubleLimits.MachineEpsilon / 2) > v * Math.Log((x / 2) * (x / 2) / v))) {
                        Y = Y_SmallArgSeries(v, x);
                        return (J, Y);
                    }
                }

                // Use Temme to find Yu where |u| <= 1/2, then use forward recurrence for Yv

                var (Yu, Yu1) = Y_Temme(u, x);
                var (Yvpn, Yvpnm1, YScale) = Recurrence.ForwardJY_B(u + 1, x, n, Yu1, Yu);
                Y = DoubleX.Ldexp(Yvpnm1, YScale);
                return (J, Y);
            }

            Debug.Assert(x > 2 && v >= 0);

            // Try asymptotics directly:

            if (x > v) {

                // x > v*v
                if (x >= HankelAsym.JYMinX(v)) {
                    var result = HankelAsym.JY(v, x);
                    if (needJ)
                        J = result.J;
                    if (needY)
                        Y = result.Y;
                    return (J, Y);
                }

                // Try Asymptotic Phase for x > 47v
                if (x >= MagnitudePhase.MinX(v)) {
                    var result = MagnitudePhase.BesselJY(v, x);
                    if (needJ)
                        J = result.J;
                    if (needY)
                        Y = result.Y;
                    return (J, Y);
                }
            }

            // fast and accurate within a limited range of v ~= x
            if (UniformAsym.IsJYPrecise(v, x)) {
                var (Jv, Yv) = UniformAsym.JY(v, x);
                if (needJ)
                    J = Jv;
                if (needY)
                    Y = Yv;
                return (J, Y);
            }

            // Try asymptotics with recurrence:
            if (x > v && x >= HankelAsym.JYMinX(v - Math.Floor(v) + 1)) {
                var (Jv, Yv) = JY_AsymRecurrence(v, x, needJ, needY);
                if (needJ)
                    J = Jv;
                if (needY)
                    Y = Yv;
                return (J, Y);
            }
            
            // Use Steed's Method
            var (SteedJv, SteedYv) = JY_Steed(v, x, needJ, needY);
            if (needJ)
                J = SteedJv;
            if (needY)
                Y = SteedYv;
            return (J, Y);

        }
    }

} // namespaces


