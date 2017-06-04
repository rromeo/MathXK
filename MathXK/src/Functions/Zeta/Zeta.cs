//  Copyright (c) 2013 Rocco Romeo
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

//  Ported to .NET based on the original:
//      Copyright (C) John Maddock 2007, Boost Software License, Version 1.0


using System;
using System.Diagnostics;

namespace MathXK
{

    internal static class _Zeta
    {

        // sc = 1-s
        public static double Imp(double s, double sc)
        {
            Debug.Assert(s != 1);

            // Use integer routines if we can. 
            if (Math.Floor(s) == s && s > int.MinValue && s <= int.MaxValue)
                return Imp((int)s);

            // Zeta(x) approx = -0.5 - 1/2 * Log(2*PI) * x
            // Expected error should be approx 3 ulps
            const double LowerLimit = DoubleLimits.RootMachineEpsilon._2;
            if (Math.Abs(s) < LowerLimit)
                return -0.5 - (Constants.Ln2PI / 2) * s;

            if (s < 0) {
                // The zeta function = 0 for -2, -4, -6...-2n
                if (Math.Floor(s / 2) == s / 2)
                    return 0;

                // http://functions.wolfram.com/ZetaFunctionsandPolylogarithms/Zeta/17/01/01/0001/
                // Use reflection: ζ(s) = Γ(1-s) * (2^s) * π^(s-1) * Sin(sπ/2) * ζ(1-s)

                if (sc <= DoubleLimits.MaxGamma )
                    return 2 * Math.Pow(2 * Math.PI, -sc) * Math2.Tgamma(sc) * Math2.SinPI(0.5 * s) * Imp(sc, s);

                // use Stirling
                double sin = Math2.SinPI(0.5 * s);
                double sign = Math.Sign(sin);
                sin = Math.Abs(sin);

                double value;
                var factor = StirlingGamma.GammaSeries(sc) * Imp(sc, s);
                var log = (sc - 0.5) * Math.Log(sc / (2 * Math.PI));
                if (log < DoubleLimits.MaxLogValue) {
                    value = factor * Math.Pow(sc / (2 * Math.PI), sc - 0.5) * Math.Exp(-sc) * (2 * sin);
                } else if (log < 2 * DoubleLimits.MaxLogValue) {
                    double e = Math.Pow(sc / (2 * Math.PI), sc / 2 - 0.25) * Math.Exp(-sc / 2);
                    value = factor * e * (sin * 2) * e;
                } else {
                    value = Math.Exp(log - sc + Math.Log(2 * sin * factor));
                }
                return sign * value;

            }

            double result;

            if (s < 1) {
                // Rational Approximation
                // Maximum Deviation Found:                     2.020e-18
                // Expected Error Term:                         -2.020e-18
                // Max error found at double precision:         3.994987e-17

                const double Y = -1.2433929443359375;

                const double p0 = 0.24339294433593750202;
                const double p1 = -0.49092470516353571651;
                const double p2 = 0.0557616214776046784287;
                const double p3 = -0.00320912498879085894856;
                const double p4 = 0.000451534528645796438704;
                const double p5 = -0.933241270357061460782e-5;

                const double q0 = 1;
                const double q1 = -0.279960334310344432495;
                const double q2 = 0.0419676223309986037706;
                const double q3 = -0.00413421406552171059003;
                const double q4 = 0.00024978985622317935355;
                const double q5 = -0.101855788418564031874e-4;

                double z = sc;
                double P = p0 + z * (p1 + z * (p2 + z * (p3 + z * (p4 + z * p5))));
                double Q = q0 + z * (q1 + z * (q2 + z * (q3 + z * (q4 + z * q5))));


                result = ((P / Q) + Y + sc) / sc;

            } else if (s <= 2) {
                // Maximum Deviation Found:        9.007e-20
                // Expected Error Term:            9.007e-20

                const double p0 = 0.577215664901532860516;
                const double p1 = 0.243210646940107164097;
                const double p2 = 0.0417364673988216497593;
                const double p3 = 0.00390252087072843288378;
                const double p4 = 0.000249606367151877175456;
                const double p5 = 0.110108440976732897969e-4;

                const double q0 = 1;
                const double q1 = 0.295201277126631761737;
                const double q2 = 0.043460910607305495864;
                const double q3 = 0.00434930582085826330659;
                const double q4 = 0.000255784226140488490982;
                const double q5 = 0.10991819782396112081e-4;


                double z = -sc;
                double P = p0 + z * (p1 + z * (p2 + z * (p3 + z * (p4 + z * p5))));
                double Q = q0 + z * (q1 + z * (q2 + z * (q3 + z * (q4 + z * q5))));

                result = P / Q + 1 / (-sc);

            } else if (s <= 4) {
                // Maximum Deviation Found:          5.946e-22
                // Expected Error Term:              -5.946e-22

                const double _Y3 = 0.6986598968505859375;

                const double p0 = -0.0537258300023595030676;
                const double p1 = 0.0445163473292365591906;
                const double p2 = 0.0128677673534519952905;
                const double p3 = 0.00097541770457391752726;
                const double p4 = 0.769875101573654070925e-4;
                const double p5 = 0.328032510000383084155e-5;

                const double q0 = 1;
                const double q1 = 0.33383194553034051422;
                const double q2 = 0.0487798431291407621462;
                const double q3 = 0.00479039708573558490716;
                const double q4 = 0.000270776703956336357707;
                const double q5 = 0.106951867532057341359e-4;
                const double q6 = 0.236276623974978646399e-7;


                double z = s - 2;
                double P = p0 + z * (p1 + z * (p2 + z * (p3 + z * (p4 + z * p5))));
                double Q = q0 + z * (q1 + z * (q2 + z * (q3 + z * (q4 + z * (q5 + z * q6)))));


                result = P / Q + _Y3 + 1 / (-sc);

            } else if (s <= 7) {
                // Maximum Deviation Found:                     2.955e-17
                // Expected Error Term:                         2.955e-17
                // Max error found at double precision:         2.009135e-16


                const double p0 = -2.49710190602259410021;
                const double p1 = -2.60013301809475665334;
                const double p2 = -0.939260435377109939261;
                const double p3 = -0.138448617995741530935;
                const double p4 = -0.00701721240549802377623;
                const double p5 = -0.229257310594893932383e-4;

                const double q0 = 1;
                const double q1 = 0.706039025937745133628;
                const double q2 = 0.15739599649558626358;
                const double q3 = 0.0106117950976845084417;
                const double q4 = -0.36910273311764618902e-4;
                const double q5 = 0.493409563927590008943e-5;
                const double q6 = -0.234055487025287216506e-6;
                const double q7 = 0.718833729365459760664e-8;
                const double q8 = -0.1129200113474947419e-9;


                double z = s - 4;
                double P = p0 + z * (p1 + z * (p2 + z * (p3 + z * (p4 + z * p5))));
                double Q = q0 + z * (q1 + z * (q2 + z * (q3 + z * (q4 + z * (q5 + z * (q6 + z * (q7 + z * q8)))))));

                result = 1 + Math.Exp(P / Q);
            } else if (s < 15) {
                // Maximum Deviation Found:                     7.117e-16
                // Expected Error Term:                         7.117e-16
                // Max error found at double precision:         9.387771e-16

                const double p0 = -4.78558028495135619286;
                const double p1 = -1.89197364881972536382;
                const double p2 = -0.211407134874412820099;
                const double p3 = -0.000189204758260076688518;
                const double p4 = 0.00115140923889178742086;
                const double p5 = 0.639949204213164496988e-4;
                const double p6 = 0.139348932445324888343e-5;

                const double q0 = 1;
                const double q1 = 0.244345337378188557777;
                const double q2 = 0.00873370754492288653669;
                const double q3 = -0.00117592765334434471562;
                const double q4 = -0.743743682899933180415e-4;
                const double q5 = -0.21750464515767984778e-5;
                const double q6 = 0.471001264003076486547e-8;
                const double q7 = -0.833378440625385520576e-10;
                const double q8 = 0.699841545204845636531e-12;


                double z = s - 7;
                double P = p0 + z * (p1 + z * (p2 + z * (p3 + z * (p4 + z * (p5 + z * p6)))));
                double Q = q0 + z * (q1 + z * (q2 + z * (q3 + z * (q4 + z * (q5 + z * (q6 + z * (q7 + z * q8)))))));

                result = 1 + Math.Exp(P/Q);
            } else if (s < 36) {
                // Max error in interpolated form:             1.668e-17
                // Max error found at long double precision:   1.669714e-17


                const double p0 = -10.3948950573308896825;
                const double p1 = -2.85827219671106697179;
                const double p2 = -0.347728266539245787271;
                const double p3 = -0.0251156064655346341766;
                const double p4 = -0.00119459173416968685689;
                const double p5 = -0.382529323507967522614e-4;
                const double p6 = -0.785523633796723466968e-6;
                const double p7 = -0.821465709095465524192e-8;

                const double q0 = 1;
                const double q1 = 0.208196333572671890965;
                const double q2 = 0.0195687657317205033485;
                const double q3 = 0.00111079638102485921877;
                const double q4 = 0.408507746266039256231e-4;
                const double q5 = 0.955561123065693483991e-6;
                const double q6 = 0.118507153474022900583e-7;
                const double q7 = 0.222609483627352615142e-14;

                double z = s - 15;
                double P = p0 + z * (p1 + z * (p2 + z * (p3 + z * (p4 + z * (p5 + z * (p6 + z * p7))))));
                double Q = q0 + z * (q1 + z * (q2 + z * (q3 + z * (q4 + z * (q5 + z * (q6 + z * q7))))));

                result = 1 + Math.Exp(P/Q);

            } else if (s < 53) {

                result = 1 + Math.Pow(2.0, -s);

            } else {

                result = 1;

            }
            return result;


        }

        internal static double Imp(int n)
        {
            // For n >= 0 this routine uses a table lookup for speed
            //
            // For n < 0, the following routines are used:
            //
            // http://functions.wolfram.com/ZetaFunctionsandPolylogarithms/Zeta/03/01/0001/
            // ζ(-n) = (-1)^n/(n+1) * B{n+1}, where B is the Bernoulli number
            //
            // http://functions.wolfram.com/ZetaFunctionsandPolylogarithms/Zeta/17/01/01/0001/
            // ζ(n) = Γ(1-n) * (2^n) * π^(n-1) * Sin(nπ/2) * ζ(1-n)

            if ( n >= 0 ) {
                if (n < Math2.ZetaTable.Length)
                    return Math2.ZetaTable[n];
                return 1.0;
            }

            // ζ(-2n) = 0
            if (!Math2.IsOdd(n))
                return 0.0;

            // underflow
            if (n < -DoubleLimits.MaxLogValue)
                return 0;

            int sc = 1 - n;
            if (sc < 2 * Math2.Bernoulli2nTable.Length)
                return -Math2.Bernoulli2nTable[sc / 2] / sc;

            double sign = Math2.SinPI(0.5 * n);


            // This won't be called because of the current length of the bernoulli table
            // if (sc <= Math2.MaxFactorialIndex)
            //    return sign * 2 * Math.Pow(2 * Math.PI, -sc) * Math2.FactorialTable[sc-1] * Imp(sc);

            // use duplication
            // Γ(2n) = 2^(2n-1)/Sqrt(Pi) * Γ(n) * Γ(n+1/2) 
            if (sc <= 2 * (Math2.MaxFactorialIndex - 1))
                return sign * Math.Pow(Math.PI, -sc - 0.5) * Math2.Tgamma(sc/2.0) * Math2.Tgamma(sc/2.0 + 0.5) * Imp(sc);

            // We're going to overflow
            return sign * Math.Exp(Math2.Lgamma(sc) - sc*Constants.Ln2PI + Constants.Ln2 + Math.Log(Imp(sc)));
        }
    };


    public partial class Math2
    {
        /// <summary>
        /// Returns the value of the Riemann Zeta function at <paramref name="n"/>
        /// <para>ζ(n) = Σ k^-n, k={1, ∞}</para>
        /// </summary>
        /// <param name="n">The function argument</param>
        public static double Zeta(int n)
        {
            if (n == 1) {
                Policies.ReportPoleError("Zeta(n: {0}): Evaluation at pole", n);
                return double.NaN;
            }

            return _Zeta.Imp(n);

        }

        /// <summary>
        /// Returns the value of the Riemann Zeta function at <paramref name="x"/>
        /// <para>ζ(x) = Σ n^-x, n={1, ∞}</para>
        /// </summary>
        /// <param name="x">The function argument</param>
        public static double Zeta(double x)
        {
            if (double.IsNaN(x)) {
                Policies.ReportDomainError("Zeta(x: {0}): NaN not allowed", x);
                return double.NaN;
            }
            if (double.IsInfinity(x)) {
                if (x > 0)
                    return 1;

                Policies.ReportDomainError("Zeta(x: {0}): -Infinity not allowed", x);
                return double.NaN;
            }
            if (x == 1) {
                Policies.ReportPoleError("Zeta(x: {0}): Evaluation at pole", x);
                return double.NaN;
            }

            return _Zeta.Imp(x, 1 - x);

        }
    }

} // namespaces



