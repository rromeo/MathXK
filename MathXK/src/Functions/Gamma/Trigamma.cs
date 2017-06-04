//  Copyright (c) 2013 Rocco Romeo
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

//  Ported to .NET based on the original:
//      Copyright (c) John Maddock 2006, Boost Software License v1.0

using System;

namespace MathXK
{

    public static partial class Math2
    {

        /// <summary>
        /// Returns the Trigamma function
        /// <para>ψ'(x) = d/dx(ψ(x)) = d''/dx(ln(Γ(x)))</para>
        /// </summary>
        /// <param name="x">Trigamma function argument</param>
        public static double Trigamma(double x)
        {
            if (double.IsNaN(x)) {
                Policies.ReportDomainError("Trigamma(x: {0}): NaN not allowed", x);
                return double.NaN;
            }
            if (double.IsInfinity(x)) {
                if (x < 0) {
                    Policies.ReportDomainError("Trigamma(x: {0}): Requires x != -Infinity", x);
                    return double.NaN;
                }
                return 0;
            }

            double result = 0;
            //
            // Check for negative arguments and use reflection:
            //
            if (x <= 0) {
                if (IsInteger(x)) {
                    Policies.ReportPoleError("Trigamma(x: {0}) Requires x is not an integer when x <= 0", x);
                    return double.NaN;
                }

                // Reflect:
                double z = 1 - x;
                double s = Math.Abs(x) < Math.Abs(z) ? Math2.SinPI(x) : Math2.SinPI(z);
                return -Trigamma(z) + (Math.PI * Math.PI) / (s * s);
            }

            if (x < 1) {
                result = 1 / (x * x);
                x += 1;
            }

            if (x <= 2) {

                // Max error in interpolated form: 3.736e-017
                const double Y_1_2 = 2.1093254089355469;

                const double p0 = -1.1093280605946045;
                const double p1 = -3.8310674472619321;
                const double p2 = -3.3703848401898283;
                const double p3 = 0.28080574467981213;
                const double p4 = 1.6638069578676164;
                const double p5 = 0.64468386819102836;

                const double q0 = 1.0;
                const double q1 = 3.4535389668541151;
                const double q2 = 4.5208926987851437;
                const double q3 = 2.7012734178351534;
                const double q4 = 0.64468798399785611;
                const double q5 = -0.20314516859987728e-6;

                double z = x;
                double P = p0 + z * (p1 + z * (p2 + z * (p3 + z * (p4 + z * p5))));
                double Q = q0 + z * (q1 + z * (q2 + z * (q3 + z * (q4 + z * q5))));

                result += (Y_1_2 + P / Q) / (x * x);

            } else if (x <= 4) {

                // Max error in interpolated form: 1.159e-017
                const double p0 = -0.13803835004508849e-7;
                const double p1 = 0.50000049158540261;
                const double p2 = 1.6077979838469348;
                const double p3 = 2.5645435828098254;
                const double p4 = 2.0534873203680393;
                const double p5 = 0.74566981111565923;

                const double q0 = 1.0;
                const double q1 = 2.8822787662376169;
                const double q2 = 4.1681660554090917;
                const double q3 = 2.7853527819234466;
                const double q4 = 0.74967671848044792;
                const double q5 = -0.00057069112416246805;

                double z = 1 / x;
                double P = p0 + z * (p1 + z * (p2 + z * (p3 + z * (p4 + z * p5))));
                double Q = q0 + z * (q1 + z * (q2 + z * (q3 + z * (q4 + z * q5))));

                result += (1 + P / Q) / x;

            } else {

                // Maximum Deviation Found:                     6.896e-018
                // Expected Error Term :                       -6.895e-018
                // Maximum Relative Change in Control Points :  8.497e-004
                const double p0 = 0.68947581948701249e-17;
                const double p1 = 0.49999999999998975;
                const double p2 = 1.0177274392923795;
                const double p3 = 2.498208511343429;
                const double p4 = 2.1921221359427595;
                const double p5 = 1.5897035272532764;
                const double p6 = 0.40154388356961734;

                const double q0 = 1.0;
                const double q1 = 1.7021215452463932;
                const double q2 = 4.4290431747556469;
                const double q3 = 2.9745631894384922;
                const double q4 = 2.3013614809773616;
                const double q5 = 0.28360399799075752;
                const double q6 = 0.022892987908906897;

                double z = 1 / x;
                double P = p0 + z * (p1 + z * (p2 + z * (p3 + z * (p4 + z * (p5 + z * p6)))));
                double Q = q0 + z * (q1 + z * (q2 + z * (q3 + z * (q4 + z * (q5 + z * q6)))));

                result += (1 + P / Q) / x;
            }


            return result;
        }
    }


}