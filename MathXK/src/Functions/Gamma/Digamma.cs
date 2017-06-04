//  Copyright (c) 2013 Rocco Romeo
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

//  Ported to .NET based on the original:
//      Copyright (c) John Maddock 2006, Boost Software License v1.0

using System;

namespace MathXK
{

    static class _Digamma
    {
        /// <summary>
        /// Digamma rational approximation for x in [1,2]
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        public static double Rational_1_2(double x)
        {
            // digamma(x) = (x - root) * (Y + R(x-1))
            //
            // Where root is the location of the positive root of digamma,
            // Y is a constant, and R is optimised for low absolute error
            // compared to Y.
            //
            // Maximum Deviation Found:               1.466e-18
            // At double precision, max error found:  2.452e-17


            const double Y = 0.99558162689208984;

            const double Root1 = 1569415565.0 / 1073741824.0;
            const double Root2 = (381566830.0 / 1073741824.0) / 1073741824.0;
            const double Root3 = 0.9016312093258695918615325266959189453125e-19;
            
            const double p0 = 0.25479851061131551;
            const double p1 = -0.32555031186804491;
            const double p2 = -0.65031853770896507;
            const double p3 = -0.28919126444774784;
            const double p4 = -0.045251321448739056;
            const double p5 = -0.0020713321167745952;

            const double q0 = 1;
            const double q1 = 2.0767117023730469;
            const double q2 = 1.4606242909763515;
            const double q3 = 0.43593529692665969;
            const double q4 = 0.054151797245674225;
            const double q5 = 0.0021284987017821144;
            const double q6 = -0.55789841321675513e-6;

            double z = x - 1;
            double P = p0 + z * (p1 + z * (p2 + z * (p3 + z * (p4 + z * p5))));
            double Q = q0 + z * (q1 + z * (q2 + z * (q3 + z * (q4 + z * (q5 + z * q6)))));

            double g = ((x - Root1) - Root2) - Root3;
            double result = g * (Y + P / Q);
            return result;
        }


    };


    public static partial class Math2
    {

        /// <summary>
        /// Returns the Digamma function
        /// <para>ψ(x) = d/dx(ln(Γ(x))) = Γ'(x)/Γ(x)</para>
        /// </summary>
        /// <param name="x">Digamma function argument</param>
        public static double Digamma(double x)
        {
            if (double.IsNaN(x)) {
                Policies.ReportDomainError("Digamma(x: {0}): NaN not allowed", x);
                return double.NaN;
            }
            if (double.IsInfinity(x)) {
                if (x < 0) {
                    Policies.ReportDomainError("Digamma(x: {0}): Requires x != -Infinity", x);
                    return double.NaN;
                }
                return double.PositiveInfinity;
            }


            if (x <= 0) {
                if (IsInteger(x)) {
                    Policies.ReportPoleError("Digamma(x: {0}) Requires x is not an integer when x <= 0", x);
                    return double.NaN;
                }


                // use the following equations to reflect
                // ψ(1-x) - ψ(x) = π * cot(π*x) 
                // ψ(-x) = 1/x + ψ(x) +  π * cot(π*x)

                if (x < -1)
                    return Digamma(1 - x) - Math.PI * Math2.CotPI(x);

                // The following is the forward recurrence:
                //      -(1/x + 1/(x+1)) + Digamma(x + 2);
                // with a little less cancellation error near Digamma root: x = -0.50408...
                return _Digamma.Rational_1_2(x + 2) - (2 * x + 1) / (x * (x + 1));

            }

            // If we're above the lower-limit for the asymptotic expansion then use it:
            if (x >= 10) {
                const double p0 = 0.083333333333333333333333333333333333333333333333333;
                const double p1 = -0.0083333333333333333333333333333333333333333333333333;
                const double p2 = 0.003968253968253968253968253968253968253968253968254;
                const double p3 = -0.0041666666666666666666666666666666666666666666666667;
                const double p4 = 0.0075757575757575757575757575757575757575757575757576;
                const double p5 = -0.021092796092796092796092796092796092796092796092796;
                const double p6 = 0.083333333333333333333333333333333333333333333333333;
                const double p7 = -0.44325980392156862745098039215686274509803921568627;

                double xm1 = x - 1;
                double z = 1 / (xm1 * xm1);
                double P = p0 + z * (p1 + z * (p2 + z * (p3 + z * (p4 + z * (p5 + z * (p6 + z * p7))))));
                return Math.Log(xm1) + 1.0 / (2 * xm1) - z * P;
            }


            double result = 0;

            //
            // If x > 2 reduce to the interval [1,2]:
            //
            while (x > 2) {
                x -= 1;
                result += 1 / x;
            }
            //
            // If x < 1 use recurrence to shift to > 1:
            //
            if (x < 1) {
                result = -1 / x;
                x += 1;
            }
            result += _Digamma.Rational_1_2(x);

            return result;
        }

    }

} // namespaces

