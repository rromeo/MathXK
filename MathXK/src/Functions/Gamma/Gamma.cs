//  Copyright (c) 2013 Rocco Romeo
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

//  Ported to .NET based on the original:
//      Copyright (c) John Maddock 2006-7, Boost Software License v1.0
//      Copyright (c) Paul A. Bristow 2007, Boost Software License v1.0

using System;
using System.Diagnostics;

namespace MathXK
{

    static class _Lgamma
    {

        /// <summary>
        /// Returns a rational approximation to Lgamma in the range [1, 3]
        /// </summary>
        /// <param name="x">z. Requires z &gt; 0</param>
        /// <param name="xm1">z-1</param>
        /// <param name="xm2">z-2</param>
        /// <returns></returns>
        public static double Rational_1_3(double x, double xm1, double xm2)
        {
            Debug.Assert(x >= 1 && x <= 3);

            // This version uses rational approximations for small
            // values of z accurate enough for 64-bit mantissas
            // (80-bit long doubles), works well for 53-bit doubles as well.

            double result = 0;
            if ((xm1 == 0) || (xm2 == 0)) {
                // nothing to do, result is zero....

            } else if (x <= 1.5) { // [1, 1.5]

                // For z in [1,1.5] use the following form:
                //
                // lgamma(z) = (z-1)(z-2)(Y + R(z-1))
                //
                // where R(z-1) is a rational approximation optimised for
                // low absolute error - as long as it's absolute error
                // is small compared to the constant Y - then any rounding
                // error in it's computation will get wiped out.
                //
                // R(z-1) has the following properties:
                //
                // At double precision: Max error found:                1.230011e-17
                // At 80-bit long double precision:   Max error found:  5.631355e-21
                // Maximum Deviation Found:                             3.139e-021
                // Expected Error Term:                                 3.139e-021

                //
                const double Y = 0.52815341949462890625;

                const double p0 = 0.490622454069039543534e-1;
                const double p1 = -0.969117530159521214579e-1;
                const double p2 = -0.414983358359495381969e0;
                const double p3 = -0.406567124211938417342e0;
                const double p4 = -0.158413586390692192217e0;
                const double p5 = -0.240149820648571559892e-1;
                const double p6 = -0.100346687696279557415e-2;

                const double q0 = 0.1e1;
                const double q1 = 0.302349829846463038743e1;
                const double q2 = 0.348739585360723852576e1;
                const double q3 = 0.191415588274426679201e1;
                const double q4 = 0.507137738614363510846e0;
                const double q5 = 0.577039722690451849648e-1;
                const double q6 = 0.195768102601107189171e-2;

                double z = xm1;
                double P = p0 + z * (p1 + z * (p2 + z * (p3 + z * (p4 + z * (p5 + z * p6)))));
                double Q = q0 + z * (q1 + z * (q2 + z * (q3 + z * (q4 + z * (q5 + z * q6)))));

                double prefix = xm1 * xm2;
                result = prefix * (Y + (P/Q));

            } else if (x <= 2) { // [1.5, 2]

                // For z in [1.5,2] use the following form:
                //
                // lgamma(z) = (2-z)(1-z)(Y + R(2-z))
                //
                // where R(2-z) is a rational approximation optimised for
                // low absolute error - as long as it's absolute error
                // is small compared to the constant Y - then any rounding
                // error in it's computation will get wiped out.
                //
                // R(2-z) has the following properties:
                //
                // At double precision, max error found:              1.797565e-17
                // At 80-bit long double precision, max error found:  9.306419e-21
                // Maximum Deviation Found:                           2.151e-021
                // Expected Error Term:                               2.150e-021
                //
                const double Y = 0.452017307281494140625;

                const double p0 = -0.292329721830270012337e-1;
                const double p1 = 0.144216267757192309184e0;
                const double p2 = -0.142440390738631274135e0;
                const double p3 = 0.542809694055053558157e-1;
                const double p4 = -0.850535976868336437746e-2;
                const double p5 = 0.431171342679297331241e-3;

                const double q0 = 0.1e1;
                const double q1 = -0.150169356054485044494e1;
                const double q2 = 0.846973248876495016101e0;
                const double q3 = -0.220095151814995745555e0;
                const double q4 = 0.25582797155975869989e-1;
                const double q5 = -0.100666795539143372762e-2;
                const double q6 = -0.827193521891290553639e-6;

                double z = -xm2;
                double P = p0 + z * (p1 + z * (p2 + z * (p3 + z * (p4 + z * p5))));
                double Q = q0 + z * (q1 + z * (q2 + z * (q3 + z * (q4 + z * (q5 + z * q6)))));

                double r = xm2 * xm1;
                result = r * (Y + (P/Q));

            } else { // [2, 3]

                // For z in [2,3) use the following form:
                // Use the following form:
                //
                // lgamma(z) = (z-2)(z+1)(Y + R(z-2))
                //
                // where R(z-2) is a rational approximation optimised for
                // low absolute error - as long as it's absolute error
                // is small compared to the constant Y - then any rounding
                // error in it's computation will get wiped out.
                //
                // R(z-2) has the following properties:
                //
                // At double: Max error found:                    4.231e-18
                // At long double: Max error found:               1.987e-21
                // Maximum Deviation Found (approximation error): 5.900e-24
                //
                const double Y = 0.158963680267333984375e0;

                const double p0 = -0.180355685678449379109e-1;
                const double p1 = 0.25126649619989678683e-1;
                const double p2 = 0.494103151567532234274e-1;
                const double p3 = 0.172491608709613993966e-1;
                const double p4 = -0.259453563205438108893e-3;
                const double p5 = -0.541009869215204396339e-3;
                const double p6 = -0.324588649825948492091e-4;

                const double q0 = 0.1e1;
                const double q1 = 0.196202987197795200688e1;
                const double q2 = 0.148019669424231326694e1;
                const double q3 = 0.541391432071720958364e0;
                const double q4 = 0.988504251128010129477e-1;
                const double q5 = 0.82130967464889339326e-2;
                const double q6 = 0.224936291922115757597e-3;
                const double q7 = -0.223352763208617092964e-6;

                double z = xm2;
                double P = p0 + z * (p1 + z * (p2 + z * (p3 + z * (p4 + z * (p5 + z * p6)))));
                double Q = q0 + z * (q1 + z * (q2 + z * (q3 + z * (q4 + z * (q5 + z * (q6 + z * q7))))));

                double r = xm2 * (x + 1);
                result = r * (Y + (P/ Q));
            }

            return result;
        }

    }

    public static partial class Math2
    {

        /// <summary>
        /// Returns the Gamma Function = Γ(x)
        /// </summary>
        /// <param name="x">Tgamma function argument</param>
        public static double Tgamma(double x)
        {
            if (double.IsNaN(x)) {
                Policies.ReportDomainError("Tgamma(x: {0}): NaN not allowed", x);
                return double.NaN;
            }
            if (double.IsInfinity(x)) {
                if (x < 0) {
                    Policies.ReportDomainError("Tgamma(x: {0}): Requires x is not negative infinity", x);
                    return double.NaN;
                }

                return double.PositiveInfinity;
            }

            double result = 1;
            if (x <= 0) {
                if (IsInteger(x)) {
                    Policies.ReportPoleError("Tgamma(x: {0}): Requires x != 0 and x be a non-negative integer", x);
                    return double.NaN;
                }

                if (x >= -GammaSmall.UpperLimit)
                    return GammaSmall.Tgamma(x);

                if (x <= -StirlingGamma.LowerLimit) {

                    const double MinX = -185;
                    if (x < MinX)
                        return 0;


                    // If we use the reflection formula directly for x < -NegMaxGamma, Tgamma will overflow; 
                    // so, use recurrence to get x > -NegMaxGamma.
                    // The result is likely to be subnormal or zero.
                    double shiftedX = x;
                    double recurrenceMult = 1.0;
                    double NegMaxGamma = -(Math2.MaxFactorialIndex + 1);
                    while (shiftedX < NegMaxGamma) {
                        recurrenceMult *= shiftedX;
                        shiftedX += 1;
                    }

                    result = -(Math.PI / StirlingGamma.Tgamma(-shiftedX)) / (shiftedX * recurrenceMult * Math2.SinPI(shiftedX));
                    return result;
                }

                // shift x to > 0:
                double product = 1;
                while (x < 0)
                    product *= x++;

                result /= product;

                // fall through
            }

            if ( x < 1 ) {

                if (x <= GammaSmall.UpperLimit) {
                    result *= GammaSmall.Tgamma(x);
                } else {
                    // Γ(x) = Exp(LogΓ(x+1))/x
                    result *= Math.Exp(_Lgamma.Rational_1_3(x + 1, x, x - 1))/x;
                }

            } else if (IsInteger(x) && (x <= Math2.FactorialTable.Length)) {
                // Γ(n) = (n-1)!
                result *= Math2.FactorialTable[(int)x - 1];

            } else if (x < StirlingGamma.LowerLimit) {

                // ensure x is in [1, 3)
                // if x >= 3, use recurrence to shift x to [2, 3)

                while (x >= 3)
                    result *= --x;

                result *= Math.Exp(_Lgamma.Rational_1_3(x, x - 1, x - 2));

            } else {

                result *= StirlingGamma.Tgamma(x);
            }

            return result;
        }

        /// <summary>
        /// Returns Γ(n)
        /// </summary>
        /// <param name="n">Argument</param>
        public static double Tgamma(int n)
        {

            if (n <= 0) {
                Policies.ReportPoleError("Tgamma(n: {0}): Requires n > 0", n);
                return double.NaN;
            }
            if (n >= FactorialTable.Length+1)
                return double.PositiveInfinity;
            return FactorialTable[n-1];
        }


        /// <summary>
        /// Returns the natural log of the Gamma Function = ln(|Γ(x)|). 
        /// Sets the sign = Sign(Γ(x)); 0 for poles.
        /// </summary>
        /// <param name="x">Lgamma function argument</param>
        /// <param name="sign">Lgamma output value = Sign(Γ(x))</param>
        public static double Lgamma(double x, out int sign)
        {
            sign = 0;
            if (double.IsNaN(x)) {
                Policies.ReportDomainError("Lgamma(x: {0}): NaN not allowed", x);
                return double.NaN;
            }
            if (double.IsInfinity(x)) {
                if (x < 0) {
                    Policies.ReportDomainError("Lgamma(x: {0}): Requires x != -Infinity", x);
                    return double.NaN;
                }
                sign = 1;
                return double.PositiveInfinity;
            }


            double result = 0;
            sign = 1;
            if (x <= 0) {
                if (IsInteger(x)) {
                    Policies.ReportPoleError("Lgamma(x: {0}): Evaluation at zero, or a negative integer", x);
                    return double.NaN;
                }

                if (x > -GammaSmall.UpperLimit)
                    return GammaSmall.Lgamma(x, out sign);

                if (x <= -StirlingGamma.LowerLimit) {
                    // Our Stirling routine does the reflection
                    // So no need to reflect here
                    result = StirlingGamma.Lgamma(x, out sign);
                    return result;

                } else {

                    double product = 1;
                    double xm2 = x - 2;
                    double xm1 = x - 1;
                    while (x < 1) {
                        product *= x;
                        xm2 = xm1;
                        xm1 = x;
                        x += 1;
                    }

                    if (product < 0) {
                        sign = -1;
                        product = -product;
                    }
                    result = -Math.Log(product) + _Lgamma.Rational_1_3(x, xm1, xm2);
                    return result;
                }

            } else if (x < 1) {

                if (x < GammaSmall.UpperLimit)
                    return GammaSmall.Lgamma(x, out sign);

                // Log(Γ(x)) = Log(Γ(x+1)/x)
                result = -Math.Log(x) + _Lgamma.Rational_1_3(x + 1, x, x - 1);

            } else if (x < 3) {

                result = _Lgamma.Rational_1_3(x, x - 1, x - 2);

            } else if (x < 8) {

                // use recurrence to shift x to [2, 3)

                double product = 1;
                while (x >= 3)
                    product *= --x;

                result = Math.Log(product) + _Lgamma.Rational_1_3(x, x - 1, x - 2);

            } else {
                // regular evaluation:
                result = StirlingGamma.Lgamma(x);
            }

            return result;

        }


        /// <summary>
        /// Returns the natural log of the Gamma Function = ln(|Γ(x)|)
        /// </summary>
        /// <param name="x">Lgamma function argument</param>
        public static double Lgamma(double x)
        {
            int sign;
            double result = Math2.Lgamma(x, out sign);
            return result;
        }

        /// <summary>
        /// Returns Γ(x + 1) - 1 with accuracy improvements in [-0.5, 2]
        /// </summary>
        /// <param name="x">Tgamma1pm1 function argument</param>
        public static double Tgamma1pm1(double x)
        {
            if (double.IsNaN(x)) {
                Policies.ReportDomainError("Tgamma1pm1(x: {0}): NaN not allowed", x);
                return double.NaN;
            }
            if (double.IsInfinity(x)) {
                if (x < 0) {
                    Policies.ReportDomainError("Tgamma1pm1(x: {0}): Requires x != -Infinity", x);
                    return double.NaN;
                }
                return double.PositiveInfinity;
            }

            double result;

            if (x < -0.5 || x >= 2)
                return Math2.Tgamma(1 + x) - 1; // Best method is simply to subtract 1 from tgamma:

            if (Math.Abs(x) < GammaSmall.UpperLimit)
                return GammaSmall.Tgamma1pm1(x);

            if (x < 0) {
                // compute exp( log( Γ(x+2)/(1+x) ) ) - 1 
                result = Math2.Expm1(-Math2.Log1p(x) + _Lgamma.Rational_1_3(x + 2, x + 1, x));
            } else {
                // compute exp( log(Γ(x+1)) ) - 1
                result = Math2.Expm1(_Lgamma.Rational_1_3(x + 1, x, x - 1));
            }

            return result;

        }
    }

} // namespaces
