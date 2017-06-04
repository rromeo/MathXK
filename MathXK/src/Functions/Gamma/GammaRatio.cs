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

    public static partial class Math2
    {

        /// <summary>
        /// Returns the ratio of gamma functions
        /// <para>TgammaDeltaRatio(x,delta) = Γ(x)/Γ(x+delta)</para>
        /// </summary>
        /// <param name="x">The argument for the gamma function in the numerator. Requires x &gt; 0</param>
        /// <param name="delta">The offset for the denominator. Requires x + delta &gt; 0</param>
        public static double TgammaDeltaRatio(double x, double delta)
        {
            double xpd = x + delta;
            if (double.IsNaN(x) || double.IsNaN(delta)) {
                Policies.ReportDomainError("TgammaDeltaRatio(x: {0}, delta: {1}): Requires x > 0; x + delta > 0", x, delta);
                return double.NaN;
            }

            // Γ(x)/Γ(x) == 1
            // this needs to be here in case x = 0
            if (delta == 0)
                return 1;

            if (!(x > 0)) {
                Policies.ReportDomainError("TgammaDeltaRatio(x: {0}, delta: {1}): Requires x > 0; negative x not implemented", x, delta);
                return double.NaN;
            }
            if (!(xpd > 0)) {
                Policies.ReportDomainError("TgammaDeltaRatio(x: {0}, delta: {1}): Requires x+delta > 0", x, delta);
                return double.NaN;
            }
            if (double.IsInfinity(x)) {
                if (double.IsInfinity(delta)) {
                    Policies.ReportDomainError("TgammaDeltaRatio(x: {0}, delta: {1}): Infinity/Infinity", x, delta);
                    return double.NaN;
                }
                return double.PositiveInfinity;
            }

            if (double.IsInfinity(delta))
                return 0.0;

            if (IsInteger(delta)) {
                if (IsInteger(x)) {
                    // Both x and delta are integers, see if we can just use table lookup
                    // of the factorials to get the result:
                    if ((x <= Math2.FactorialTable.Length) && (x + delta <= Math2.FactorialTable.Length))
                        return Math2.FactorialTable[(int)x - 1] / Math2.FactorialTable[(int)(x + delta) - 1];
                }


                // if delta is a small integer, we can use a finite product
                if (Math.Abs(delta) < 20) {
                    if (delta == 0)
                        return 1;
                    if (delta < 0) {
                        x -= 1;
                        double result = x;
                        while ( ++delta < 0 ) {
                            x -= 1;
                            result *= x;
                        }
                        return result;
                    } else {
                        double result = x;
                        while ( --delta > 0 ) {
                            x += 1;
                            result *= x;
                        }
                        return 1/result;
                    }
                }
            }

            // our Tgamma ratio already handles cases where 
            // one or both arguments are small 

            const double LowerLimit = DoubleLimits.MachineEpsilon / Constants.EulerMascheroni; // 1.73eps
            if (x < LowerLimit || xpd < LowerLimit)
                return Math2.TgammaRatio(x, xpd);

            if ((x < 1 && delta > Math2.FactorialTable.Length) || (xpd < 1 && x > Math2.FactorialTable.Length))
                return Math2.TgammaRatio(x, xpd);

            // Use the Lanczos approximation to compute the delta ratio 
            return Lanczos.TgammaDeltaRatio(x, delta);
        }


        /// <summary>
        /// Returns the ratio of gamma functions
        /// <para>TgammaRatio(x,y) = Γ(x)/Γ(y)</para>
        /// </summary>   
        /// <param name="x">The argument for the gamma function in the numerator. Requires x &gt; 0</param>
        /// <param name="y">The argument for the gamma function in the denominator. Requires y &gt; 0</param>
        public static double TgammaRatio(double x, double y)
        {
            if (!(x > 0) || !(y > 0)) {
                Policies.ReportDomainError("TgammaRatio(x: {0}, y: {1}): Requires x > 0, y > 0; negative arguments not implemented", x, y);
                return double.NaN;
            }
            if (double.IsInfinity(x)) {
                if (double.IsInfinity(y)) {
                    Policies.ReportDomainError("TgammaRatio(x: {0}, y: {1}): Infinity/Infinity", x, y);
                    return double.NaN;
                }
                return double.PositiveInfinity;
            }
            if (double.IsInfinity(y)) 
                return 0.0;


            // if we're not going to overflow, take the direct approach
            if (x <= DoubleLimits.MaxGamma && y <= DoubleLimits.MaxGamma) {
                // need to protect against denormalized numbers, so
                // use the first term of the power series: Gamma[x] ~= 1/x - Euler ...
                const double LowerLimit = DoubleLimits.MachineEpsilon / Constants.EulerMascheroni;
                if (x < LowerLimit) {
                    if (y < LowerLimit)
                        return y / x;

                    // careful: Tgamma(y) can be < 1, x can be denorm or Min Normal
                    double result = 1 / (x * Tgamma(y));
                    return result;
                }
                if (y < LowerLimit) {
                    double result = y * Tgamma(x);
                    return result;
                }

                return Tgamma(x) / Tgamma(y);
            }

            // Handle the case where we have one large argument and one small argument

            Debug.Assert(x > DoubleLimits.MaxGamma || y > DoubleLimits.MaxGamma);

            if (x <= DoubleLimits.MaxGamma) {
                // At max values
                // min denorm:  Γ(2^-1074)/Γ(170) ~= 4.7e18
                // min norm:    Γ(2^-1022)/Γ(170) ~= 1052.74
                //              Γ(1)/Γ(170)       ~= 2.3e-305

                // use duplication
                // Γ(2z) = 2^(2z-1)/Sqrt(Pi) * Γ(z) * Γ(z+1/2) 
                if ( y <= 2 * DoubleLimits.MaxGamma - 1) {
                    var t1 = Tgamma(y / 2);
                    var t2 = Math.Pow(2, y) / (2 * Constants.SqrtPI);
                    var t3 = Tgamma(y / 2 + 0.5);
                    if (x < DoubleLimits.MachineEpsilon/Constants.EulerMascheroni ) {
                        return 1.0/(x * t1 * t2 * t3);
                    } else {
                        return Tgamma(x)/ t1 / t2 / t3;
                    }
                }

                // this will underflow
                return Math.Exp(Lgamma(x) - Lgamma(y));

            }
            if (y <= DoubleLimits.MaxGamma) {
                // At max values
                // min denorm:  Γ(170)/Γ(2^-1074)    ~=    2.1e-19
                // min norm     Γ(170)/Γ(2^-1022)    ~=    9.5e-4
                //              Γ(170)/Γ(1)          ~=    4.3e304

                // use duplication
                // Γ(2z) = 2^(2z-1)/Sqrt(Pi) * Γ(z) * Γ(z+1/2) 
                if (y <= 2 * DoubleLimits.MaxGamma - 1) {
                    var t1 = Tgamma(x / 2);
                    var t2 = Math.Pow(2, x) / (2 * Constants.SqrtPI);
                    var t3 = Tgamma(x / 2 + 0.5);
                    if (y < DoubleLimits.MachineEpsilon / Constants.EulerMascheroni) {
                        return (y * t1 * t2 * t3);
                    } else {
                        return t1/Tgamma(y) * t2 * t3;
                    }
                }

                // this will overflow
                return Math.Exp(Lgamma(x) - Lgamma(y));

            }

            // Regular case, x and y both large and similar in magnitude:

            return TgammaDeltaRatio(x, y - x);

        }
    }

} //namespaces