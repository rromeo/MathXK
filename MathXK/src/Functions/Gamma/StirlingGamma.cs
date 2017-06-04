//  Copyright (c) 2013 Rocco Romeo
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

using System;
using System.Diagnostics;
using MathXK.Numerics;

namespace MathXK
{

    /// <summary>
    /// Stirling asymptotic approximation to the Gamma Function
    /// </summary>
    internal static class StirlingGamma
    {
        /// <summary>
        /// The lower absolute limit on the Stirling series routines = 8
        /// </summary>
        public const double LowerLimit = 8;


        /// <summary>
        /// Returns Ln(|Γ(x)|) for |x| &gt;= <c>LowerLimit</c> using the Stirling asymptotic series
        /// </summary>
        /// <param name="x">The argument</param>
        /// <param name="sign">Sets sign = Sign(Γ(x))</param>
        /// <returns></returns>
        public static double Lgamma(double x, out int sign)
        {
            Debug.Assert(Math.Abs(x) >= LowerLimit);

            // (Ln(2*Pi)-1)/2 = Ln(Sqrt(2*Pi/E))
            const double Half_Ln2PIm1 = 0.418938533204672741780329736405617639861397473637783412817151;

            sign = 0;
            if (x > 0) {
                //return (x - 0.5) * Math.Log(x) - x + Constants.Ln2PI / 2 + LgammaSeries(x);
                sign = 1;
                return (x - 0.5) * (Math.Log(x) - 1) + (Half_Ln2PIm1 + LgammaSeries(x));
            }

            // Use the reflection formula: Γ(-x) = -π/(x * Γ(x) * Sin(π * x)) 
            x = -x;
            double sin = Math2.SinPI(x);
            if (sin == 0)
                return double.NaN;
            if (sin < 0) {
                sin = -sin;
                sign = 1;
            } else
                sign = -1;

            sin *= (x / Math.PI);
            return -Math.Log(sin) - Lgamma(x);
        }

        /// <summary>
        /// Returns Ln(|Γ(x)|) for |x| &gt;= <c>LowerLimit</c> using the Stirling asymptotic series
        /// </summary>
        /// <param name="x">The argument</param>
        /// <returns></returns>
        public static double Lgamma(double x)
        {
            Debug.Assert(Math.Abs(x) >= LowerLimit);

            int sign = 0;
            return Lgamma(x, out sign);
        }

        /// <summary>
        /// Returns the sum of the Stirling series fractional 
        /// <para>LgammaSeries(x) = 1/(12 x) - 1/(360 x^3) + 1/(1260 x^5) - 1/(1680 x^7) + ... </para>
        /// <para>Ln(|Γ(x)|) = ((x - 0.5) * Math.Log(x) - x + Ln(2π)/2 + LgammaSeries(x)</para>
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        public static double LgammaSeries(double x)
        {
            Debug.Assert(x >= LowerLimit);

            // 1/(12 x) - 1/(360 x^3) + 1/(1260 x^5) - 1/(1680 x^7) + 1/(1188 x^9) - 691/(360360 x^11) + 1/(156 x^13) - 3617/(122400 x^15) + 43867/(244188 x^17) - 174611/(125400 x^19) + 77683/(5796 x^21)

            const double c0 = 1.0 / 12;
            const double c1 = -1.0 / 360;
            const double c2 = 1.0 / 1260;
            const double c3 = -1.0 / 1680;
            const double c4 = 1.0 / 1188;
            const double c5 = -691.0 / 360360;

            const double c6 = 1.0 / 156;
            const double c7 = -3617.0 / 122400;
            const double c8 = 43867.0 / 244188;
            const double c9 = -174611.0 / 125400;

            double r = 1 / (x * x);
            double s5 = 0;
            if (x < 27)
                s5 = c5 + r * (c6 + r * (c7 + r * (c8 + r * c9)));
            double result = (c0 + r * (c1 + r * (c2 + r * (c3 + r * (c4 + r * s5)))));
            return result / x;

        }


        /// <summary>
        /// Returns the gamma series, aka the scaled gamma function Γ<sup>*</sup>(x)
        /// <para>GammaSeries(x) = 1 + 1/(12x) + 1/(288x^2) - 139/(51840x^3) ...</para>
        /// <para>Γ(x) = (Sqrt(2*π) * Pow(x, x-1/2) * Exp(-x)) * GammaSeries(x)</para>
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        public static double GammaSeries(double x)
        {
            Debug.Assert(x >= LowerLimit);

            // At x == 26, 1-(c8/x^8)/Series error =  1.5e-16
            // At x == 16, 1-(c10/x^10)/Series error = 1.1e-16
            // At x == 8, 1-(c16/x^16)/Series error = 7.1e-17

            // N[Series[Gamma[x]/(Sqrt[2*Pi]*Power[x, x - 1/2]*Exp[-x]), {x, Infinity, 16}], 50]

            const double c0  = 1;
            const double c1 = 0.083333333333333333333333333333333333333333333333333;
            const double c2 = 0.0034722222222222222222222222222222222222222222222222;
            const double c3 = -0.0026813271604938271604938271604938271604938271604938;
            const double c4 = -0.00022947209362139917695473251028806584362139917695473;
            const double c5 = 0.00078403922172006662747403488144228884969625710366451;
            const double c6 = 0.000069728137583658577742939882857578330829359635943998;
            const double c7 = -0.00059216643735369388286483622560440118739158519679782;
            const double c8 = -0.000051717909082605921933705784300205882281785345342730;
            const double c9 = 0.00083949872067208727999335751676498344519818211159301;
            const double c10 = 0.000072048954160200105590857193022501505206345173797548;
            const double c11 = -0.0019144384985654775265008988583285225448768935789535;
            const double c12 = -0.00016251626278391581689863512398027099810587259193225;
            const double c13 = 0.0064033628338080697948236380902657958304018940939629;
            const double c14 = 0.00054016476789260451518046750857024173554725441597922;
            const double c15 = -0.029527880945699120505440651054693824446565482825438;
            const double c16 = -0.0024817436002649977309156583687434643239751680472290;

            double z = 1/x;
            double s9 = 0;
            if ( x < 26 )
                s9 = c9 + z * (c10 + z * (c11 + z * (c12 + z * (c13 + z * (c14 + z * (c15 + z * c16))))));
            double result = c0 + z * (c1 + z * (c2 + z * (c3 + z * (c4 + z * (c5 + z * (c6 + z * (c7 + z * (c8 + z * s9 ))))))));
            return result;
        }

        ///<summary>
        /// Returns Γ(x) for |x| &gt;= <c>LowerLimit</c> using the Stirling asymptotic series
        /// <para>Γ(x) = x^x * e^-x * Sqrt(2π/x) * GammaSeries(x)</para>
        /// </summary>  
        public static double Tgamma(double x)
        {
            Debug.Assert(x >= LowerLimit);

            // check for large x to avoid potential Inf/Inf below
            if (x > DoubleLimits.MaxLogValue )
                return double.PositiveInfinity;

            double result = Constants.Sqrt2PI * GammaSeries(x);
            if ( x > DoubleLimits.MaxArgXPowX ) {
                // we're going to overflow unless this is done with care:
                double hp = Math.Pow(x, (x / 2) - 0.25);
                result *= (hp / Math.Exp(x));
                result *= hp;
            } else {
                result *= (Math.Pow(x, x - 0.5) / Math.Exp(x));
            }

            return result;
        }

        /// <summary>
        /// Returns Γ(z)/Γ(z+delta) using Stirling approximations 
        /// </summary>
        public static double TgammaDeltaRatio(double x, double delta)
        {
            Debug.Assert(x >= LowerLimit && delta + x >= LowerLimit);

            // Compute Γ(x)/Γ(x+Δ), which reduces to:
            //
            // (x/(x+Δ))^(x-0.5) * (e/(x+Δ))^Δ * (GammaSeries(x)/GammaSeries(x+Δ))
            //  OR
            // (1+Δ/x)^-(x+Δ-0.5) * (e/x)^Δ * (GammaSeries(x)/GammaSeries(x+Δ))

            double h = delta / x;

            double factor = (Math.Abs(delta) < 1) ? Math.Exp(LgammaSeries(x) - LgammaSeries(x + delta)) : GammaSeries(x) / GammaSeries(x + delta);
            factor *= Math.Sqrt(1 + h);

            double d = delta / (x + delta);
            if (Math.Abs(d) <= 0.5) {
                double l = (x + delta) * Math2.Log1pmx(-d);
                double deltalogx = -delta * Math.Log(x);
                if (deltalogx > DoubleLimits.MinLogValue && deltalogx < DoubleLimits.MaxLogValue)
                    return factor * Math.Pow(x, -delta) * Math.Exp(l);
                return factor * Math.Exp(l + deltalogx);
            }


            double result;
            if (Math.Abs(h) <= 0.5) {
                result = Math.Exp(-(x + delta) * Math2.Log1p(h));
            } else {
                result = Math.Pow(x/(x+delta), (x + delta));
            }

            result *= factor;
            result *= Math.Pow(Math.E/x, delta);

            return result;
        }

        /// <summary>
        /// Returns Γ(x) - Γ(x + Δ), for large x, small Δ
        /// </summary>
        /// <param name="x"></param>
        /// <param name="delta"></param>
        /// <returns></returns>
        public static double LgammaDelta(double x, double delta)
        {
            Debug.Assert(x >= LowerLimit && delta + x >= LowerLimit);

            // Compute Log(Γ(x)/Γ(x+Δ)):
            
            // Γ(x)/Γ(x+Δ) reduces to:
            // (x/(x+Δ))^(x-0.5) * (e/(x+Δ))^Δ * (GammaSeries(x)/GammaSeries(x+Δ))
            //  OR
            // (1+Δ/x)^-(x+Δ-0.5) * (e/x)^Δ * (GammaSeries(x)/GammaSeries(x+Δ))

            double lf;
            if (Math.Abs(delta) < 1)
                lf = LgammaSeries(x) - LgammaSeries(x + delta);
            else
                lf = Math.Log(GammaSeries(x) / GammaSeries(x + delta));

            double h = delta / x;
            double result;
            if (Math.Abs(h) <= 0.5) {
                result = -(x - 0.5 + delta) * Math2.Log1p(h)  + delta * (1 - Math.Log(x)) + lf;
            } else {
                result = (x - 0.5) * Math.Log(x / (x + delta)) + delta * (1 - Math.Log(x + delta)) + lf;
            }

            return result;


        }

        /// <summary>
        /// Returns Beta(a, b) using Stirling approximations
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <returns></returns>
        public static double Beta(double a, double b)
        {
            Debug.Assert(a >= LowerLimit && b >= LowerLimit);

            // Beta(a, b) = 
            // (a/c)^a * (b/c)^b * Sqrt(2πc/(a*b)) * (GammaSeries(a)*GammaSeries(b)/GammaSeries(c)) 

            if (a < b)
                Utility.Swap(ref a, ref b);

            double c = a + b;

            double result = (GammaSeries(a) * GammaSeries(b) / GammaSeries(c));
            result *= Constants.Sqrt2PI * Math.Sqrt(1 / a + 1 / b);

            // calculate (a/c)^a or the equivalent (1 - b/c)^a 
            // calculate (b/c)^b or the equivalent (1 - a/c)^b 

            double t1 = b / c;
            result *= (t1 <= 0.5) ? Math.Exp(a * Math2.Log1p(-t1)) : Math.Pow(a / c, a);

            double t2 = a / c;
            result *= (t2 <= 0.5) ? Math.Exp(b * Math2.Log1p(-t2)) : Math.Pow(b / c, b);

            return result;
        }

        /// <summary>
        /// Returns Log(Beta(a,b)) using Stirling approximations
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <returns></returns>
        public static double LogBeta(double a, double b) 
        {
            Debug.Assert(a >= LowerLimit && b >= LowerLimit);

            double c = a + b;

            // Beta(a, b) = 
            // (a/c)^a * (b/c)^b * Sqrt(2πc/(a*b)) * (GammaSeries(a)*GammaSeries(b)/GammaSeries(c)) 

            double logValue = 0;
            double factor = GammaSeries(a) * GammaSeries(b) / GammaSeries(c);
            factor *= Constants.Sqrt2PI * Math.Sqrt(1/a + 1/b);


            // calculate (a/c)^a or the equivalent(1-b/c)^a 
            // choose the latter method if agh/cgh is sufficiently close to 1

            double t1 = b / c;
            if (t1 < 0.5)
                logValue += a * Math2.Log1p(-t1);
            else
                logValue += a * Math.Log(a / c);

            // calculate (b/c)^b or the equivalent (1 - a/c)^b 
            // choose the latter method if b/c is sufficiently close to 1

            double t2 = a / c;
            if (t2 < 0.5)
                logValue += b * Math2.Log1p(-t2);
            else
                logValue += b * Math.Log(b / c);

            logValue += Math.Log(factor);

            return logValue;

        
        }


        /// <summary>
        /// Returns x^a/Γ(a) using the Stirling approximation 
        /// </summary>
        /// <param name="a"></param>
        /// <param name="x"></param>
        public static double PowDivGamma(double x, double a)
        {
            Debug.Assert(a >= LowerLimit);

            // Stirling calculation. Compute:
            // (x/a)^a * e^a * Sqrt(a/(2pi)) / GammaSeries(a)


            double factor = Constants.RecipSqrt2PI * Math.Sqrt(a) / GammaSeries(a);
            Debug.Assert(factor >= 1);

            double t = x / a;
            if (t >= 1)
                return Math.Pow(t, a) * Math.Exp(a) * factor;

            // at this point a > 8 && 0 < t < 1
            double alogt = a * Math.Log(t);
            if (a <= DoubleLimits.MaxLogValue && alogt >= DoubleLimits.MinLogValue) {
                return Math.Pow(t, a) * Math.Exp(a) * factor;
            }
            if (a <= 2 * DoubleLimits.MaxLogValue && alogt >= 2 * DoubleLimits.MinLogValue) {
                double hp = Math.Pow(t, a / 2) * Math.Exp(a / 2);
                return hp * factor * hp;
            }

            return Math.Exp(alogt + a + Math.Log(factor));

        }

        /// <summary>
        /// Returns m*e^exp
        /// </summary>
        /// <param name="exp"></param>
        /// <param name="m"></param>
        /// <returns></returns>
        private static double MultExp(double m, double exp)
        {
            if (exp >= DoubleLimits.MinLogValue && exp <= DoubleLimits.MaxLogValue)
                return m * Math.Exp(exp);

            double exph = Math.Exp(exp / 2);
            return m * exph * exph;
        }


        /// <summary>
        /// Compute (x^a)(e^-x)/Γ(a) using the Stirling approximation
        /// </summary>
        /// <param name="a"></param>
        /// <param name="x"></param>
        /// <returns></returns>
        public static double PowExpDivGamma(double a, double x)
        {
            Debug.Assert(a >= LowerLimit && x >= 0);

            // Compute: (x/a)^a * e^(a-x) * Sqrt(a/(2pi))  / GammaSeries(a)

            double factor = Constants.RecipSqrt2PI * Math.Sqrt(a) / GammaSeries(a);
            Debug.Assert(factor >= 1);

            // for a ~ x
            // a * (ln(1+(x-a)/a) - ((x-a)/a)) === a*ln(x/a) + a - x 
            var d = (x - a) / a;
            if (Math.Abs(d) < 0.5) {
                double u = a * Math2.Log1pmx(d);
                return MultExp(factor, u);
            }

            double t = x / a;
            double alogt = a * Math.Log(t);
            double amx = a - x;
            if ( t > 1 ) {
                // x > a

                if (amx >= DoubleLimits.MinLogValue && alogt <= DoubleLimits.MaxLogValue)
                    return Math.Pow(t, a) * Math.Exp(amx) * factor;
                if (amx >= 2*DoubleLimits.MinLogValue && alogt <= 2 * DoubleLimits.MaxLogValue) {
                    double h = Math.Pow(t, a / 2) * Math.Exp(amx / 2);
                    return h * factor * h;
                }

            } else {
                // x < a
                if (amx <= DoubleLimits.MaxLogValue && alogt >= DoubleLimits.MinLogValue)
                    return Math.Pow(t, a) * Math.Exp(amx) * factor;
                if (amx <= 2 * DoubleLimits.MaxLogValue && alogt >= 2 * DoubleLimits.MinLogValue) {
                    double h = Math.Pow(t, a / 2) * Math.Exp(amx / 2);
                    return h * factor * h;
                }

            }

            // okay, use logs
            return Math.Exp(alogt + amx + Math.Log(factor));

        }



    }

} // namespace
