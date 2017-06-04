//  Copyright (c) 2013 Rocco Romeo
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

//  Ported to .NET based on the original:
//      Copyright (c) 2006 Xiaogang Zhang, Boost Software License v1.0
//  History:
//      XZ wrote original code.
//      RR added full Hankel series for Bessel K, I

using System;
using System.Diagnostics;

namespace MathXK
{

    /// <summary>
    /// Hankel's asymptotic expansions for the Bessel Functions
    /// </summary>
    internal static class HankelAsym
    {
        #region Bessel J, Y


        /// <summary>
        /// Returns the minimum x value required to use the Hankel asymptotic for J{v) and Y{v)
        /// </summary>
        /// <param name="v"></param>
        /// <returns></returns>
        public static double JYMinX(double v)
        {
            // http://dlmf.nist.gov/10.74.i
            // If |x| is large compared with v^2 then the Hankel asymptotic expansion
            // is available

            // Note: 32 works consistently for double;
            // for half-integers, Hankel converges very quickly

            double v2 = v * v;
            double minX = ((v - Math.Floor(v)) == 0.5) ? 1 : 32;
            return Math.Max(v2, minX);
        }

        /// <summary>
        /// Computes J{v}(x) and Y{v}(x) using Hankel's asymptotic expansion.
        /// Generally available if |x| &gt; v^2 
        /// </summary>
        /// <param name="v"></param>
        /// <param name="x"></param>
        /// <returns>(J, Y) if the series converged; otherwise (NaN, NaN)</returns>
        public static (double J, double Y) JY(double v, double x)
        {
            // Preconditions
            Debug.Assert(x >= 0 && v >= 0);


            // Hankel Asymtotic expansion for large argument
            // valid for v fixed and x->inf
            // http://dlmf.nist.gov/10.17.i

            double z8 = 8 * x;
            double sq = 1;
            double mu = 4 * v * v;
            double term = 1;
            bool ok = false;


            double pSum = 1;
            double pPrev;
            double qSum = 0;
            double qPrev;

            double k = 1;
            while (k < Policies.MaxSeriesIterations) {

                // check to see when the series starts to diverge        
                qPrev = qSum;
                double mult = (2 * v - sq) * (2 * v + sq) / (k * z8);
                if (Math.Abs(mult) >= 1.0) 
                    break;

                term *= mult;
                qSum += term;
                k += 1;
                sq += 2;


                // check to see when the series starts to diverge
                pPrev = pSum;
                mult = -(2 * v - sq) * (2 * v + sq) / (k * z8);
                if (Math.Abs(mult) >= 1.0) 
                    break;

                term *= mult;
                pSum += term;
                k += 1;
                sq += 2;

                // if there are no changes, exit
                if (pPrev == pSum && qPrev == qSum) {
                    ok = true;
                    break;
                }

            }


            if (!ok) {
#if EXTRA_DEBUG
                Debug.WriteLine("Hankel did not converge: v = {0}, x = {1}", v, x);
#endif

                return (double.NaN, double.NaN);
            } 


            // The following is:
            //double ang = x - (v / 2 + 0.25) * Math.PI;
            //double sc = sin(ang);
            //double cc = cos(ang);
            //chi = sqrt(2 / (Math.PI * x));
            //Y = chi * (p * sc + q * cc);
            //J = chi * (p * cc - q * sc);
            //
            // Using angle addition rules
            // Note: sin(PI/4) = cos(PI/4) = 1/sqrt(2)

            var (sin, cos) = Math2.SinCos(x, -(Math2.Mod(v / 2, 2) + 0.25));
            double chi = Constants.RecipSqrtHalfPI / Math.Sqrt(x);

            double Y = chi * (pSum * sin + qSum * cos);
            double J = chi * (pSum * cos - qSum * sin);

            return (J,Y);
        }

        #endregion

        #region Bessel I, K

        /// <summary>
        /// Returns the minimum x value required to use the Hankel asymptotic for I{v) and K{v)
        /// </summary>
        /// <param name="v"></param>
        /// <returns></returns>
        public static double IKMinX(double v)
        {
            // http://dlmf.nist.gov/10.74.i
            // If |x| is large compared with v^2 then the Hankel asymptotic expansion
            // is available

            // Note: 32 seems to work consistently;

            return Math.Max(v*v, 32);
        }



        /// <summary>
        /// Returns I{v}(x) using Hankel’s Expansion for large argument. 
        /// <para>If <paramref name="expScale"/> == true, returns e^-x * I{v}(x)</para>
        /// </summary>
        /// <param name="v"></param>
        /// <param name="x"></param>
        /// <param name="expScale"></param>
        /// <returns></returns>
        /// <seealso href="http://dlmf.nist.gov/10.40#E1"/>
        public static double I(double v, double x, bool expScale = false)
        {
            // See: http://dlmf.nist.gov/10.40#i

            double mu = 4 * v * v;
            double ex = 8 * x;

            bool ok = false;
            double sum = 1, lastSum = sum;
            double sq = 1;
            double term = 1;
            int n = 1;
            for (; n < Policies.MaxSeriesIterations; n++) {
                // check to see if the series diverges
                double mult = -((mu - sq * sq) / (n * ex));
                if (Math.Abs(mult) >= 1.0)
                    break;


                lastSum = sum;
                term *= mult;
                sum += term;
                sq += 2;

                if (Math.Abs(term) <= Math.Abs(lastSum) * Policies.SeriesTolerance) { 
                    ok = true;
                    break;
                }

            }

            if (!ok) {
                Policies.ReportConvergenceError("HankelAsym.I(v: {0}, x: {1}): No convergence after {2} iterations", v, x, n);
                return double.NaN;
            }

            double result;
            if (expScale) {
                result = (sum / Math.Sqrt(x)) * Constants.RecipSqrt2PI;
            } else {
                // Try to avoid overflow:
                double e = Math.Exp(x / 2);
                result = e / Math.Sqrt(x) * Constants.RecipSqrt2PI * sum * e;
            }

            return result;
        }


        /// <summary>
        /// Returns Kv(x) using Hankel’s Expansion for large argument.
        /// <para>If <paramref name="expScale"/> == true, returns e^x * K{v}(x)</para>
        /// </summary>
        /// <param name="v"></param>
        /// <param name="x"></param>
        /// <param name="expScale"></param>
        /// <returns></returns>
        /// <seealso href="http://dlmf.nist.gov/10.40#E1"/>
        public static double K(double v, double x, bool expScale = false)
        {
            // See: http://dlmf.nist.gov/10.40#i

            double mu = 4 * v * v;
            double ex = 8 * x;

            bool ok = false;
            double sum = 1, lastSum = sum;
            double sq = 1;
            double term = 1;
            int n = 1;
            for (; n < Policies.MaxSeriesIterations; n++) {
                // check to see if the series diverges
                double mult = ((mu - sq * sq) / (n * ex));
                if (Math.Abs(mult) >= 1.0)
                    break;


                lastSum = sum;
                term *= mult;
                sum += term;
                sq += 2;

                if (Math.Abs(term) <= Math.Abs(lastSum) * Policies.SeriesTolerance) {
                    ok = true;
                    break;
                }

            }

            if (!ok) {
                Policies.ReportConvergenceError("HankelAsym.K(v: {0}, x: {1}): No convergence after {2} iterations", v, x, n);
                return double.NaN;
            }

            double result;
            if (expScale) {
                result = (sum / Math.Sqrt(x)) * Constants.SqrtHalfPI;
            } else {
                result = (Math.Exp(-x) / Math.Sqrt(x)) * sum * Constants.SqrtHalfPI;
            }

            return result;
        }

        #endregion



    }


}