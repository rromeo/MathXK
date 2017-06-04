//  Copyright (c) 2013 Rocco Romeo
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

//  Ported to .NET based on the original:
//      Copyright (c) John Maddock 2006, Boost Software License v1.0
//  History:
//      JM wrote original C++ code.
//      RR improved edge cases and added asymptotics (based on DiDonato and Morris).

#if DEBUG
//#define EXTRA_DEBUG
#endif

using System;
using System.Diagnostics;

namespace MathXK
{

    internal static class _Ibeta
    {

        /// <summary>
        /// Compute the leading power terms in the incomplete Beta:
        /// <para>(x^a)(y^b)/Beta(a,b) when regularized</para>
        /// <para>(x^a)(y^b) otherwise</para>
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="x"></param>
        /// <param name="y">1-x</param>
        /// <param name="normalised"></param>
        /// <param name="multiplier"></param>
        /// <returns></returns>
        /// <remarks>
        /// Almost all of the error in the incomplete beta comes from this
        /// function: particularly when a and b are large. Computing large
        /// powers are *hard* though, and using logarithms just leads to
        /// horrendous cancellation errors.
        /// </remarks>
        public static double PowerTerms(double a, double b, double x, double y, bool normalised, double multiplier)
        {
            if (multiplier == 0)
                return 0;

            if (!normalised)
                return multiplier * (Math.Pow(x, a) * Math.Pow(y, b));

            double c = a + b;

            // combine power terms with Lanczos approximation:
            double agh = a + (Lanczos.G - 0.5);
            double bgh = b + (Lanczos.G - 0.5);
            double cgh = c + (Lanczos.G - 0.5);

            // compute:
            // mutiplier * (x*cgh/agh)^a * (y*cgh/bgh)^b * Sqrt((agh*bgh)/(e*cgh)) * (Sum_expG_scaled(c) / (Sum_expG_scaled(a) * Sum_expG_scaled(b)))

            double factor = (Lanczos.SeriesExpGScaled(c) / Lanczos.SeriesExpGScaled(a)) / Lanczos.SeriesExpGScaled(b);
            if (a > b)
                factor *= Math.Sqrt((bgh / Math.E) * (agh / cgh));
            else
                factor *= Math.Sqrt((agh / Math.E) * (bgh / cgh));

            if (double.IsInfinity(factor * multiplier)) {

#if EXTRA_DEBUG
                Debug.WriteLine("Using Logs in _Beta.PowerTerms: a = {0}; b = {1}; x = {2}; y= {3}; mult = {4}", a, b, x, y, multiplier);
#endif
                // this will probably fail but...
                double logValue;
                if (x <= y)
                    logValue = a * Math.Log(x) + b * Math2.Log1p(-x);
                else
                    logValue = a * Math2.Log1p(-y) + b * Math.Log(y);
                logValue += Math.Log(multiplier);
                logValue -= Math2.LogBeta(a, b);
                return Math.Exp(logValue);
            } else {
                factor *= multiplier;
            }

            double result = 1;

            //Debug.WriteLine("IbetaPowerTerms(a: {0}, b: {1}, x: {2}, y: {3}) - Prefix = {4}", a, b, x, y, result);


            // l1 and l2 are the base of the exponents minus one:
            double l1 = (x * b - y * agh) / agh;
            double l2 = (y * a - x * bgh) / bgh;

            // t1 and t2 are the exponent terms
            double t1 = (x * cgh) / agh;
            double t2 = (y * cgh) / bgh;


            // if possible, set the result to partially cancel out with the first term
            bool sameDirection = false;
            if (t1 >= 1 && t2 >= 1) {
                sameDirection = true;
                if (factor < 1) {
                    result = factor;
                    factor = 1;
                }
            } else if (t1 <= 1 && t2 <= 1) {
                sameDirection = true;
                if (factor > 1) {
                    result = factor;
                    factor = 1;
                }
            }

            if (sameDirection) {
                // This first branch handles the simple case where the two power terms 
                // both go in the same direction (towards zero or towards infinity). 
                // In this case if either term overflows or underflows, 
                // then the product of the two must do so also.  

                if (Math.Abs(l1) < 0.5)
                    result *= Math.Exp(a * Math2.Log1p(l1));
                else
                    result *= Math.Pow(t1, a);

                if (Math.Abs(l2) < 0.5)
                    result *= Math.Exp(b * Math2.Log1p(l2));
                else
                    result *= Math.Pow(t2, b);

                result *= factor;

            } else {
                // This second branch handles the case where the two power terms 
                // go in opposite directions (towards zero or towards infinity). 

                Debug.Assert(result == 1);

                bool useExpT1 = false;
                bool useExpT2 = false;

                double logt1;
                if (Math.Abs(l1) < 0.5) {
                    useExpT1 = true;
                    logt1 = a * Math2.Log1p(l1);
                } else
                    logt1 = a * Math.Log(t1);

                double logt2;
                if (Math.Abs(l2) < 0.5) {
                    useExpT2 = true;
                    logt2 = b * Math2.Log1p(l2);
                } else
                    logt2 = b * Math.Log(t2);

                if ((logt1 >= DoubleLimits.MaxLogValue) || (logt1 <= DoubleLimits.MinLogValue)
                || (logt2 >= DoubleLimits.MaxLogValue) || (logt2 <= DoubleLimits.MinLogValue)
                ) {
                    double logSum = logt1 + logt2;
                    if ((logSum >= DoubleLimits.MaxLogValue) || (logSum <= DoubleLimits.MinLogValue))
                        result = Math.Exp(logSum + Math.Log(factor));
                    else
                        result = factor * Math.Exp(logSum);
                } else {

                    // ensure that t1 and result will partially cancel
                    if (t1 >= 1) {
                        if (factor < 1) {
                            result = factor;
                            factor = 1;
                        }
                    } else {
                        if (factor > 1) {
                            result = factor;
                            factor = 1;
                        }
                    }

                    if (useExpT1)
                        result *= Math.Exp(logt1);
                    else
                        result *= Math.Pow(t1, a);

                    if (useExpT2)
                        result *= Math.Exp(logt2);
                    else
                        result *= Math.Pow(t2, b);

                    result *= factor;

                }


            }
            return result;
        }


        public static double PowerTerms(double a, double b, double x, double y, bool normalised)
        {
            return PowerTerms(a, b, x, y, normalised, 1);
        }

        /// <summary>
        /// Compute: multiplier * x^a/Beta(a,b) while trying to avoid overflows/underflows
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="multiplier"></param>
        /// <returns></returns>
        static double SeriesRegularizedPrefix(double a, double b, double x, double y, double multiplier)
        {
            double c = a + b;

            // incomplete beta power term, combined with the Lanczos approximation:
            double agh = a + (Lanczos.G - 0.5);
            double bgh = b + (Lanczos.G - 0.5);
            double cgh = c + (Lanczos.G - 0.5);


            // compute:
            // mutiplier * (x*cgh/agh)^a * (cgh/bgh)^b * Sqrt((agh*bgh)/(e*cgh)) * (SumGScaled(c)/(SumGScaled(a)*SumGScaled(b)) 

            double factor = (Lanczos.SeriesExpGScaled(c) / Lanczos.SeriesExpGScaled(a)) / Lanczos.SeriesExpGScaled(b);
            if (a > b)
                factor *= Math.Sqrt((bgh / Math.E) * (agh / cgh));
            else
                factor *= Math.Sqrt((agh / Math.E) * (bgh / cgh));
            factor *= multiplier;

            double result = 1;


            // l1 and l2 are the base of the exponents minus one:
            double l1 = (x * b - y * agh) / agh;
            double l2 = a / bgh;

            // t1 and t2 are the exponent terms
            double t1 = (x * cgh) / agh;
            double t2 = cgh / bgh;


            // if possible, set the result to partially cancel out with the first term
            Debug.Assert(t2 >= 1);

            if (t1 >= 1) {
                // This first branch handles the simple case where the two power terms 
                // both go in the same direction (towards zero or towards infinity). 
                // In this case if either term overflows or underflows, 
                // then the product of the two must do so also.  

                if (factor < 1) {
                    result = factor;
                    factor = 1;
                }

                if (Math.Abs(l1) < 0.5)
                    result *= Math.Exp(a * Math2.Log1p(l1));
                else
                    result *= Math.Pow(t1, a);

                if (Math.Abs(l2) < 0.5)
                    result *= Math.Exp(b * Math2.Log1p(l2));
                else
                    result *= Math.Pow(t2, b);

                result *= factor;

            } else {
                // This second branch handles the case where the two power terms 
                // go in opposite directions (towards zero or towards infinity). 

                Debug.Assert(result == 1);
                Debug.Assert(t1 < 1);

                bool useExpT1 = false;
                bool useExpT2 = false;

                double logt1;
                if (Math.Abs(l1) < 0.5) {
                    useExpT1 = true;
                    logt1 = a * Math2.Log1p(l1);
                } else
                    logt1 = a * Math.Log(t1);

                double logt2;
                if (Math.Abs(l2) < 0.5) {
                    useExpT2 = true;
                    logt2 = b * Math2.Log1p(l2);
                } else
                    logt2 = b * Math.Log(t2);

                if ((logt1 >= DoubleLimits.MaxLogValue) || (logt1 <= DoubleLimits.MinLogValue)
                || (logt2 >= DoubleLimits.MaxLogValue) || (logt2 <= DoubleLimits.MinLogValue)
                ) {
                    double logSum = logt1 + logt2;
                    if ((logSum >= DoubleLimits.MaxLogValue) || (logSum <= DoubleLimits.MinLogValue))
                        result = Math.Exp(logSum + Math.Log(factor));
                    else
                        result = factor * Math.Exp(logSum);
                } else {

                    // Set result so that t1 and result will partially cancel
                    if (factor > 1) {
                        result = factor;
                        factor = 1;
                    }

                    if (useExpT1)
                        result *= Math.Exp(logt1);
                    else
                        result *= Math.Pow(t1, a);

                    if (useExpT2)
                        result *= Math.Exp(logt2);
                    else
                        result *= Math.Pow(t2, b);

                    result *= factor;

                }
            }

            return result;

        }


        /// <summary>
        /// Series approximation to the incomplete beta:
        /// <para>Σ( (1-b)_{k} * x^k/( k! * (a+k) ), k={0, Inf}</para>
        /// <para>Bx(a, b) = z^a/a * SeriesSum</para>
        /// <para>Ix(a, b) = z^a/(a*B(a,b)) * SeriesSum</para>
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="x"></param>
        /// <returns></returns>
        static double SeriesSum(double a, double b, double x)
        {
            const double DefaultTolerance = 2 * DoubleLimits.MachineEpsilon;

            double sum = 1 / a;
            double term = 1.0;
            for (int n = 1; n < Policies.MaxSeriesIterations; n++) {
                double prevSum = sum;
                term *= (1.0 - b / n) * x;
                double delta = term / (a + n);
                sum += delta;

                if (Math.Abs(delta) <= Math.Abs(prevSum) * DefaultTolerance)
                    return sum;
            }

            Policies.ReportConvergenceError("Series did not converge after {0} iterations", Policies.MaxSeriesIterations);
            return double.NaN;
        }



        /// <summary>
        /// Evaluate the incomplete beta function using a series expansion
        /// <para>Bx(a, b) = z^a/a * Σ( (1-b)_{k} * x^k/( k! * (a+k) ), k={0, Inf}</para>
        /// <para>Ix(a, b) = z^a/(a*B(a,b)) * Σ( (1-b)_{k} * x^k/( k! * (a+k) ), k={0, Inf}</para>
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="regularized"></param>
        /// <returns></returns>
        public static double Series1(double a, double b, double x, double y, bool regularized)
        {
            Debug.Assert(b <= 1 || b * x <= 0.7);

            // See: http://functions.wolfram.com/GammaBetaErf/Beta3/06/01/03/01/01/0003/

            if (regularized) {
                double seriesSum = SeriesSum(a, b, x);
                double result = SeriesRegularizedPrefix(a, b, x, y, seriesSum);
                return result;
            }

            return Math.Pow(x, a) * SeriesSum(a, b, x);

        }



        /// <summary>
        /// Evaluate the incomplete beta via the continued fraction representation 
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="normalised"></param>
        /// <returns></returns>
        /// <remarks>
        /// Continued fraction for Ibeta- See: http://dlmf.nist.gov/8.17.E22
        /// Converges rapidly when x &lt; (a+1)/(a+b+2); for x &gt; (a+1)/(a+b+2) || 1-x &lt; (b+1)/(a+b+2), use I{1-x}(a, b)
        /// Continued fraction for the incomplete beta used when a &gt; 1 and b &gt; 1
        /// <para>According to NR this is O(sqrt(max(a,b))</para>
        /// </remarks>
        public static double Fraction2(double a, double b, double x, double y, bool normalised)
        {
            Debug.Assert(x <= a / (a + b) || y <= b / (a + b));

            double result = PowerTerms(a, b, x, y, normalised);
            if (result == 0)
                return result;

            // Define the continued fraction function

            int m = 0;
            Func<(double an, double bn)> fracF = () => {
                double aN = (a + m - 1) * (a + b + m - 1) * m * (b - m) * x * x;
                double denom = (a + 2 * m - 1);
                aN /= denom * denom;

                double bN = m;
                bN += (m * (b - m) * x) / (a + 2 * m - 1);
                bN += ((a + m) * (a - (a + b) * x + 1 + m * (2 - x))) / (a + 2 * m + 1);

                ++m;

                return (aN, bN);
            };


            var (_, b0) = fracF(); // get b0 - skip a0;
            double fract = ContinuedFraction.Eval(b0, fracF);

#if EXTRA_DEBUG
            double cvg = (a + 1) / (a + b + 2);
            Debug.WriteLine("CF B({0}, {1}, {2}): cvg = {3}; iterations = {4}", a, b, x, cvg, m);
#endif

            return result / fract;
        }


        /// <summary>
        /// Implementation for the non-regularized cases: Bx(a,b) and B(a,b)-Bx(a,b)
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="x"></param>
        /// <param name="invert"></param>
        /// <returns></returns>
        public static double BetaImp(double a, double b, double x, bool invert)
        {
            Debug.Assert(a > 0);
            Debug.Assert(b > 0 && b < int.MaxValue);
            Debug.Assert(x > 0 && x < 1);

            const bool normalised = false;
            double fract;
            double y = 1 - x;


            if (a <= 1 && b <= 1) {
                if (x > 0.5) {
                    Utility.Swap(ref a, ref b);
                    Utility.Swap(ref x, ref y);
                    invert = !invert;
                }

                // Both a,b < 1:
                if ((a >= Math.Min(0.2, b)) || (Math.Pow(x, a) <= 0.9)) {
                    fract = Series1(a, b, x, y, false);
                    return (invert) ? Math2.Beta(a, b) - fract : fract;
                }


                Utility.Swap(ref a, ref b);
                Utility.Swap(ref x, ref y);
                invert = !invert;

                if (y >= 0.3) {
                    fract = Series1(a, b, x, y, false);
                    return (invert) ? Math2.Beta(a, b) - fract : fract;
                }

                // Sidestep on a, and then use the series representation:
                double prefix = RisingFactorialRatio(a + b, a, 20);
                fract = IbetaAStep(a, b, x, y, 20, normalised);
                fract += IbetaLargeASmallBSeries(a + 20, b, x, y, normalised, prefix);
                return (invert) ? Math2.Beta(a, b) - fract : fract;
            }


            if (a <= 1 || b <= 1) {
                if (x > 0.5) {
                    Utility.Swap(ref a, ref b);
                    Utility.Swap(ref x, ref y);
                    invert = !invert;
                }

                // One of a, b < 1 only:
                if ((b <= 1) || ((x < 0.1) && (Math.Pow(b * x, a) <= 0.7))) {
                    fract = Series1(a, b, x, y, false);
                    return (invert) ? Math2.Beta(a, b) - fract : fract;
                }

                Utility.Swap(ref a, ref b);
                Utility.Swap(ref x, ref y);
                invert = !invert;

                if (y >= 0.3) {
                    fract = Series1(a, b, x, y, false);
                    return (invert) ? Math2.Beta(a, b) - fract : fract;
                }

                if (a >= 15) {
                    fract = IbetaLargeASmallBSeries(a, b, x, y, normalised, 1);
                    return (invert) ? Math2.Beta(a, b) - fract : fract;
                }

                // Sidestep to improve errors:
                double prefix = RisingFactorialRatio(a + b, a, 20);
                fract = IbetaAStep(a, b, x, y, 20, normalised);
                fract += IbetaLargeASmallBSeries(a + 20, b, x, y, normalised, prefix);
                return (invert) ? Math2.Beta(a, b) - fract : fract;
            }


            // Both a,b >= 1:
            double lambda = (a < b) ? a - (a + b) * x : (a + b) * y - b;
            if (lambda < 0) {
                Utility.Swap(ref a, ref b);
                Utility.Swap(ref x, ref y);
                invert = !invert;
            }

            if (b >= 40) {
                fract = Fraction2(a, b, x, y, normalised);
                return invert ? Math2.Beta(a, b) - fract : fract;
            }

            /*
            if( IsInteger(a) && IsInteger(b)) {
                // relate to the binomial distribution and use a finite sum:
                double k = a - 1;
                double n = b + k;
                fract = binomial_ccdf(n, k, x, y);
                if ( !invert)
                    return Math2.Beta(a, b) * fract;
                return Math2.Beta(a, b) * (1.0 - fract); 
            }
            */

            if (b * x <= 0.7) {
                fract = Series1(a, b, x, y, false);
                return (invert) ? Math2.Beta(a, b) - fract : fract;
            }

            if (a > 15) {
                // sidestep so we can use the series representation:
                int n = (int)Math.Floor(b);
                if (n == b)
                    --n;
                double bbar = b - n;
                double prefix = RisingFactorialRatio(a + bbar, bbar, n);
                fract = IbetaAStep(bbar, a, y, x, n, normalised);
                fract += IbetaLargeASmallBSeries(a, bbar, x, y, normalised, 1);
                fract /= prefix;
                return invert ? Math2.Beta(a, b) - fract : fract;
            }


            fract = Fraction2(a, b, x, y, normalised);
            return invert ? Math2.Beta(a, b) - fract : fract;

        }


        // Computes the difference between ibeta(a,b,x) and ibeta(a+k,b,x):
        static double IbetaAStep(double a, double b, double x, double y, int k, bool normalised)
        {

            double sum = 0;
            double term = 1;

            // series summation from 0 to k-1:
            for (int i = 0; i < k - 1; ++i) {
                term *= (a + b + i) * x / (a + i + 1);
                sum += term;
            }

            sum += 1; // add sum afterwards to maintain precision

            Debug.Assert(!double.IsInfinity(sum), "Overflow: a = " + a + "; b = " + b);


            double result = PowerTerms(a, b, x, y, normalised, sum) / a;
            return result;
        }

        /// <summary>
        /// Calculates the ratio of two rising factorial for small k
        /// <para>(a)_k/(b)_k</para> 
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="k"></param>
        /// <returns></returns>
        static double RisingFactorialRatio(double a, double b, int k)
        {
            //
            // This function is only needed for the non-regular incomplete beta,
            // it computes the delta in:
            // beta(a,b,x) = prefix + delta * beta(a+k,b,x)
            // it is currently only called for small k.
            //

            // calculate:
            // (a)(a+1)(a+2)...(a+k-1)
            // _______________________
            // (b)(b+1)(b+2)...(b+k-1)

            // This is only called with small k, for large k
            // it is grossly inefficient, do not use outside it's
            // intended purpose!!!
            if (k == 0)
                return 1;
            double result = 1;
            for (int i = 0; i < k; ++i)
                result *= (a + i) / (b + i);
            return result;
        }

        static double IbetaLargeASmallBSeries(double a, double b, double x, double y, bool normalised, double mult)
        {
            Debug.Assert(a >= 15, "Requires a >= 15: a= " + a);
            Debug.Assert(b >= 0 && b <= 1, "Requires b >= 0 && b <= 1: b= " + b);

            Debug.Assert(mult >= 0);

            //
            // Routine for a > 15, b < 1
            //
            // Begin by figuring out how large our table of Pn's should be,
            // quoted accuracies are "guestimates" based on empiracal observation.
            // Note that the table size should never exceed the size of our
            // tables of factorials.
            //

            const int Pn_size = 30; // 16-20 digit accuracy  


            //
            // This is DiDonato and Morris's BGRAT routine, see Eq's 9 through 9.6.
            //
            // Some values we'll need later, these are Eq 9.1:
            //
            double bm1 = b - 1;
            double t = a + bm1 / 2;
            double lx = (y < 0.35) ? Math2.Log1p(-y) : Math.Log(x);
            double u = -t * lx;
            // and from from 9.2:

            double prefix;
            double j; // = Math2.GammaQ(b, u) / IgammaPrefixRegularized(b, u);

            if (-u > DoubleLimits.MinLogValue) {
                double xpt = (y < 0.35) ? Math.Exp(-u) : Math.Pow(x, t);

                // h = u^b * e^-u
                double h = Math.Pow(u, b) * xpt; //IgammaPrefix(b, u);
                j = Math2.Tgamma(b, u) / h;

                prefix = Math.Pow(-lx, b) * xpt;
                if (normalised)
                    prefix /= Math2.Beta(a, b);
                prefix *= mult;


            } else {
                // if we use the previous equation, we'll end up with denorm/denorm or 0/0
                // So, use the asymptotic expansion of Γ(a, x) and cancel terms
                // Γ(b, u) = e^-u * u^(a-1) * TgammaAsymSeries(b, u)  
                j = _Igamma.Asym_SeriesLargeZ(b, u) / u;

                // use logs, though we'll probably get zero
                double lVal = -u + Math.Log(Math.Pow(-lx, b) * mult);
                if (normalised)
                    lVal -= Math2.LogBeta(a, b);
                prefix = Math.Exp(lVal);

            }

            if (prefix == 0)
                return 0;

            //
            // now we need the quantity Pn, unfortunatately this is computed
            // recursively, and requires a full history of all the previous values
            // so no choice but to declare a big table and hope it's big enough...
            //
            double[] p = new double[Pn_size];
            p[0] = 1;  // see 9.3.


            //
            // Now we can start to pull things together and evaluate the sum in Eq 9:
            //
            double sum = j;  // Value at N = 0
            // some variables we'll need:
            uint tnp1 = 1; // 2*N+1
            double lx2 = lx / 2;
            lx2 *= lx2;
            double lxp = 1;
            double t4 = 4 * t * t;
            double b2n = b;

            for (int n = 1; n < p.Length; ++n) {
                /*
                // debugging code, enable this if you want to determine whether
                // the table of Pn's is large enough...
                //
                static int max_count = 2;
                if(n > max_count)
                {
                    max_count = n;
                    Debug.WriteLine("Max iterations in BGRAT was {0}", n);
                }
                */
                //
                // begin by evaluating the next Pn from Eq 9.4:
                //
                tnp1 += 2;
                p[n] = 0;
                double mbn = b - n;
                int tmp1 = 3;
                for (int m = 1; m < n; ++m) {
                    mbn = m * b - n;
                    p[n] += mbn * p[n - m] / Math2.FactorialTable[tmp1];
                    tmp1 += 2;
                }
                p[n] /= n;
                p[n] += bm1 / Math2.FactorialTable[tnp1];
                //
                // Now we want Jn from Jn-1 using Eq 9.6:
                //
                j = (b2n * (b2n + 1) * j + (u + b2n + 1) * lxp) / t4;
                lxp *= lx2;
                b2n += 2;
                //
                // pull it together with Eq 9:
                //
                double r = p[n] * j;

                double prevSum = sum;
                sum += r;
                if (prevSum == sum)
                    break;
            }

            return prefix * sum;

        }

#if false
        // unused function left for future reference

        //
        // For integer arguments we can relate the incomplete beta to the
        // complement of the binomial distribution cdf and use this finite sum.
        //
        static double binomial_ccdf(double n, double k, double x, double y)
        {
            Debug.Assert( IsInteger(n) && ( (n >= int.MinValue && n <= int.MaxValue )) ); 
    
           double result = Math.Pow(x, n);
           double term = result;
           for(int i = (int)(n - 1); i > k; --i)
           {
              term *= ((i + 1) * y) / ((n - i) * x);
              result += term;
           }

           return result;
        }
#endif

        // 
        /// <summary>
        /// Calculates Ix(a,b) for large a and b (i.e &gt;= 15) using an asymptotic expansion
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="lambda"></param>
        /// <returns></returns>
        static double LargeABAsymptotic(double a, double b, double lambda)
        {
            const double eps = 10 * DoubleLimits.MachineEpsilon; // relative tolerance

            Debug.Assert(a >= 15 && b >= 15, "Requires a >= 15 && b >= 15");
            Debug.Assert(lambda >= 0); // lambda  = (a + b)*(1-x) - b

            // This is BASYM from TOMS708. 
            // See: Armido Didonato, Alfred Morris, Algorithm 708: 
            //      Significant Digit Computation of the Incomplete Beta Function Ratios,
            //      ACM Transactions on Mathematical Software,
            //      Volume 18, Number 3, 1992, pages 360-373.


            //
            //  num is the maximum value that n can take in the do loop
            //  ending at statement 50. it is required that num be even.
            //  the arrays a0, b0, c, d have dimension num + 2.
            //

            const double e0 = 2.0 * Constants.RecipSqrtPI; // 2/Sqrt(PI)
            const double e1 = 0.5 * Constants.RecipSqrt2; // 0.5/Sqrt(2)
            const int num = 20;

            double[] a0 = new double[num + 2];
            double[] b0 = new double[num + 2];
            double[] c = new double[num + 2];
            double[] d = new double[num + 2];

            double h, r0, r1, w0;

            if (a < b) {
                h = a / b;
                r0 = 1.0 / (1.0 + h);
                r1 = (b - a) / b;
                w0 = 1.0 / Math.Sqrt(a * (1.0 + h));
            } else {
                h = b / a;
                r0 = 1.0 / (1.0 + h);
                r1 = (b - a) / a;
                w0 = 1.0 / Math.Sqrt(b * (1.0 + h));
            }

            //double f = -(a * Math2.Log1pmx( -lambda/a ) + b * Math2.Log1pmx(lambda/b));
            double f = -(a * Math2.Log1p(-lambda / a) + b * Math2.Log1p(lambda / b));
            double t = Math.Exp(-f);

            if (t == 0)
                return 0;

            double z0 = Math.Sqrt(f);
            double z = 0.5 * (z0 / e1);
            double z2 = f + f;

            a0[1] = (2.0 / 3.0) * r1;
            c[1] = -0.5 * a0[1];
            d[1] = -c[1];

            // double j0 = ( 0.5 / e0 ) * Math.Exp(z0*z0)*Math2.Erfc(z0);
            double j0 = (0.5 / e0) * Math2.Erfcx(z0);
            double j1 = e1;
            double sum2 = j0 + d[1] * w0 * j1;

            double s = 1.0;
            double h2 = h * h;
            double hn = 1.0;
            double w = w0;
            double znm1 = z;
            double zn = z2;

            for (int n = 2; n <= num; n += 2) {

                hn = h2 * hn;
                a0[n] = 2.0 * r0 * (1.0 + h * hn) / (n + 2.0);
                int np1 = n + 1;
                s = s + hn;
                a0[np1] = 2.0 * r1 * s / (n + 3.0);

                for (int i = n; i <= np1; i++) {

                    double r = -0.5 * (i + 1.0);
                    b0[1] = r * a0[1];

                    for (int m = 2; m <= i; m++) {

                        double bsum = 0.0;
                        int mm1 = m - 1;

                        for (int j = 1; j <= mm1; j++) {
                            int mmj = m - j;
                            bsum += (j * r - mmj) * a0[j] * b0[mmj];
                        }

                        b0[m] = r * a0[m] + bsum / m;

                    }

                    c[i] = b0[i] / (i + 1.0);

                    double dsum = 0.0;
                    int im1 = i - 1;

                    for (int j = 1; j <= im1; j++)
                        dsum += d[i - j] * c[j];


                    d[i] = -(dsum + c[i]);

                }

                j0 = e1 * znm1 + (n - 1.0) * j0;
                j1 = e1 * zn + n * j1;
                znm1 = z2 * znm1;
                zn = z2 * zn;
                w = w0 * w;
                double t0 = d[n] * w * j0;
                w = w0 * w;
                double t1 = d[np1] * w * j1;
                sum2 += (t0 + t1);

                if ((Math.Abs(t0) + Math.Abs(t1)) <= eps * sum2)
                    break;

            }


            double u = StirlingGamma.GammaSeries(a + b)/(StirlingGamma.GammaSeries(a) * StirlingGamma.GammaSeries(b)); //Math.Exp(-Bcorr(a, b));
            // below is the equivalent:
            // double p = a/(a+b);
            // double q = b/(a+b);
            // return Math.Sqrt(2*Math.PI*(a+b)/(a*b))*_Beta.PowerTerms(a, b, p, q, true);

            return e0 * t * u * sum2;

        }

        /// <summary>
        /// Implementation for the Regularized Beta cases Ix(a,b) and 1-Ix(a,b)
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="x"></param>
        /// <param name="invert"></param>
        /// <returns></returns>
        public static double BetaRegularizedImp(double a, double b, double x, bool invert)
        {
            Debug.Assert(a >= 0);
            Debug.Assert(b >= 0);
            Debug.Assert(!(a == 0 && b == 0));
            Debug.Assert(x >= 0 && x <= 1);


            const bool normalised = true;
            double fract;
            double y = 1 - x;

            if (a <= 1 && b <= 1) {
                // Both a,b <= 1:

                if (x > 0.5) {
                    Utility.Swap(ref a, ref b);
                    Utility.Swap(ref x, ref y);
                    invert = !invert;
                }

                if (a >= Math.Min(0.2, b) || (Math.Pow(x, a) <= 0.9)) {
                    fract = Series1(a, b, x, y, true);
                    return (invert) ? 1 - fract : fract;
                }

                Utility.Swap(ref a, ref b);
                Utility.Swap(ref x, ref y);
                invert = !invert;
                if (y >= 0.3) {
                    fract = Series1(a, b, x, y, true);
                    return (invert) ? 1 - fract : fract;
                }

                // Sidestep on a, and then use the series representation:
                fract = IbetaAStep(a, b, x, y, 20, normalised);
                fract += IbetaLargeASmallBSeries(a + 20, b, x, y, normalised, 1);
                return (invert) ? 1 - fract : fract;
            }

            if (a <= 1 || b <= 1) {


                if (x > 0.5) {
                    Utility.Swap(ref a, ref b);
                    Utility.Swap(ref x, ref y);
                    invert = !invert;
                }

                // One of a, b < 1 only:
                if ((b <= 1) || ((x < 0.1) && (Math.Pow(b * x, a) <= 0.7))) {
                    fract = Series1(a, b, x, y, true);
                    return (invert) ? 1 - fract : fract;
                }

                Utility.Swap(ref a, ref b);
                Utility.Swap(ref x, ref y);
                invert = !invert;

                if (y >= 0.3) {
                    fract = Series1(a, b, x, y, true);
                    return (invert) ? 1 - fract : fract;
                }

                if (a >= 15) {
                    fract = IbetaLargeASmallBSeries(a, b, x, y, normalised, 1);
                    return (invert) ? 1 - fract : fract;
                }

                // Sidestep to improve errors:
                fract = IbetaAStep(a, b, x, y, 20, normalised);
                fract += IbetaLargeASmallBSeries(a + 20, b, x, y, normalised, 1);
                return (invert) ? 1 - fract : fract;
            }


            // Both a,b >= 1:
            // lambda < 0 
            //      => x > a/(a+b) for a <= b 
            //      or y < b/(a+b) for a > b
            double lambda = (a <= b) ? a - (a + b) * x : (a + b) * y - b;
            if (lambda < 0) {
                Utility.Swap(ref a, ref b);
                Utility.Swap(ref x, ref y);
                lambda = Math.Abs(lambda);
                invert = !invert;
            }


            if (b < 40) {

                if (b * x <= 0.7) {
                    fract = Series1(a, b, x, y, true);
                    return (invert) ? 1 - fract : fract;
                }

                /*
                if(IsInteger(a) && IsInteger(b)) {
                    // relate to the binomial distribution and use a finite sum:
                    double k = a - 1;
                    double n = b + k;
                    fract = binomial_ccdf(n, k, x, y);
                    return invert ? 1 - fract : fract;
                }
                */

                int n = (int)b;
                double bbar = b - n;

                if (bbar == 0.0) {
                    n = n - 1;
                    bbar = 1.0;
                }

                fract = IbetaAStep(bbar, a, y, x, n, normalised);
                if (x <= 0.7) {
                    fract += Series1(a, bbar, x, y, true);
                    return (invert) ? 1 - fract : fract;
                }

                if (a <= 15) {
                    n = 20;
                    fract += IbetaAStep(a, bbar, x, y, n, normalised);
                    a += n;
                }

                fract += IbetaLargeASmallBSeries(a, bbar, x, y, normalised, 1);
                return invert ? 1 - fract : fract;
            }

            // below this point,  a >= 1 && b >= 40 

            if (a <= b) {
                if (a <= 100.0 || lambda > 0.03 * a) { // x < 0.97 * a/(a+b)
                    fract = Fraction2(a, b, x, y, normalised);
                    return invert ? 1 - fract : fract;
                }
            } else {
                if (b <= 100.0 || lambda > 0.03 * b) { // y > 1.03 * b/(a+b)
                    fract = Fraction2(a, b, x, y, normalised);
                    return invert ? 1 - fract : fract;
                }
            }

#if EXTRA_DEBUG
            Debug.WriteLine("Using Asymp: a = {0}; b = {1}; lambda = {2}", a, b, lambda);
#endif

            fract = LargeABAsymptotic(a, b, lambda);
            return invert ? 1 - fract : fract;

        }

    }





    public static partial class Math2
    {

        /// <summary>
        /// Returns the partial derivative of IBeta(a,b,x) with respect to <paramref name="x"/>
        /// <para>IbetaDerivative(a, b, x) = ∂( I<sub>x</sub>(a,b) )/∂x = (1-x)^(b-1)*x^(a-1)/B(a,b)</para>
        /// </summary>
        /// <param name="a">Requires a ≥ 0 and not both a and b are zero</param>
        /// <param name="b">Requires b ≥ 0 and not both a and b are zero</param>
        /// <param name="x">Defined for 0 ≤ x ≤ 1</param>
        public static double IbetaDerivative(double a, double b, double x)
        {

            if (!(a > 0)
            || !(b > 0)
            || !(x >= 0 && x <= 1)) {
                Policies.ReportDomainError("IbetaDerivative(a: {0}, b: {1}, x: {2}): Requires a,b >= 0; x in [0,1]", a, b, x);
                return double.NaN;
            }

            //
            // Now the corner cases:
            //
            if (x == 0) {
                if (a > 1)
                    return 0;
                if (a == 1)
                    return 1 / Math2.Beta(a, b);

                return double.PositiveInfinity;
            }
            if (x == 1) {
                if (b > 1)
                    return 0;
                if (b == 1)
                    return 1 / Math2.Beta(a, b);

                return double.PositiveInfinity;
            }

            // handle denorm x
            if (x < DoubleLimits.MinNormalValue) {
                double logValue = a * Math.Log(x) + b * Math2.Log1p(-x) - Math.Log(x * (1 - x));
                logValue -= Math2.LogBeta(a, b);
                double result = Math.Exp(logValue);
#if EXTRA_DEBUG
                Debug.WriteLine("IbetaDerivative(a: {0}, b: {1}, x: {2}): Denorm = {3}", a, b, x, result);
#endif
                return result;
            }

            //
            // Now the regular cases:
            //
            double mult = (1 - x) * x;
            double f1 = _Ibeta.PowerTerms(a, b, x, 1 - x, true, 1 / mult);

            return f1;
        }



        /// <summary>
        /// Returns the Incomplete Beta function
        /// <para>B<sub>x</sub>(a, b) = ∫ (t^(a-1) * (1-t)^(b-1)) dt, t = {0,x}</para>
        /// </summary>
        /// <param name="a">Requires a ≥ 0 and not both a and b are zero</param>
        /// <param name="b">Requires b ≥ 0 and not both a and b are zero</param>
        /// <param name="x">0 ≤ x ≤ 1</param>
        public static double Beta(double a, double b, double x)
        {
            if ((!(a > 0) || double.IsInfinity(a))
            || (!(b > 0) || double.IsInfinity(b))
            || !(x >= 0 && x <= 1)) {
                Policies.ReportDomainError("Beta(a: {0}, b: {1}, x: {2}): Requires finite a,b > 0; x in [0,1]", a, b, x);
                return double.NaN;
            }

            // special values
            if (x == 0)
                return 0;
            if (x == 1)
                return Math2.Beta(a, b);

            return _Ibeta.BetaImp(a, b, x, false);

        }

        /// <summary>
        /// Returns the Incomplete Beta Complement
        /// <para>Betac(a, b, x) = B(a,b) - B<sub>x</sub>(a,b)</para>
        /// </summary>
        /// <param name="a">Requires a ≥ 0 and not (a == b == 0)</param>
        /// <param name="b">Requires b ≥ 0 and not (a == b == 0)</param>
        /// <param name="x">0 ≤ x ≤ 1</param>
        public static double Betac(double a, double b, double x)
        {
            if ((!(a > 0) || double.IsInfinity(a))
            || (!(b > 0) || double.IsInfinity(b))
            || !(x >= 0 && x <= 1)) {
                Policies.ReportDomainError("Betac(a: {0}, b: {1}, x: {2}): Requires finite a,b > 0; x in [0,1]", a, b, x);
                return double.NaN;
            }

            // special values
            if (x == 0)
                return Math2.Beta(a, b);
            if (x == 1)
                return 0;

            return _Ibeta.BetaImp(a, b, x, true);
        }


        /// <summary>
        /// Returns the Regularized Incomplete Beta function 
        /// <para>I<sub>x</sub>(a,b) = B<sub>x</sub>(a,b)/B(a,b)</para>
        /// </summary>
        /// <param name="a">Requires a ≥ 0 and not both a and b are zero</param>
        /// <param name="b">Requires b ≥ 0 and not both a and b are zero</param>
        /// <param name="x">0 ≤ x ≤ 1</param>
        public static double Ibeta(double a, double b, double x)
        {
            if ((!(a >= 0) || double.IsInfinity(a))
            || (!(b >= 0) || double.IsInfinity(b))
            || !(x >= 0 && x <= 1)) {
                Policies.ReportDomainError("Ibeta(a: {0}, b: {1}, x: {2}): Requires finite a,b >= 0; x in [0,1]", a, b, x);
                return double.NaN;
            }
            if (a == 0 && b == 0) {
                Policies.ReportDomainError("Ibeta(a: {0}, b: {1}, x: {2}): Requires that a and b cannot both be zero", a, b, x);
                return double.NaN;
            }

            // special cases:
            if (a == 0)
                return 1;

            if (b == 0)
                return 0;

            if (x == 0)
                return 0;

            if (x == 1)
                return 1;

            return _Ibeta.BetaRegularizedImp(a, b, x, false);
        }

        /// <summary>
        /// Returns the Regularized Incomplete Beta Complement = 1 - I<sub>x</sub>(a,b)
        /// </summary>
        /// <param name="a">Requires a ≥ 0 and not both a and b are zero</param>
        /// <param name="b">Requires b ≥ 0 and not both a and b are zero</param>
        /// <param name="x">Defined for 0 ≤ x ≤ 1</param>
        public static double Ibetac(double a, double b, double x)
        {
            if ((!(a >= 0) || double.IsInfinity(a))
            || (!(b >= 0) || double.IsInfinity(b))
            || !(x >= 0 && x <= 1)) {
                Policies.ReportDomainError("Ibetac(a: {0}, b: {1}, x: {2}): Requires finite a,b >= 0; x in [0,1]", a, b, x);
                return double.NaN;
            }
            if (a == 0 && b == 0) {
                Policies.ReportDomainError("Ibetac(a: {0}, b: {1}, x: {2}): Requires that a and b cannot both be zero", a, b, x);
                return double.NaN;
            }

            // extend to a few very special cases:
            if (a == 0)
                return 0;

            if (b == 0)
                return 1;

            if (x == 0)
                return 1;

            if (x == 1)
                return 0;


            return _Ibeta.BetaRegularizedImp(a, b, x, true);
        }

    }



} // namespace 
