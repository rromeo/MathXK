//  Copyright (c) 2014 Rocco Romeo
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

//  Ported to .NET based on the original:
//      Copyright (C) 2013 Nikhar Agrawal, Christopher Kormanyos, John Maddock, Paul Bristow, Boost Software License, Version 1.0

//#define EXTRA_DEBUG

using System;
using System.Collections.Generic;
using System.Diagnostics;

namespace MathXK
{
    internal static class _Polygamma
    {


        /// <summary>
        /// Returns polygamma using the asymptotic series
        /// <para>Ψ(n)(z) = ((-1)^(n-1)*(n-1)!*(n+2*z))/(2*z^(n+1)) - (-1)^n * Σ[((2k+n-1)!/((2k)!*z^(2k+n)))*BernoulliB[2k], {k, 1, Infinity}]</para>
        /// </summary>
        /// <param name="n"></param>
        /// <param name="x"></param>
        /// <returns></returns>
        static double AtInfinityPlus(int n, double x) 
        {
            Debug.Assert(n > 0 && x >= 0, $"AtInfinityPlus(n: {n}, x:{x}): Requires n > 0 && x >= 0");

            // See http://functions.wolfram.com/GammaBetaErf/PolyGamma2/06/02/0001/
            //
            // sum       == current value of accumulated sum.
            // term      == value of current term to be added to sum.
            // part_term == value of current term excluding the Bernoulli number part
            //
            if (n + x == x) {
                // x is crazy large, just concentrate on the first part of the expression and use logs:
                if (n == 1) return 1 / x;
                double nlx = n * Math.Log(x);
                double sign = Math2.IsOdd(n) ? 1 : -1;
                if ((nlx < DoubleLimits.MaxLogValue) && (n < Math2.MaxFactorialIndex))
                    return sign * Math2.Factorial(n - 1) * Math.Pow(x, -n);
                else
                    return sign * Math.Exp(Math2.Lgamma(n) - n * Math.Log(x));
            }
            double term, sum, part_term;
            double x_squared = x * x;
            //
            // Start by setting part_term to:
            //
            // (n-1)! / x^(n+1)
            //
            // which is common to both the first term of the series (with k = 1)
            // and to the leading part.  
            // We can then get to the leading term by:
            //
            // part_term * (n + 2 * x) / 2
            //
            // and to the first term in the series 
            // (excluding the Bernoulli number) by:
            //
            // part_term n * (n + 1) / (2x)
            //
            // If either the factorial would overflow,
            // or the power term underflows, this just gets set to 0 and then we
            // know that we have to use logs for the initial terms:
            //
            double np1_lnx = (n + 1) * Math.Log(x);
            if ((n > Math2.MaxFactorialIndex) || (np1_lnx > DoubleLimits.MaxLogValue)) {
                // Either n is very large, or the power term underflows,
                // set the initial values of part_term, term and sum via logs:
                part_term = Math2.Lgamma(n) - np1_lnx;
                sum = Math.Exp(part_term + Math.Log(n + 2 * x) - Constants.Ln2);
                part_term += Math.Log(((double)n) * (n + 1)) - Constants.Ln2 - Math.Log(x);
                part_term = Math.Exp(part_term);
            } else {
                part_term = Math2.Factorial(n - 1) * Math.Pow(x, -n - 1);
                sum = part_term * (n + 2 * x) / 2;
                part_term *= ((double)n) * (n + 1) / 2.0;
                part_term /= x;
            }
            //
            // If the leading term is 0, so is the result:
            //
            if (sum == 0)
                return sum;

            for (int k = 1; k < Math2.Bernoulli2nTable.Length; ) {
                term = part_term * Math2.Bernoulli2nTable[k];
                sum += term;
                //
                // Normal termination condition:
                //
                if (Math.Abs(term / sum) < DoubleLimits.MachineEpsilon)
                    return Math2.IsOdd(n - 1) ? -sum : sum;

                //
                // Increment our counter, and move part_term on to the next value:
                //
                ++k;
                part_term *= ((double)(n + 2 * k - 2)) * (n - 1 + 2 * k);
                part_term /= (2 * k - 1) * 2 * k;
                part_term /= x_squared;
            }

            Policies.ReportConvergenceError("Polygamma(n: {0}, x: {1}): series did not converge after {2} iterations", n, x, Policies.MaxSeriesIterations);
            return double.NaN;
        }

        /// <summary>
        /// Returns the Polygamma function using recurrence
        /// <para>Ψ(n)(z-m) = Ψ(n)(z) - (-1)^n * n! * Σ[1/(z - k)^(n + 1), {k, 1, m}]</para>
        /// </summary>
        /// <param name="n"></param>
        /// <param name="x"></param>
        /// <param name="nIterations"></param>
        /// <returns></returns>
        /// <seealso href="http://functions.wolfram.com/GammaBetaErf/PolyGamma2/16/01/01/0017/"/>
        static double ForwardRecurrence(int n, double x, int nIterations)
        {
            Debug.Assert(n > 0 && x >= 0 && nIterations >= 0 && nIterations < Policies.MaxSeriesIterations, $"ForwardRecurrence(n: {n}, x:{x}, nIterations: {nIterations}): Requires n > 0 && x >= 0 && nIterations in [0, Policies.MaxSeriesIterations]");

            double sum = 0;
            int np1 = n + 1;

            // Forward recursion to larger x, need to check for overflow first though:
            if (Math.Log(x + nIterations) * np1 < DoubleLimits.MaxLogValue) {
                for (int k = 0; k < nIterations; k++) {
                    sum += 1 / Math.Pow(x + k, np1);
                }
                sum *= Math2.Factorial(n);
            } else {
                double lgamma = Math2.Lgamma(np1);
                for (int k = 0; k < nIterations; k++) {
                    double log_term = lgamma - Math.Log(x + k) * np1;
                    sum += Math.Exp(log_term);
                }
            }
            if (Math2.IsOdd(n - 1))
                sum = -sum;

            return sum;
        }


        /// <summary>
        /// Returns Polygamma using the small series when x is near zero
        /// </summary>
        /// <param name="n"></param>
        /// <param name="x"></param>
        /// <returns></returns>
        static double SmallSeriesNearZero(int n, double x)
        {
            Debug.Assert(n >= 1, "Series not for digamma");

            // Use a series expansion for x near zero which uses poly_gamma(m, 1) which,
            // in turn, uses the Riemann zeta function for integer arguments.
            // http://functions.wolfram.com/GammaBetaErf/PolyGamma2/06/01/03/01/02/

            double n_fact = Math2.Factorial(n);
            if (double.IsInfinity(n_fact)) {
                Debug.Assert(Math.Abs(x) < 1, "Assuming x < 1 when n > MaxFactorial");
#if EXTRA_DEBUG
                Debug.WriteLine("Polygamma({0}, {1}): {2} iterations", n, x, 0);
#endif

                return (Math2.IsOdd(n) || (x < 0)) ? double.PositiveInfinity : double.NegativeInfinity;
            }

            // psi(n, z) = ((-1)^(n-1) n!)/z^(n + 1) + Sum[(PolyGamma[n + k, 1]/k!) z^k, {k, 0, Infinity}]
            // psi(n, 1) = (-1)^(n+1) * n! * Zeta(n+1) 
            // so,
            // psi(n, z) = ((-1)^(n-1) n!)/z^(n + 1) + Sum[(((-1)^(n+k+1) * (n+k)! * Zeta(n+k+1))/k!) z^k, {k, 0, Infinity}]

            // First term in the series is necessarily < zeta(2) < 2, so
            // ignore the sum if it will have no effect on the result anyway:
            double mult = Math2.IsOdd(n - 1) ? -1 : 1;
            double sum = mult / Math.Pow(x, n + 1);
            if (Math.Abs(sum) > 2 / DoubleLimits.MachineEpsilon) {
#if EXTRA_DEBUG
                Debug.WriteLine("Polygamma({0}, {1}): {2} iterations", n, x, 0);
#endif
                return sum * n_fact;
            }

            double term = mult * Math2.Zeta(n + 1);
            if (Math.Abs(term / sum) < DoubleLimits.MachineEpsilon) {
#if EXTRA_DEBUG
                Debug.WriteLine("Polygamma({0}, {1}): {2} iterations", n, x, 1);
#endif
                return n_fact * sum;
            }
            sum += term;

            for (int k = 1; k < Policies.MaxSeriesIterations; k++) {

                mult *= -(n + k) * x / k;
                term = mult * Math2.Zeta(k + n + 1);

                if (Math.Abs(term) <= DoubleLimits.MachineEpsilon * Math.Abs(sum) ) {
#if EXTRA_DEBUG
                    Debug.WriteLine("Polygamma({0}, {1}): {2} iterations", n, x, k + 1);
#endif
                    return n_fact * sum;
                }

                sum += term;
            }

            Policies.ReportConvergenceError("Polygamma(n: {0}, x: {1}): series did not converge after {2} iterations", n, x, Policies.MaxSeriesIterations);
            return double.NaN;
        }


        static readonly double[][] _P = { 
            new double[] { -1 },
            new double[] { 2 },
            new double[] { -2, -4 },
            new double[] { 16, 8 },
            new double[] { -16, -88, -16 },
            new double[] { 272, 416, 32 },
            new double[] { -272, -2880, -1824, -64 },
            new double[] { 7936, 24576, 7680, 128 },
            new double[] { -7936, -137216, -185856, -31616, -256 },
            new double[] { 353792, 1841152, 1304832, 128512, 512 },
            new double[] { -353792, -9061376, -21253376, -8728576, -518656, -1024},
            new double[] { 22368256, 175627264, 222398464, 56520704, 2084864, 2048 },
            new double[] { -22368256, -795300864, -2868264960, -2174832640, -357888000, -8361984, -4096 },
            new double[] { 1903757312, 21016670208, 41731645440, 20261765120, 2230947840, 33497088, 8192 },
            new double[] { -1903757312, -89702612992, -460858269696, -559148810240, -182172651520, -13754155008, -134094848, -16384 },
            new double[] { 209865342976, 3099269660672, 8885192097792, 7048869314560, 1594922762240, 84134068224, 536608768, 32768 },
            new double[] { -209865342976, -12655654469632, -87815735738368, -155964390375424, -84842998005760, -13684856848384, -511780323328, -2146926592, -65536 },
            new double[] { 29088885112832, 553753414467584, 2165206642589696, 2550316668551168, 985278548541440, 115620218667008, 3100738912256, 8588754944, 131072 },
            new double[] { -29088885112832, -2184860175433728, -19686087844429824, -48165109676113920, -39471306959486976, -11124607890751488, -965271355195392, -18733264797696, -34357248000, -262144 },
            new double[] { 4951498053124096, 118071834535526400, 603968063567560704, 990081991141490688, 584901762421358592, 122829335169859584, 7984436548730880, 112949304754176, 137433710592, 524288 },
        };

        /// <summary>
        /// Cache of the derivatives of Cot
        /// </summary>
        static List<double[]> _CotTableCache = new List<double[]>();


        static double PolyCotPi(int n, double x, double xc)
        {
            // Return n'th derivative of cot(pi*x) at x, these are simply
            // tabulated for up to n = 20, beyond that it is possible to
            // calculate coefficients as follows:
            //
            // The general form of each derivative is:
            //
            // pi^n * SUM{k=0, n} C[k,n] * cos^k(pi * x) * csc^(n+1)(pi * x)
            //
            // With constant C[0,1] = -1 and all other C[k,n] = 0;
            // Then for each k < n+1:
            // C[k-1, n+1]  -= k * C[k, n];
            // C[k+1, n+1]  += (k-n-1) * C[k, n];
            //
            // Note that there are many different ways of representing this derivative thanks to
            // the many trigomonetric identies available.  In particular, the sum of powers of
            // cosines could be replaced by a sum of cosine multiple angles, and indeed if you
            // plug the derivative into Mathematica this is the form it will give.  The two
            // forms are related via the Chebeshev polynomials of the first kind and
            // T_n(cos(x)) = cos(n x).  The polynomial form has the great advantage that
            // all the cosine terms are zero at half integer arguments - right where this
            // function has it's minumum - thus avoiding cancellation error in this region.
            //
            // And finally, since every other term in the polynomials is zero, we can save
            // space by only storing the non-zero terms.  This greatly complexifies
            // subscripting the tables in the calculation, but halves the storage space
            // (and complexity for that matter).
            //
            double s = Math.Abs(x) < Math.Abs(xc) ? Math2.SinPI(x) : Math2.SinPI(xc);
            double c = Math2.CosPI(x);
            double prefix = Math.Pow(Math.PI / s, n) / s;
            int index = n - 1;
            if (index < _P.Length) {
                double result = Polynomial.EvalEven(_P[index], c);
                if (Math2.IsOdd(index))
                    result *= c;  // First coeffient is order 1, and really an odd polynomial.

                return prefix * result;
            }

            if (index >= _CotTableCache.Count) {
                // We'll have to compute the coefficients up to n:
                if (_CotTableCache.Count == 0) {
                    double[] a = new double[1];
                    a[0] = -1;
                    _CotTableCache.Add(a);
                }

                for (int i = _CotTableCache.Count - 1; i < index; ++i) {
                    int offset = i & 1; // 1 if the first cos power is 0, otherwise 0.
                    int sin_order = i + 2;  // order of the sin term
                    int max_cos_order = sin_order - 1;  // largest order of the polynomial of cos terms
                    int max_columns = (max_cos_order - offset) / 2;  // How many entries there are in the current row.
                    int next_offset = (offset != 0) ? 0 : 1;
                    int next_max_columns = (max_cos_order + 1 - next_offset) / 2;  // How many entries there will be in the next row
                    double[] newRow = new double[next_max_columns + 1];

                    for (int column = 0; column <= max_columns; ++column) {
                        int cos_order = 2 * column + offset;  // order of the cosine term in entry "column"
                        Debug.Assert(column < _CotTableCache[i].Length);
                        Debug.Assert((cos_order + 1) / 2 < newRow.Length);
                        newRow[(cos_order + 1) / 2] += ((cos_order - sin_order) * _CotTableCache[i][column]) / (sin_order - 1);
                        if (cos_order != 0)
                            newRow[(cos_order - 1) / 2] += (-cos_order * _CotTableCache[i][column]) / (sin_order - 1);
                    }

                    _CotTableCache.Add(newRow);
                }

            }

            double sum = Polynomial.EvalEven(_CotTableCache[index], c);
            if (Math2.IsOdd(index))
                sum *= c;  // First coeffient is order 1, and really an odd polynomial.
            if (sum == 0)
                return sum;
            //
            // the remaining terms are computed using logs since the powers and factorials
            // get real large real quick:
            //
            if (s == 0)
                return sum * double.PositiveInfinity;
            double power_terms = n * Constants.LnPI;
            power_terms -= Math.Log(Math.Abs(s)) * (n + 1);
            power_terms += Math2.Lgamma(n);
            power_terms += Math.Log(Math.Abs(sum));

            return Math.Exp(power_terms) * ((s < 0) && Math2.IsOdd(n + 1) ? -1 : 1) * Math.Sign(sum);
        }



        public static double Imp(int n, double x)
        {
            Debug.Assert(n > 0);

            // Limit for use of small-x-series is chosen 
            // so that the series doesn't go too divergent in the first few terms

            double small_x_limit = Math.Min(5.0 / n, 0.25);
            if (x < 0) {
                if (Math2.IsInteger(x)) {
                    // Result is infinity if x is odd, and a pole error if x is even.
                    if (!Math2.IsOdd(x)) {
                        Policies.ReportPoleError("Polygamma(n: {0}, x: {1}): Expecting non-negative integer x", n, x);
                        return double.NaN;
                    }
                    return double.PositiveInfinity;
                }

                // if |x| < small_x_limit, compute the series directly
                if (x >= -small_x_limit)
                    return SmallSeriesNearZero(n, x);

                // We have tabulated the derivatives of cot(x) up to the 21st derivative, which
                // allows us to use: http://functions.wolfram.com/06.15.16.0001.01
                double z = 1 - x;
                double result = Imp(n, z) + Math.PI * PolyCotPi(n, z, x);
                return Math2.IsOdd(n) ? -result : result;

#if false
                // Recurrence currently not used but saved
                // Try http://functions.wolfram.com/06.15.16.0007.01

                Debug.Assert(x > -int.MaxValue);

                int m = (int)Math.Ceiling(-x);
                double zpm = x + m;
                double sum = 0;
                for (int k = 1; k <= m; ++k) {
                    sum += Math.Pow(zpm - k, -n - 1);
                }
                sum *= Math2.Factorial(n);
                if (Math2.IsOdd(n))
                    sum = -sum;
                return Imp(n, zpm) - sum;
#endif
            }

            if (x < small_x_limit) 
                return SmallSeriesNearZero(n, x);

            if (x == 0.5) {
                double result = Math2.Factorial(n) * Math2.Zeta(n + 1);
                if (!Math2.IsOdd(n))
                    result = -result;
                result *= Math2.Ldexp(1, n + 1) - 1;
                return result;
            }

            if (x == 1) {
                double result = Math2.Factorial(n) * Math2.Zeta(n + 1);
                if (!Math2.IsOdd(n))
                    result = -result;
                return result;
            }


            const int Digits = 16;

            // Minimum x for to use the asymptotic series
            double AsymMinX = Math.Floor(0.4 * Digits + 4.0 * n);
            if (x >= AsymMinX) 
                return AtInfinityPlus(n, x);

            int iterations = (int)(AsymMinX - x);

            double sum = ForwardRecurrence(n, x, iterations);

            return sum + AtInfinityPlus(n, x + iterations);
        }
    }

    public static partial class Math2
    {


        /// <summary>
        /// Returns the Polygamma function = d(n)/dx Lgamma(x)
        /// </summary>
        /// <param name="n">Order. Requires <paramref name="n"/> ≥ 0. Limited to <paramref name="n"/> ≤ 170</param>
        /// <param name="x">Argument</param>
        /// <returns></returns>
        public static double Polygamma(int n, double x)
        {
            if (n < 0) {
                Policies.ReportDomainError("Polygamma(n: {0}, x: {1}): Expecting n >= 0", n, x);
                return double.NaN;
            }
            if (n > Math2.MaxFactorialIndex ) {
                Policies.ReportNotImplementedError("Polygamma(n: {0}, x: {1}): Not implemented for n > {2}", n, x, Math2.MaxFactorialIndex);
                return double.NaN;
            }
            if (double.IsNaN(x)) {
                Policies.ReportDomainError("Polygamma(n: {0}, x: {1}): NaN not allowed", n, x);
                return double.NaN;
            }

            // Handle Digamma and Trigamma using polynomial approximations
            if (n == 0)
                return Math2.Digamma(x);
            if (n == 1)
                return Math2.Trigamma(x);

            if (double.IsInfinity(x)) {
                if (x < 0) {
                    Policies.ReportDomainError("Polygamma(n: {0}, x: {1}): Requires x != -Infinity", n, x);
                    return double.NaN;
                }
                return 0;
            }

            return _Polygamma.Imp(n, x);
        }


    }


}

