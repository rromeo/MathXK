//  Copyright (c) 2013 Rocco Romeo
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

//  Ported to .NET based on the original:
//      Copyright (c) John Maddock 2006, Boost Software License v1.0


// This inverts the incomplete
// beta functions Ibeta and Ibetac on the first parameters "a"
// and "b" using a generic root finding algorithm (TOMS Algorithm 748).

#if DEBUG
//#define EXTRA_DEBUG
#endif


using System;
using System.Diagnostics;
using MathXK.Roots;

namespace MathXK
{

    internal static class _IbetaInvAB
    {
        private const int UlpsTolerance = 64;
        private static readonly RootFinder.ToleranceFunc _tol = RootFinder.GetToleranceFunc(UlpsTolerance * DoubleLimits.MachineEpsilon, 0);

        static double SolveA(double b, double x, double p, double q, double guess, double step, double lower, double upper)
        {
            Func<double, double> f;
            if (p <= q) {
                f = a => Math2.Ibeta(a, b, x) - p;
            } else {
                f = a => Math2.Ibetac(a, b, x) - q;
            }

#if EXTRA_DEBUG
            Debug.WriteLine("IBetaInvA(guess =? {0}, b={1}, x={2}) == p({3}), q({4}); range = [{5}, {6}]", guess, b, x, p, q, lower, upper);
#endif

            RootResults r = RootFinder.Toms748Bracket(f, guess, step, FunctionShape.Unknown, lower, upper, _tol, Policies.MaxRootIterations);
            if (r == null) {
                Policies.ReportRootNotFoundError("Invalid parameter in root solver");
                return double.NaN;
            }
            if (!r.Success) {
                Policies.ReportRootNotFoundError("Root not found after {0} iterations", r.Iterations);
                return double.NaN;
            }

#if EXTRA_DEBUG
            Debug.WriteLine("Result = {0};  Toms748 iterations: {1}",r.SolutionX, r.Iterations);
#endif

            return r.SolutionX;
        }

        public static double SolveB(double a, double x, double p, double q, double guess, double step, double lower, double upper)
        {
            Func<double, double> f;
            if (p <= q) {
                f = b => Math2.Ibeta(a, b, x) - p;
            } else {
                f = b => Math2.Ibetac(a, b, x) - q;
            }

#if EXTRA_DEBUG
            Debug.WriteLine("IBetaInvB(a={0}, guess?={1}, x={2}) == p({3}), q({4}); range = [{5}, {6}]", a, guess, x, p, q, lower, upper);
#endif

            RootResults r = RootFinder.Toms748Bracket(f, guess, step, FunctionShape.Unknown, lower, upper, _tol, Policies.MaxRootIterations);
            if (r == null) {
                Policies.ReportRootNotFoundError("Invalid parameter in root solver");
                return double.NaN;
            }
            if (!r.Success) {
                Policies.ReportRootNotFoundError("Root not found after {0} iterations", r.Iterations);
                return double.NaN;
            }

#if EXTRA_DEBUG
            Debug.WriteLine("Result = {0};  Toms748 iterations: {1}",r.SolutionX, r.Iterations);
#endif

            return r.SolutionX;
        }


        static double InverseNegativeBinomialCornishFisher(double n, double sf, double sfc, double p, double q)
        {

            // mean:
            double m = n * (sfc) / sf;
            double t = Math.Sqrt(n * (sfc));
            // standard deviation:
            double sigma = t / sf;
            // skewness
            double sk = (1 + sfc) / t;
            // kurtosis:
            double k = (6 - sf * (5 + sfc)) / (n * (sfc));
            // Get the inverse of a std normal distribution:
            double x = Math2.ErfcInv(p > q ? 2 * q : 2 * p) * Constants.Sqrt2;
            // Set the sign:
            if (p < 0.5)
                x = -x;
            double x2 = x * x;
            // w is correction term due to skewness
            double w = x + sk * (x2 - 1) / 6;
            //
            // Add on correction due to kurtosis.
            //
            if (n >= 10)
                w += k * x * (x2 - 3) / 24 + sk * sk * x * (2 * x2 - 5) / -36;

            w = m + sigma * w;
            if (w < DoubleLimits.MinNormalValue)
                return DoubleLimits.MinNormalValue;
            return w;
        }

        public static double IbetaInvAImp(double b, double x, double xc, double p, double q)
        {
            //
            // Special cases first:
            //
#if EXTRA_DEBUG
            Debug.WriteLine("b = {0}; x = {1};  p = {2}; q = {3}", b, x, p, q);
#endif
            if (p == 0)
                return double.PositiveInfinity; 
            if (q == 0)
                return 0; 


            // Set the upper and lower limits for the root finder
            const double lower = 0;
            const double upper = double.MaxValue;


            //
            // Now figure out a starting guess for what a may be, 
            // we'll start out with a value that'll put p or q
            // right bang in the middle of their range, the functions
            // are quite sensitive so we should need too many steps
            // to bracket the root from there:
            //
            double guess = 0;
            double factor = 5;
            //
            // Convert variables to parameters of a negative binomial distribution:
            //
            double n = b;

            double sf = xc;
            double sfc = x;
            double u = q;
            double v = p;

            // I_x(0, b) = 1
            // I_x(1, b) = 1-(1-x)^b

            if (u <= Math.Pow(sf, n)) {
                //
                // Result is less than 1, negative binomial approximation
                // is useless....
                //
                guess = Math.Min(((p < q) ? b * 2 : b / 2), 1.0);
            }
            if (n * n * n * u * sf > 0.005)
                guess = 1 + InverseNegativeBinomialCornishFisher(n, sf, sfc, u, v);

            if (guess < 10) {
                //
                // Negative binomial approximation not accurate in this area:
                //
                guess = Math.Min(((p < q) ? b * 2 : b / 2), 10.0);
            } else
                factor = (v < DoubleLimits.RootMachineEpsilon._2) ? 2 : (guess < 20 ? 1.2 : 1.1);


            double step = guess * (-factor);

            return SolveA(b, x, p, q, guess, step, lower, upper);
        }


        public static double IbetaInvBImp(double a, double x, double xc, double p, double q)
        {
            //
            // Special cases first:
            //
#if EXTRA_DEBUG
            Debug.WriteLine("a = {0}; x = {1};  p = {2}; q = {3}", a, x, p, q);
#endif

            if (p == 0)
                return 0; 
            if (q == 0)
                return double.PositiveInfinity; 


            // Set the upper and lower limits for the root finder
            const double lower = 0;
            const double upper = double.MaxValue;

            double guess = 0;
            double factor = 5;
            //
            // Convert variables to parameters of a negative binomial distribution:
            //
            double n = a;
            double sf = x;
            double sfc = xc;
            double u = p;
            double v = q;

            // I_x(a, 0) = 0
            // I_x(a, 1) = x^a
            if (u <= Math.Pow(sf, n)) {
                //
                // Result is less than 1, negative binomial approximation
                // is useless....
                //
                guess = Math.Min(((p >= q) ? a * 2 : a / 2), 1.0);

            }

            if (n * n * n * u * sf > 0.005)
                guess = 1 + InverseNegativeBinomialCornishFisher(n, sf, sfc, u, v);

            if (guess < 10) {
                //
                // Negative binomial approximation not accurate in this area:
                //
                guess = Math.Min(((p >= q) ? a * 2 : a / 2), 10.0);
            } else
                factor = (v < DoubleLimits.RootMachineEpsilon._2) ? 2 : (guess < 20 ? 1.2 : 1.1);

#if EXTRA_DEBUG
            Debug.WriteLine("guess = {0}", guess);
#endif

            double step = -factor * guess;
            return SolveB(a, x, p, q, guess, step, lower, upper);
        }
    }

    public static partial class Math2
    {

        /// <summary>
        /// Returns "a" such that: I<sub>x</sub>(a, b) == p
        /// </summary>
        /// <param name="b">Requires b &gt; 0</param>
        /// <param name="x">Requires 0 ≤ x ≤ 1</param>
        /// <param name="p">Requires 0 ≤ p ≤ 1</param>
        public static double IbetaInvA(double b, double x, double p)
        {
            if ((!(b > 0) || double.IsInfinity(b))
            || (!(x >= 0 && x <= 1))
            || (!(p >= 0 && p <= 1))) {
                Policies.ReportDomainError("IbetaInvA(b: {0}, x: {1}, p: {2}): Requires finite b > 0; x,p in [0,1]", b, x, p);
                return double.NaN;
            }

            if (p == 1)
                return 0; 
            if (p == 0)
                return double.PositiveInfinity; 

            // use I_x(a,b) = 1 - I_{1-x}(b, a)
            if (p > 0.5 && x > 0.5)
                return _IbetaInvAB.IbetaInvBImp(b, 1 - x, x, 1 - p, p);


            return _IbetaInvAB.IbetaInvAImp(b, x, 1 - x, p, 1 - p);
        }

        /// <summary>
        /// Returns "a" such that: 1 - I<sub>x</sub>(a, b) == q
        /// </summary>
        /// <param name="b">Requires b &gt; 0</param>
        /// <param name="x">Requires 0 ≤ x ≤ 1</param>
        /// <param name="q">Requires 0 ≤ q ≤ 1</param>
        public static double IbetacInvA(double b, double x, double q)
        {
            if ((!(b > 0) || double.IsInfinity(b))
            || (!(x >= 0 && x <= 1))
            || (!(q >= 0 && q <= 1))) {
                Policies.ReportDomainError("IbetacInvA(b: {0}, x: {1}, q: {2}): Requires finite b > 0; x,q in [0,1]", b, x, q);
                return double.NaN;
            }

            if (q == 0)
                return 0; 
            if (q == 1)
                return double.PositiveInfinity; 

            // use I_x(a,b) = 1 - I_{1-x}(b, a)
            if (q > 0.5 && x > 0.5)
                return _IbetaInvAB.IbetaInvBImp(b, 1 - x, x, q, 1 - q);

            return _IbetaInvAB.IbetaInvAImp(b, x, 1 - x, 1 - q, q);
        }

        /// <summary>
        /// Returns "b" such that: I<sub>x</sub>(a, b) == p
        /// </summary>
        /// <param name="a">Requires a &gt; 0</param>
        /// <param name="x">Requires 0 ≤ x ≤ 1</param>
        /// <param name="p">Requires 0 ≤ p ≤ 1</param>
        public static double IbetaInvB(double a, double x, double p)
        {
            if ((!(a > 0) || double.IsInfinity(a))
            || (!(x >= 0 && x <= 1))
            || (!(p >= 0 && p <= 1))) {
                Policies.ReportDomainError("IbetaInvB(a: {0}, x: {1}, p: {2}): Requires finite a > 0; x,p in [0,1]", a, x, p);
                return double.NaN;
            }

            if (p == 0)
                return 0; 
            if (p == 1)
                return Double.PositiveInfinity; 

            // use I_x(a,b) = 1 - I_{1-x}(b, a)
            if (p > 0.5 && x > 0.5)
                return _IbetaInvAB.IbetaInvAImp(a, 1 - x, x, 1 - p, p);

            return _IbetaInvAB.IbetaInvBImp(a, x, 1 - x, p, 1 - p);

        }

        /// <summary>
        /// Returns "b" such that: 1 - I<sub>x</sub>(a, b) == q
        /// </summary>
        /// <param name="a">Requires a ≥ 0</param>
        /// <param name="x">Requires 0 ≤ x ≤ 1</param>
        /// <param name="q">Requires 0 ≤ q ≤ 1</param>
        public static double IbetacInvB(double a, double x, double q)
        {
            if ((!(a > 0) || double.IsInfinity(a))
            || (!(x >= 0 && x <= 1))
            || (!(q >= 0 && q <= 1))) {
                Policies.ReportDomainError("IbetaInvB(a: {0}, x: {1}, q: {2}): Requires finite a > 0; x,q in [0,1]", a, x, q);
                return double.NaN;
            }

            if (q == 0)
                return double.PositiveInfinity;
            if (q == 1)
                return 0; 

            // use I_x(a,b) = 1 - I_{1-x}(b, a)
            if (q > 0.5 && x > 0.5)
                return _IbetaInvAB.IbetaInvAImp(a, 1 - x, x, q, 1 - q);

            return _IbetaInvAB.IbetaInvBImp(a, x, 1 - x, 1 - q, q);
        }

    }

} // namespaces



