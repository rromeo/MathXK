//  Copyright (c) 2013 Rocco Romeo
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

//  Ported to .NET based on the original:
//      Copyright (c) John Maddock 2006, Boost Software License v1.0
//
//  History:
//      JM created original C++ code
//      RR ported to C# and added limits and guess refinements

// This inverts the incomplete gamma functions P and Q on the first parameter "a" using a generic
// root finding algorithm (TOMS Algorithm 748).

#if DEBUG
//#define EXTRA_DEBUG
#endif


using System;
using System.Collections.Generic;
using System.Diagnostics;
using MathXK.Roots;

namespace MathXK
{

    static class _GammaInvA {

        [Conditional("EXTRA_DEBUG"), Conditional("DEBUG")]
        static void ExtraDebugLog(double z, double p, double q, double aResult, double aGuess, double step, double minA, double maxA, int iterations, string notes)
        {
            if (iterations > 6) {
                if (p <= q) {
                    Debug.WriteLine($"PInvA(z: {z}, p: {p}) = {aResult}; g: {aGuess}; step={step}; Iter = {iterations}; Range: [{minA}, {maxA}]; {notes}");
                } else {
                    Debug.WriteLine($"QInvA(z: {z}, q: {q}) = {aResult}; g: {aGuess}; step={step}; Iter = {iterations}; Range: [{minA}, {maxA}]; {notes}");
                }
            }

        }

        /// <summary>
        /// Simple, stable, but not very accurate LogGammaP for small z &lt; 0.7*a
        /// <para>Log((z^a)(e^-z)/(a*Γ(a)) * 1F1(1, a+1, z))</para>
        /// </summary>
        /// <param name="a"></param>
        /// <param name="z"></param>
        /// <returns></returns>
        static double LogGammaPSmallZ(double a, double z)
        {
            // See: http://dlmf.nist.gov/8.5#E1
            double t = HypergeometricSeries.Sum1F1(1, a + 1, z, -1);
            double result = a * Math.Log(z) - z - Math2.Lgamma(1 + a) + Math2.Log1p(t);
            return result;
        }

        public static (double Result, int Iterations) SolveA(double z, double p, double q, double aGuess, double step, double min, double max)
        {
            Debug.Assert(aGuess >= min && aGuess <= max, "Guess out of range");

            Func<double, double> f;
            if (p < q) {
                f = a => Math2.GammaP(a, z) - p;
            } else {
                f = a => Math2.GammaQ(a, z) - q;
            }

            //
            // Use our generic derivative-free root finding procedure.
            // We could use Newton steps here, taking the PDF of the
            // Poisson distribution as our derivative, but that's
            // even worse performance-wise than the generic method :-(
            //

            // dQ/da is increasing, dP/da is decreasing
            FunctionShape fShape = (p < q) ? FunctionShape.Decreasing : FunctionShape.Increasing;
            RootResults rr = RootFinder.Toms748Bracket(f, aGuess, step, fShape, min, max);
            if (rr == null) {
                Policies.ReportRootNotFoundError($"Invalid Parameter: PQInvA(z: {z}, p: {p}, q: {q}) [{min}, {max}] g: {aGuess}");
                return (double.NaN, 0);
            }

            if (!rr.Success) {
                Policies.ReportRootNotFoundError("Root not found after {0} iterations", rr.Iterations);
                return (double.NaN, rr.Iterations);
            }

            return (rr.SolutionX, rr.Iterations);
        }

        static double InversePoissonCornishFisher(double lambda, double p, double q)
        {
            // mean:
            double m = lambda;
            // standard deviation:
            double sigma = Math.Sqrt(lambda);
            // skewness
            double sk = 1 / sigma;
            // kurtosis:
            // double k = 1/lambda;
            // Get the inverse of a std normal distribution:
            double x = Math2.ErfcInv(p > q ? 2 * q : 2 * p) * Constants.Sqrt2;
            // Set the sign:
            if (p < 0.5)
                x = -x;
            double x2 = x * x;
            // w is correction term due to skewness
            // double w = x + sk * (x2 - 1) / 6;

            double w = x + (x2 - 1) / (6 * sigma);

            //if(lambda >= 10)
            //   w += k * x * (x2 - 3) / 24 + sk * sk * x * (2 * x2 - 5) / -36;

            if (lambda >= 10)
                w += x / 72 * (1 - x2) / lambda;

            w = m + sigma * w;

            return w;
        }

        /// <summary>
        /// Returns a new estimate of "a" based on "aGuess" such that P(a, z) == p
        /// <para>Expecting a >> z </para>
        /// </summary>
        /// <param name="z"></param>
        /// <param name="p"></param>
        /// <param name="q"></param>
        /// <param name="aGuess"></param>
        /// <param name="aMax"></param>
        /// <param name="aMin"></param>
        /// <returns></returns>
        public static (double Result, int Iterations) RefineP_LargeA_SmallZ(double z, double p, double q, double aGuess, ref double aMin, ref double aMax)
        {
            Debug.Assert(z >= 0, "Requires z >= 0");
            Debug.Assert(aGuess >= 0 && aGuess <= double.MaxValue, "Requires finite aGuess >= 0");
            Debug.Assert(p <= 0.5, "Requires p <= 0.5");
            Debug.Assert(aGuess > z);

            // Refinement approach
            // Given p = GammaP(a, z)
            // Assume our guess for a, ag = a + delta, where delta is small
            // Create a first order model for GammaP(a+delta, z)/GammaP(a, z) then solve for delta
            // 
            // As a -> Infinity, GammaP(a+delta, z)/GammaP(a, z) ~ (z/a)^delta *(1 - 6*delta*(1+delta)/(-1 + 12a + 12z) ...)
            // Mathematica Eqn:
            // Simplify[Normal[ Series[GammaRegularized[a + delta, 0, z]/ GammaRegularized[a, 0, z], {a, Infinity, 2}]], {a > 0, z > 0, a + delta > 0}]

            // Let Pa = GammaP(a, z) = p
            // Pg = GammaP(a+delta, z)
            //
            // Using linear extrapolation, Pa = Pg + (Pa-Pg)/delta * delta
            // => p = Pg + Pg*(Pa/Pg-1) ~= Pg * (z/(a+delta))^-delta ~= (z/ag)^-delta (when delta << a )
            // So, delta ~= Log(Pg/p) / Log(z/ag)
            // And our next guess, aNext = ag - delta

            int i = 0;
            for (; i < 5; i++) {
                const double eps = 1.0 / 16384; // through logs we can safely achieve 1e-6
                // GammaP(a, z) can underflow when dealing with small p, so use logs for stability
                // Since we are in p < 0.5, aGuess > z
                double num = LogGammaPSmallZ(aGuess, z) - Math.Log(p);
                double d = (z - aGuess) / aGuess;
                double den = (Math.Abs(d) < 0.5) ? Math2.Log1p(d) : Math.Log(z) - Math.Log(aGuess);
                double delta = num / den;

                double aLast = aGuess;
                aGuess -= delta;

                if (delta == 0)
                    break;

                if (delta > 0)
                    aMax = aLast;
                else
                    aMin = aLast;

                if (Math.Abs(delta) <= eps * Math.Abs(aLast))
                    break;
            }

            return (aGuess, i);
        }

        public static (double Result, int Iterations) RefineQ_SmallA_LargeZ(double z, double p, double q, double aGuess)
        {
            Debug.Assert(z >= 0, "Requires z >= 0");
            Debug.Assert(aGuess >= 0 && aGuess <= double.MaxValue, "Requires finite aGuess >= 0");
            Debug.Assert(q <= 0.5, "Requires q <= 0.5");
            Debug.Assert(aGuess < z+1);

            // f = Log(z^a/Gamma(a))-Log(q);
            // f/f' = (Log(z^a/Gamma(a))-Log(q))/ (Log(z) - Polygamma(0, a))
            //
            // http://functions.wolfram.com/GammaBetaErf/GammaRegularized/06/02/02/0001/
            // For z->Infinity, Q(a,z) ~= z^(a-1)*e^-z/Gamma(a) * (1 - (1-a)/z ...)

            // Find a quick estimate of a 
            // where Log[q] = Log[z^(a-1)*e^-z/Gamma(a)] to get into the ballpark

            double lz = Math.Log(z);
            double lq = Math.Log(q);

            int i = 0;
            for (; i < 1; i++) {
                const double eps = 1.0 / 2; 

                double f = (aGuess - 1) * lz - z - Math2.Lgamma(aGuess) - lq;
                if (f == 0)
                    break;

                var delta = f/(lz - Math2.Digamma(aGuess));
                double aLast = aGuess;
                aGuess -= delta;

                if (Math.Abs(delta) <= eps * Math.Abs(aLast))
                    break;
            }

            return (aGuess, i);
        }

        public static double ImpQ(double z, double p, double q)
        {
            Debug.Assert(q <= 0.5);

            if (q == 0)
                return DoubleLimits.MinNormalValue;

#if EXTRA_DEBUG
            var gs = new Stack<double>();
#endif

            // Set outside limits
            // Since Q(a, z) > 1/2 for a >= z+1, so maxA < z+1
            double minA = DoubleLimits.MinNormalValue;
            double maxA = z + 1;
            if (maxA == z)
                maxA = Math2.FloatNext(z);


            double step = 0.5;
            double aGuess = 0;

            // We can use the relationship between the incomplete 
            // gamma function and the poisson distribution to
            // calculate an approximate inverse. 
            // Mostly it is pretty accurate, except when a is small or q is tiny.  

            // Case 1: a <= 1 : Inverse Cornish Fisher is unreliable in this area

            var e = Math.Exp(-z);
            if (q <= e) {
                if (q == e)
                    return 1;

                maxA = 1;

                // If a <= 1, P(a, x) = 1 - Q(a, x) >= (1-e^-x)^a
                var den = (z > 0.5) ? Math2.Log1p(-Math.Exp(-z)) : Math.Log(-Math2.Expm1(-z));
                var limit = Math2.Log1p(-q) / den;
                if ( limit < 1)
                    minA = Math.Max(minA, limit);

                step = (1 - minA) * 0.125;
                aGuess = minA * (1 + step);

                if (minA < 0.20) {
                    // Use a first order approximation Q(a, z) at a=0
                    // Mathematica: Assuming[ a >= 0 && z > 0, FunctionExpand[Normal[Series[GammaRegularized[a, z], {a, 0, 1}]]]] 
                    // q = -a * ExpIntegralEi[-z]
                    // The smaller a is the better the guess
                    // as z->0 term2/term1->1/2*Log[z], which is adequate for double precision

                    aGuess = -q / Math2.Expint(-z);
                    if (aGuess <= DoubleLimits.MachineEpsilon) {
                        ExtraDebugLog(z, p, q, aGuess, aGuess, step, minA, maxA, 0, "");
                        return aGuess;
                    }
                    step = aGuess * 0.125;

                } else if (minA >= 0.40) {
                    // Use a first order approximation Q(a, z) at a=1
                    // Mathematica: Assuming[ a >= 0 && z > 0, FunctionExpand[Normal[Series[GammaRegularized[a, z], {a, 1, 1}]]]] 

                    var d = Math.Exp(-z) * (Math.Log(z) + Constants.EulerMascheroni) + Math2.Expint(1, z);
                    aGuess = (q - Math.Exp(-z)) / d + 1;

                }

                aGuess = Math.Max(aGuess, minA);
                aGuess = Math.Min(aGuess, maxA);

                var r = SolveA(z, p, q, aGuess, step, minA, maxA);

                ExtraDebugLog(z, p, q, r.Result, aGuess, step, minA, maxA, r.Iterations, "");

                return r.Result;
            }


            // Case 2: a > 1 : Inverse Poisson Corning Fisher approximation improves as a->Infinity

            minA = 1;
            aGuess = 0.5 + InversePoissonCornishFisher(z, q, p);
            if (aGuess > 100) {
                step = 0.01;
            } else if (aGuess > 10) {
                step = 0.01;
            } else {
                // our poisson approximation is weak.
                step = 0.05;
            }


#if EXTRA_DEBUG
            double og = aGuess;
#endif

            aGuess = Math.Max(aGuess, minA);
            aGuess = Math.Min(aGuess, maxA);

#if EXTRA_DEBUG
            if (og != aGuess)
                gs.Push(og);
#endif

            // For small values of q, our Inverse Poisson Cornish Fisher guess can be way off
            // Try to come up with a better guess
            if (q < DoubleLimits.MachineEpsilon ) {
                step *= 10;
                if (Math.Abs((aGuess - 1) / z) < 1.0 / 8) {
#if EXTRA_DEBUG
                    double old = aGuess;
#endif
                    (aGuess, _) = RefineQ_SmallA_LargeZ(z, p, q, aGuess);
#if EXTRA_DEBUG
                    if (old != aGuess)
                        gs.Push(old);
#endif
                }
            }

            var rr = SolveA(z, p, q, aGuess, step, minA, maxA);

#if EXTRA_DEBUG
            if (rr.Iterations > 6) {
                string guesses = (gs.Count == 0) ? string.Empty : ", og: (" + string.Join(", ", gs) + ")";
                ExtraDebugLog(z, p, q, rr.Result, aGuess, step, minA, maxA, rr.Iterations, guesses);
            }
#endif

            return rr.Result;

        }

        public static double ImpP(double z, double p, double q)
        {
            Debug.Assert(p > 0 && p < 0.5);

            // Set outside limits
            // Since P(a, z) >= 1/2 for a >= z, so minA < z

            double minA = z;
            double maxA = double.MaxValue;

            double aGuess = 0;
            double step = 0.5;

            if (z <= 1) {
                var e = -Math2.Expm1(-z);
                if (p >= e) {
                    if (p == e)
                        return 1;

                    // This is one of the outside limits of the GammaP/GammaQ
                    // From: http://dlmf.nist.gov/8.10#E11
                    // If a >= 1, P(a, x) <= (1-e^-x)^a
                    // If a <= 1, P(a, x) >= (1-e^-x)^a
                    // not a bad guess, for small z and small a
                    // but gets much worse as x increases or a >> 1      

                    maxA = 1;
                    var limit = Math.Log(p) / Math.Log(e);
                    if (limit < 1) 
                        minA = Math.Max(minA, limit);

                    // could use bisection, but the root often lies nearer to minA
                    step = (1 - minA) * 0.125;
                    aGuess = minA * (1 + step);

                    var r = SolveA(z, p, q, aGuess, step, minA, maxA);

                    ExtraDebugLog(z, p, q, r.Result, aGuess, step, minA, maxA, r.Iterations,"");

                    return r.Result;

                } else {
                    minA = 1;
                    var limit = Math.Log(p) / Math.Log(e);
                    if (limit > 1)
                        maxA = Math.Min(maxA, limit);
                }
            }

            aGuess = 0.5 + InversePoissonCornishFisher(z, q, p);
            if (aGuess > 100) {
                step = 0.01;
            } else if (aGuess > 10) {
                step = 0.01;
            } else {
                step = 0.05;
            }

#if EXTRA_DEBUG
            double og = aGuess;
#endif
            aGuess = Math.Max(aGuess, minA);
            aGuess = Math.Min(aGuess, maxA);

            if (p < DoubleLimits.MachineEpsilon)
                step *= 10;

#if EXTRA_DEBUG
            var gs = new Stack<double>();
            if (og != aGuess)
                gs.Push(og);
#endif

            // When p is very small, our Inverse Poisson Cornish Fisher can be significantly off.
            // Refining the guess here can significantly reduce iterations later
            int refineIter = 0;
            if (aGuess >= 5 && z / (1 + aGuess) <= 0.7 ) {

#if EXTRA_DEBUG
                gs.Push(aGuess);
#endif
                (aGuess, refineIter) = RefineP_LargeA_SmallZ(z, p, q, aGuess, ref minA, ref maxA);
                step = aGuess * 1e-6;

                aGuess = Math.Max(aGuess, minA);
                aGuess = Math.Min(aGuess, maxA);
            }

            var rr = SolveA(z, p, q, aGuess, step, minA, maxA);

#if EXTRA_DEBUG
            if (rr.Iterations + refineIter > 6) {

                string refString = (refineIter == 0) ? string.Empty : "(" + rr.Iterations + "+" + refineIter + ")";
                string guesses = (gs.Count == 0) ? string.Empty : ", og: (" + string.Join(", ", gs) + ")";

                Debug.WriteLine($"PInvA(z: {z}, p: {p}) = {rr.Result}, g: {aGuess}, TIter = {rr.Iterations + refineIter}{refString}{guesses}");
            }
#endif
            return rr.Result;
        }
    }



    public static partial class Math2
    {

        /// <summary>
        /// Returns the value "a" such that: p == GammaP(a, x);
        /// </summary>
        /// <param name="x">Requires x ≥ 0</param>
        /// <param name="p">Requires 0 ≤ p ≤ 1</param>
        public static double GammaPInvA(double x, double p)
        {
            if ((!(x >= 0) || double.IsInfinity(x))
            || (!(p >= 0 && p <= 1))) {
                Policies.ReportDomainError("GammaPInvA(x: {0}, p: {1}): Requires finite x >= 0; p in [0,1]", x, p);
                return double.NaN;
            }

            if (p == 0)
                return double.MaxValue;
            if (p == 1)
                return DoubleLimits.MinNormalValue;

            if (p < 0.5)
                return _GammaInvA.ImpP(x, p, 1 - p);

            return _GammaInvA.ImpQ(x, p, 1 - p);
        }

        /// <summary>
        /// Returns the value "a" such that q == GammaQ(a, x);
        /// </summary>
        /// <param name="x">Requires x ≥ 0</param>
        /// <param name="q">Requires 0 ≤ q ≤ 1</param>
        public static double GammaQInvA(double x, double q)
        {
            if ((!(x >= 0) || double.IsInfinity(x))
            || (!(q >= 0 && q <= 1))) {
                Policies.ReportDomainError("GammaPInvA(x: {0}, q: {1}): Requires finite x >= 0; q in [0,1]", x, q);
                return double.NaN;
            }

            if (q == 0)
                return DoubleLimits.MinNormalValue;
            if (q == 1)
                return double.MaxValue;

            var p = 1 - q;
            if (q <= 0.5)
                return _GammaInvA.ImpQ(x, p, q);

            return _GammaInvA.ImpP(x, p, q);
        }

    }



} // namespaces




