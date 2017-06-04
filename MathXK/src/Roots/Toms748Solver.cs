//  Copyright (c) 2013 Rocco Romeo
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

//  Ported to .NET based on the original:
//      Copyright (C) John Maddock 2006, Boost Software License, Version 1.0

//  TOMS algorithm 748 is described in detail in:
//  Algorithm 748: Enclosing Zeros of Continuous Functions, G. E. Alefeld, F. A. Potra and Yixun Shi, ACM Transactions on Mathematica1 Software, Vol. 21. No. 3. September 1995. Pages 327-344.
//  and http://www.netlib.org/toms/748

#if DEBUG
//#define EXTRA_DEBUG
#endif

using System;
using System.Diagnostics;

namespace MathXK.Roots
{


    public static partial class RootFinder
    {


        /// <summary>
        /// Performs standard secant interpolation of [a,b] given
        /// function evaluations f(a) and f(b).  
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="fa"></param>
        /// <param name="fb"></param>
        /// <returns></returns>
        static double SecantStep(double a, double b, double fa, double fb)
        {
            double c = a - (fa / (fb - fa)) * (b - a);
            return c;
        }

        static double QuadraticStep(double a, double b, double d, double fa, double fb, double fd, int count)
        {
            //
            // Performs quadratic interpolation to determine the next point,
            // takes count Newton steps to find the location of the
            // quadratic polynomial.
            //
            // Point d must lie outside of the interval [a,b], it is the third
            // best approximation to the root, after a and b.
            //
            // Note: this does not guarentee to find a root
            // inside [a, b], so we fall back to a secant step should
            // the result be out of range.
            //
            // Start by obtaining the coefficients of the quadratic polynomial:
            //

            double C = fa;
            double B = (fb - fa)/(b - a);
            double A = ((fd - fb) / (d - b) - B) / (d - a);

            // if a == 0, we do not have an quadratic 
            // if a is not finite, we have an error in our calculations
            // try a secant step, instead:
            if (A == 0 || double.IsNaN(A) || double.IsInfinity(A)) 
                return SecantStep(a, b, fa, fb);
            
            //
            // Determine the starting point of the Newton steps:
            //
            double c = (Math.Sign(A) * Math.Sign(fa) > 0) ? a : b;

            //
            // Take the Newton steps:
            //
            bool error = false;
            for (int i = 0; i < count; i++) {
                double fi = C + (B + A * (c - b)) * (c - a); // f(x_i)
                double fpi = B + A * (2 * c - a - b);        // f'(x_i)
                double delta = fi/fpi;
                if (double.IsInfinity(delta) || double.IsNaN(delta)) {
                    error = true;
                    break;
                }
                c -= delta;
            }
            if ( error || (c <= a) || (c >= b)) {
                // Oops, failure, try a secant step:
                c = SecantStep(a, b, fa, fb);
            }
            return c;
        }

        static double CubicStep(double a, double b, double d, double e, double fa, double fb, double fd, double fe)
        {
            //
            // Uses inverse cubic interpolation of f(x) at points 
            // [a,b,d,e] to obtain an approximate root of f(x).
            // Points d and e lie outside the interval [a,b]
            // and are the third and forth best approximations
            // to the root that we have found so far.
            //
            // Note: this does not guarentee to find a root
            // inside [a, b], so we fall back to quadratic
            // interpolation in case of an erroneous result.
            //

            double q11 = (d - e) * (fd / (fe - fd));
            double q21 = (b - d) * (fb / (fd - fb));
            double q31 = (a - b) * (fa / (fb - fa));
            double d21 = (b - d) * (fd / (fd - fb));
            double d31 = (a - b) * (fb / (fb - fa));

            double q22 = (d21 - q11) * (fb / (fe - fb));
            double q32 = (d31 - q21) * (fa / (fd - fa));
            double d32 = (d31 - q21) * (fd / (fd - fa));
            double q33 = (d32 - q22) * (fa / (fe - fa));
            double c = q31 + q32 + q33 + a;

#if EXTRA_DEBUG
            Debug.WriteLine("a = {0} b = {1} d = {2} e = {3} fa = {4} fb = {5} fd = {6} fe = {7}", a, b, d, e, fa, fb, fd, fe);
            Debug.WriteLine("q11 = {0} q21 = {1} q31 = {2} d21 = {3} d31 = {4}", q11, q21, q31, d21, d31);
            Debug.WriteLine("q22 = {0} q32 = {1} d32 = {2} q33 = {3} c = {4}", q22, q32, d32, q33, c);
#endif

            return c;
        }

        static double DoubleSecantStep(double a, double b, double fa, double fb)
        {
            double u;
            double fu;
            if (Math.Abs(fa) < Math.Abs(fb)) {
                u = a;
                fu = fa;
            } else {
                u = b;
                fu = fb;
            }
            double c = u - 2 * (fu / (fb - fa)) * (b - a);
            if (Math.Abs(c - u) > (b - a) / 2) {
                c = a + (b - a) / 2;
            }


            return c;
        }


        private static RootResults Toms748Int(Func<double, double> f, double min, double fmin, double max, double fmax, ToleranceFunc areNear, int max_iter)
        {
            // private routine so these conditions should already be true
            Debug.Assert(f != null, "Function cannot be null");
            Debug.Assert(min <= max, $"Requires min <= max, min: {min}, max: {max}");
            Debug.Assert(IsBracket(fmin, fmax), "Root is not bracketed");

            //
            // Main entry point and logic for Toms Algorithm 748
            // root finder.
            //

            double xSolution = double.NaN;
            int iterations = 0;


            double a = min;
            double fa = fmin;

            double b = max;
            double fb = fmax;

            if (fa == 0) {
                xSolution = a;
                goto exit;
            }
            if (fb == 0) {
                xSolution = b;
                goto exit;
            }
            if (areNear(a, b)) {
                xSolution = a + (b - a) / 2.0; // don't have an exact solution so take a midpoint
                goto exit;
            }


            double d, fd; // d is the third best guess to the root
            double e, fe; // e is the fourth best guess to the root


            // initialize
            fe = e = d = fd = double.NaN;
            double c;

            int stepNumber = 0;
            while (iterations < max_iter) {
                // the root is in [a, b]; save it 
                double a0 = a;
                double b0 = b;

                switch(stepNumber) {
                case 0:
                    // *Only* on the first iteration we take a secant step:
                    c = SecantStep(a, b, fa, fb);
                    break;
                case 1:
                    // *Only* on the second iteration, we take a quadratic step:
                    c = QuadraticStep(a, b, d, fa, fb, fd, 2);
                    break;
                case 2:
                case 3:
                    //
                    // Starting with the third step taken
                    // we can use either quadratic or cubic interpolation.
                    // Cubic interpolation requires that all four function values
                    // fa, fb, fd, and fe are distinct, should that not be the case
                    // then variable prof will get set to true, and we'll end up
                    // taking a quadratic step instead.
                    //
                    const double min_diff = DoubleLimits.MinNormalValue;
                    bool prof = (Math.Abs(fa - fb) < min_diff) 
                        || (Math.Abs(fa - fd) < min_diff) 
                        || (Math.Abs(fa - fe) < min_diff) 
                        || (Math.Abs(fb - fd) < min_diff) 
                        || (Math.Abs(fb - fe) < min_diff) 
                        || (Math.Abs(fd - fe) < min_diff);

                    // the first time use 2 newton steps; the second time use three
                    int newtonSteps = (stepNumber == 2) ? 2 : 3;
                    if (prof) {
                        c = QuadraticStep(a, b, d, fa, fb, fd, newtonSteps);
                    } else {
                        c = CubicStep(a, b, d, e, fa, fb, fd, fe);
                        if (!(c > a && c < b))
                            c = QuadraticStep(a, b, d, fa, fb, fd, newtonSteps);
                    }
                    break;
                case 4:
                    // Now we take a double-length secant step:
                    c = DoubleSecantStep(a, b, fa, fb);
                    break;
                case 5:
                    //
                    // And finally... check to see if an additional bisection step is 
                    // to be taken, we do this if we're not converging fast enough:
                    //
                    const double mu = 0.5;
                    if ((b - a) < mu * (b0 - a0)) {
                        stepNumber = 2;
                        continue;
                    }
                    c = a + (b - a) / 2;
                    break;
                default:
                    Policies.ReportEvaluationError("Switch out of bounds");
                    return null;
                
                }

                // save our next best bracket
                e = d;
                fe = fd;


                // We require c in the interval [ a + |a|*minTol, b - |b|*minTol ]
                // If the interval is not valid: (b-a) < (|a|+|b|)*minTol, then take a midpoint

                const double minTol = DoubleLimits.MachineEpsilon * 10;

                double minA = (a < 0) ? a * (1 - minTol) : a * (1 + minTol);
                double minB = (b < 0) ? b * (1 + minTol) : b * (1 - minTol);
                if ((b - a) < minTol * (Math.Abs(a) + Math.Abs(b))) {
                    c = a + (b - a) / 2;
                } else if (c < minA) {
                    c = minA;
                } else if (c > minB) {
                    c = minB;
                }

                // OK, lets invoke f(c):
                double fc = f(c);
                ++iterations;

                // if we have a zero then we have an exact solution to the root:
                if (fc == 0) {
                    xSolution = c;
                    break;
                }

                if (double.IsNaN(fc)) {
                    Policies.ReportDomainError("Requires a continuous function: f({0}) = {1}", c, fc);
                    break;
                }

                //
                // Non-zero fc, update the interval:
                //
                if ( IsBracket(fa, fc) ) {
                    d = b; fd = fb;
                    b = c; fb = fc;
                } else {
                    d = a; d = fa;
                    a = c; fa = fc;
                }


                if (areNear(a, b)) {
                    // we're close enough, but we don't have an exact solution
                    // so take a midpoint
                    xSolution = a + (b - a) / 2.0; 
                    break;
                }

                // reset our step if necessary;
                if ( ++stepNumber > 5 )
                    stepNumber = 2;

            }

        exit:
            RootResults result = new RootResults();
            result.SolutionX = xSolution;
            result.Bracket = new RootBracketResult(a, fa, b, fb);
            result.Iterations = iterations;
            return result;

        }


        /// <summary>
        /// TOMS Algorithm 748 uses a mixture of cubic, quadratic and linear (secant) interpolation to locate the root of f(x)
        /// Note: min and max must bracket the root.
        /// </summary>
        /// <param name="f">Function to solve</param>
        /// <param name="min">The root must lie within [min, max]</param>
        /// <param name="max">The root must lie within [min, max]</param>
        /// <param name="areNear">A function that returns true when the required tolerance is reached, or null to use the default tolerance</param>
        /// <param name="maxIterations">The maximum number of iterations to keep refining the root estimate</param>
        /// <returns>
        /// Null, if any of the parameters were incorrect.
        /// The RootResults which contains the root if one was found, or the final state before failing.
        /// </returns>
        public static RootResults Toms748(Func<double, double> f, double min, double max, ToleranceFunc areNear, int maxIterations)
        {
            const RootResults nullResult = null;

            if (f == null) {
                Policies.ReportDomainError("Function cannot be null");
                return nullResult;
            }
            if (!(min <= max)) {
                Policies.ReportDomainError("Requires min <= max : min = {0}; max = {1}", min, max);
                return nullResult;
            }

            // use the default if max iterations is negative
            if (maxIterations < 3) {
                Policies.ReportDomainError("Requires maxIterations >= 3: maxIterations = {0}", maxIterations);
                return nullResult;
            }

            // use the default if tolerance is null
            if (areNear == null)
                areNear = GetToleranceFunc();


            double fmin = f(min);
            double fmax = f(max);

            if (!IsBracket(fmin, fmax)) {
                RootResults results = new RootResults()
                {
                    Iterations = 2,
                    Bracket = new RootBracketResult(min, fmin, max, fmax)
                };
                Policies.ReportDomainError("The root is not bracketed: f({0}) = {1}; f({2}) = {3}", min, fmin, max, fmax);
                return results;
            }

            var r = Toms748Int(f, min, fmin, max, fmax, areNear, maxIterations - 2);
            r.Iterations += 2;
            return r;
        }

        /// <summary>
        /// TOMS Algorithm 748 uses a mixture of cubic, quadratic and linear (secant) interpolation to locate the root of f(x)
        /// Note: min and max must bracket the root.
        /// </summary>
        /// <param name="f">Function to solve</param>
        /// <param name="min">The root must lie within [min, max]</param>
        /// <param name="max">The root must lie within [min, max]</param>
        /// <param name="relTolerance">The maximum relative tolerance to be reached before stopping</param>
        /// <param name="absTolerance">The maximum absolute tolerance to be reached before stopping</param>
        /// <param name="maxIterations">The maximum number of iterations to keep refining the root estimate</param>
        /// <returns>
        /// Null, if any of the parameters were incorrect.
        /// The RootResults which contains the root if one was found, or the final state before failing.
        /// </returns>
        public static RootResults Toms748(Func<double, double> f, double min, double max, double relTolerance, double absTolerance, int maxIterations)
        {
            const RootResults nullResult = null;

            if (!(relTolerance >= 0) || double.IsInfinity(relTolerance)) {
                Policies.ReportDomainError("Requires finite relative tolerance >= 0: relTolerance = {0}", relTolerance);
                return nullResult;
            }

            if (!(absTolerance >= 0) || double.IsInfinity(absTolerance)) {
                Policies.ReportDomainError("Requires finite absolute tolerance >= 0: absTolerance = {0}", absTolerance);
                return nullResult;
            }

            var areNear = GetToleranceFunc(relTolerance, absTolerance);
            return Toms748(f, min, max, areNear, maxIterations);
        }

        /// <summary>
        /// TOMS Algorithm 748 uses a mixture of cubic, quadratic and linear (secant) interpolation to locate the root of f(x)
        /// Note: min and max must bracket the root.
        /// </summary>
        /// <param name="f">Function to solve</param>
        /// <param name="min">The root must lie within [min, max]</param>
        /// <param name="max">The root must lie within [min, max]</param>
        /// <returns>
        /// Null, if any of the parameters were incorrect.
        /// The RootResults which contains the root if one was found, or the final state before failing.
        /// </returns>
        public static RootResults Toms748(Func<double, double> f, double min, double max)
        {
            return Toms748(f, min, max, null, Policies.MaxRootIterations);
        }


        /// <summary>
        /// TOMS Algorithm 748 uses a mixture of cubic, quadratic and linear (secant) interpolation to locate the root of f(x)
        /// Note: TOMS748 requires the root to be bracketed, so this routine will attempt to bracket the root.
        /// </summary>
        /// <param name="f">Function to solve</param>
        /// <param name="guess">The initial guess for the root.</param>
        /// <param name="step">The starting step size</param>
        /// <param name="fType">Unknown, increasing, or decreasing</param>
        /// <param name="min">Optional. The minimum location for the root, or NaN</param>
        /// <param name="max">Optional. The maximum location for the root, or NaN</param>
        /// <param name="areNear">A function that returns true when the required tolerance is reached, or null to use the default tolerance</param>
        /// <param name="maxIterations">The maximum number of iterations to keep refining the root estimate</param>
        /// <returns>
        /// Null, if any of the parameters were incorrect.
        /// The RootResults which contains the root if one was found, or the final state before failing.
        /// </returns>
        public static RootResults Toms748Bracket(Func<double, double> f, double guess, double step, FunctionShape fType, double min, double max, ToleranceFunc areNear, int maxIterations)
        {
            const RootResults nullResult = null;

            if (f == null) {
                Policies.ReportDomainError("Function cannot be null");
                return nullResult;
            }
            // NaN is valid here for min or max, so only fail if both have a value and min > max
            if (min > max) {
                Policies.ReportDomainError("Requires min <= max : min = {0}; max = {1}", min, max);
                return nullResult;
            }
            // use the default if max iterations is negative
            if (maxIterations < 3) {
                Policies.ReportDomainError("Requires maxIterations >= 3: maxIterations = {0}", maxIterations);
                return nullResult;
            }

            // use the default if tolerance is null
            if (areNear == null)
                areNear = GetToleranceFunc();


            int iterations;
            RootBracketResult b = FindBracket(f, guess, step, fType, min, max, maxIterations, out iterations);
            if (b == null) {
                Policies.ReportDomainError("Invalid parameter in bracket");
                return nullResult;
            }

            if (!b.IsValid) {
                RootResults result = new RootResults()
                {
                    Bracket = b,
                    Iterations = iterations
                };
                return result;
            }

            var r = Toms748Int(f, b.XMin, b.FxMin, b.XMax, b.FxMax, areNear, maxIterations - iterations);
            r.Iterations += iterations;
            return r;
        }

        /// <summary>
        /// TOMS Algorithm 748 uses a mixture of cubic, quadratic and linear (secant) interpolation to locate the root of f(x)
        /// <para>Note: TOMS748 requires the root to be bracketed, so this routine will attempt to bracket the root before calling TOMS748.</para>
        /// </summary>
        /// <param name="f">Function to solve</param>
        /// <param name="guess">The initial guess for the root.</param>
        /// <param name="step">The starting step size</param>
        /// <param name="fType">Unknown, increasing, or decreasing</param>
        /// <param name="min">Optional. The minimum location for the root, or NaN</param>
        /// <param name="max">Optional. The maximum location for the root, or NaN</param>
        /// <returns>
        /// Null, if any of the parameters were incorrect.
        /// The RootResults which contains the root if one was found, or the final state before failing.
        /// </returns>
        public static RootResults Toms748Bracket(Func<double, double> f, double guess, double step, FunctionShape fType, double min, double max)
        {
            return Toms748Bracket(f, guess, step, fType, min, max, null, Policies.MaxRootIterations);
        }


        /// <summary>
        /// TOMS Algorithm 748 uses a mixture of cubic, quadratic and linear (secant) interpolation to locate the root of f(x)
        /// Note: TOMS748 requires the root to be bracketed, so this routine will attempt to bracket the root.
        /// </summary>
        /// <param name="f">Function to solve</param>
        /// <param name="guess">The initial guess for the root.</param>
        /// <param name="step">The starting step size</param>
        /// <param name="fType">Unknown, increasing, or decreasing</param>
        /// <param name="min">Optional. The minimum location for the root, or NaN</param>
        /// <param name="max">Optional. The maximum location for the root, or NaN</param>
        /// <param name="relTolerance">The maximum relative tolerance to be reached before stopping</param>
        /// <param name="absTolerance">The maximum absolute tolerance to be reached before stopping</param>
        /// <param name="maxIterations">The maximum number of iterations to keep refining the root estimate</param>
        /// <returns>
        /// Null, if any of the parameters were incorrect.
        /// The RootResults which contains the root if one was found, or the final state before failing.
        /// </returns>
        public static RootResults Toms748Bracket(Func<double, double> f, double guess, double step, FunctionShape fType, double min, double max, double relTolerance, double absTolerance, int maxIterations)
        {
            const RootResults nullResult = null;

            if (!(relTolerance >= 0) || double.IsInfinity(relTolerance)) {
                Policies.ReportDomainError("Requires finite relative tolerance >= 0: relTolerance = {0}", relTolerance);
                return nullResult;
            }

            if (!(absTolerance >= 0) || double.IsInfinity(absTolerance)) {
                Policies.ReportDomainError("Requires finite absolute tolerance >= 0: absTolerance = {0}", absTolerance);
                return nullResult;
            }

            var tol = GetToleranceFunc(relTolerance, absTolerance);
            return Toms748Bracket(f, guess, step, fType, min, max, tol, maxIterations);
        }

    }

}