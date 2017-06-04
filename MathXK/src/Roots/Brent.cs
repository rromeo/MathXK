//  Copyright (c) 2013 Rocco Romeo
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)


using System;
using System.Diagnostics;

namespace MathXK.Roots
{

    public static partial class RootFinder
    {

        // For more information, see: http://en.wikipedia.org/wiki/Brent's_method


        private static RootResults BrentInt(Func<double, double> f, double min, double fmin, double max, double fmax, ToleranceFunc areNear, int maxIterations)
        {

            // private routine so these conditions should already be true
            Debug.Assert(f != null, "Function cannot be null");
            Debug.Assert(min < max, "Parameter min < max");
            Debug.Assert(IsBracket(fmin, fmax), "Root is not bracketed");

            // code based on pseudocode in http://en.wikipedia.org/wiki/Brent's_method

            double xSolution = double.NaN;
            int iterations = 0;
            double x = double.NaN;

            if (fmin == 0) {
                xSolution = min;
                goto exit;
            }
            if (fmax == 0) {
                xSolution = max;
                goto exit;
            }
            if (areNear(min, max)) {
                xSolution = min + (max - min) / 2.0; // don't have an exact solution so take a midpoint
                goto exit;
            }

            // b is chosen to be the closer point to the root (i.e. Abs(fb) < Abs(fa) )
            // a is always of the opposite sign as b so that [a,b] or [b,a] contains the root
            // c = the previous value of b except initially when c = a
            // d = the previous value of c

            double a, b, c, fa, fb, fc, d = 0, df = 0;
            if (Math.Abs(fmin) < Math.Abs(fmax)) {
                a = max; fa = fmax;
                b = min; fb = fmin;
            } else {
                b = max; fb = fmax;
                a = min; fa = fmin;
            }

            c = a; fc = fa; // first step will be a secant step


            bool usedBisection = true; // last step was bisection -- is d set
            while (iterations < maxIterations) {


                if (fa != fc && fb != fc) {
                    // use inverse quadratic interpolation step
                    // xn = b; xnm1=a; xnm2 = c;

                    x = (a * fb * fc) / ((fa - fb) * (fa - fc)) +
                        (b * fa * fc) / ((fb - fa) * (fb - fc)) +
                        (c * fa * fb) / ((fc - fa) * (fc - fb));

                } else { // fa != fb here because Brent requires fa,fb to have opposite signs
                    // secant step
                    // xn = b; xnm1=a; 
                    x = b - fb * (b - a) / (fb - fa);
                }

                if (((3 * a + b) / 4 <= x && x <= b) ||
                        (usedBisection && (Math.Abs(x - b) >= Math.Abs(b - c) / 2 || areNear(b, c))) ||  //Math.Abs(b-c) < tol.AbsTolerance 
                        (!usedBisection && (Math.Abs(x - b) >= Math.Abs(c - d) / 2) || areNear(c, d))) {  //Math.Abs(c-d) < tol.AbsTolerance
                    // bisection step
                    x = (a + b) / 2;
                    usedBisection = true;
                } else {
                    usedBisection = false;
                }

                double fx = f(x);
                iterations++;
                if (fx == 0) {
                    xSolution = x;
                    break;
                }

                if (double.IsNaN(fx)) {
                    Policies.ReportDomainError("Requires that function returns a value: f({0}) = {1}", x, fx);
                    break;
                }


                d = c; df = fc;
                c = b; fc = fb;

                if (IsBracket(fa, fx)) {
                    b = x; fb = fx;
                } else {
                    a = x; fa = fx;
                }

                // a is always further away from the root
                if (Math.Abs(fa) < Math.Abs(fb)) {
                    double t, ft;
                    t = b; ft = fb;
                    b = a; fb = fa;
                    a = t; fa = ft;
                }

                if (areNear(a, b)) {
                    xSolution = x = a + 0.5 * (b - a); // don't have an exact solution, so return the midpoint
                    break;
                }


            }

            min = a; fmin = fa;
            max = b; fmax = fb;
        exit:
            RootResults result = new RootResults();
            if (double.IsNaN(xSolution))
                result.LastX = x;
            else
                result.SolutionX = xSolution;

            result.Bracket = new RootBracketResult(min, fmin, max, fmax);
            result.Iterations = iterations;
            return result;
        }

        /// <summary>
        /// Brent-Dekker uses a combination of inverse quadratic interpolation, secant interpolation, and bisection to locate the root of f(x).
        /// Note: min and max must bracket the root.
        /// </summary>
        /// <param name="f">function to solve</param>
        /// <param name="min">root must lie within [min, max]</param>
        /// <param name="max">root must lie within [min, max]</param>
        /// <param name="areNear">A function that returns true when the required tolerance is reached, or null to use the default tolerance</param>
        /// <param name="maxIterations">The maximum number of iterations to keep refining the root estimate</param>
        /// <returns>
        /// Null, if any of the parameters were incorrect.
        /// The RootResults which contains the root if one was found, or the final state before failing.
        /// </returns>
        public static RootResults Brent(Func<double, double> f, double min, double max, ToleranceFunc areNear, int maxIterations)
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

            // use the default if tolerance is null
            if (areNear == null)
                areNear = GetToleranceFunc();

            if (maxIterations < 3) {
                Policies.ReportDomainError("Requires maxIterations >= 3: maxIterations = {0}", maxIterations);
                return nullResult;
            }

            double fmin = f(min);
            double fmax = f(max);

            if (!IsBracket(fmin, fmax)) {
                RootResults results = new RootResults() {
                    Iterations = 2,
                    Bracket = new RootBracketResult(min, fmin, max, fmax)
                };
                Policies.ReportDomainError("The root is not bracketed: f({0}) = {1}; f({2}) = {3}", min, fmin, max, fmax);
                return results;
            }


            var r = BrentInt(f, min, fmin, max, fmax, areNear, maxIterations - 2);
            r.Iterations += 2;
            return r;

        }

        /// <summary>
        /// Brent-Dekker uses a combination of inverse quadratic interpolation, secant interpolation, and bisection to locate the root of f(x).
        /// Note: min and max must bracket the root.
        /// </summary>
        /// <param name="f">function to solve</param>
        /// <param name="min">root must lie within [min, max]</param>
        /// <param name="max">root must lie within [min, max]</param>
        /// <returns>
        /// Null, if any of the parameters were incorrect.
        /// The RootResults which contains the root if one was found, or the final state before failing.
        /// </returns>
        public static RootResults Brent(Func<double, double> f, double min, double max)
        {
            return Brent(f, min, max, null, Policies.MaxRootIterations);
        }

        /// <summary>
        /// Brent-Dekker uses a combination of inverse quadratic interpolation, secant interpolation, and bisection to locate the root of f(x).
        /// Note: min and max must bracket the root.
        /// </summary>
        /// <param name="f">function to solve</param>
        /// <param name="min">root must lie within [min, max]</param>
        /// <param name="max">root must lie within [min, max]</param>
        /// <param name="relTolerance">The maximum relative tolerance</param>
        /// <param name="absTolerance">The maximum absolute tolerance</param>
        /// <param name="maxIterations">The maximum number of iterations to keep refining the root estimate</param>
        /// <returns>
        /// Null, if any of the parameters were incorrect.
        /// The RootResults which contains the root if one was found, or the final state before failing.
        /// </returns>
        public static RootResults Brent(Func<double, double> f, double min, double max, double relTolerance, double absTolerance, int maxIterations)
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
            return Brent(f, min, max, tol, maxIterations);
        }


        /// <summary>
        /// Brent-Dekker method for finding roots. Robust and relatively quick.
        /// Uses a combination of inverse quadratic interpolation, secant interpolation, and bisection.
        /// Note: Brent-Dekker requires that the root be bracketed. This routine will attempt to 
        /// bracket the root before using Brent-Dekker.
        /// </summary>
        /// <param name="f">function to solve</param>
        /// <param name="guess">the initial guess for the root.</param>
        /// <param name="step">the starting step size</param>
        /// <param name="fType">Unknown, increasing, or decreasing</param>
        /// <param name="min">Optional. The minimum location for the root, or NaN</param>
        /// <param name="max">Optional. The maximum location for the root, or NaN</param>
        /// <param name="areNear">Optional. The maximum tolerance, or null to use the default tolerance</param>
        /// <param name="maxIterations">The maximum number of iterations to keep refining the root estimate</param>
        /// <returns>
        /// Null, if any of the parameters were incorrect.
        /// The RootResults which contains the root if one was found, or the final state before failing.
        /// </returns>
        public static RootResults BrentBracket(Func<double, double> f, double guess, double step, FunctionShape fType, double min, double max, ToleranceFunc areNear, int maxIterations)
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
                RootResults result = new RootResults();
                result.Bracket = b;
                result.Iterations = iterations;
                return result;
            }

            return BrentInt(f, b.XMin, b.FxMin, b.XMax, b.FxMax, areNear, maxIterations - iterations);
        }

        /// <summary>
        /// Brent-Dekker method for finding roots. Robust and relatively quick.
        /// Uses a combination of inverse quadratic interpolation, secant interpolation, and bisection.
        /// Note: Brent-Dekker requires that the root be bracketed. This routine will attempt to 
        /// bracket the root before using Brent-Dekker.
        /// </summary>
        /// <param name="f">function to solve</param>
        /// <param name="guess">the initial guess for the root.</param>
        /// <param name="step">the starting step size</param>
        /// <param name="fType">Unknown, increasing, or decreasing</param>
        /// <param name="min">Optional. The minimum location for the root, or NaN</param>
        /// <param name="max">Optional. The maximum location for the root, or NaN</param>
        /// <returns>
        /// Null, if any of the parameters were incorrect.
        /// The RootResults which contains the root if one was found, or the final state before failing.
        /// </returns>
        public static RootResults BrentBracket(Func<double, double> f, double guess, double step, FunctionShape fType, double min, double max)
        {
            return BrentBracket(f, guess, step, fType, min, max, null, Policies.MaxRootIterations);

        }


        /// <summary>
        /// Brent-Dekker method for finding roots. Robust and relatively quick.
        /// Uses a combination of inverse quadratic interpolation, secant interpolation, and bisection.
        /// <para>Note: Brent-Dekker requires that the root be bracketed. This routine will attempt to 
        /// bracket the root before using Brent-Dekker.</para>
        /// </summary>
        /// <param name="f">function to solve</param>
        /// <param name="guess">the initial guess for the root.</param>
        /// <param name="step">the starting step size</param>
        /// <param name="fType">Unknown, increasing, or decreasing</param>
        /// <param name="min">Optional. The minimum location for the root, or NaN</param>
        /// <param name="max">Optional. The maximum location for the root, or NaN</param>
        /// <param name="relTolerance">The maximum relative tolerance</param>
        /// <param name="absTolerance">The maximum absolute tolerance</param>
        /// <param name="maxIterations">The maximum number of iterations to keep refining the root estimate</param>
        /// <returns>
        /// Null, if any of the parameters were incorrect.
        /// The RootResults which contains the root if one was found, or the final state before failing.
        /// </returns>
        public static RootResults BrentBracket(Func<double, double> f, double guess, double step, FunctionShape fType, double min, double max, double relTolerance, double absTolerance, int maxIterations)
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
            return BrentBracket(f, guess, step, fType, min, max, areNear, maxIterations);

        }
    }

}
