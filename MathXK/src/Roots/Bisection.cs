//  Copyright (c) 2013 Rocco Romeo
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

//  Ported to .NET based on the original:
//      Copyright (C) John Maddock 2006, Boost Software License, Version 1.0


using System;
using System.Diagnostics;

namespace MathXK.Roots
{

    public static partial class RootFinder
    {

        // See: http://en.wikipedia.org/wiki/Bisection_method

        private static RootResults BisectionInt(Func<double, double> f, double min, double fmin, double max, double fmax, ToleranceFunc areNear, int maxIterations)
        {
            // private routine so these conditions should already be true
            Debug.Assert(f != null);
            Debug.Assert(areNear != null);
            Debug.Assert(min < max);
            Debug.Assert(IsBracket(fmin, fmax));
            Debug.Assert(maxIterations > 0);

            // check the initial conditions

            double xSolution = double.NaN;
            int iterations = 0;
            double mid = double.NaN;

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

            while (iterations < maxIterations) {

                // evaluate the function at the midpoint
                mid = min + (max - min) / 2;
                double fmid = f(mid);
                iterations++;

                if (fmid == 0) {
                    xSolution = mid;
                    break;
                }



                if (IsBracket(fmin, fmid)) {
                    max = mid;
                    fmax = fmid;
                } else {
                    min = mid;
                    fmin = fmid;
                }

                if (areNear(max, min)) {
                    xSolution = min + (max - min) / 2.0; // don't have an exact solution so take a midpoint
                    break;
                }


            }

        exit:
            RootResults result = new RootResults();
            if (double.IsNaN(xSolution))
                result.LastX = mid;
            else
                result.SolutionX = xSolution;

            result.Bracket = new RootBracketResult(min, fmin, max, fmax);
            result.Iterations = iterations;
            return result;
        }


        /// <summary>
        /// Uses the Bisection method to locate the root of f(x).
        /// Note: min and max must bracket the root.
        /// </summary>
        /// <param name="f">Function to solve</param>
        /// <param name="min">Root must lie within [min, max]</param>
        /// <param name="max">Root must lie within [min, max]</param>
        /// <param name="areNear">A function that returns true when the required tolerance is reached, or null to use the default tolerance</param>
        /// <param name="maxIterations">The maximum number of iterations to keep refining the root estimate</param>
        /// <returns>
        /// Null, if any of the parameters were incorrect.
        /// The RootResults which contains the root if one was found, or the final state before failing.
        /// </returns>
        public static RootResults Bisection(Func<double, double> f, double min, double max, ToleranceFunc areNear, int maxIterations)
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
            if (maxIterations < 3) {
                Policies.ReportDomainError("Requires maxIterations >= 3: maxIterations = {0}", maxIterations);
                return nullResult;
            }

            // use the default if tolerance is null
            if (areNear == null)
                areNear = GetToleranceFunc();


            // the root must be bracketed
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


            var r = BisectionInt(f, min, fmin, max, fmax, areNear, maxIterations - 2);
            r.Iterations += 2;
            return r;
        }

        /// <summary>
        /// Uses the Bisection method to locate the root of f(x).
        /// Note: min and max must bracket the root.
        /// </summary>
        /// <param name="f">Function to solve</param>
        /// <param name="min">The minimum x-value. Root must lie within [min, max]</param>
        /// <param name="max">The maximum x-value. Root must lie within [min, max]</param>
        /// <param name="relTolerance">The maximum relative tolerance</param>
        /// <param name="absTolerance">The maximum absolute tolerance</param>
        /// <param name="maxIterations">The maximum number of iterations to keep refining the root estimate</param>
        /// <returns>
        /// Null, if any of the parameters were incorrect.
        /// The RootResults which contains the root if one was found, or the final state before failing.
        /// </returns>
        public static RootResults Bisection(Func<double, double> f, double min, double max, double relTolerance, double absTolerance, int maxIterations)
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
            return Bisection(f, min, max, areNear, maxIterations);
        }



        /// <summary>
        /// Uses the Bisection method to locate the root of f(x).
        /// Note: min and max must bracket the root.
        /// </summary>
        /// <param name="f">Function to solve</param>
        /// <param name="min">root must lie within [min, max]</param>
        /// <param name="max">root must lie within [min, max]</param>
        /// <returns>
        /// Null, if any of the parameters were incorrect.
        /// The RootResults which contains the root if one was found, or the final state before failing.
        /// </returns>
        public static RootResults Bisection(Func<double, double> f, double min, double max)
        {
            return Bisection(f, min, max, null, Policies.MaxRootIterations);
        }
    }


}
