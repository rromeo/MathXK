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

        // For more information see: http://en.wikipedia.org/wiki/Newton's_method

        /// <summary>
        /// The Newton Safe method combines the Newton–Raphson method with the bisection method.
        /// Unlike the Newton method, NewtonSafe requires that the root be bracketed.
        /// </summary>
        /// <param name="f">function to solve along with first derivative at point x</param>
        /// <param name="guess">the initial guess for the root</param>
        /// <param name="min">root must lie within [min, max]</param>
        /// <param name="max">root must lie within [min, max]</param>
        /// <param name="areNear">A function that returns true when the required tolerance is reached, or null to use the default tolerance</param>
        /// <param name="maxIterations">The maximum number of iterations to keep refining the root estimate</param>
        /// <returns>
        /// Null, if any of the parameters were incorrect.
        /// The RootResults which contains the root if one was found, or the final state before failing.
        /// </returns>
        public static RootResults NewtonSafe(Func<double, (double, double)> f, double guess, double min, double max, ToleranceFunc areNear, int maxIterations)
        {
            const RootResults nullResult = null;

            if (f == null) {
                Policies.ReportDomainError("Function cannot be null");
                return nullResult;
            }

            // use the default if tolerance is null
            if (areNear == null)
                areNear = GetToleranceFunc();

            if (maxIterations < 3) {
                Policies.ReportDomainError("Requires maxIterations >= 3: maxIterations = {0}", maxIterations);
                return nullResult;
            }

            if (!(min <= max)) {
                Policies.ReportDomainError("Requires min <= max : min = {0}; max = {1}", min, max);
                return nullResult;
            }
            if (!(guess >= min && guess <= max)) {
                Policies.ReportDomainError("Guess must be between min and max: guess = {0}; min = {1}; max = {2}", guess, min, max);
                return nullResult;
            }


            double xSolution = double.NaN;
            int iterations = 0;

            // check to make sure the root is bracketed
            double fmin = f(min).Item1;
            double fmax = f(max).Item1;
            iterations = 2;
            if (!IsBracket(fmin, fmax)) {
                RootResults results = new RootResults() {
                    Iterations = 2,
                    Bracket = new RootBracketResult(min, fmin, max, fmax)
                };
                Policies.ReportDomainError("The root is not bracketed: f({0}) = {1}; f({2}) = {3}", min, fmin, max, fmax);
                return results;
            }


            if (fmin == 0) {
                xSolution = min;
                goto exit;
            }
            if (fmax == 0) {
                xSolution = max;
                goto exit;
            }

            // check if the bracket size is within the tolerance
            // if it is, we don't have an exact solution so take a midpoint
            if (areNear(min, max)) {
                xSolution = min + 0.5 * (max - min);
                goto exit;
            }



            double delta = max - min;
            double delta1 = double.MaxValue;
            double delta2 = double.MaxValue;


            double x = guess;
            bool usedBisection = true; // in this case delta, delta1, delta2 are invalid

            while (iterations < maxIterations) {

                delta2 = delta1;
                delta1 = delta;

                // get the value and first derivative
                var (f0, f1) = f(x);
                iterations++;


                if (f0 == 0) {
                    xSolution = x;
                    break;
                }

                if (double.IsNaN(f0)) {
                    Policies.ReportDomainError("Requires that function returns a value: f({0}) = {1}", x, f0);
                    break;
                }


                // update our bracket
                if (IsBracket(f0, fmax)) {
                    min = x; fmin = f0;
                } else {
                    Debug.Assert(IsBracket(fmin, f0));
                    max = x; fmax = f0;
                }


                double midpoint = min + 0.5 * (max - min);
                if (areNear(min, max)) {
                    xSolution = midpoint; // don't have an exact solution so take a midpoint
                    break;
                }


                if (f1 == 0) {
                    // derivative == zero
                    x = midpoint;
                    usedBisection = true;
                    continue;
                }


                // Newton's method: x_{n+1} = x_n - f(x_n)/f'(x_n)
                delta = f0 / f1;
                x -= delta;


                // if the next x is outside of (min,max), use bisection
                if (!(x > min && x < max)) {
                    x = midpoint;
                    usedBisection = true;
                    continue;
                }

                // or if the last two steps haven't converged fast enough
                // then try bisection:

                if (!usedBisection && (Math.Abs(delta * 2) > Math.Abs(delta2))) {

                    x = midpoint;
                    usedBisection = true;
                } else
                    usedBisection = false;
            }

        exit:

            RootResults result = new RootResults();
            result.SolutionX = xSolution;
            result.Bracket = new RootBracketResult(min, fmin, max, fmax);
            result.Iterations = iterations;
            return result;

        }


        /// <summary>
        /// Newton's method for finding roots.
        /// Note: min and max must bracket the root.
        /// </summary>
        /// <param name="f">function to solve along with first derivative at point x</param>
        /// <param name="guess">the initial guess for the root</param>
        /// <param name="min">root must lie within [min, max]</param>
        /// <param name="max">root must lie within [min, max]</param>
        /// <returns>
        /// Null, if any of the parameters were incorrect.
        /// The RootResults which contains the root if one was found, or the final state before failing.
        /// </returns>
        public static RootResults NewtonSafe(Func<double, (double, double)> f, double guess, double min, double max)
        {
            return NewtonSafe(f, guess, min, max, null, Policies.MaxRootIterations);

        }

        /// <summary>
        /// Newton's method for finding roots.
        /// Note: min and max must bracket the root.
        /// </summary>
        /// <param name="f">function to solve along with first derivative at point x</param>
        /// <param name="guess">the initial guess for the root</param>
        /// <param name="min">root must lie within [min, max]</param>
        /// <param name="max">root must lie within [min, max]</param>
        /// <param name="relTolerance">The maximum relative tolerance</param>
        /// <param name="absTolerance">The maximum absolute tolerance</param>
        /// <param name="maxIterations">The maximum number of iterations to keep refining the root estimate</param>
        /// <returns>
        /// Null, if any of the parameters were incorrect.
        /// The RootResults which contains the root if one was found, or the final state before failing.
        /// </returns>
        public static RootResults NewtonSafe(Func<double, (double, double)> f, double guess, double min, double max, double relTolerance, double absTolerance, int maxIterations)
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
            return NewtonSafe(f, guess, min, max, areNear, maxIterations);

        }

    }

}
