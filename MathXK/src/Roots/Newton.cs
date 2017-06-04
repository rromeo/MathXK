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
        /// Newton method for finding roots. Does not require the root to be bracketed.
        /// Note: Only use this routine for well behaved functions within the given limits; NewtonSafe is more robust.
        /// </summary>
        /// <param name="f">function to solve along with first derivative at point x</param>
        /// <param name="guess">the initial guess for the root.</param>
        /// <param name="min">Optional. The minimum location for the root, or NaN</param>
        /// <param name="max">Optional. The maximum location for the root, or NaN</param>
        /// <param name="areNear">A function that returns true when the required tolerance is reached, or null to use the default tolerance</param>
        /// <param name="maxIterations">The maximum number of iterations to keep refining the root estimate</param>
        /// <returns>
        /// Null, if any of the parameters were incorrect.
        /// The RootResults which contains the root if one was found, or the final state before failing.
        /// </returns>
        public static RootResults Newton(Func<double, (double, double)> f, double guess, double min, double max, ToleranceFunc areNear, int maxIterations)
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

            // NaN is valid here for min or max, so only fail if both have a value and min > max
            if (min > max) {
                Policies.ReportDomainError("Requires min <= max : min = {0}; max = {1}", min, max);
                return nullResult;
            }

            if (double.IsNaN(guess) || guess < min || guess > max) {
                Policies.ReportDomainError("Requires a valid initial guess between min and max: guess = {0}; min = {1}; max = {2}", guess, min, max);
                return nullResult;
            }

            double xSolution = double.NaN;
            int iterations = 0;

            // We have not tried f(fmin) or f(fmax), so we do not know if the root is truly bracketed.
            double x = guess;
            double fmin = double.NaN; 
            double fmax = double.NaN;

            while (iterations < maxIterations) {

                iterations++;

                // get the value and first derivative
                var (f0, f1) = f(x);

                if (f0 == 0) {
                    xSolution = x;
                    break;
                }

                if (double.IsNaN(f0) || double.IsNaN(f1)) {
                    Policies.ReportDomainError("Requires that function returns a value: f({0}) = ({1}, {2})", x, f0, f1);
                    break;
                }

                // check for zero derivative
                if (f1 == 0)
                    break;

                // Newton's method: x_{n+1} = x_n - f(x_n)/f'(x_n)
                double last_x = x;
                double delta = (f0 / f1);
                x -= delta;

                // update the local minimum, maximum
                if (delta > 0) {
                    max = last_x; fmax = f0;
                } else {
                    min = last_x; fmin = f0;
                }

                // if we have outer limits on the root, enforce them 
                if (x < min) {
                    if (double.IsNaN(fmin))
                        x = min;
                    else {
                        // make sure we have a bracket
                        if ( !IsBracket(f0, fmin))
                            break;
                        // we've already used the min, so let's try a bisection
                        x = 0.5 * (last_x + min); 
                    }
                } else if (x > max) {
                    if (double.IsNaN(fmax))
                        x = max;
                    else {
                        // make sure we have a bracket
                        if ( !IsBracket(f0, fmax))
                            break;
                        // we've already used the max so, let's try a bisection
                        x = 0.5*(last_x + max); 
                    }
                }

                if (areNear(last_x, x)) {
                    // We can't be sure we have a root here because the root isn't bracketed.
                    // Assume the function is well behaved within the given (min, max) limits.
                    xSolution = x;
                    break;
                }

            }

            RootResults result = new RootResults();
            if (double.IsNaN(xSolution))
                result.LastX = x;
            else
                result.SolutionX = x;
            result.Bracket = new RootBracketResult(min, fmin, max, fmax);
            result.Iterations = iterations;
            return result;

        }

        /// <summary>
        /// Newton method for finding roots.
        /// Does not require the root to be bracketed.
        /// </summary>
        /// <param name="f">function to solve along with first derivative at point x</param>
        /// <param name="guess">the initial guess for the root.</param>
        /// <param name="min">Optional. The minimum location for the root, or NaN</param>
        /// <param name="max">Optional. The maximum location for the root, or NaN</param>
        /// <returns>
        /// Null, if any of the parameters were incorrect.
        /// The RootResults which contains the root if one was found, or the final state before failing.
        /// </returns>
        public static RootResults Newton(Func<double, (double, double)> f, double guess, double min, double max)
        {
            return Newton(f, guess, min, max, null, Policies.MaxRootIterations);
        }

        /// <summary>
        /// Newton method for finding roots.
        /// Does not require the root to be bracketed.
        /// </summary>
        /// <param name="f">function to solve along with first derivative at point x</param>
        /// <param name="guess">the initial guess for the root.</param>
        /// <param name="min">Optional. The minimum location for the root, or NaN</param>
        /// <param name="max">Optional. The maximum location for the root, or NaN</param>
        /// <param name="relTolerance">The maximum relative tolerance</param>
        /// <param name="absTolerance">The maximum absolute tolerance</param>
        /// <param name="maxIterations">The maximum number of iterations to keep refining the root estimate</param>
        /// <returns>
        /// Null, if any of the parameters were incorrect.
        /// The RootResults which contains the root if one was found, or the final state before failing.
        /// </returns>
        public static RootResults Newton(Func<double, (double, double)> f, double guess, double min, double max, double relTolerance, double absTolerance, int maxIterations)
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
            return Newton(f, guess, min, max, areNear, maxIterations);
        }

        /// <summary>
        /// Newton method for finding roots.
        /// Does not require the root to be bracketed.
        /// </summary>
        /// <param name="f">function to solve along with first derivative at point x</param>
        /// <param name="guess">the initial guess for the root</param>
        /// <returns>
        /// Null, if any of the parameters were incorrect.
        /// The RootResults which contains the root if one was found, or the final state before failing.
        /// </returns>
        public static RootResults Newton(Func<double, (double, double)> f, double guess)
        {
            return Newton(f, guess, double.NaN, double.NaN);
        }


    }


}