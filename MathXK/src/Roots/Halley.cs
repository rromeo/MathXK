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

        /// <summary>
        /// Halley's method for finding roots.
        /// Note: min and max must bracket the root.
        /// </summary>
        /// <param name="f">function to solve along with first and second derivative at point x</param>
        /// <param name="guess">the initial guess for the root</param>
        /// <param name="min">root must lie within [min, max]</param>
        /// <param name="max">root must lie within [min, max]</param>
        /// <param name="areNear">A function that returns true when the required tolerance is reached, or null to use the default tolerance</param>
        /// <param name="maxIterations">The maximum number of iterations to keep refining the root estimate</param>
        /// <returns>
        /// Null, if any of the parameters were incorrect.
        /// The RootResults which contains the root if one was found, or the final state before failing.
        /// </returns>
        public static RootResults Halley(Func<double, ValueTuple<double, double, double>> f, double guess, double min, double max, ToleranceFunc areNear, int maxIterations)
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
                Policies.ReportDomainError("Requires a valid initial guess between min and max: guess = {0}; min = {1}; max = {2}", guess, min, max);
                return nullResult;
            }


            double xSolution = double.NaN;
            int iterations = 0;



            double f0, f1, f2;
            f0 = double.NaN;

            double x = guess;

            double delta = max - min;
            double delta1 = double.MaxValue;
            double delta2 = double.MaxValue;

            bool out_of_bounds_sentry = false;


            while (iterations < maxIterations) {
                double last_f0 = f0;
                delta2 = delta1;
                delta1 = delta;

                // get the value and first derivative
                iterations++;
                (f0, f1, f2) = f(x);

                if (f0 == 0) {
                    xSolution = x;
                    break;
                }

                if (double.IsNaN(f0)) {
                    Policies.ReportDomainError("Requires that function returns a value: f({0}) = {1}", x, f0);
                    break;
                }


                if ( f1 == 0 ) {
                    // Oops zero derivative!!!

                    // this is initially null
                    if (double.IsNaN(last_f0)) {
                        // this must be the first iteration, pretend that we had a
                        // previous one at either min or max:
                        x = (guess == min) ? max : min;

                        iterations++;
                        (f0, f1, f2) = f(x);

                        last_f0 = f0;

                        delta = x - guess;
                    }
                    if (IsBracket(last_f0, f0)) {
                        // we've crossed over so move in opposite direction to last step:
                        if (delta < 0) {
                            delta = (x - min) / 2;
                        } else {
                            delta = (x - max) / 2;
                        }

                    } else {
                        // move in same direction as last step:
                        if (delta < 0) {
                            delta = (x - max) / 2;
                        } else {
                            delta = (x - min) / 2;
                        }
                    }

                } else {

                    // Use Halley's method:
                    // Let f = f(x{n}); f' = f'(x{n}); f'' = f''(x{n});
                    //
                    // x{n+1} = x{n} - (2 * f * f')/(2 * f'^2 -  f * f'' ) 
                    // OR
                    // x{n+1} = x{n} - (2 * f)/(2 * f' -  f * f''/ f' ) 
                    // OR
                    // x{n+1} = x{n} - (f/f')/(1 - (1/2)*(f/f')(f''/f')) 


                    if (f2 != 0 && !double.IsInfinity(f2)) {
                        double num = 2 * f0;
                        double denom = 2 * f1 - f0 * (f2 / f1);

                        delta = num / denom;
                        if (double.IsInfinity(delta) || double.IsNaN(delta)) {
                            // possible overflow, use Newton step:
                            delta = f0 / f1;
                        }

                        if (delta * f1 / f0 < 0) {
                            // probably cancellation error, try a Newton step instead:
                            delta = f0 / f1;
                        }
                    } else
                        delta = f0 / f1;
                }

                double convergence = Math.Abs(delta / delta2);
                if ((convergence > 0.8) && (convergence < 2)) {
                    // last two steps haven't converged, try bisection:
                    delta = (delta > 0) ? (x - min) / 2 : (x - max) / 2;
                    if (Math.Abs(delta) > x)
                        delta = Math.Sign(delta) * x; // protect against huge jumps!
                    // reset delta2 so that this branch will *not* be taken on the
                    // next iteration:
                    delta2 = delta * 3;

                }



                guess = x;
                x -= delta;

                // check for out of bounds step:
                if (x < min) {
                    double diff = ((Math.Abs(min) < 1) && (Math.Abs(x) > 1) && (double.MaxValue / Math.Abs(x) < Math.Abs(min))) ? 1000.0 : x / min;
                    if (Math.Abs(diff) < 1)
                        diff = 1 / diff;
                    if (!out_of_bounds_sentry && (diff > 0) && (diff < 3)) {
                        // Only a small out of bounds step, lets assume that the result
                        // is probably approximately at min:
                        delta = 0.99 * (guess - min);
                        x = guess - delta;
                        out_of_bounds_sentry = true; // only take this branch once!
                    } else {
                        delta = (guess - min) / 2;
                        x = guess - delta;
                    }
                } else if (x > max) {
                    double diff = ((Math.Abs(max) < 1) && (Math.Abs(x) > 1) && (double.MaxValue / Math.Abs(x) < Math.Abs(max))) ? 1000.0 : x / max;
                    if (Math.Abs(diff) < 1)
                        diff = 1 / diff;
                    if (!out_of_bounds_sentry && (diff > 0) && (diff < 3)) {
                        // Only a small out of bounds step, lets assume that the result
                        // is probably approximately at min:
                        delta = 0.99 * (guess - max);
                        x = guess - delta;
                        out_of_bounds_sentry = true; // only take this branch once!
                    } else {
                        delta = (guess - max) / 2;
                        x = guess - delta;
                    }
                }
                // update brackets:
                if (delta > 0)
                    max = guess;
                else
                    min = guess;


                if (areNear(x, guess)) {
                    xSolution = x;
                    break;
                }

            }

            RootResults result = new RootResults()
            {
                SolutionX = xSolution,
                // Note: this halley implementation does not keep track of fmin, fmax
                // for the brackets
                Bracket = new RootBracketResult(min, double.NaN, max, double.NaN),
                Iterations = iterations
            };
            return result;

        }


        /// <summary>
        /// Halley's method for finding roots.
        /// Note: min and max must bracket the root.
        /// </summary>
        /// <param name="f">function to solve along with first and second derivative at point x</param>
        /// <param name="guess">the initial guess for the root</param>
        /// <param name="min">root must lie within [min, max]</param>
        /// <param name="max">root must lie within [min, max]</param>
        /// <returns>
        /// Null, if any of the parameters were incorrect.
        /// The RootResults which contains the root if one was found, or the final state before failing.
        /// </returns>
        public static RootResults Halley(Func<double, ValueTuple<double, double, double>> f, double guess, double min, double max)
        {
            return Halley(f, guess, min, max, null, Policies.MaxRootIterations);
        }

        /// <summary>
        /// Halley's method for finding roots.
        /// Note: min and max must bracket the root.
        /// </summary>
        /// <param name="f">function to solve along with first and second derivative at point x</param>
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
        public static RootResults Halley(Func<double, ValueTuple<double, double, double>> f, double guess, double min, double max, double relTolerance, double absTolerance, int maxIterations)
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
            return Halley(f, guess, min, max, areNear, maxIterations);
        }

    }

}
