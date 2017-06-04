//  Copyright (c) 2013 Rocco Romeo
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)


using System;
using System.Diagnostics;

namespace MathXK.Roots
{

    /// <summary>
    /// Represents the bracket (possibly) containing the root
    /// </summary>
    public class RootBracketResult
    {
        /// <summary>
        /// Construct a bracket such that the root lies between the points [xMin, xMax]
        /// Note: fxMin = f(xMin) and fxMax = f(xMax).
        /// </summary>
        /// <param name="xMin">The minimum x position of the bracket</param>
        /// <param name="fxMin">The corresponding y position, f(xMin)</param>
        /// <param name="xMax">The maximum x position of the bracket</param>
        /// <param name="fxMax">The corresponding y position, f(xMax)</param>
        public RootBracketResult(double xMin, double fxMin, double xMax, double fxMax)
        {
            XMin = xMin;
            FxMin = fxMin;

            XMax = xMax;
            FxMax = fxMax;
        }


        /// <summary>
        /// Gets the lower limit for x.
        /// </summary>
        public double XMin { get; internal set; }

        /// <summary>
        /// Gets f(xMin).
        /// </summary>
        public double FxMin { get; internal set; }

        /// <summary>
        /// Gets the upper limit for x.
        /// </summary>
        public double XMax { get; internal set; }

        /// <summary>
        /// Gets f(xMax).
        /// </summary>
        public double FxMax { get; internal set; }

        /// <summary>
        /// Returns true if [xmin, xmax] contains a root(zero).
        /// </summary>
        public bool IsValid { get { return IsBracket(FxMin, FxMax); } }

        /// <summary>
        /// Returns true if the interval [fa, fb] contains a root(zero)
        /// </summary>
        /// <param name="fa">f(a)</param>
        /// <param name="fb">f(b)</param>
        /// <returns></returns>
        private static bool IsBracket(double fa, double fb)
        {
            // note: Math.Sign(x) throws an exception if x == NaN. 
            // This routine will return false.

            return ((fa >= 0 && fb <= 0) || (fa <= 0 && fb >= 0));
        }

    };


    /// <summary>
    /// Provides a hint to the root finding functions
    /// <para>Is it increasing or decreasing? Or Unknown</para>
    /// </summary>
    public enum FunctionShape {
        /// <summary>
        /// Function behavior is unknown
        /// </summary>
        Unknown,

        /// <summary>
        /// Increasing function
        /// </summary>
        Increasing,

        /// <summary>
        /// Decreasing Function
        /// </summary>
        Decreasing
    };

    public static partial class RootFinder
    {

        /// <summary>
        /// Tries to locate a bracket given an initial guess and a step. 
        /// Routine walks downhill/uphill taking steps of increasing size until the function changes sign.
        /// </summary>
        /// <param name="f">function to solve</param>
        /// <param name="guess">initial guess of x value of the root</param>
        /// <param name="step">the starting step size</param>
        /// <param name="fType">unknown, increasing, or decreasing</param>
        /// <param name="min">Optional. The minimum limit on x if it exists, otherwise NaN</param>
        /// <param name="max">Optional. The maximum limit on x if it exists, otherwise NaN</param>
        /// <param name="maxIterations">The maximum number of iterations to find the bracket</param>
        /// <param name="nIterations">Output. Return the number of iterations used</param>
        /// <returns>
        /// Null, if any of the parameters were incorrect.
        /// The bracket if one was found, otherwise the last attempt.
        /// </returns>
        public static RootBracketResult FindBracket(Func<double, double> f, double guess, double step, FunctionShape fType, double min, double max, int maxIterations, out int nIterations)
        {
            const RootBracketResult nullResult = null;
            nIterations = 0;

            if (f == null) {
                Policies.ReportDomainError("Function cannot be null");
                return nullResult;
            }
            if (min > max) {
                Policies.ReportDomainError("Requires min <= max : min = {0}; max = {1}", min, max);
                return nullResult;
            }
            if (maxIterations < 3) {
                Policies.ReportDomainError("Iterations must be > 3: maxIterations = {0}", maxIterations);
                return nullResult;
            }
            if (step == 0 || double.IsNaN(step) ) {
                Policies.ReportDomainError("Requires |Step| > 0 : step = {0}", step);
                return nullResult;
            }


            Interval bounds = new Interval(min, max);
            if (guess != bounds.Clip(guess)) {
                Policies.ReportDomainError("Guess must be between min and max: guess = {0}; min = {1}; max = {2}", guess, min, max);
                return nullResult;
            }


            // choose a step multiplier > 1 so that we don't go back and forth to
            // the same locations
            const double stepMultiplier = 1.5;

            double x = guess;
            double fx = f(guess);
            int iterations = 1;
            if (fx == 0 || double.IsNaN(fx)) {
                nIterations = iterations;
                return new RootBracketResult(x, fx, x, fx);
            }


            // the bracket is contained within [a, b], where a <= b
            // fa = f(a), fb = f(b)
            double a, fa, b, fb;
            bool force = false;

            // set initial values at guess and guess+step
            // If we don't know what the function looks like
            // count on the user to tell us with the sign of step
            if (fType == FunctionShape.Increasing) {
                // if ( fx < 0 ) positive step, else negative step
                force = true;
                step = -Math.Sign(fx)*Math.Abs(step); 
            } else if (fType == FunctionShape.Decreasing ) {
                // if ( fx < 0 ) negative step, else positive step
                force = true;
                step = Math.Sign(fx) * Math.Abs(step);
            }

            a = x; fa = fx;
            b = x; fb = fx;

            iterations = 1;
            while (iterations < maxIterations) {

                Debug.Assert(step != 0, "FindBracket: Step == 0");

                if ( step > 0  ) {
                    if ( x == max ) {
                        // We've gone as far as we can in this direction
                        // Maybe the root is in the other direction?
                        if ( !force ) {
                            step = -step;
                            force = true;
                            continue;
                        }

                        // we've hit both limits, 
                        // return the last (unsuccessful) bracket
                        nIterations = iterations;
                        return new RootBracketResult(a, fa, b, fb);
                    }

                    // is our step size big enough
                    double nextx = bounds.Clip(x + step);
                    if ( nextx == x ) {
                        step = Math.Abs(x * DoubleLimits.MachineEpsilon);
                        nextx = bounds.Clip(x + step);
                    }

                    x = nextx;
                    fx = f(x);
                    iterations++;

                    a = b; fa = fb;
                    b = x; fb = fx;

                } else {
                    // step < 0

                    if ( x == min ) {
                        // We've gone as far as we can in this direction
                        // Maybe the root is in the other direction?
                        if (!force) {
                            step = -step;
                            force = true;
                            continue;
                        }

                        // We've hit both limits, 
                        // return the last (unsuccessful) bracket
                        nIterations = iterations;
                        return new RootBracketResult(a, fa, b, fb);
                    }

                    // is our step size big enough
                    double nextx = bounds.Clip(x + step);
                    if ( nextx == x ) {
                        step = -Math.Abs(x * DoubleLimits.MachineEpsilon);
                        nextx = bounds.Clip(x + step);
                    }

                    x = nextx;
                    fx = f(x);
                    iterations++;

                    b = a; fb = fa;
                    a = x; fa = fx;

                }

                // check to see if the root is bracketed
                if (IsBracket(fa, fb)) {
                    nIterations = iterations;
                    return new RootBracketResult(a, fa, b, fb);
                }

                // at this point a and b have the same sign
                // otherwise the root would have been bracketed


                // if it is a rising function (i.e. fa < fb )
                //      if both are in negative territory then root > b (move b right)
                //      if both are in positive territory then root < a (move a left)
                // else
                //      if both are in negative territory then root < a (move a left)
                //      if both are in positive territory then root > b (move b right)

                // if fa == fb, continue with the same direction
                if (!force) {
                    double astep = Math.Abs(step);
                    if (fa < fb) {
                        if (fa < 0) {
                            step = +astep;
                        } else {
                            step = -astep;
                        }
                    } else if (fa > fb) {
                        if (fa < 0) {
                            step = -astep;
                        } else {
                            step = +astep;
                        }
                    }
                }

                // set our new step size
                step *= stepMultiplier;

            }

            // no bracket found -- return the (unsuccessful) bracket
            nIterations = iterations;
            return new RootBracketResult(a, fa, b, fb);

        }
    }

}
