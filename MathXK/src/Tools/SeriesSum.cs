//  Copyright (c) 2013 Rocco Romeo
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

//  Ported to .NET based on the original:
//      Copyright (c) John Maddock 2005-2006, Boost Software License, Version 1.0


using System;

namespace MathXK
{
    /// <summary>
    /// Routines for evaluating a series sum
    /// </summary>
    public static class SeriesSum
    {
        /// <summary>
        /// Returns the sum of a series with each term given by the function f(). 
        /// </summary>
        /// <param name="f">Delegate returning the next term</param>
        /// <param name="initialValue">initial value</param>
        /// <param name="tolerance">Stops when the relative tolerance is reached</param>
        /// <param name="maxIterations">The maximum number of terms</param>
        /// <returns>
        /// The sum of the series or NaN if the series does not converge within the maximum
        /// number of iterations specified in the policy.
        /// </returns>
        public static double Eval(Func<double> f, double initialValue, double tolerance, uint maxIterations)
        {
            if (f == null) {
                Policies.ReportDomainError("Requires f != null");
                return double.NaN;
            }
            if (!(tolerance >= 0)) {
                Policies.ReportDomainError("Requires relative tolerance >= 0");
                return double.NaN;
            }


            double sum = initialValue;
            for (UInt32 i = 0; i < maxIterations; i++) {
                double prevSum = sum;
                double delta = f();
                sum += delta;

                if (Math.Abs(delta) <= Math.Abs(prevSum) * tolerance) {
                    return sum;
                }

            }

            Policies.ReportConvergenceError("Series did not converge in {0} iterations", maxIterations);
            return double.NaN;
        }


        /// <summary>
        /// Returns the sum of a series with each term given by the function f(). 
        /// </summary>
        /// <param name="f">Delegate returning the next term in the series</param>
        /// <param name="initialValue">initial value</param>
        /// <param name="tolerance">Stops when the relative tolerance is reached</param>
        /// <returns>
        /// The sum of the series or NaN if the series does not converge within the maximum
        /// number of iterations specified in the policy.
        /// </returns>
        public static double Eval(Func<double> f, double initialValue, double tolerance)
        {
            return Eval(f, initialValue, tolerance, Policies.MaxSeriesIterations);
        }


        /// <summary>
        /// Returns the sum of a series with each term given by the function f(). 
        /// </summary>
        /// <param name="f">Delegate returning the next term in the series</param>
        /// <param name="initialValue">initial value</param>
        /// <returns>
        /// The sum of the series or NaN if the series does not converge within the maximum
        /// number of iterations specified in the policy.
        /// </returns>
        public static double Eval(Func<double> f, double initialValue = 0)
        {
            return Eval(f, initialValue, Policies.SeriesTolerance, Policies.MaxSeriesIterations);
        }

    }

}
