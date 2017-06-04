//  (C) Copyright Rocco Romeo 2013.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)


using System;
using System.Diagnostics;
using System.Runtime.InteropServices;

namespace MathXK.Roots
{

    /// <summary>
    /// The results of the root finding operations
    /// </summary>
    public class RootResults
    {
        private bool _isSuccess;
        private double _xSolutionOrLast;
        private RootBracketResult _bracket;
        private int _iterations;


        /// <summary>
        /// Returns true if the root was located within the specified tolerances
        /// </summary>
        public bool Success { get { return _isSuccess; } }

        /// <summary>
        /// Returns the root or NaN if the root was not found
        /// </summary>
        /// <value>Requires that the solution value is not NaN. If there is no solution set LastX to the last value</value>
        public double SolutionX
        {
            get { return _isSuccess ? _xSolutionOrLast : double.NaN; }
            set
            {
                _isSuccess = !double.IsNaN(value);
                _xSolutionOrLast = value;
            }
        }

        /// <summary>
        /// Returns the last (x-value) attempt at finding the root
        /// </summary>
        /// <value>The value of the last attempt</value>
        public double LastX
        {
            get { return _xSolutionOrLast; }
            set { _isSuccess = false; _xSolutionOrLast = value; }
        }

        /// <summary>
        /// Returns the bracket, the last attempt at finding a bracket, or null if not supported
        /// </summary>
        /// <value>The bracket</value>
        public RootBracketResult Bracket
        {
            get { return _bracket; }
            set { _bracket = value; }
        }


        /// <summary>
        /// The number of iterations before the routine was terminated
        /// </summary>
        /// <value>Requires number of iterations &gt;=0</value>
        public int Iterations
        {
            get { return _iterations; }
            set
            {
                if (value < 0)
                    throw new Exceptions.DomainException("Requires Iterations >= 0: Iterations = {0}", value);
                _iterations = value;
            }
        }


    };


    /// <summary>
    /// Contains all the root finders
    /// </summary>
    public static partial class RootFinder
    {

        /// <summary>
        /// Gets the default tolerance for the root finder in ulps.
        /// </summary>
        public const int DefaultUlpsTolerance = 10;

        /// <summary>
        /// Gets the default relative tolerance for the root finder.
        /// </summary>
        public const double DefaultRelTolerance = DefaultUlpsTolerance * DoubleLimits.MachineEpsilon;

        /// <summary>
        /// Gets the default absolute tolerance for the root finder.
        /// </summary>
        public const double DefaultAbsTolerance = 0;

        /// <summary>
        /// Gets the default maximum number of iterations for the root finder.
        /// </summary>
        public const int DefaultIterations = Policies.MaxRootIterations;


        /// <summary>
        /// Returns true if the interval [fa, fb] contains the root(zero)
        /// </summary>
        /// <param name="fa">f(a)</param>
        /// <param name="fb">f(b)</param>
        /// <returns></returns>
        private static bool IsBracket(double fa, double fb)
        {
            return ((fa >= 0 && fb <= 0) || (fa <= 0 && fb >= 0));
        }


    };

}