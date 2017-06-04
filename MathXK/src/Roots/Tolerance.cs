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
        /// <summary>
        /// A function (delegate) that returns true if x is near y; otherwise false
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <returns></returns>
        public delegate bool ToleranceFunc(double x, double y);

        /// <summary>
        /// A convenience function that returns a ToleranceFunc for use in the Root Finder routines.
        /// <para>If relTolerance > 0, returns (x, y) => |x-y| &lt;= relTolerance * Min(|x|,|y|)</para> 
        /// <para>If absTolerance > 0, returns (x, y) => |x-y| &lt;= absTolerance</para> 
        /// <para>If both are greater than 0, returns a combination</para>
        /// <para>If both are 0, returns (x, y) => (x == y) </para>
        /// </summary>
        /// <param name="relTolerance">The relative tolerance. Requires relTolerance >= 0</param>
        /// <param name="absTolerance">The absolute tolerance. Requires absTolerance >= 0</param>
        /// <returns></returns>
        public static ToleranceFunc GetToleranceFunc(double relTolerance = DefaultRelTolerance, double absTolerance = DefaultAbsTolerance)
        {
            if (!(relTolerance >= 0) || double.IsInfinity(relTolerance))
                throw new Exceptions.DomainException("Requires finite relative tolerance >= 0: relTolerance = {0}", relTolerance);
            if (!(absTolerance >= 0) || double.IsInfinity(absTolerance))
                throw new Exceptions.DomainException("Requires finite absolute tolerance >= 0: absTolerance = {0}", absTolerance);

            bool AreNearRelAndAbs(double x, double y) {
                if (x == y)
                    return true;

                double rdiff = relTolerance * Math.Min(Math.Abs(x), Math.Abs(y));
                double diff = Math.Abs(x - y);
                return (diff <= rdiff || diff <= absTolerance);
            };

            bool AreNearRel(double x, double y) {
                if (x == y)
                    return true;

                return (Math.Abs(x - y) <= relTolerance * Math.Min(Math.Abs(x), Math.Abs(y)));
            };

            bool AreNearAbs(double x, double y) {
                if (x == y)
                    return true;

                return (Math.Abs(x - y) <= absTolerance);
            };

            bool AreEqual(double x, double y) {
                return (x == y);
            };

            // given the relative and absolute tolerances, return the appropriate
            // comparison function

            ToleranceFunc f;

            if (relTolerance > 0) {
                if (absTolerance > 0)
                    f = AreNearRelAndAbs;
                else
                    f = AreNearRel;
            } else {
                if (absTolerance > 0) 
                    f = AreNearAbs;
                else 
                    f = AreEqual;
            }

            return f;
        }

        /// <summary>
        /// A convenience function that returns a ToleranceFunc for use in the Root Finder routines.
        /// <para>If ulpsTolerance > 0, returns (x, y) => AreNearUlps(x, y, absTolerance)</para> 
        /// <para>Otherwise, returns (x, y) => (x == y)</para>
        /// </summary>
        /// <param name="ulpsTolerance"></param>
        /// <returns></returns>
        public static ToleranceFunc AreNearUlps(int ulpsTolerance = DefaultUlpsTolerance)
        {
            if (ulpsTolerance < 0)
                throw new Exceptions.DomainException("Requires ulpsTolerance >= 0: ulpsTolerance = {0}", ulpsTolerance);

            bool AreNearUlps(double x, double y) {
                return Math2.AreNearUlps(x, y, ulpsTolerance);
            };

            bool AreEqual(double x, double y) {
                return (x == y);
            };

            ToleranceFunc f;
            if (ulpsTolerance > 0)
                f = AreNearUlps;
            else 
                f = AreEqual;

            return f;
        }

    }




}