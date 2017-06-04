//  (C) Copyright Rocco Romeo 2013.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

using System;
using System.Diagnostics;
using MathXK;

namespace MathXK
{

    public static partial class Math2
    {

        /// <summary>
        /// Returns true if | x - y | ≤ tolerance.
        /// </summary>
        /// <param name="x">Argument 1</param>
        /// <param name="y">Argument 2</param>
        /// <param name="tolerance">The absolute tolerance. Requires tolerance ≥ 0</param>
        /// <returns></returns>
        public static bool AreNearAbs(double x, double y, double tolerance)
        {
            if (!(tolerance >= 0))
                throw new Exceptions.DomainException("Requires Relative Error(tolerance) >= 0: tolerance = " + tolerance);

            return (Math.Abs(x - y) <= tolerance);
        }



        /// <summary>
        /// Returns true if |x-y| ≤ tolerance*max(|x|,|y|).
        /// </summary>
        /// <param name="x">Argument 1</param>
        /// <param name="y">Argument 2</param>
        /// <param name="tolerance">The relative tolerance. Requires tolerance ≥ 0</param>
        /// <returns></returns>
        public static bool AreNearRel(double x, double y, double tolerance)
        {
            if (!(tolerance >= 0))
                throw new Exceptions.DomainException("Requires Relative Error(tolerance) >= 0: tolerance = " + tolerance);

            double absX = Math.Abs(x);
            double absY = Math.Abs(y);
            double diff = Math.Abs(x - y);

            return (diff <= tolerance * Math.Max(absX, absY));
        }





        /// <summary>
        /// Using integer representations, returns true if |x-y| ≤ maxUlps 
        /// </summary>
        /// <param name="x">Argument 1</param>
        /// <param name="y">Argument 2</param>
        /// <param name="maxUlps">Maximum number of representable floating point numbers by which x and y may differ</param>
        /// <returns></returns>
        public static bool AreNearUlps(double x, double y, int maxUlps)
        {
            // Uses Bruce Dawson's technique for comparing IEEE doubles
            // See http://www.cygnus-software.com/papers/comparingfloats/comparingfloats.htm

            // There are several optional checks that you can do, depending
            // on what behavior you want from your floating point comparisons.
            // These checks should not be necessary and they are included
            // mainly for completeness.

            // If A or B are a NAN, return false. NANs are equal to nothing,
            // not even themselves.
            // This check is only needed if you will be generating NANs
            // and you use a maxUlps greater than 4 million or you want to
            // ensure that a NAN does not equal itself.

            if (maxUlps < 0)
                throw new Exceptions.DomainException("Requires Relative Error(maxUlps) >= 0: maxUlps = " + maxUlps);

            if (double.IsNaN(x) || double.IsNaN(y))
                return false;

            // if identical or if both are equal to infinity
            if (x == y)
                return true;

            // this check needs to be here so that inf is not close to double.MaxValue
            if (double.IsInfinity(x) || double.IsInfinity(y))
                return false;



#if false
            // After adjusting floats so their representations are lexicographically
            // ordered as twos-complement integers a very small positive number
            // will compare as 'close' to a very small negative number. If this is
            // not desireable, and if you are on a platform that supports
            // subnormals (which is the only place the problem can show up) then
            // you need this check.
            // The check for x == y is because zero and negative zero have different
            // signs but are equal to each other.
            if (Math.Sign(x) != Math.Sign(y))
                return A == B;
#endif

            // Make aInt lexicographically ordered as a twos-complement int
            // int.MinValue = 0x80000000
            Int64 aInt = BitConverter.DoubleToInt64Bits(x);
            if (aInt < 0)
                aInt = ~aInt + 1; // Int64.MinValue - aInt;   


            // Make bInt lexicographically ordered as a twos-complement int
            Int64 bInt = BitConverter.DoubleToInt64Bits(y);
            if (bInt < 0)
                bInt = ~bInt + 1; // Int64.MinValue - bInt;   


            // Now we can compare aInt and bInt to find out how far apart A and B are.
            return (Math.Abs(aInt - bInt) <= (Int64)maxUlps);
        }

    }





} // namespace

