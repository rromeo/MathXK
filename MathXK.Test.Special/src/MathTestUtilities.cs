//  (C) Copyright Rocco Romeo 2013.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

using Microsoft.VisualStudio.TestTools.UnitTesting;
using System;
using System.Numerics;


namespace MathXK.Test
{
    static class Assert2
    {
        public static void AreNear(double expected, double actual, int maxUlps)
        {
            double ulps = Math2.FloatDistance(expected, actual);
            if (!(expected == actual || Math.Abs(ulps) <= maxUlps)) 
                throw new AssertFailedException($"Failed: actual = {actual}; expected = {expected}; ulps = {ulps}; MaxUlps = {maxUlps}");
        }

        public static void AreNear(double expected, double actual, int maxUlps, Func<string> onFalureString)
        {
            if (!Math2.AreNearUlps(expected, actual, maxUlps))
                throw new AssertFailedException(onFalureString());
        }


#if false
        // Not supporting complex yet, but left here for future reference

        public static void AreNear(Complex expected, Complex actual, int maxUlps)
        {

            double ulpsRe = Math2.FloatDistance(expected.Real, actual.Real);
            double ulpsIm = Math2.FloatDistance(expected.Imaginary, expected.Imaginary);

            bool nearRe = (expected.Real == actual.Real || Math.Abs(ulpsRe) <= maxUlps);
            bool nearIm = (expected.Imaginary == actual.Imaginary || Math.Abs(ulpsIm) <= maxUlps);

            if ( nearRe && nearIm )
                return;

            throw new AssertFailedException($"Failed: actual = {actual}; expected = {expected}; ulps = ({ulpsRe}, {ulpsIm}); MaxUlps = {maxUlps}");
        }

        public static void AreNear(Complex expected, Complex actual, int maxUlps, Func<string> onFalureString) {

            bool nearRe = Math2.AreNearUlps(expected.Real, actual.Real, maxUlps);
            bool nearIm = Math2.AreNearUlps(expected.Imaginary, actual.Imaginary, maxUlps);

            if (nearRe && nearIm)
                return;

            throw new AssertFailedException(onFalureString());
        }
#endif

        public static void Throws<T>(Action func) where T: Exception
        {
            var exceptionThrown = false;
            try {
                func.Invoke();
            } catch (T) {
                exceptionThrown = true;
            }

            if (!exceptionThrown) {
                throw new AssertFailedException($"An exception of type {typeof(T)} was expected, but not thrown");
            }
        }


        /// <summary>
        /// Returns true if expected == actualFunc() or if the expected exception was thrown.
        /// Most commonly used to test for NaN or Domain Exceptions
        /// </summary>
        /// <typeparam name="T"></typeparam>
        /// <param name="expected"></param>
        /// <param name="actualFunc"></param>
        public static void AreEqual<T>(double expected, Func<double> actualFunc) where T: Exception
        {
            double actual;
            try {
                actual = actualFunc.Invoke();
            } catch (T) {
                return;
            }

            // watch for NaN
            if (!actual.Equals(expected))
                throw new AssertFailedException($"Expected a result of {expected}, got {actual}");

        }

    }
}