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
                throw new AssertFailedException("Failed: actual = " + actual + "; expected = " + expected + "; ulps = " + ulps + "; MaxUlps = " + maxUlps);
        }

        public static void AreNear(double expected, double actual, int maxUlps, Func<string> onFalureString)
        {
            if (!Math2.AreNearUlps(expected, actual, maxUlps))
                throw new AssertFailedException(onFalureString());
        }


#if false
        // Not supporting complex yet, but left here for future reference

        public static void AreNear(Complex expected, Complex actual, int maxUlps) {

            bool nearRe = Math2.AreNearUlps(expected.Real, actual.Real, maxUlps);
            bool nearIm = Math2.AreNearUlps(expected.Imaginary, actual.Imaginary, maxUlps);

            if ( nearRe && nearIm )
                return;

            if ( !nearRe && nearIm )
                throw new AssertFailedException("Failure: actual.Real not close to expected.Real: actual = " + actual + "; expected = " + expected + "; ulps = " + maxUlps);

            if ( nearRe && !nearIm )
                throw new AssertFailedException("Failure: actual.Imaginary not close to expected.Imaginary: actual = " + actual + "; expected = " + expected + "; ulps = " + maxUlps);

            throw new AssertFailedException("Failure: actual not close to expected: actual = " + actual + "; expected = " + expected + "; ulps = " + maxUlps);
        }

        public static void AreNear(Complex expected, Complex actual, int maxUlps, Func<string> onFalureString) {

            bool nearRe = Math2.AreNearUlps(expected.Real, actual.Real, maxUlps);
            bool nearIm = Math2.AreNearUlps(expected.Imaginary, actual.Imaginary, maxUlps);

            if (nearRe && nearIm)
                return;

            throw new AssertFailedException(onFalureString());
        }
#endif


        // See:
        // http://stackoverflow.com/questions/933613/c-how-do-i-use-assert-unit-testing-to-verify-that-an-exception-has-been-thro

        public static void Throws<T>(Action func) where T: Exception
        {
            var exceptionThrown = false;
            try {
                func.Invoke();
            } catch (T) {
                exceptionThrown = true;
            }

            if (!exceptionThrown) {
                throw new AssertFailedException(
                    String.Format("An exception of type {0} was expected, but not thrown", typeof(T))
                    );
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
                throw new AssertFailedException(
                    String.Format("Expected a result of {0}, got {1}", expected, actual)
                    );

        }

    }
}