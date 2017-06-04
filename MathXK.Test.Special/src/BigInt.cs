//  (C) Copyright Rocco Romeo 2013.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

using System.Numerics;

namespace MathXK.Test
{
    public class BigInt
    {

        /// <summary>
        /// Compute n!! using BigInteger multiplication
        /// </summary>
        /// <param name="n"></param>
        /// <returns></returns>
        public static double Factorial2(int n)
        {
            if (n < -1)
                return double.NaN;

            if (n == 0 || n == -1)
                return 1.0;

            BigInteger b = 1;
            for (int i = n; i > 1; i -= 2)
                b *= i;
            return (double)b;

        }

        /// <summary>
        /// Compute n! using BigInteger multiplication
        /// </summary>
        /// <param name="n"></param>
        /// <returns></returns>
        public static double Factorial(int n)
        {
            if (n < 0)
                return double.NaN;

            if (n == 0)
                return 1.0;

            BigInteger b = 1;
            for (int i = 2; i <= n; i++)
                b *= i;
            return (double)b;

        }


        /// <summary>
        /// Returns the rising factorial: 
        /// <para> if n > 0, x*(x+1)*(x+2)...*(x+(n-1))</para>
        /// <para> if n == 0, 1</para>
        /// <para> if n &lt; 0, 1/((x-1)*(x-2)*...*(x-n))</para>
        /// </summary>
        public static double FactorialRising(BigInteger x, int n)
        {

            if (n == 0)
                return 1;

            // x*(x+1)*(x+2)...*(x+(n-1))
            if (n > 0) {
                BigInteger result1 = x;
                for (int i = 1; i < n; i++)
                    result1 *= (x + i);
                return (double)(result1);
            }

            BigInteger result = 1;
            for (int i = n; i < 0; i++)
                result *= (x + i);

            return 1.0 / (double)result;
        }

        /// <summary>
        /// Returns the falling factorial = x*(x-1)*(x-2)...*(x-(n-1))
        /// </summary>
        /// <param name="x"></param>
        /// <param name="n">Requires n ≥ 0</param>
        public static double FactorialFalling(BigInteger x, int n)
        {

            if (n == 0)
                return 1;
            if (n < 0)
                return double.NaN;

            BigInteger result = x;
            for (int i = 1; i < n; i++)
                result *= (x - i);

            return (double)result;
        }


    }
}
