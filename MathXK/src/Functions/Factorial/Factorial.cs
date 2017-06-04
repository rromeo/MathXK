//  Copyright (c) 2013 Rocco Romeo
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

//  Ported to .NET based on the original:
//      Copyright (C) John Maddock 2006, 2010, Boost Software License, Version 1.0

using System;
using System.Diagnostics;

namespace MathXK
{

    public static partial class Math2
    {

        /// <summary>
        /// Returns n!
        /// </summary>
        /// <param name="n">Requires n ≥ 0</param>
        public static double Factorial(int n)
        {

            if (n < 0) {
                Policies.ReportDomainError("Factorial(n: {0}): Requires n >= 0", n);
                return double.NaN;
            }
            if (n >= FactorialTable.Length) 
                return double.PositiveInfinity;
            return FactorialTable[n];
        }


        /// <summary>
        /// Returns the double factorial = n!!
        /// <para>If n is odd, returns: n * (n-2) * ... * 5 * 3 * 1</para>
        /// <para>If n is even, returns: n * (n-2) * ... * 6 * 4 * 2</para>
        /// <para>If n = 0, -1, returns 1</para>
        /// </summary>
        /// <param name="n">Requires n ≥ -1</param>
        public static double Factorial2(int n)
        {
            if (n < -1) {
                Policies.ReportDomainError("Factorial2(n: {0}): Requires n >= -1", n);
                return double.NaN;
            }

            // common values
            if (n == -1 || n == 0)
                return 1;

            double result;
            if (IsOdd(n)) {
                // odd n:
                if (n < FactorialTable.Length) {
                    int k = (n - 1) / 2;
                    return Math.Ceiling(Math2.Ldexp(FactorialTable[n] / FactorialTable[k], -k) - 0.5);
                }
                //
                // Fallthrough: n is too large to use table lookup, try the
                // gamma function instead.
                //
                result = Constants.RecipSqrtPI * Math2.Tgamma(0.5 * n + 1.0);
                result = Math.Ceiling(Math2.Ldexp(result, (n + 1) / 2) - 0.5);
            } else {
                // even n:
                int k = n / 2;
                result = Math2.Ldexp(Factorial(k), k);
            }

            return result;
        }

        /// <summary>
        /// Returns the rising factorial:
        /// <para>(x)_{n} = Γ(x+n)/Γ(x)</para>
        /// </summary>
        /// <returns>
        /// <para> if n > 0, x*(x+1)*(x+2)...*(x+(n-1))</para>
        /// <para> if n == 0, 1</para>
        /// <para> if n &lt; 0, 1/((x-1)*(x-2)*...*(x-n))</para>
        /// </returns>
        public static double FactorialRising(double x, int n)
        {
            // This section uses the following notation: 
            // falling factorial: (x)_{_n}  
            // rising factorial: (x)_{n}  

            // Standard eqn for the rising factorial:
            // (x)_{n} = Γ(x+n)/Γ(x)

            double result = double.NaN;

            // (x)_{0} = 1
            if (n == 0)
                return 1;

            if (double.IsNaN(x)) {
                Policies.ReportDomainError("FactorialRising(x: {0}, n: {1}): Requires x not NaN", x, n);
                return double.NaN;
            }


            if (n == int.MinValue) {
                // avoid -int.MinValue error
                // this will overflow/underflow, but here it is
                if (x > 0 && IsInteger(x) && (x + n) <= 0) {
                    Policies.ReportPoleError("FactorialRising(x: {0}, n: {1}): Divide by Zero", x, n);
                    return double.NaN;
                }

                // (x)_{n} = (x)_{n+1}/(x+n)
                return FactorialRising(x, n + 1) / (x + n);
            }


            if (x == 0) {

                // (0)_{n} = 0, n > 0
                if (n > 0)
                    return 0;

                // (0)_{-n} = (-1)^n / n!
                result = Math2.Factorial(-n);
                if (IsOdd(n))
                    result = -result;
                return 1 / result;
            }

            if (x < 0) {

                // For x less than zero, we really have a falling
                // factorial, modulo a possible change of sign.
                // Note: the falling factorial isn't defined for negative n
                // so cover that first

                if (n < -x) {
                    // handle (n < 0 || ( n > 0 && n < -x)) directly
                    // (x)_{n} = Γ(1-x)/Γ(1-x-n)
                    result = Math2.TgammaDeltaRatio(1 - x, -n);
                } else {
                    Debug.Assert(n > 0);
                    result = FactorialFalling(-x, n);
                }

                if (IsOdd(n))
                    result = -result;

                return result;
            }


            // x > 0

            // standard case: (n > 0 || (n < 0 && x > -n))
            // (x)_{n} = Γ(x+n)/Γ(x)
            if (-x < n)
                return 1 / Math2.TgammaDeltaRatio(x, n);

            Debug.Assert((x + n) <= 0);

            // (x)_{n} = (-1)^n *  Γ(1-x)/Γ(1-x-n)
            if (x < 1) {
                result = TgammaDeltaRatio(1 - x, -n);
                if (IsOdd(n))
                    result = -result;
                return result;
            }

            // x >= 1
            result = FactorialFalling(x - 1, -n);
            if (result == 0) {
                if (IsInteger(x)) {
                    Policies.ReportPoleError("FactorialRising(x: {0}, n: {1}): Divide by Zero", x, n);
                    return double.NaN;
                }

                return double.PositiveInfinity;
            }
            return 1 / result;
        }

        /// <summary>
        /// Returns the falling factorial:
        /// <para>(x)_{_n} = x*(x-1)*(x-2)...*(x-(n-1)) = Γ(x+1)/Γ(x-n+1)</para>
        /// </summary>
        /// <param name="x"></param>
        /// <param name="n">Requires n ≥ 0</param>
        /// <returns>
        /// <para>if n > 0, x*(x-1)*(x-2)...*(x-(n-1))</para>
        /// <para>if n == 0, 1</para>
        /// </returns>
        public static double FactorialFalling(double x, int n)
        {
            // This section uses the following notation: 
            // falling factorial: (x)_{_n}  
            // rising factorial: (x)_{n}  

            // Standard eqn for the falling factorial:
            // (x)_{_n} = Γ(x+1)/Γ(x+1-n)

            double result = double.NaN;

            // (x)_{_0} = 1
            if (n == 0)
                return 1;

            if ((double.IsNaN(x))
            || (n < 0)) {
                Policies.ReportDomainError("FactorialFalling(x: {0}, n: {1}): Requires finite x; n >= 0", x, n);
                return double.NaN;
            }

            if (n == 1)
                return x;

            if (x == 0)
                return 0;
            if (x < 0) {
                //
                // For x < 0 we really have a rising factorial
                // modulo a possible change of sign:
                // (x)_{_n} = (-1)^n * (-x)_{n}
                result = FactorialRising(-x, n);
                return IsOdd(n) ? -result : result;
            }

            // standard case: Γ(x+1)/Γ(x+1-n)
            if (x > n - 1) {
                // (n)_{_n} = n!
                if (x == n)
                    return Math2.Factorial(n);

                // (n)_{_n+1} = n!, n>1
                if (x == n + 1) {
                    Debug.Assert(n > 1);
                    return Math2.Factorial(n + 1);
                }

                return Math2.TgammaDeltaRatio(x + 1, -n);
            }


            Debug.Assert(x <= n - 1);


            if (x < 1) {
                int m = n - 1;

                double MaxFactorial = Math2.MaxFactorialIndex;
                if (m < MaxFactorial) {
                    result = x / TgammaDeltaRatio(1 - x, m);
                } else {
                    // try not to overflow/underflow too soon

                    double m2 = m - MaxFactorial;
                    result = (x / TgammaDeltaRatio(1 - x, MaxFactorial)) / TgammaDeltaRatio(1 - x + MaxFactorial, m2);
                }


                return IsOdd(m) ? -result : result;
            }

            //
            // x+1-n will be negative and TgammaDeltaRatio won't
            // handle it, split the product up into three parts:
            //
            double xp1 = x + 1;
            int n2 = (int)xp1;
            if (n2 == xp1)
                return 0;
            result = Math2.TgammaDeltaRatio(xp1, -n2);
            x -= n2;
            result *= x;
            ++n2;
            if (n2 < n)
                result *= FactorialFalling(x - 1, n - n2);
            return result;
        }

    }

} // namespaces



