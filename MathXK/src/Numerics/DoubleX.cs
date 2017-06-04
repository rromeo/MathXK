//  (C) Copyright Rocco Romeo 2013.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)


/*
 * Exponent reduction code from fdlibm
 * http://www.netlib.org/fdlibm/e_exp.c
 * ====================================================
 * Copyright (C) 2004 by Sun Microsystems, Inc. All rights reserved.
 *
 * Permission to use, copy, modify, and distribute this
 * software is freely granted, provided that this notice 
 * is preserved.
 * ====================================================
 */

using System;
using System.Diagnostics;

namespace MathXK.Numerics
{

    ///<summary>
    /// A double with 16 bit exponent.
    /// Used to facilitate double precision calculations where overflow/underflow is likely.
    ///</summary>
    internal struct DoubleX : IEquatable<DoubleX>
    {

        //Data

        ///<summary>mantissa</summary>
        private readonly double _m;

        ///<summary>exponent</summary>
        private readonly int _e;


        // Constants for Log(2) and 1/Log(2)

        private const double Ln2 = 6.93147180559945286227e-01; /* 0x3FE62E42, 0xFEFA39EF */
        private const double Ln2_Hi = 6.93147182464599609375e-01; /* 0x3FE62E43, 0x00000000 */
        private const double Ln2_Lo = -1.90465429995776804525e-09; /* 0xBE205C61, 0x0CA86C39 */

        private const double RLn2 = 1.44269504088896338700e+00; /* 0x3FF71547, 0x652B82FE =1/ln2 */
        private const double RLn2_Hi = 1.44269502162933349609e+00; /* 0x3FF71547, 0x60000000 =24b 1/ln2*/
        private const double RLn2_Lo = 1.92596299112661746887e-08; /* 0x3E54AE0B, 0xF85DDF44 =1/ln2 tail*/


        private DoubleX(double m, int e)
        {
            Debug.Assert(e >= MinExponent && e <= MaxExponent);

            _m = m;
            _e = e;
        }



        public static readonly DoubleX Zero = new DoubleX(0, 0);
        public static readonly DoubleX One = new DoubleX(1);
        public static readonly DoubleX Two = new DoubleX(2);
        public static readonly DoubleX NaN = new DoubleX(double.NaN, 0);
        public static readonly DoubleX PositiveInfinity = new DoubleX(double.PositiveInfinity, 0);
        public static readonly DoubleX NegativeInfinity = new DoubleX(double.NegativeInfinity, 0);

        // MaxExponent and MinExponent are defined similarly to the C++ std::numeric_limits
        // 2^exp = 0.5^(exp+1), so they need to be offset by one

        /// <summary>
        /// The maximum exponent such that 2^(MaxExponent-1) does not overflow. 
        /// MaxExponent = 32767
        /// </summary>
        public const int MaxExponent = short.MaxValue;

        /// <summary>
        /// The minimum exponent such that 2^(MinExponent-1) does not underflow.
        /// MinExponent = -32768
        /// </summary>
        public const int MinExponent = short.MinValue;


        /// <summary>
        /// Maximum value of x before e^x overflows
        /// <para>Value = 22712</para>
        /// </summary>
        public const double MaxLogValue = 22712;

        /// <summary>
        /// Minimum value of x before e^x underflows
        /// <para>Value = -22712</para>
        /// </summary>
        public const double MinLogValue = -MaxLogValue;


        /// <summary>
        /// Maximum value of x such that x^x will not overflow
        /// <para>Value = 2854</para>
        /// </summary>
        public const double MaxPowValue = 2854;

        /// <summary>
        /// Static constructor used for debug
        /// </summary>
        static DoubleX()
        {
            Debug.Assert(MaxLogValue == Math.Floor(MaxExponent * Ln2), "Wrong Value for MaxLogValue");

        }


        public DoubleX(double x)
        {
            _m = Math2.Frexp(x, out _e);
        }

        /// <summary>
        /// Gets the mantissa
        /// </summary>
        public double Mantissa { get { return _m; } }

        /// <summary>
        /// Gets the binary exponent
        /// </summary>
        public int Exponent { get { return _e; } }

        /// <summary>
        /// Returns true if x is NaN
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        public static bool IsNaN(DoubleX x) { return double.IsNaN(x._m); }

        /// <summary>
        /// Returns true if |x| is Infinity
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        public static bool IsInfinity(DoubleX x) { return double.IsInfinity(x._m); }


        /// <summary>
        /// Unary minus
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        public static DoubleX operator -(DoubleX x)
        {
            return new DoubleX(-x._m, x._e);
        }

        /// <summary>
        /// Returns x*y
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <returns></returns>
        public static DoubleX operator *(DoubleX x, DoubleX y)
        {
            return Ldexp(x._m * y._m, x._e + y._e);
        }

        /// <summary>
        /// Returns x/y
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <returns></returns>
        public static DoubleX operator /(DoubleX x, DoubleX y)
        {
            return Ldexp(x._m / y._m, x._e - y._e);
        }

        /// <summary>
        /// Returns x*y
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <returns></returns>
        public static DoubleX operator *(DoubleX x, double y)
        {
            return x * (new DoubleX(y));
        }

        /// <summary>
        /// Return x*y
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <returns></returns>
        public static DoubleX operator *(double x, DoubleX y)
        {
            return (new DoubleX(x)) * y;
        }

        /// <summary>
        /// Returns x/y
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <returns></returns>
        public static DoubleX operator /(DoubleX x, double y)
        {
            return x / (new DoubleX(y));
        }


        /// <summary>
        /// Returns x/y
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <returns></returns>
        public static DoubleX operator /(double x, DoubleX y)
        {
            return (new DoubleX(x)) / y;
        }


        /// <summary>
        /// Returns x+y
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <returns></returns>
        public static DoubleX operator +(DoubleX x, DoubleX y)
        {
            double a = x._m + y._m;
            if (double.IsNaN(a) || double.IsInfinity(a) )
                return new DoubleX(a, 0);


            if ( x._m == 0 )
                return y;
            if ( y._m == 0 )
                return x;

            // if xm > (1+eps)*ym 
            int le = x._e - y._e;
            if ( Math.Abs(le) > DoubleLimits.MantissaBits )
                return (le > 0) ? x : y;

            double m;
            int e;

            int ediff = le;
            if ( ediff == 0 ) {
                e = x._e;
                m = x._m+ y._m;
            } else if ( ediff > 0 ) {
                // |x| > |y|
                e = x._e;
                m = x._m + Math2.Ldexp(y._m, -ediff);
            } else {
                // |x| < |y|
                e = y._e;
                m = y._m + Math2.Ldexp(x._m, ediff);
            }

            // now renormalize
            return Ldexp(m, e);
        }

        /// <summary>
        /// Returns x-y
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <returns></returns>
        public static DoubleX operator -(DoubleX x, DoubleX y)
        {
            return x + (-y);
        }

        /// <summary>
        /// returns true if x == y, otherwise false
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <returns></returns>
        public static bool operator ==(DoubleX x, DoubleX y)
        {
            return (x._m == y._m && x._e == y._e);
        }

        /// <summary>
        /// returns true if x != y, otherwise false
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <returns></returns>
        public static bool operator !=(DoubleX x, DoubleX y)
        {
            return !(x == y);
        }

        /// <summary>
        /// returns true if x &gt; y, otherwise false
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <returns></returns>
        public static bool operator >(DoubleX x, DoubleX y)
        {
            if (double.IsNaN(x._m) || double.IsNaN(y._m))
                return false;

            if (Math.Sign(x._m) <= Math.Sign(y._m))
                return false;

            return ((x._e > y._e) || ((x._e == y._e) && (x._m > y._m)));
        }

        /// <summary>
        /// returns true if x &gt;= y, otherwise false
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <returns></returns>
        public static bool operator >=(DoubleX x, DoubleX y)
        {
            if (double.IsNaN(x._m) || double.IsNaN(y._m))
                return false;

            if (Math.Sign(x._m) < Math.Sign(y._m))
                return false;

            return ((x._e > y._e) || ((x._e == y._e) && (x._m >= y._m)));
        }

        /// <summary>
        /// returns true if x &lt; y, otherwise false
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <returns></returns>
        public static bool operator <(DoubleX x, DoubleX y)
        {
            if (double.IsNaN(x._m) || double.IsNaN(y._m))
                return false;

            if (Math.Sign(x._m) >= Math.Sign(y._m))
                return false;

            return ((x._e < y._e) || ((x._e == y._e) && (x._m < y._m)));
        }

        /// <summary>
        /// returns true if x &lt;= y, otherwise false
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <returns></returns>
        public static bool operator <=(DoubleX x, DoubleX y)
        {
            if (double.IsNaN(x._m) || double.IsNaN(y._m))
                return false;

            if (Math.Sign(x._m) > Math.Sign(y._m))
                return false;

            return ((x._e < y._e) || ((x._e == y._e) && (x._m <= y._m)));
        }

        /// <summary>
        /// Returns 2^n for integer n
        /// </summary>
        /// <param name="n"></param>
        /// <returns></returns>
        public static DoubleX Exp2(int n)
        {
            // Since normalized 2^n = 0.5 * 2^(n+1)
            // we will only overflow if n > MaxBinaryExponent-1 
            // ( or n == MaxBinaryExponent )

            if (n >= MaxExponent)
                return PositiveInfinity;
            if (n < MinExponent-1)
                return Zero;

            return new DoubleX(0.5, n+1);
        }

        /// <summary>
        /// Returns 2^x
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        public static DoubleX Exp2(double x)
        {

            if (double.IsNaN(x)) {
                Policies.ReportDomainError("Exp2(x: {0}): x must have a value", x);
                return NaN;
            }

            if (x >= MaxExponent)
                return PositiveInfinity;
            if (x < MinExponent)
                return Zero;


            // 2^0 = 1
            if (x == 0)
                return One;
            if (x == 1)
                return Two;

            // truncate so there's no chance of overflow in the exponent
            int exp_i = (int)x;
            double exp_f = x - exp_i;

            Debug.Assert(Math.Abs(exp_f) < 1);

            if (exp_f == 0)
                return new DoubleX(0.5, exp_i + 1);

            return Ldexp(Math.Pow(2, exp_f), exp_i);
            //return new DoubleX(0.5, exp_i + 1) * Math.Pow(2, exp_f);
        }

        /// <summary>
        /// Returns e^x
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        public static DoubleX Exp(double x)
        {
            if (double.IsNaN(x)) {
                Policies.ReportDomainError("Exp(x: {0}): x must have a value", x);
                return NaN;
            }


            // definite overflow/underflow
            // e^(MaxLogValue) > 2^(MaxBinaryExponent)
            // e^(MinLogValue) < 2^(MinBinaryExponent)
            if (x > MaxLogValue)
                return PositiveInfinity;
            if (x < MinLogValue)
                return Zero;

            bool invert = false;
            double absX = x;
            if (x < 0) {
                absX = -x;
                invert = true;
            }

            // Reduce x to an r so that |r| <= 0.5*ln2 ~ 0.34658.
            // Given x, find r and integer k such that x = k*ln2 + r,  |r| <= 0.5*ln2.  
            double t = Math.Floor(RLn2 * absX + 0.5);
            double r = (absX - t * Ln2_Hi) - t * Ln2_Lo;

            if ( invert ) {
                t = -t;
                r = -r;
            }

            return Ldexp(Math.Exp(r), (int)t);
        }

        /// <summary>
        /// Return x^y
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <returns></returns>
        public static DoubleX Pow(double x, double y)
        {
            if (double.IsNaN(x) || double.IsNaN(y)) {
                Policies.ReportDomainError("Pow(x: {0}, y: {1}): x,y must have a value", x, y);
                return NaN;
            }

            // x^0 = 1, including NaN
            // But, sometimes NaN can communicate an error condition, so assume the latter
            if (y == 0)
                return One;

            if (y < 0)
                return 1 / Pow(x, -y);

            // (0)^y = 0, 1^y = 1, x^1 = 1
            if (x == 0 || x == 1 || y == 1)
                return new DoubleX(x);

            if (x < 0) {
                if (!Math2.IsInteger(y)) {
                    Policies.ReportDomainError("Pow(x: {0}, y: {1}): with x<0, y must be an integer", x, y);
                    return NaN;
                }
                return (Math2.IsOdd(y)) ? -Pow(-x, y) : Pow(-x, y);
            }

            // when 0 < y < 1, there's no chance of overflow/underflow
            if (y < 1)
                return new DoubleX(Math.Pow(x, y));

            if (x == 2)
                return Exp2(y);

            double ylx = y * Math.Log(x);
            if (ylx >= DoubleLimits.MinLogValue && ylx <= DoubleLimits.MaxLogValue)
                return new DoubleX(Math.Pow(x, y));

            // TODO: Take a less granular approach here
            if (ylx >= 2 * DoubleLimits.MinLogValue && ylx <= 2 * DoubleLimits.MaxLogValue) {
                DoubleX result = new DoubleX(Math.Pow(x, y / 2));
                return result * result;
            }
            if (ylx >= 4 * DoubleLimits.MinLogValue && ylx <= 4 * DoubleLimits.MaxLogValue) {
                DoubleX result = new DoubleX(Math.Pow(x, y / 4));
                DoubleX result2 = result * result;
                return result2 * result2;
            }
           
            return Exp(ylx);
        }


        /// <summary>
        /// Convert to double
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        public static explicit operator double(DoubleX x)
        {
            return Math2.Ldexp(x._m, x._e);
        }


        /// <summary>
        /// Convert from double
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        public static implicit operator DoubleX(double x)
        {
            return new DoubleX(x);
        }

        /// <summary>
        /// Returns x*2^exp
        /// </summary>
        /// <param name="x"></param>
        /// <param name="exp"></param>
        /// <returns></returns>
        public static DoubleX Ldexp(double x, int exp)
        {
            if (x == 0 || double.IsNaN(x) || double.IsInfinity(x))
                return new DoubleX(x, 0);

            // normalize the result
            int xExp;
            double xMant = Math2.Frexp(x, out xExp);

            // watch for overflows/underflows on the exponent
            Int64 esum = exp;
            esum += xExp;

            if (esum < MinExponent)
                return Zero; // underflow
            if (esum > MaxExponent)
                return (x >= 0) ? PositiveInfinity : NegativeInfinity;

            return new DoubleX(xMant, (int)esum);

        }

        /// <summary>
        /// Returns x*2^exp
        /// </summary>
        /// <param name="x"></param>
        /// <param name="exp"></param>
        /// <returns></returns>
        public static DoubleX Ldexp(DoubleX x, int exp)
        {
            if (x._m == 0 || double.IsNaN(x._m) || double.IsInfinity(x._m))
                return x;

            Int64 nexp = x._e; 
            nexp += exp;

            if (nexp > MaxExponent) 
                return (x._m > 0) ? PositiveInfinity : NegativeInfinity;
            if (nexp < MinExponent)
                return Zero;


            return new DoubleX(x._m, (int)nexp);
        }

        public bool Equals(DoubleX other)
        {
            // this == x, except that NaNs are equal
            return (_e == other._e && _m.Equals(other));
        }

        public override bool Equals(object obj)
        {
            if ( obj == null )
                return false;

            if (obj is DoubleX o)
                return this.Equals(o);

            return false;
        }

        public override int GetHashCode()
        {
            return _m.GetHashCode() ^ _e.GetHashCode();
        }


    };


}