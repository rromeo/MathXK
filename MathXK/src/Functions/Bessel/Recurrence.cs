//  Copyright (c) 2013 Rocco Romeo
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

//  Ported to .NET based on the original:
//      Copyright (c) 2006 Xiaogang Zhang, Boost Software License v1.0
//  History:
//      Source originally by XZ. 
//      RR restructured code and made extensive changes to improve accuracy.

#if DEBUG
// uncomment to get more detailed debug messages
//#define EXTRA_DEBUG

#endif

using System;
using System.Diagnostics;
using System.Runtime.InteropServices;
using MathXK.Numerics;

namespace MathXK {

    static partial class _Bessel {


        static class Recurrence {


            /// <summary>
            /// Returns a tuple (J{v+n}(x), J{v+n-1}(x), BScale) OR (Y{v+n}(x), Y{v+n-1}(x), BScale) after n forward recurrences.
            /// Note: to get actual values multiply by 2^BScale.
            /// </summary>
            /// <param name="v">Order of the current parameter</param>
            /// <param name="x">Argument</param>
            /// <param name="n">Number of forward recurrence steps. 0 returns (current,prev)</param>
            /// <param name="current">J{v}(x) OR Y{v}(x)</param>
            /// <param name="prev">J{v-1}(x) OR Y{v-1}(x)</param>
            /// <returns></returns>
            public static (double JYvpn, double JYvpnm1, int BScale) ForwardJY_B(double v, double x, int n, double current, double prev)
            {
                Debug.Assert(n >= 0);
                Debug.Assert(!double.IsInfinity(2 * (v + n) / x));

                int binaryScale = 0;
                for (int k = 0; k < n; k++) {
                    double fact = 2 * (v + k) / x;
                    double value = fact * current - prev;

                    if (double.IsInfinity(value)) {
                        int exp = 0;
                        current = Math2.Frexp(current, out exp);
                        prev = Math2.Ldexp(prev, -exp);
                        binaryScale += exp;
                        value = fact * current - prev;
                    }
                    prev = current;
                    current = value;
                }

                {
                    int exp = 0;
                    current = Math2.Frexp(current, out exp);
                    prev = Math2.Ldexp(prev, -exp);
                    binaryScale += exp;
                }

                return (current, prev, binaryScale);
            }


            /// <summary>
            /// Returns (J{v+n}(x), J{v+n-1}(x)) OR (Y{v+n}(x), Y{v+n-1}(x)) after n forward recurrences
            /// </summary>
            /// <param name="v">Order of the current parameter</param>
            /// <param name="x">Argument</param>
            /// <param name="n">Number of forward recurrence steps. 0 returns (current,prev)</param>
            /// <param name="current">J{v}(x) OR Y{v}(x)</param>
            /// <param name="prev">J{v-1}(x) OR Y{v-1}(x)</param>
            /// <returns></returns>
            public static (double JYvpn, double JYvpnm1) ForwardJY(double v, double x, int n, double current, double prev)
            {
                int scale = 0;
                for (int k = 0; k < n; k++) {
                    double fact = 2 * (v + k) / x;
                    double value = fact * current - prev;

                    if (double.IsInfinity(value)) {
                        int exp = 0;
                        current = Math2.Frexp(current, out exp);
                        prev = Math2.Ldexp(prev, -exp);
                        scale += exp;
                        value = fact * current - prev;
                    }
                    prev = current;
                    current = value;
                }

                if (scale != 0) {
                    current = Math2.Ldexp(current, scale);
                    prev = Math2.Ldexp(prev, scale);
                }

                return (current, prev);
            }


            /// <summary>
            /// Returns (J{v-n}(x), J{v-n+1}(x), scale) OR (Y{v-n}(x), Y{v-n+1}(x), scale) after n backward recurrences
            /// <para>Note: Multiply by 2^BScale for the true value</para>
            /// </summary>
            /// <param name="v">Order of the current parameter</param>
            /// <param name="x">Argument</param>
            /// <param name="n">Number of backward recurrence steps. 0 returns (current,prev)</param>
            /// <param name="current">J{v}(x) OR Y{v}(x)</param>
            /// <param name="prev">J{v+1}(x) OR Y{v+1}(x)</param>
            /// <returns></returns>
            public static (double JYvmn, double JYvmnp1, int BScale) BackwardJY_B(double v, double x, int n, double current, double prev)
            {
                Debug.Assert(n >= 0);
                Debug.Assert(!double.IsInfinity(2 * v / x));

                int binaryScale = 0;
                for (int k = 0; k < n; k++) {
                    double fact = 2 * (v - k) / x;
                    double next = fact * current - prev;

                    if (double.IsInfinity(next)) {
                        int exp = 0;
                        current = Math2.Frexp(current, out exp);
                        prev = Math2.Ldexp(prev, -exp);
                        binaryScale += exp;

                        next = fact * current - prev;

                    }

                    prev = current;
                    current = next;
                }

                {
                    int exp = 0;
                    current = Math2.Frexp(current, out exp);
                    prev = Math2.Ldexp(prev, -exp);
                    binaryScale += exp;
                }

                return (current, prev, binaryScale);
            }




            /// <summary>
            /// Returns (J{v-n}(x), J{v-n+1}(x)) OR (Y{v-n}(x), Y{v-n+1}(x)) after n backward recurrences
            /// </summary>
            /// <param name="v">Order of the current parameter</param>
            /// <param name="x">Argument</param>
            /// <param name="n">Number of backward recurrence steps. 0 returns (current,prev)</param>
            /// <param name="current">J{v}(x) OR Y{v}(x)</param>
            /// <param name="prev">J{v+1}(x) OR Y{v+1}(x)</param>
            /// <returns></returns>
            public static (double JYvmn, double JYvmnp1) BackwardJY(double v, double x, int n, double current, double prev)
            {
                Debug.Assert(n >= 0);
                Debug.Assert(!double.IsInfinity(2 * v / x));

                int scale = 0;
                for (int k = 0; k < n; k++) {
                    double fact = 2 * (v - k) / x;
                    double next = fact * current - prev;

                    if (double.IsInfinity(next)) {
                        int exp = 0;
                        current = Math2.Frexp(current, out exp);
                        prev = Math2.Ldexp(prev, -exp);
                        scale += exp;

                        next = fact * current - prev;

                    }

                    prev = current;
                    current = next;
                }

                if (scale != 0) {
                    current = Math2.Ldexp(current, scale);
                    prev = Math2.Ldexp(prev, scale);
                }

                return (current, prev);
            }


            /// <summary>
            /// Returns (K{v+n}, K{v+n-1}) after n forward recurrences
            /// </summary>
            /// <param name="v">Order of the current parameter</param>
            /// <param name="x">Argument</param>
            /// <param name="n">Number of forward recurrence steps. 0 returns (current,prev)</param>
            /// <param name="current">K{v}(x)</param>
            /// <param name="prev">K{v-1}(x)</param>
            /// <returns></returns>
            public static (double Kvpn, double Kvpnm1) ForwardK(double v, double x, int n, double current, double prev)
            {
                Debug.Assert(n >= 0);

                // binary scale
                int binaryScale = 0;
                for (int k = 0; k < n; k++) {
                    double factor = 2 * (v + k) / x;
                    double value = factor * current + prev;

                    if (double.IsInfinity(value)) {
                        int exp = 0;
                        current = Math2.Frexp(current, out exp);
                        prev = Math2.Ldexp(prev, -exp);
                        binaryScale += exp;
                        value = factor * current + prev;
                    }

                    prev = current;
                    current = value;
                }

                if (binaryScale != 0) {
                    current = Math2.Ldexp(current, binaryScale);
                    prev = Math2.Ldexp(prev, binaryScale);

                }

                return (current, prev);
            }

            /// <summary>
            /// Returns binary scaled (K{v+n}, K{v+n-1}, scale) after n forward recurrences.
            /// Note: to get actual values, multiply by 2^BScale.
            /// </summary>
            /// <param name="v">Order of the current parameter</param>
            /// <param name="x">Argument</param>
            /// <param name="n">Number of forward recurrence steps. 0 returns (current,prev)</param>
            /// <param name="current">K{v}(x)</param>
            /// <param name="prev">K{v-1}(x)</param>
            /// <returns></returns>
            public static (double Kvpn, double Kvpnm1, int BScale) ForwardK_B(double v, double x, int n, double current, double prev)
            {
                Debug.Assert(n >= 0);

                // binary scale
                int binaryScale = 0;
                for (int k = 0; k < n; k++) {
                    double factor = 2 * (v + k) / x;
                    double value = factor * current + prev;

                    if (double.IsInfinity(value)) {
                        int exp = 0;
                        current = Math2.Frexp(current, out exp);
                        prev = Math2.Ldexp(prev, -exp);
                        binaryScale += exp;
                        value = factor * current + prev;
                    }

                    prev = current;
                    current = value;
                }

                {
                    int exp = 0;
                    current = Math2.Frexp(current, out exp);
                    prev = Math2.Ldexp(prev, -exp);
                    binaryScale += exp;
                }

                return (current, prev, binaryScale);
            }




        }
    }
}