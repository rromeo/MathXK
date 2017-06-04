//  Copyright (c) 2013 Rocco Romeo
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

//  Ported to .NET based on the original:
//      Copyright (c) 2006 Xiaogang Zhang, Boost Software License v1.0

using System;
using System.Diagnostics;
using MathXK.Numerics;

namespace MathXK {

    internal static partial class _Bessel {

        /// <summary>
        /// Returns K{0}(x), or if expScaled == true, e^x * K{0}(x)
        /// </summary>
        /// <param name="x"></param>
        /// <param name="expScaled"></param>
        /// <returns></returns>
        public static double K0(double x, bool expScaled = false)
        {
            if (!(x >= 0)) {
                Policies.ReportDomainError("BesselK(v: 0, x: {0}) : Requires x >= 0; complex number result not supported", x);
                return double.NaN;
            }
            if (x == 0)
                return double.PositiveInfinity;

            // Modified Bessel function of the second kind of order zero
            // minimax rational approximations on intervals, see
            // Russon and Blair, Chalk River Report AECL-3461, 1969,
            // as revised by Pavel Holoborodko in "Rational Approximations 
            // for the Modified Bessel Function of the Second Kind - K0(x) 
            // for Computations with Double Precision", see 
            // http://www.advanpix.com/2015/11/25/rational-approximations-for-the-modified-bessel-function-of-the-second-kind-k0-for-computations-with-double-precision/
            //
            // The actual coefficients used are boost's (derived by JM)
            // since we extend to both greater and lesser precision than the
            // references above.  We can also improve performance WRT to
            // Holoborodko without loss of precision.

            double value;
            if (x <= 1) {

                // Maximum Deviation Found:                     6.077e-17
                // Expected Error Term : -6.077e-17
                // Maximum Relative Change in Control Points : 7.797e-02
                // Max Error found at double precision = Poly : 1.003156e-16
                const double Y = 1.137250900268554688;
                const double p0 = -1.372509002685546267e-01;
                const double p1 = 2.574916117833312855e-01;
                const double p2 = 1.395474602146869316e-02;
                const double p3 = 5.445476986653926759e-04;
                const double p4 = 7.125159422136622118e-06;

                const double q0 = 1.000000000000000000e+00;
                const double q1 = -5.458333438017788530e-02;
                const double q2 = 1.291052816975251298e-03;
                const double q3 = -1.367653946978586591e-05;

                var a = x * x / 4;
                double P = p0 + a * (p1 + a * (p2 + a * (p3 + a * p4)));
                double Q = q0 + a * (q1 + a * (q2 + a * q3));
                var r1 = (P / Q + Y) * a + 1;

                // Maximum Deviation Found:                     3.429e-18
                // Expected Error Term : 3.392e-18
                // Maximum Relative Change in Control Points : 2.041e-02
                // Max Error found at double precision = Poly : 2.513112e-16

                const double c0 = 1.159315156584124484e-01;
                const double c1 = 2.789828789146031732e-01;
                const double c2 = 2.524892993216121934e-02;
                const double c3 = 8.460350907213637784e-04;
                const double c4 = 1.491471924309617534e-05;
                const double c5 = 1.627106892422088488e-07;
                const double c6 = 1.208266102392756055e-09;
                const double c7 = 6.611686391749704310e-12;

                a = x * x;
                var r2 = c0 + a * (c1 + a * (c2 + a * (c3 + a * (c4 + a * (c5 + a * (c6 + a * c7))))));

                value = r2 - Math.Log(x) * r1;
                if (expScaled)
                    value *= Math.Exp(x);

            } else {
                // Maximum Deviation Found:                     4.316e-17
                // Expected Error Term : 9.570e-18
                // Maximum Relative Change in Control Points : 2.757e-01
                // Max Error found at double precision = Poly : 1.001560e-16

                const double Y = 1;

                const double p0 = 2.533141373155002416e-01;
                const double p1 = 3.628342133984595192e+00;
                const double p2 = 1.868441889406606057e+01;
                const double p3 = 4.306243981063412784e+01;
                const double p4 = 4.424116209627428189e+01;
                const double p5 = 1.562095339356220468e+01;
                const double p6 = -1.810138978229410898e+00;
                const double p7 = -1.414237994269995877e+00;
                const double p8 = -9.369168119754924625e-02;

                const double q0 = 1.000000000000000000e+00;
                const double q1 = 1.494194694879908328e+01;
                const double q2 = 8.265296455388554217e+01;
                const double q3 = 2.162779506621866970e+02;
                const double q4 = 2.845145155184222157e+02;
                const double q5 = 1.851714491916334995e+02;
                const double q6 = 5.486540717439723515e+01;
                const double q7 = 6.118075837628957015e+00;
                const double q8 = 1.586261269326235053e-01;

                var a = 1 / x;
                double P = p0 + a * (p1 + a * (p2 + a * (p3 + a * (p4 + a * (p5 + a * (p6 + a * (p7 + a * p8)))))));
                double Q = q0 + a * (q1 + a * (q2 + a * (q3 + a * (q4 + a * (q5 + a * (q6 + a * (q7 + a * q8)))))));
                value = (P / Q + Y);

                value /= Math.Sqrt(x);
                if (!expScaled) {
                    // The following is: value * Exp(-x)
                    // when x is large, don't underflow too soon
                    var ex = Math.Exp(-x / 2);
                    value = (value * ex) * ex;
                }
            }

            return value;
        }

        /// <summary>
        /// Returns K{1}(x), or if expScaled == true, e^x * K{1}(x)
        /// </summary>
        /// <param name="x"></param>
        /// <param name="expScaled"></param>
        /// <returns></returns>
        public static double K1(double x, bool expScaled = false)
        {
            if (x < 0) {
                Policies.ReportDomainError("BesselK(v: 1, x: {0}): Requires x >= 0; complex number result not supported", x);
                return double.NaN;
            }
            if (x == 0)
                return double.PositiveInfinity;

            // Modified Bessel function of the second kind of order zero
            // minimax rational approximations on intervals, see
            // Russon and Blair, Chalk River Report AECL-3461, 1969,
            // as revised by Pavel Holoborodko in "Rational Approximations 
            // for the Modified Bessel Function of the Second Kind - K0(x) 
            // for Computations with Double Precision", see 
            // http://www.advanpix.com/2016/01/05/rational-approximations-for-the-modified-bessel-function-of-the-second-kind-k1-for-computations-with-double-precision/
            //
            // The actual coefficients used are boost's (derived by JM)
            // since we extend to both greater and lesser precision than the
            // references above.  We can also improve performance WRT to
            // Holoborodko without loss of precision.

            double value;
            if (x <= 1) {
                double r1, r2;
                {
                    // Maximum Deviation Found:                     1.922e-17
                    // Expected Error Term : 1.921e-17
                    // Maximum Relative Change in Control Points : 5.287e-03
                    // Max Error found at double precision = Poly : 2.004747e-17
                    const double Y = 8.69547128677368164e-02;
                    const double p0 = -3.62137953440350228e-03;
                    const double p1 = 7.11842087490330300e-03;
                    const double p2 = 1.00302560256614306e-05;
                    const double p3 = 1.77231085381040811e-06;


                    const double q0 = 1.00000000000000000e+00;
                    const double q1 = -4.80414794429043831e-02;
                    const double q2 = 9.85972641934416525e-04;
                    const double q3 = -8.91196859397070326e-06;

                    var a = x * x / 4;
                    double P = p0 + a * (p1 + a * (p2 + a * p3));
                    double Q = q0 + a * (q1 + a * (q2 + a * q3));

                    r1 = (P / Q + Y);
                    r1 = (1 + a * (0.5 + r1 * a)) * (0.5 * x);
                }
                {
                    // Maximum Deviation Found:                     4.053e-17
                    // Expected Error Term : -4.053e-17
                    // Maximum Relative Change in Control Points : 3.103e-04
                    // Max Error found at double precision = Poly : 1.246698e-16

                    const double p0 = -3.07965757829206184e-01;
                    const double p1 = -7.80929703673074907e-02;
                    const double p2 = -2.70619343754051620e-03;
                    const double p3 = -2.49549522229072008e-05;

                    const double q0 = 1.00000000000000000e+00;
                    const double q1 = -2.36316836412163098e-02;
                    const double q2 = 2.64524577525962719e-04;
                    const double q3 = -1.49749618004162787e-06;

                    var a = x * x;
                    double P = p0 + a * (p1 + a * (p2 + a * p3));
                    double Q = q0 + a * (q1 + a * (q2 + a * q3));
                    r2 = P / Q;
                }

                value = r2 * x + 1 / x + Math.Log(x) * r1;
                if (expScaled)
                    value *= Math.Exp(x);


            } else {
                // Maximum Deviation Found:                     8.883e-17
                // Expected Error Term : -1.641e-17
                // Maximum Relative Change in Control Points : 2.786e-01
                // Max Error found at double precision = Poly : 1.258798e-16

                const double Y = 1.45034217834472656;

                const double p0 = -1.97028041029226295e-01;
                const double p1 = -2.32408961548087617e+00;
                const double p2 = -7.98269784507699938e+00;
                const double p3 = -2.39968410774221632e+00;
                const double p4 = 3.28314043780858713e+01;
                const double p5 = 5.67713761158496058e+01;
                const double p6 = 3.30907788466509823e+01;
                const double p7 = 6.62582288933739787e+00;
                const double p8 = 3.08851840645286691e-01;


                const double q0 = 1.00000000000000000e+00;
                const double q1 = 1.41811409298826118e+01;
                const double q2 = 7.35979466317556420e+01;
                const double q3 = 1.77821793937080859e+02;
                const double q4 = 2.11014501598705982e+02;
                const double q5 = 1.19425262951064454e+02;
                const double q6 = 2.88448064302447607e+01;
                const double q7 = 2.27912927104139732e+00;
                const double q8 = 2.50358186953478678e-02;

                var a = 1 / x;
                double P = p0 + a * (p1 + a * (p2 + a * (p3 + a * (p4 + a * (p5 + a * (p6 + a * (p7 + a * p8)))))));
                double Q = q0 + a * (q1 + a * (q2 + a * (q3 + a * (q4 + a * (q5 + a * (q6 + a * (q7 + a * q8)))))));
                value = (P / Q + Y);
                value /= Math.Sqrt(x);

                if (!expScaled) {
                    // The following is: value * Exp(-x)
                    // when x is large, don't underflow too soon
                    var ex = Math.Exp(-x / 2);
                    value = (value * ex) * ex;
                }
            }

            return value;

        }

        /// <summary>
        /// Returns K{n}(x) for integer order
        /// </summary>
        /// <param name="n"></param>
        /// <param name="x"></param>
        /// <returns></returns>
        public static double KN(int n, double x)
        {
            if (x < 0) {
                Policies.ReportDomainError("BesselK(v: {0}, x: {1}): Requires x >= 0 for real result", n, x);
                return double.NaN;
            }
            if (x == 0)
                return double.PositiveInfinity;


            // even function
            // K{-n}(z) = K{n}(z)
            if (n < 0)
                n = -n;

            if (n == 0)
                return K0(x);

            if (n == 1)
                return K1(x);

            double v = n;

            // Hankel is fast, and reasonably accurate, saving us from many recurrences.
            if (x >= HankelAsym.IKMinX(v))
                return HankelAsym.K(v, x);

            // the uniform expansion is here as a last resort 
            // to limit the number of recurrences, but it is less accurate.
            if (UniformAsym.IsIKAvailable(v, x))
                return UniformAsym.K(v, x);

            // Since K{v}(x) has a (e^-x)/sqrt(x) multiplier
            // using recurrence can underflow too quickly for large x,
            // so, use a scaled version
            double result;
            if (x > 1) {

                double prev = K0(x, true);
                double current = K1(x, true);

                // for large v and x this number can get very large
                // maximum observed K(1000,10) = 2^6211 

                var (Kv, _, binaryScale) = Recurrence.ForwardK_B(1, x, n - 1, current, prev);


                // Compute: value * 2^(binaryScale) * e^-x

                if (x < -DoubleX.MinLogValue) {
                    DoubleX exs = DoubleX.Ldexp(DoubleX.Exp(-x), binaryScale);
                    result = Math2.Ldexp(Kv * exs.Mantissa, exs.Exponent);
                } else {
                    result = Math.Exp(-x + Math.Log(Kv) + binaryScale * Constants.Ln2);
                }

            } else {


                double prev = K0(x);
                double current = K1(x);

                result = Recurrence.ForwardK(1, x, n - 1, current, prev).Kvpn;
            }


            return result;
        }

    }

} // namespaces


