//  Copyright (c) 2017 Rocco Romeo
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

//  Ported to .NET based on the original:
//      Copyright (c) 2006 Xiaogang Zhang, Boost Software License v1.0
//      Copyright (c) 2017 John Maddock, Boost Software License v1.0


using System;

namespace MathXK {

    static partial class _Bessel {

        /// <summary>
        /// Returns I{0}(x), or if expScale == true, e^-|x|*I{0}(x)
        /// </summary>
        /// <param name="x"></param>
        /// <param name="expScale"></param>
        /// <returns></returns>
        public static double I0(double x, bool expScale = false)
        {
            // Modified Bessel function of the first kind of order zero
            // we use the approximating forms derived in:
            // "Rational Approximations for the Modified Bessel Function of the First Kind – I0(x) for Computations with Double Precision"
            // by Pavel Holoborodko, 
            // see http://www.advanpix.com/2015/11/11/rational-approximations-for-the-modified-bessel-function-of-the-first-kind-i0-computations-double-precision
            // The actual coefficients used are boost's, and extend Pavel's work to precision's other than double.

            if (x < 0)
                x = -x;                         // even function

            if (x == 0)
                return 1;

            if (x < 7.75) {

                // Bessel I0 over[10 ^ -16, 7.75]
                // Max error in interpolated form : 3.042e-18
                // Max Error found at double precision = Poly : 5.106609e-16 Cheb : 5.239199e-16

                const double c0 = 1.00000000000000000e+00;
                const double c1 = 2.49999999999999909e-01;
                const double c2 = 2.77777777777782257e-02;
                const double c3 = 1.73611111111023792e-03;
                const double c4 = 6.94444444453352521e-05;
                const double c5 = 1.92901234513219920e-06;
                const double c6 = 3.93675991102510739e-08;
                const double c7 = 6.15118672704439289e-10;
                const double c8 = 7.59407002058973446e-12;
                const double c9 = 7.59389793369836367e-14;
                const double c10 = 6.27767773636292611e-16;
                const double c11 = 4.34709704153272287e-18;
                const double c12 = 2.63417742690109154e-20;
                const double c13 = 1.13943037744822825e-22;
                const double c14 = 9.07926920085624812e-25;

                var a = x * x / 4;
                var r = c0 + a * (c1 + a * (c2 + a * (c3 + a * (c4 + a * (c5 + a * (c6 + a * (c7 + a * (c8 + a * (c9 + a * (c10 + a * (c11 + a * (c12 + a * (c13 + a * c14)))))))))))));
                var result = a * r + 1;
                if (expScale)
                    result *= Math.Exp(-x);
                return result;

            } else if (x < 500) {

                // Max error in interpolated form : 1.685e-16
                // Max Error found at double precision = Poly : 2.575063e-16 Cheb : 2.247615e+00

                const double c0 = 3.98942280401425088e-01;
                const double c1 = 4.98677850604961985e-02;
                const double c2 = 2.80506233928312623e-02;
                const double c3 = 2.92211225166047873e-02;
                const double c4 = 4.44207299493659561e-02;
                const double c5 = 1.30970574605856719e-01;
                const double c6 = -3.35052280231727022e+00;
                const double c7 = 2.33025711583514727e+02;
                const double c8 = -1.13366350697172355e+04;
                const double c9 = 4.24057674317867331e+05;
                const double c10 = -1.23157028595698731e+07;
                const double c11 = 2.80231938155267516e+08;
                const double c12 = -5.01883999713777929e+09;
                const double c13 = 7.08029243015109113e+10;
                const double c14 = -7.84261082124811106e+11;
                const double c15 = 6.76825737854096565e+12;
                const double c16 = -4.49034849696138065e+13;
                const double c17 = 2.24155239966958995e+14;
                const double c18 = -8.13426467865659318e+14;
                const double c19 = 2.02391097391687777e+15;
                const double c20 = -3.08675715295370878e+15;
                const double c21 = 2.17587543863819074e+15;

                var a = 1 / x;
                var s11 = c11 + a * (c12 + a * (c13 + a * (c14 + a * (c15 + a * (c16 + a * (c17 + a * (c18 + a * (c19 + a * (c20 + a * c21)))))))));
                var r = c0 + a * (c1 + a * (c2 + a * (c3 + a * (c4 + a * (c5 + a * (c6 + a * (c7 + a * (c8 + a * (c9 + a * (c10 + a * s11))))))))));

                if (expScale)
                    return r / Math.Sqrt(x);
                return Math.Exp(x) / Math.Sqrt(x) * r;

            } else {

                // Max error in interpolated form : 2.437e-18
                // Max Error found at double precision = Poly : 1.216719e-16

                const double c0 = 3.98942280401432905e-01;
                const double c1 = 4.98677850491434560e-02;
                const double c2 = 2.80506308916506102e-02;
                const double c3 = 2.92179096853915176e-02;
                const double c4 = 4.53371208762579442e-02;

                var a = 1 / x;
                var r = c0 + a * (c1 + a * (c2 + a * (c3 + a * c4)));

                // Unscaled result = Exp(x)/Sqrt(x) * r;

                if (expScale)
                    return r / Math.Sqrt(x);

                // careful not to overflow too soon
                var ex = Math.Exp(x / 2);
                return (ex / Math.Sqrt(x)) * r * ex;
            }
        }

        /// <summary>
        /// Returns I{1}(x), or if expScale == true, e^-|x|*I{1}(x)
        /// </summary>
        /// <param name="x"></param>
        /// <param name="expScale"></param>
        /// <returns></returns>
        public static double I1(double x, bool expScale = false)
        {
            // Modified Bessel function of the first kind of order zero
            // we use the approximating forms derived in:
            // "Rational Approximations for the Modified Bessel Function of the First Kind – I1(x) for Computations with Double Precision"
            // by Pavel Holoborodko, 
            // see http://www.advanpix.com/2015/11/12/rational-approximations-for-the-modified-bessel-function-of-the-first-kind-i1-for-computations-with-double-precision/
            // The actual coefficients used are boost's, and extend Pavel's work to precision's other than double.

            double value;

            if (x == 0)
                return 0;

            if (x < 7.75) {

                // Bessel I0 over[10 ^ -16, 7.75]
                // Max error in interpolated form: 5.639e-17
                // Max Error found at double precision = Poly: 1.795559e-16

                const double c0 = 8.333333333333333803e-02;
                const double c1 = 6.944444444444341983e-03;
                const double c2 = 3.472222222225921045e-04;
                const double c3 = 1.157407407354987232e-05;
                const double c4 = 2.755731926254790268e-07;
                const double c5 = 4.920949692800671435e-09;
                const double c6 = 6.834657311305621830e-11;
                const double c7 = 7.593969849687574339e-13;
                const double c8 = 6.904822652741917551e-15;
                const double c9 = 5.220157095351373194e-17;
                const double c10 = 3.410720494727771276e-19;
                const double c11 = 1.625212890947171108e-21;
                const double c12 = 1.332898928162290861e-23;

                var a = x * x / 4;
                var r = c0 + a * (c1 + a * (c2 + a * (c3 + a * (c4 + a * (c5 + a * (c6 + a * (c7 + a * (c8 + a * (c9 + a * (c10 + a * (c11 + a * c12)))))))))));

                value = (1 + a * (0.5 + a * r)) * (0.5 * x);
                if (expScale)
                    value *= Math.Exp(-x);

            } else if (x < 500) {

                // Max error in interpolated form: 1.796e-16
                // Max Error found at double precision = Poly: 2.898731e-16

                const double c0 = 3.989422804014406054e-01;
                const double c1 = -1.496033551613111533e-01;
                const double c2 = -4.675104253598537322e-02;
                const double c3 = -4.090895951581637791e-02;
                const double c4 = -5.719036414430205390e-02;
                const double c5 = -1.528189554374492735e-01;
                const double c6 = 3.458284470977172076e+00;
                const double c7 = -2.426181371595021021e+02;
                const double c8 = 1.178785865993440669e+04;
                const double c9 = -4.404655582443487334e+05;
                const double c10 = 1.277677779341446497e+07;
                const double c11 = -2.903390398236656519e+08;
                const double c12 = 5.192386898222206474e+09;
                const double c13 = -7.313784438967834057e+10;
                const double c14 = 8.087824484994859552e+11;
                const double c15 = -6.967602516005787001e+12;
                const double c16 = 4.614040809616582764e+13;
                const double c17 = -2.298849639457172489e+14;
                const double c18 = 8.325554073334618015e+14;
                const double c19 = -2.067285045778906105e+15;
                const double c20 = 3.146401654361325073e+15;
                const double c21 = -2.213318202179221945e+15;

                var a = 1 / x;
                var s11 = c11 + a * (c12 + a * (c13 + a * (c14 + a * (c15 + a * (c16 + a * (c17 + a * (c18 + a * (c19 + a * (c20 + a * c21)))))))));
                var r = c0 + a * (c1 + a * (c2 + a * (c3 + a * (c4 + a * (c5 + a * (c6 + a * (c7 + a * (c8 + a * (c9 + a * (c10 + a * s11))))))))));

                if (expScale)
                    value = r / Math.Sqrt(x);
                else
                    value = Math.Exp(x) / Math.Sqrt(x) * r;

            } else {

                // Max error in interpolated form: 1.320e-19
                // Max Error found at double precision = Poly: 7.065357e-17

                const double c0 = 3.989422804014314820e-01;
                const double c1 = -1.496033551467584157e-01;
                const double c2 = -4.675105322571775911e-02;
                const double c3 = -4.090421597376992892e-02;
                const double c4 = -5.843630344778927582e-02;

                var a = 1 / x;
                var r = c0 + a * (c1 + a * (c2 + a * (c3 + a * c4)));

                if (expScale)
                    value = r / Math.Sqrt(x);
                else {
                    // careful not to overflow too soon
                    var ex = Math.Exp(x / 2);
                    return (ex / Math.Sqrt(x)) * r * ex;

                }
            }

            if (x < 0)
                value *= -1; // odd function

            return value;

        }

    }

} // namespaces


