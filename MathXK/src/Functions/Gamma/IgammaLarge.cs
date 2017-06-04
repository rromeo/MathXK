//  Copyright (c) 2013 Rocco Romeo
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

//  Ported to .NET based on the original:
//      Copyright (c) John Maddock 2006, Boost Software License v1.0

using System;

namespace MathXK {

    static partial class _Igamma {

        /// <summary>
        /// Asymptotic expansions of the incomplete gamma functions 
        /// P(a, x) and Q(a, x), used when a is large and x ~ a.
        /// </summary>
        static class TemmeSymmetricAsym {

            //
            // The primary reference is:
            //
            // "The Asymptotic Expansion of the Incomplete Gamma Functions"
            // N. M. Temme.
            // Siam J. Math Anal. Vol 10 No 4, July 1979, p757.
            //
            // A different way of evaluating these expansions,
            // plus a lot of very useful background information is in:
            // 
            // "A Set of Algorithms For the Incomplete Gamma Functions."
            // N. M. Temme.
            // Probability in the Engineering and Informational Sciences,
            // 8, 1994, 291.
            //
            // An alternative implementation is in:
            //
            // "Computation of the Incomplete Gamma Function Ratios and their Inverse."
            // A. R. Didonato and A. H. Morris.
            // ACM TOMS, Vol 12, No 4, Dec 1986, p377.
            //
            // There are various versions of the same code below, each accurate
            // to a different precision.  To understand the code, refer to Didonato
            // and Morris, from Eq 17 and 18 onwards.
            //
            // The coefficients used here are not taken from Didonato and Morris:
            // the domain over which these expansions are used is slightly different
            // to theirs, and their constants are not quite accurate enough for
            // 128-bit long double's.  Instead the coefficients were calculated
            // using the methods described by Temme p762 from Eq 3.8 onwards.
            // The values obtained agree with those obtained by Didonato and Morris
            // (at least to the first 30 digits that they provide).
            // At double precision the degrees of polynomial required for full
            // machine precision are close to those recomended to Didonato and Morris,
            // but of course many more terms are needed for larger types.
            //

            private static readonly double[][] _C = {

                new double[] {
                    -0.33333333333333333,
                    0.083333333333333333,
                    -0.014814814814814815,
                    0.0011574074074074074,
                    0.0003527336860670194,
                    -0.00017875514403292181,
                    0.39192631785224378e-4,
                    -0.21854485106799922e-5,
                    -0.185406221071516e-5,
                    0.8296711340953086e-6,
                    -0.17665952736826079e-6,
                    0.67078535434014986e-8,
                    0.10261809784240308e-7,
                    -0.43820360184533532e-8,
                    0.91476995822367902e-9,
                },

                new double[] {
                    -0.0018518518518518519,
                    -0.0034722222222222222,
                    0.0026455026455026455,
                    -0.00099022633744855967,
                    0.00020576131687242798,
                    -0.40187757201646091e-6,
                    -0.18098550334489978e-4,
                    0.76491609160811101e-5,
                    -0.16120900894563446e-5,
                    0.46471278028074343e-8,
                    0.1378633446915721e-6,
                    -0.5752545603517705e-7,
                    0.11951628599778147e-7,
                },


                new double[] {
                    0.0041335978835978836,
                    -0.0026813271604938272,
                    0.00077160493827160494,
                    0.20093878600823045e-5,
                    -0.00010736653226365161,
                    0.52923448829120125e-4,
                    -0.12760635188618728e-4,
                    0.34235787340961381e-7,
                    0.13721957309062933e-5,
                    -0.6298992138380055e-6,
                    0.14280614206064242e-6,
                },


                new double[] {
                    0.00064943415637860082,
                    0.00022947209362139918,
                    -0.00046918949439525571,
                    0.00026772063206283885,
                    -0.75618016718839764e-4,
                    -0.23965051138672967e-6,
                    0.11082654115347302e-4,
                    -0.56749528269915966e-5,
                    0.14230900732435884e-5,
                },


                new double[] {
                    -0.0008618882909167117,
                    0.00078403922172006663,
                    -0.00029907248030319018,
                    -0.14638452578843418e-5,
                    0.66414982154651222e-4,
                    -0.39683650471794347e-4,
                    0.11375726970678419e-4,
                },

                new double[] {
                    -0.00033679855336635815,
                    -0.69728137583658578e-4,
                    0.00027727532449593921,
                    -0.00019932570516188848,
                    0.67977804779372078e-4,
                    0.1419062920643967e-6,
                    -0.13594048189768693e-4,
                    0.80184702563342015e-5,
                    -0.22914811765080952e-5,
                },

                new double[] {
                    0.00053130793646399222,
                    -0.00059216643735369388,
                    0.00027087820967180448,
                    0.79023532326603279e-6,
                    -0.81539693675619688e-4,
                    0.56116827531062497e-4,
                    -0.18329116582843376e-4,
                },

                new double[] {
                    0.00034436760689237767,
                    0.51717909082605922e-4,
                    -0.00033493161081142236,
                    0.0002812695154763237,
                    -0.00010976582244684731,
                },


                new double[] {
                    -0.00065262391859530942,
                    0.00083949872067208728,
                    -0.00043829709854172101,
                },

                new double[] {  -0.00059676129019274625 },


            };

            /// <summary>
            /// Asymptotic expansions of the incomplete gamma function P(a, x), 
            /// used when a is large and x ~= a.
            /// </summary>
            /// <param name="a"></param>
            /// <param name="x"></param>
            /// <returns></returns>
            public static double GammaP(double a, double x)
            {
                double sigma = (x - a) / a;
                double phi = -Math2.Log1pmx(sigma);
                double y = a * phi;
                double z = Math.Sqrt(2 * phi);
                if (x < a)
                    z = -z;

                double p = Polynomial.Eval(_C, z, 1 / a);
                double result = p * (Math.Exp(-y) * (Constants.RecipSqrt2PI / Math.Sqrt(a)));
                if (x < a)
                    result = -result;

                result += Math2.Erfc(Math.Sqrt(y)) / 2;

                return result;
            }
        }
    }

}  // namespace 



