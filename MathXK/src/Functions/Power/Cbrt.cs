//  Copyright (c) 2013 Rocco Romeo
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

//  Ported to .NET based on the original:
//      Copyright (c) John Maddock 2006, Boost Software License, Version 1.0


namespace MathXK
{

    public partial class Math2
    {
        
        /// <summary>
        /// Table of 2^((n-2)/3) used in Cbrt
        /// </summary>
        private static readonly double[] _TwoPowThirds = {
            0.62996052494743658238360530363911,  // 2^-2/3
            0.79370052598409973737585281963615,  // 2^-1/3
            1,
            1.2599210498948731647672106072782,   // 2^1/3
            1.5874010519681994747517056392723,   // 2^2/3
        };

        /// <summary>
        /// Returns the cube root of <paramref name="x"/>
        /// </summary>
        /// <param name="x">Cbrt function argument</param>
        public static double Cbrt(double x)
        {
            if (x < 0)
                return -Cbrt(-x);

            if (double.IsNaN(x)) {
                Policies.ReportDomainError("Cbrt(x: {0}): NaN not allowed", x);
                return x;
            }

            if (double.IsInfinity(x))
                return x;

            if (x == 0)
                return 0;

            // let x = y* 2^n
            // Cbrt(y*(2^n)) = Cbrt(2^n)*Cbrt(y)

            int n;
            double y = Math2.Frexp(x, out n);

            // cbrt approximation for z in the range [0.5,1]
            // It's hard to say what number of terms gives the optimum
            // trade off between precision and performance, this seems
            // to be about the best for double precision.
            //
            // Maximum Deviation Found:                     1.231e-006
            // Expected Error Term:                         -1.231e-006
            // Maximum Relative Change in Control Points:   5.982e-004

            const double c0 = 0.37568269008611818;
            const double c1 = 1.3304968705558024;
            const double c2 = -1.4897101632445036;
            const double c3 = 1.2875573098219835;
            const double c4 = -0.6398703759826468;
            const double c5 = 0.13584489959258635;

            // Compute the polynomial
            double guess = c0 + y * (c1 + y * (c2 + y * (c3 + y * (c4 + y * c5))));

            //
            // Now Halley iteration.
            // We do this here rather than calling the root finding routine since we can
            // simplify the expressions algebraically, and don't need most of the error
            // checking of the boilerplate version as we know in advance that the function
            // is well behaved. 
            //
            // Since our appoximation function's error is O(10^-6), 
            // one Halley iteration should suffice.

            double g3 = guess * guess * guess;
            guess -= guess * (g3 - y) / (2 * g3 + y);

            int exp = n / 3;
            guess *= _TwoPowThirds[(n % 3) + 2];
            guess  = Math2.Ldexp(guess, exp);

            return guess;
        }
    }


} // namespaces





