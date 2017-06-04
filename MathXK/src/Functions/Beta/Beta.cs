//  Copyright (c) 2017 Rocco Romeo
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

//  Ported to .NET based on the original:
//      Copyright (c) John Maddock 2006, Boost Software License v1.0
//  History:
//      JM wrote original C++ code.
//      RR ported to C#, improved Beta accuracy, and added LogBeta.

//#define EXTRA_DEBUG

using System;
using System.Diagnostics;

namespace MathXK {


    public partial class Math2
    {

        /// <summary>
        /// Returns the Beta function
        /// <para>B(a,b) = Γ(a)*Γ(b)/Γ(a+b)</para>
        /// </summary>
        /// <param name="a">Requires finite a > 0</param>
        /// <param name="b">Requires finite b > 0</param>
        public static double Beta(double a, double b)
        {
            double c = a + b;

            if ((!(a > 0) || double.IsInfinity(a))
                || (!(b > 0) || double.IsInfinity(b))) {
                Policies.ReportDomainError("Beta(a: {0}, b: {1}): Requires finite a,b > 0", a, b);
                return double.NaN;
            }
            if (double.IsInfinity(c)) {
                Policies.ReportDomainError("Beta(a: {0}, b: {1}): Requires finite c == a+b: c = {2}", a, b, c);
                return double.NaN;
            }


            // Special cases:
            if ((c == a) && (b < DoubleLimits.MachineEpsilon))
                return Math2.Tgamma(b);
            if ((c == b) && (a < DoubleLimits.MachineEpsilon))
                return Math2.Tgamma(a);
            if (b == 1)
                return 1 / a;
            if (a == 1)
                return 1 / b;

            // B(a,b) == B(b, a)
            if (a < b)
                Utility.Swap(ref a, ref b);

            // from this point a >= b

            if (a < 1) {
                // When x < TinyX, Γ(x) ~ 1/x 
                const double SmallX = DoubleLimits.MachineEpsilon / Constants.EulerMascheroni;
                if ( a < SmallX ) {
                    if (c < SmallX)
                        return (c / a) / b;
                    Debug.Assert(c <= GammaSmall.UpperLimit);
                    return 1/(GammaSmall.Tgamma(c) * a * b);
                }
                if ( b < SmallX ) 
                    return Tgamma(a) / (b * Tgamma(c));
                return Tgamma(a) * ( Tgamma(b)/Tgamma(c) );
            }

            // Our Stirling series is more accurate than Lanczos
            if (a >= StirlingGamma.LowerLimit ) {
                if (b >= StirlingGamma.LowerLimit)
                    return StirlingGamma.Beta(a, b);

                // for large a, Γ(a)/Γ(a + b) ~ a^-b, so don't underflow too soon
                if (b < 1 || b * Math.Log(a) <= -DoubleLimits.MinLogValue)
                    return Tgamma(b) * StirlingGamma.TgammaDeltaRatio(a, b);

                // fall through for very large a small 1 < b < LowerLimit

            }


            return Lanczos.Beta(a, b);
        }

        /// <summary>
        /// Returns Log(Beta(a,b))
        /// </summary>
        /// <param name="a">Requires finite a > 0</param>
        /// <param name="b">Requires finite b > 0</param>
        public static double LogBeta(double a, double b)
        {
            double c = a + b;

            if ((!(a > 0) || double.IsInfinity(a))
                || (!(b > 0) || double.IsInfinity(b))) {
                Policies.ReportDomainError("Beta(a: {0}, b: {1}): Requires finite a,b > 0", a, b);
                return double.NaN;
            }
            if (double.IsInfinity(c)) {
                Policies.ReportDomainError("Beta(a: {0}, b: {1}): Requires finite c == a+b: c = {2}", a, b, c);
                return double.NaN;
            }

            // Special cases:
            if ((c == a) && (b < DoubleLimits.MachineEpsilon))
                return Math2.Lgamma(b);
            if ((c == b) && (a < DoubleLimits.MachineEpsilon))
                return Math2.Lgamma(a);
            if (b == 1)
                return -Math.Log(a);
            if (a == 1)
                return -Math.Log(b);

            // B(a,b) == B(b, a)
            if (a < b)
                Utility.Swap(ref a, ref b);

            // from this point a >= b

            if (a < 1) {
                // When x < TinyX, Γ(x) ~ 1/x 
                const double SmallX = DoubleLimits.MachineEpsilon / Constants.EulerMascheroni;
                if (a < SmallX) {
                    if (c < SmallX) {
                        // In this range, Beta(a, b) ~= 1/a + 1/b 
                        // so, we won't have an argument overflow as long as 
                        // the min(a, b) >= 2*DoubleLimits.MinNormalValue

                        if (b < 2 * DoubleLimits.MinNormalValue)
                            return Math.Log(c / a) - Math.Log(b);

                        return Math.Log((c / a) / b);
                    }
                    Debug.Assert(c <= GammaSmall.UpperLimit);
                    return -Math.Log((GammaSmall.Tgamma(c) * a * b));
                }
                if (b < SmallX)
                    return Math.Log(Tgamma(a) / (b * Tgamma(c)));

                return Math.Log((Tgamma(a) / Tgamma(c)) * Tgamma(b));

            } else if (a >= StirlingGamma.LowerLimit) {

                if ( b >= StirlingGamma.LowerLimit )
                    return StirlingGamma.LogBeta(a, b);

                return StirlingGamma.LgammaDelta(a, b) + Lgamma(b);
            }

            return Lanczos.LogBeta(a, b);
        }

    }



} // namespace