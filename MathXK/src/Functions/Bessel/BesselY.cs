//  Copyright (c) 2013 Rocco Romeo
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

//  Ported to .NET based on the original:
//      Copyright (c) 2007 John Maddock, Boost Software License v1.0


using System;
using System.Diagnostics;

using MathXK.Numerics;

namespace MathXK
{

    internal static partial class _Bessel
    {

        //
        // Series form for BesselY as z -> 0, 
        // see: http://functions.wolfram.com/Bessel-TypeFunctions/BesselY/06/01/04/01/01/0003/
        // This series is only useful when the second term is small compared to the first
        // otherwise we get catestrophic cancellation errors.
        //
        // Approximating tgamma(v) by v^v, and assuming |tgamma(-z)| < eps we end up requiring:
        // eps/2 * v^v(x/2)^-v > (x/2)^v or Math.Log(eps/2) > v Math.Log((x/2)^2/v)
        //


        /// <summary>
        /// Series form for Y{v}(x) as z -&gt; 0 and v is not an integer
        /// <para>Y{v}(z) = (-((z/2)^v Cos(πv) Γ(-v))/π)) * 0F1(; 1+v; -z^2/4) - ((z/2)^-v Γ(v))/π) * 0F1(; 1-v; -z^2/4)</para>
        /// </summary>
        /// <param name="v"></param>
        /// <param name="x"></param>
        /// <returns></returns>
        /// <seealso href="http://functions.wolfram.com/Bessel-TypeFunctions/BesselY/26/01/02/0001/"/>
        public static DoubleX Y_SmallArgSeries(double v, double x)
        {
            Debug.Assert(v >= 0 && x >= 0);
            Debug.Assert( !Math2.IsInteger(v), "v cannot be an integer");


            DoubleX powDivGamma;
            double lnp = v * Math.Log(x / 2);
            if (v > DoubleLimits.MaxGamma || lnp < DoubleLimits.MinLogValue || lnp > DoubleLimits.MaxLogValue) {
                powDivGamma = DoubleX.Exp(lnp - Math2.Lgamma(v));
            } else {
                DoubleX p = Math.Pow(x / 2, v);
                DoubleX gam = Math2.Tgamma(v);
                powDivGamma = p / gam;
            }


            // Series 1: -((z/2)^-v Γ(v))/π) * 0F1(; 1-v; -z^2/4)

            double s1 = HypergeometricSeries.Sum0F1(1 - v, -x * x / 4);
            DoubleX result1 = -(s1/Math.PI)/ powDivGamma;

            // Series2: -((z/2)^v Cos(πv) Γ(-v))/π)) * 0F1(; 1+v; -z^2/4)
            // using the reflection formula: Γ(-v) = -π/(Γ(v) * v * Sin(π * v)) 
            // prefix = (z/2)^v * Cos(π * v)/(Γ(v) * v * Sin(π * v)) 
            // prefix = (z/2)^v * Cot(π * v) / (Γ(v) * v) 

            DoubleX result2 = DoubleX.Zero;
            double cot = Math2.CotPI(v);
            if (cot != 0) {
                double s2 = HypergeometricSeries.Sum0F1(v + 1, -x * x / 4);
                result2 = powDivGamma * cot * s2 / v; // Do this all in DoubleX
            }

            return (result1 - result2);
        }



    }

    public partial class Math2
    {

        /// <summary>
        /// Returns the Cylindrical Bessel Y function: Y<sub>v</sub>(x)
        /// </summary>
        /// <param name="v">Order</param>
        /// <param name="x">Argument</param>
        /// <returns></returns>
        public static double BesselY(double v, double x)
        {
            if (!(x >= 0)) {
                Policies.ReportDomainError("BesselY(v: {0}, x: {1}): Requires x >= 0 for real result", v, x);
                return double.NaN;
            }
            if (x == 0) {
                if (v == 0) {
                    Policies.ReportDomainError("BesselY(v: {0}, x: {1}): Requires x > 0 when v != 0 for real result", v, x);
                    return double.NaN;
                }
                return double.NegativeInfinity;
            }

            if (double.IsInfinity(x))
                return 0;

            if (IsInteger(v))
                return (double)_Bessel.YN(v, x);

            var Y = _Bessel.JY(v, x, false, true).Y;
            return (double)Y;
        }
    }


} // namespaces