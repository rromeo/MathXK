//  Copyright (c) 2013 Rocco Romeo
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

using System;
using System.Diagnostics;

namespace MathXK
{

    /// <summary>
    /// Helper functions for Gamma
    /// </summary>
    static class GammaSmall
    {

        /// <summary>
        /// The upper absolute limit on the small gamma series routines ~= 0.0625
        /// </summary>
        public const double UpperLimit = DoubleLimits.RootMachineEpsilon._13;

        ///<summary>
        /// Returns Γ(x) for |x| &lt; UpperLimit using a Taylor series
        /// </summary>  
        public static double Tgamma(double x)
        {
            Debug.Assert(Math.Abs(x) <= UpperLimit && x !=0, "Requires |x| < UpperLimit and x!= 0");

            double absX = Math.Abs(x);
            double tgamma = 1 / x - Constants.EulerMascheroni;
            if (absX >= DoubleLimits.RootMachineEpsilon._2)
                tgamma += Series(x);

            return tgamma;

        }

        ///<summary>
        /// Returns Log(|Γ(x)|) for |x| &lt;= UpperLimit using a Taylor series
        /// </summary>  
        public static double Lgamma(double x, out int sign)
        {
            Debug.Assert(Math.Abs(x) <= UpperLimit && x != 0, "Requires |x| < UpperLimit and x!= 0");

            sign = Math.Sign(x);

            // Log(|Γ(x)|) = -Log(|x|), |x| < eps/Euler
            // Watch for denormalized numbers
            const double TinyX = DoubleLimits.MachineEpsilon / Constants.EulerMascheroni;
            double absX = Math.Abs(x);
            if (absX < TinyX)
                return -Math.Log(absX);

            double tgamma = 1 / x - Constants.EulerMascheroni;
            if (absX >= DoubleLimits.RootMachineEpsilon._2)
                tgamma += Series(x);

            return Math.Log(Math.Abs(tgamma));

        }

        ///<summary>
        /// Returns Γ(x+1)-1 for |x| &lt; UpperLimit using a Taylor series
        /// </summary>  
        public static double Tgamma1pm1(double x)
        {
            // x*Γ(x) = Γ(x+1)
            // So, Γ(x+1)-1 = x(1/x - γ + SmallSeries(x))-1 = x*(-γ+SmallSeries(x))
            Debug.Assert(Math.Abs(x) <= UpperLimit);

            return x * (-Constants.EulerMascheroni + Series(x));
        }

        ///<summary>
        /// Returns Γ(x) - (1/x - γ) for |x| &lt; Limit using a Taylor series
        ///</summary>  
        private static double Series(double x)
        {
            Debug.Assert(Math.Abs(x) <= UpperLimit);

            // see: http://functions.wolfram.com/GammaBetaErf/Gamma/06/01/01/01/

            const double c1 = 0.98905599532797255539539565150063470793918352072821;
            const double c2 = -0.90747907608088628901656016735627511492861144907256;
            const double c3 = 0.98172808683440018733638029402185085036057367972347;
            const double c4 = -0.98199506890314520210470141379137467551742650714720;
            const double c5 = 0.99314911462127619315386725332865849803749075523943;
            const double c6 = -0.99600176044243153397007841966456668673529880955458;
            const double c7 = 0.99810569378312892197857540308836723752396852479018;
            const double c8 = -0.99902526762195486779467805964888808853230396352566;
            const double c9 = 0.99951565607277744106705087759437019443450329799460;
            const double c10 = -0.99975659750860128702584244914060923599695138562883;
            const double c11 = 0.99987827131513327572617164259000321938762910895432;
            const double c12 = -0.99993906420644431683585223136895513185794350282804;

            return x * (c1 + x * (c2 + x * (c3 + x * (c4 + x * (c5 + x * (c6 + x * (c7 + x * (c8 + x * (c9 + x * (c10 + x * (c11 + x * c12)))))))))));

        }

    };


} // namespace
