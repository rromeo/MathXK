//  Copyright (c) 2013 Rocco Romeo
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

//  Ported to .NET based on the original:
//      Copyright (c) 2012 John Maddock, Boost Software License v1.0
// History:
//      Source originally by JM. 
//      RR made changes to use small series and added asymptotic functions.

using System;
using System.Diagnostics;


namespace MathXK
{


    /// <summary>
    /// Airy Asymptotic expansions
    /// </summary>
    /// <remarks>
    /// These are equivalent in accuracy but better on performance. Also, some of the asymptotics make it easier
    /// to avoid early overflow/underflow
    /// </remarks>
    internal class AiryAsym
    {
        // Asymptotic series for positive arguments >=16
        // Mathematica eqns to check "convergence" of asymptotic series: 
        // AiTerm[z_, k_] :=  Pochhammer[1/6, k]*Pochhammer[5/6, k]/k! * (-3/(4*z^(3/2)))^k
        // BiTerm[z_, k_] :=  Pochhammer[1/6, k]*Pochhammer[5/6, k]/k! * (3/(4*z^(3/2)))^k
        // Table[N[AiTerm[16, n], 30], {n, 1, 30}]
        // Table[N[BiTerm[16, n], 30], {n, 1, 30}]


        /// <summary>
        /// Returns the AiryAi for positive x using the asymptotic series
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        /// <seealso href="http://functions.wolfram.com/Bessel-TypeFunctions/AiryAi/06/02/01/01/"/>
        public static double Ai(double x) 
        {
            Debug.Assert(x > 0);

            double z = x;
            double sqrtz = Math.Sqrt(z);
            double z1_5 = z * Math.Sqrt(z);

            double hyp_arg = (-3.0/4)/z1_5;
            double s1 = HypergeometricSeries.Sum2F0(1.0/6, 5.0/6, hyp_arg);

            double prefix = ((Constants.RecipSqrtPI / 2) / Math.Sqrt(sqrtz)) * Math.Exp(-2 * z1_5 / 3);
            double ai = prefix * s1; 

            return ai;
        }

        /// <summary>
        /// Returns the AiryAiPrime for positive x using the asymptotic series
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        /// <seealso href="http://functions.wolfram.com/Bessel-TypeFunctions/AiryAiPrime/06/02/01/01/0005/"/>
        public static double AiPrime(double x)
        {
            Debug.Assert(x > 0);

            double z = x;
            double sqrtz = Math.Sqrt(z);
            double z1_5 = z * sqrtz;

            double hyp_arg = (-3.0 / 4) / z1_5;
            double s1 = HypergeometricSeries.Sum2F0(-1.0 / 6, 7.0 / 6, hyp_arg);

            double aip;
            double prefix = ((Constants.RecipSqrtPI / 2) * Math.Sqrt(sqrtz));
            double expArg = -2 * z1_5 / 3;
            if (expArg > DoubleLimits.MinLogValue) {
                aip = prefix * s1 * Math.Exp(expArg);
            } else if (expArg > 2 * DoubleLimits.MinLogValue) {
                double e = Math.Exp(expArg / 2);
                aip = (prefix * s1 * e) * e;
            } else {
                aip = Math.Exp(expArg - Math.Log(prefix * s1));
            }            

            return -aip;
        }


        /// <summary>
        /// Returns the AiryAi for positive x using the asymptotic series
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        /// <seealso href="http://functions.wolfram.com/Bessel-TypeFunctions/AiryBi/06/02/01/01/0005/"/>
        public static double Bi(double x)
        {
            Debug.Assert(x > 0);

            double z = x;
            double sqrtz = Math.Sqrt(z);
            double z1_5 = z * sqrtz;

            double hyp_arg = (3.0 / 4) / z1_5;
            double s1 = HypergeometricSeries.Sum2F0(1.0 / 6, 5.0 / 6, hyp_arg);

            // watch for overflows 
            double bi;
            double prefix = (Constants.RecipSqrtPI / Math.Sqrt(sqrtz));
            double expArg = 2 * z1_5 / 3;
            if ( expArg < DoubleLimits.MaxLogValue ) {
                bi = prefix * s1 * Math.Exp(expArg); 
            } else if ( expArg < 2*DoubleLimits.MaxLogValue ) {
                double e = Math.Exp(expArg / 2);
                bi = (prefix * s1 * e) * e;
            } else {
                bi = Math.Exp(expArg - Math.Log(prefix * s1));
            }            

            return bi;
        }

        /// <summary>
        /// Returns the AiryBiPrime for positive x using the asymptotic series
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        /// <seealso href="http://functions.wolfram.com/Bessel-TypeFunctions/AiryBiPrime/06/02/01/01/0005/"/>
        public static double BiPrime(double x)
        {
            Debug.Assert(x > 0);

            double z = x;
            double sqrtz = Math.Sqrt(z);
            double z1_5 = z * sqrtz;

            double hyp_arg = (3.0 / 4) / z1_5;
            double s1 = HypergeometricSeries.Sum2F0(-1.0 / 6, 7.0 / 6, hyp_arg);

            double bip;
            double prefix = (Constants.RecipSqrtPI * Math.Sqrt(sqrtz));
            double expArg = 2 * z1_5 / 3;
            if (expArg < DoubleLimits.MaxLogValue) {
                bip = prefix * s1 * Math.Exp(expArg);
            } else if (expArg < 2 * DoubleLimits.MaxLogValue) {
                double e = Math.Exp(expArg / 2);
                bip = (prefix * s1 * e) * e;
            } else {
                bip = Math.Exp(expArg - Math.Log(prefix * s1));
            }

            return bip;
        }


        /// <summary>
        /// Asymptotic AiryAi and AiryBi series for large -x
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        /// <seealso href="http://functions.wolfram.com/Bessel-TypeFunctions/AiryAi/06/02/01/02/"/>
        /// <seealso href="http://functions.wolfram.com/Bessel-TypeFunctions/AiryBi/06/02/01/02/0004/"/>
        public static (double Ai, double Bi) AiBiNeg(double x) 
        {
            Debug.Assert(x < 0);

            double z = -x;
            double sqrtz = Math.Sqrt(z);
            double z1_5 = z * sqrtz;
            double hyp_arg = (-9.0/4) / (z * z * z);
            double s1 = HypergeometricSeries.Sum4F1(1.0 / 12, 5.0 / 12, 7.0 / 12, 11.0 / 12, 1.0 / 2, hyp_arg);
            double s2 = HypergeometricSeries.Sum4F1(7.0 / 12, 11.0 / 12, 13.0 / 12, 17.0 / 12, 3.0 / 2, hyp_arg);

            double zeta = (2 * z1_5) / 3;
            double p = (5.0 / (48 * z1_5));

            // The following is:
            // ai = (1/(sqrt(pi) * z^(1/4))) * (sin(zeta + pi/4) * s1 - p * cos(zeta + pi/4) * s2)
            // bi = (1/(sqrt(pi) * z^(1/4))) * (cos(zeta + pi/4) * s1 + p * sin(zeta + pi/4) * s2)
            // Using trig expansion.

            double sin = Math.Sin(zeta);
            double cos = Math.Cos(zeta);

            double prefix = Constants.RecipSqrt2PI / Math.Sqrt(sqrtz);
            double ai = prefix * ((s1 - p * s2) * cos + (s1 + p * s2) * sin);
            double bi = prefix * ((s1 + p * s2) * cos - (s1 - p * s2) * sin);


            return (ai, bi);
        }

        /// <summary>
        /// Asymptotic AiryAiPrime and AiryBiPrime series for large -x
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        /// <seealso href="http://functions.wolfram.com/Bessel-TypeFunctions/AiryAiPrime/06/02/01/02/"/>
        /// <seealso href="http://functions.wolfram.com/Bessel-TypeFunctions/AiryBiPrime/06/02/01/02/"/>
        public static (double Ai, double Bi) AiBiPrimeNeg(double x)
        {
            Debug.Assert(x < 0);

            double z = -x;
            double sqrtz = Math.Sqrt(z);
            double z1_5 = z * sqrtz;
            double hyp_arg = (-9.0/4) / (z * z * z);
            double s1 = HypergeometricSeries.Sum4F1(-1.0 / 12, 5.0 / 12, 7.0 / 12, 13.0 / 12, 1.0 / 2, hyp_arg);
            double s2 = HypergeometricSeries.Sum4F1(5.0 / 12, 11.0 / 12, 13.0 / 12, 19.0 / 12, 3.0 / 2, hyp_arg);

            double zeta = (2 * z1_5) / 3;
            double p = (7.0 / (48 * z1_5));

            // The following is:
            // ai = -(1/(sqrt(pi) * z^(1/4))) * (cos(zeta + pi/4) * s1 - p * sin(zeta + pi/4) * s2)
            // bi = (1/(sqrt(pi) * z^(1/4))) * (sin(zeta + pi/4) * s1 + p * cos(zeta + pi/4) * s2)
            // Using trig expansion.

            double sin = Math.Sin(zeta);
            double cos = Math.Cos(zeta);

            double prefix = Constants.RecipSqrt2PI * Math.Sqrt(sqrtz);
            double ai = -prefix * ((s1 - p * s2) * cos - (s1 + p * s2) * sin);
            double bi = prefix * ((s1 + p * s2) * cos + (s1 - p * s2) * sin);

            return (ai, bi);
        }

    }

    public static partial class Math2
    {
        /// <summary>
        /// Returns the value of the Airy Ai function at <paramref name="x"/>
        /// </summary>
        /// <param name="x">The function argument</param>
        /// <returns></returns>
        public static double AiryAi(double x)
        {
            if (double.IsNaN(x)) {
                Policies.ReportDomainError("AiryAi(x: {0}): NaN not allowed", x);
                return double.NaN;
            }
            if (double.IsInfinity(x))
                return 0;

            double absX = Math.Abs(x);
            if ( x >= -3 && x <= 2 ) {
                // Use small series. See:
                // http://functions.wolfram.com/Bessel-TypeFunctions/AiryAi/06/01/02/01/0004/
                // http://dlmf.nist.gov/9.4

                double z = x * x * x / 9;
                double s1 = Constants.AiryAi0 * HypergeometricSeries.Sum0F1(2 / 3.0, z);
                double s2 = x * Constants.AiryAiPrime0 * HypergeometricSeries.Sum0F1(4 / 3.0, z);
                return s1 + s2; 
            }

            const double v = 1.0 / 3;
            double zeta = 2 * absX * Math.Sqrt(absX) / 3;
            if (x < 0) {

                if (x < -32)
                    return AiryAsym.AiBiNeg(x).Ai;

                // The following is
                //double j1 = Math2.BesselJ(v, zeta);
                //double j2 = Math2.BesselJ(-v, zeta);
                //double ai = Math.Sqrt(-x) * (j1 + j2) / 3;

                var (J, Y) = _Bessel.JY(v, zeta, true, true);
                double s = 0.5*(J - ((double)Y)/Constants.Sqrt3);
                double ai = Math.Sqrt(-x)* s;

                return ai;
            } else {
                // Use the K relationship: http://dlmf.nist.gov/9.6

                if (x >= 16)
                    return AiryAsym.Ai(x);

                double ai = Math2.BesselK(v, zeta) * Math.Sqrt(x / 3) / Math.PI;
                return ai;
            }
        }


        /// <summary>
        /// Returns the value of the Airy Bi function at <paramref name="x"/>
        /// </summary>
        /// <param name="x">The function argument</param>
        /// <returns></returns>
        public static double AiryBi(double x)
        {
            if (double.IsNaN(x)) {
                Policies.ReportDomainError("AiryBi(x: {0}): NaN not allowed", x);
                return double.NaN;
            }
            if (double.IsInfinity(x)) {
                if (x < 0)
                    return 0;
                return double.PositiveInfinity;
            }


            double absX = Math.Abs(x);
            if (x >= -3 && x <= 3) {
                // Use small series. See:
                // http://functions.wolfram.com/Bessel-TypeFunctions/AiryBi/26/01/01/
                // http://dlmf.nist.gov/9.4

                double z = x * x * x / 9;
                double s1 = Constants.AiryBi0 * HypergeometricSeries.Sum0F1(2 / 3.0, z);
                double s2 = x * Constants.AiryBiPrime0 * HypergeometricSeries.Sum0F1(4 / 3.0, z);

                return s1 + s2;
            }


            const double v = 1.0 / 3;
            double zeta = (absX * Math.Sqrt(absX) * 2) / 3;
            if (x < 0) {

                if (x < -32)
                    return AiryAsym.AiBiNeg(x).Bi;

                // The following is
                //double j1 = Math2.BesselJ(v, zeta);
                //double j2 = Math2.BesselJ(-v, zeta);
                //double bi = Math.Sqrt(-x / 3) * (j2 - j1);

                var (J, Y) = _Bessel.JY(v, zeta, true, true);
                double s = -0.5 * (J/ Constants.Sqrt3 + ((double)Y) );
                double bi = Math.Sqrt(-x) * s;
                return bi;

            } else {

                if (x >= 16)
                    return AiryAsym.Bi(x);

                // The following is
                //double j1 = Math2.BesselI(v, zeta);
                //double j2 = Math2.BesselI(-v, zeta);
                //double bi = Math.Sqrt(x / 3) * (j1 + j2);

                var (I, K) = _Bessel.IK(v, zeta, true, true);
                double s = (2/ Constants.Sqrt3 * I + K/Math.PI);
                var bi = Math.Sqrt(x) * s;

                return bi;

            }
        }


        /// <summary>
        /// Returns the value of the first derivative to the Airy Ai function at <paramref name="x"/>
        /// </summary>
        /// <param name="x">The function argument</param>
        /// <returns></returns>
        public static double AiryAiPrime(double x)
        {
            if (double.IsNaN(x)) {
                Policies.ReportDomainError("AiryAiPrime(x: {0}): NaN not allowed", x);
                return double.NaN;
            }
            if (double.IsInfinity(x)) {
                if (x > 0)
                    return 0;

                Policies.ReportDomainError("AiryAiPrime(x: {0}): Requires x != -Infinity", x);
                return double.NaN;
            }

            double absX = Math.Abs(x);
            if (x >= -3 && x <= 2) {
                // Use small series. See:
                // http://functions.wolfram.com/Bessel-TypeFunctions/AiryAiPrime/26/01/01/
                // http://dlmf.nist.gov/9.4

                double z = x * x * x / 9;
                double s1 = x * x * (Constants.AiryAi0 / 2) * HypergeometricSeries.Sum0F1(5 / 3.0, z);
                double s2 = Constants.AiryAiPrime0 * HypergeometricSeries.Sum0F1(1 / 3.0, z);

                return s1 + s2; 
            }

            const double v = 2.0 / 3;
            double zeta = 2 * absX * Math.Sqrt(absX) / 3;
            if (x < 0) {

                if (x < -32)
                    return AiryAsym.AiBiPrimeNeg(x).Ai;

                // The following is
                //double j1 = Math2.BesselJ(v, zeta);
                //double j2 = Math2.BesselJ(-v, zeta);
                //double aip = -x * (j1 - j2) / 3;

                var (J, Y) = _Bessel.JY(v, zeta, true, true);
                double s = 0.5 * ( J + ((double)Y) / Constants.Sqrt3);
                double aip = -x  * s;
                return aip;

            } else {

                if (x >= 16)
                    return AiryAsym.AiPrime(x);

                double aip = -Math2.BesselK(v, zeta) * x / (Constants.Sqrt3 * Math.PI);
                return aip;
            }
        }


        /// <summary>
        /// Returns the value of first derivative to the Airy Bi function at <paramref name="x"/>
        /// </summary>
        /// <param name="x">AiryBiPrime function argument</param>
        /// <returns></returns>
        public static double AiryBiPrime(double x)
        {
            if (double.IsNaN(x)) {
                Policies.ReportDomainError("AiryBiPrime(x: {0}): NaN not allowed", x);
                return double.NaN;
            }
            if (double.IsInfinity(x)) {
                if (x < 0)
                    return 0;
                return double.PositiveInfinity;
            }

            double absX = Math.Abs(x);
            if (x >= -3 && x <= 3) {
                // Use small series. See:
                // http://functions.wolfram.com/Bessel-TypeFunctions/AiryBiPrime/26/01/01/
                // http://dlmf.nist.gov/9.4

                double z = x * x * x / 9;
                double s1 = x * x * (Constants.AiryBi0 / 2) * HypergeometricSeries.Sum0F1(5 / 3.0, z);
                double s2 = Constants.AiryBiPrime0 * HypergeometricSeries.Sum0F1(1 / 3.0, z);

                return s1 + s2;
            }


            const double v = 2.0 / 3;
            double zeta = 2 * absX * Math.Sqrt(absX) / 3;
            if (x < 0) {

                if (x < -32)
                    return AiryAsym.AiBiPrimeNeg(x).Bi;

                // The following is
                //double j1 = Math2.BesselJ(v, zeta);
                //double j2 = Math2.BesselJ(-v, zeta);
                //double bip = -x * (j1 + j2) / Constants.Sqrt3;
                //return bip;

                var (J, Y) = _Bessel.JY(v, zeta, true, true);
                double s = 0.5 * (J/ Constants.Sqrt3 - ((double)Y));
                double bip = -x * s;
                return bip;

            } else {

                if (x >= 16)
                    return AiryAsym.BiPrime(x);

                // The following is
                //double j1 = Math2.BesselI(v, zeta);
                //double j2 = Math2.BesselI(-v, zeta);
                //double bip = (x / Constants.Sqrt3) * (j1 + j2);

                var (I, K) = _Bessel.IK(v, zeta, true, true);
                double s = (2 / Constants.Sqrt3 * I + K / Math.PI);
                var bip = x * s;
                return bip;

            }
        }



    }



} // namespaces


