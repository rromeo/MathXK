//  Copyright (c) 2013 Rocco Romeo
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

//  Ported to .NET based on the original:
//      Copyright (C) John Maddock 2006, Boost Software License, Version 1.0
//  History:
//      JM wrote original C++ code. 
//      RR ported to C# and added Erfcx.

using System;
using System.Diagnostics;

namespace MathXK
{

    public static partial class Math2
    {

        /// <summary>
        /// Computes Erf, Erfc, or Erfcx
        /// </summary>
        /// <param name="x"></param>
        /// <param name="complement">Compute 1-Erf</param>
        /// <param name="scaled">Compute Erfcx</param>
        /// <returns></returns>
        private static double Erf_Imp(double x, bool complement, bool scaled)
        {
            const double MaxErfArg = 5.81;
            const double MaxErfcArg = 27.25;

            if (x < 0) {
                // Transformations:
                // Erf(-x) = -Erf(x);
                // Erfc(-x) = 1+Erf(x)
                // Erfcx(-x) = Math.Exp(x * x) * (1 + Erf(-x))

                if (!complement)
                    return -Erf_Imp(-x, false, false);
                if (!scaled)
                    return 1 + Erf_Imp(-x, false, false);
                return Math.Exp(x * x) * (1 + Erf_Imp(-x, false, false));
            }

            double result = 0;
            if (x <= 0.5) {
                // In this region compute Erf

                if (x < DoubleLimits.RootMachineEpsilon._2) {
                    // _Y0 = 2/Sqrt(PI)
                    const double Y = 1.128379167095512573896158903121545171688101258657997713688171;

                    result = x * Y;
                } else {
                    // Maximum Deviation Found:                     1.561e-17
                    // Expected Error Term:                         1.561e-17
                    // Maximum Relative Change in Control Points:   1.155e-04
                    // Max Error found at double precision =        2.961182e-17

                    const double Y = 1.044948577880859375;

                    const double p0 = 0.0834305892146531832907;
                    const double p1 = -0.338165134459360935041;
                    const double p2 = -0.0509990735146777432841;
                    const double p3 = -0.00772758345802133288487;
                    const double p4 = -0.000322780120964605683831;

                    const double q0 = 1;
                    const double q1 = 0.455004033050794024546;
                    const double q2 = 0.0875222600142252549554;
                    const double q3 = 0.00858571925074406212772;
                    const double q4 = 0.000370900071787748000569;

                    double z = x * x;
                    double P = p0 + z * (p1 + z * (p2 + z * (p3 + z * p4)));
                    double Q = q0 + z * (q1 + z * (q2 + z * (q3 + z * q4)));

                    result = x * (Y + P / Q);
                }

                // Erfc = 1-Erf
                // Erfcx = exp(z*z)*(1-Erf)
                if (complement)
                    result = 1 - result;
                if (scaled)
                    result *= Math.Exp(x * x);

                return result;
            }

            // For double precision,
            // Erf(z) ~= 1, when z >= 5.81 
            // Erfc(x) ~= 0, when x >= 27.25  
            if (!complement) {
                if (x >= MaxErfArg)
                    return 1;
            } else if (!scaled) {
                if (x >= MaxErfcArg)
                    return 0;
            }

            // compute Erfcx = exp(z*z) * erfc(z)

            if (x < 1.5) {

                // Maximum Deviation Found:                     3.702e-17
                // Expected Error Term:                         3.702e-17
                // Maximum Relative Change in Control Points:   2.845e-04
                // Max Error found at double precision =        4.841816e-17
                const double Y = 0.405935764312744140625;

                const double p0 = -0.098090592216281240205;
                const double p1 = 0.178114665841120341155;
                const double p2 = 0.191003695796775433986;
                const double p3 = 0.0888900368967884466578;
                const double p4 = 0.0195049001251218801359;
                const double p5 = 0.00180424538297014223957;

                const double q0 = 1;
                const double q1 = 1.84759070983002217845;
                const double q2 = 1.42628004845511324508;
                const double q3 = 0.578052804889902404909;
                const double q4 = 0.12385097467900864233;
                const double q5 = 0.0113385233577001411017;
                const double q6 = 0.337511472483094676155e-5;

                double z = x - 0.5;

                double P = p0 + z * (p1 + z * (p2 + z * (p3 + z * (p4 + z * p5))));
                double Q = q0 + z * (q1 + z * (q2 + z * (q3 + z * (q4 + z * (q5 + z * q6)))));

                result = (Y + P / Q)/x;

            } else if (x < 2.5) {

                // Max Error found at double precision =        6.599585e-18
                // Maximum Deviation Found:                     3.909e-18
                // Expected Error Term:                         3.909e-18
                // Maximum Relative Change in Control Points:   9.886e-05
                const double Y = 0.50672817230224609375;

                const double p0 = -0.0243500476207698441272;
                const double p1 = 0.0386540375035707201728;
                const double p2 = 0.04394818964209516296;
                const double p3 = 0.0175679436311802092299;
                const double p4 = 0.00323962406290842133584;
                const double p5 = 0.000235839115596880717416;

                const double q0 = 1;
                const double q1 = 1.53991494948552447182;
                const double q2 = 0.982403709157920235114;
                const double q3 = 0.325732924782444448493;
                const double q4 = 0.0563921837420478160373;
                const double q5 = 0.00410369723978904575884;

                double z = x - 1.5;
                double P = p0 + z * (p1 + z * (p2 + z * (p3 + z * (p4 + z * p5))));
                double Q = q0 + z * (q1 + z * (q2 + z * (q3 + z * (q4 + z * q5))));

                result = (Y + P / Q) / x;

            } else if (x < 4.5) {

                // Maximum Deviation Found:                     1.512e-17
                // Expected Error Term:                         1.512e-17
                // Maximum Relative Change in Control Points:   2.222e-04
                // Max Error found at double precision =        2.062515e-17
                const double Y = 0.5405750274658203125;

                const double p0 = 0.00295276716530971662634;
                const double p1 = 0.0137384425896355332126;
                const double p2 = 0.00840807615555585383007;
                const double p3 = 0.00212825620914618649141;
                const double p4 = 0.000250269961544794627958;
                const double p5 = 0.113212406648847561139e-4;

                const double q0 = 1;
                const double q1 = 1.04217814166938418171;
                const double q2 = 0.442597659481563127003;
                const double q3 = 0.0958492726301061423444;
                const double q4 = 0.0105982906484876531489;
                const double q5 = 0.000479411269521714493907;


                double z = x - 3.5;
                double P = p0 + z * (p1 + z * (p2 + z * (p3 + z * (p4 + z * p5))));
                double Q = q0 + z * (q1 + z * (q2 + z * (q3 + z * (q4 + z * q5))));

                result = (Y + P / Q) / x;

            } else if (x < 28) {

                // Max Error found at double precision =        2.997958e-17
                // Maximum Deviation Found:                     2.860e-17
                // Expected Error Term:                         2.859e-17
                // Maximum Relative Change in Control Points:   1.357e-05
                const double Y = 0.5579090118408203125;

                const double p0 = 0.00628057170626964891937;
                const double p1 = 0.0175389834052493308818;
                const double p2 = -0.212652252872804219852;
                const double p3 = -0.687717681153649930619;
                const double p4 = -2.5518551727311523996;
                const double p5 = -3.22729451764143718517;
                const double p6 = -2.8175401114513378771;

                const double q0 = 1;
                const double q1 = 2.79257750980575282228;
                const double q2 = 11.0567237927800161565;
                const double q3 = 15.930646027911794143;
                const double q4 = 22.9367376522880577224;
                const double q5 = 13.5064170191802889145;
                const double q6 = 5.48409182238641741584;


                double z = 1 / x;
                double P = p0 + z * (p1 + z * (p2 + z * (p3 + z * (p4 + z * (p5 + z * p6)))));
                double Q = q0 + z * (q1 + z * (q2 + z * (q3 + z * (q4 + z * (q5 + z * q6)))));

                result = (Y + P / Q) / x;

            } else {

                // otherwise use the asymptotic series
                result = Constants.RecipSqrtPI / x;
                if (x < Constants.RecipSqrt2 / DoubleLimits.RootMachineEpsilon._2)
                    result *= HypergeometricSeries.Sum2F0(1, 0.5, -1 / (x * x));

            }

            // Erfc = exp(-z*z) * Erfcx
            // Erf = 1-Erfc
            if (!scaled)
                result *= Math.Exp(-x * x);
            if (!complement)
                result = 1 - result;

            return result;


        }

        /// <summary>
        /// Returns the Error Function
        /// <para>Erf(x) = 2/Sqrt(π) * ∫ (e^(-t^2) dt, t = {0,x}</para>
        /// </summary>
        /// <param name="x">The argument</param>
        public static double Erf(double x)
        {
            if (double.IsNaN(x)) {
                Policies.ReportDomainError("Erf(x: {0}): NaN not allowed", x);
                return double.NaN;
            }

            return Erf_Imp(x, false, false);
        }

        /// <summary>
        /// Returns the Complementary Error Function
        /// <para>Erfc(x) = 1 - Erf(x)</para>
        /// </summary>
        /// <param name="x">The argument</param>
        public static double Erfc(double x)
        {
            if (double.IsNaN(x)) {
                Policies.ReportDomainError("Erfc(x: {0}): NaN not allowed", x);
                return double.NaN;
            }

            return Erf_Imp(x, true, false);
        }


        /// <summary>
        /// Returns the Scaled Complementary Error Function
        /// <para>Erfcx(x) = e^(x^2)*Erfc(x)</para>
        /// </summary>
        /// <param name="x">The Argument</param>
        /// <returns></returns>
        public static double Erfcx(double x)
        {
            if (double.IsNaN(x)) {
                Policies.ReportDomainError("Erfcx(x: {0}): NaN not allowed", x);
                return double.NaN;
            }

            return Erf_Imp(x, true, true);
        }


    }

} // namespace 





