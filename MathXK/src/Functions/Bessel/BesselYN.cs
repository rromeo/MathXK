//  Copyright (c) 2013 Rocco Romeo
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

//  Ported to .NET based on the original:
//      Copyright (c) 2006 Xiaogang Zhang, Boost Software License v1.0
//  History:
//      XZ wrote original code.
//      RR improved results in for large arguments in Y0, Y1  

using System;
using System.Diagnostics;
using MathXK.Numerics;

namespace MathXK
{

    internal static partial class _Bessel
    {

        /// <summary>
        /// Returns Y{0}(x)
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        public static double Y0(double x)
        {

            if (x < 0) {
                Policies.ReportDomainError("BesselY(v: 0, x: {0}): Requires x >= 0; complex number result not supported", x);
                return double.NaN;
            }
            if (x == 0) {
                Policies.ReportPoleError("BesselY(v: 0, x: {0}): Overflow", x);
                return double.NegativeInfinity;
            }

            if (double.IsInfinity(x))
                return 0;


            double value, factor;

            if (x < DoubleLimits.MachineEpsilon)
                return (2 / Math.PI) * (Math.Log(x / 2) + Constants.EulerMascheroni);

            // Bessel function of the second kind of order zero
            // x <= 8, minimax rational approximations on root-bracketing intervals
            // x > 8, Hankel asymptotic expansion in Hart, Computer Approximations, 1968


            if (x <= 3) {                      // x in [eps, 3]

                const double Zero1 = 8.9357696627916752158e-01;
                const double Zero1a = 3837883847.0 / 4294967296;
                const double Zero1b = -8.66317862382940502854117587748531380699855129307710548e-11;

                const double p0 = 1.0723538782003176831e+11;
                const double p1 = -8.3716255451260504098e+09;
                const double p2 = 2.0422274357376619816e+08;
                const double p3 = -2.1287548474401797963e+06;
                const double p4 = 1.0102532948020907590e+04;
                const double p5 = -1.8402381979244993524e+01;

                const double q0 = 5.8873865738997033405e+11;
                const double q1 = 8.1617187777290363573e+09;
                const double q2 = 5.5662956624278251596e+07;
                const double q3 = 2.3889393209447253406e+05;
                const double q4 = 6.6475986689240190091e+02;
                const double q5 = 1.0;


                double z = x * x;
                double P = p0 + z * (p1 + z * (p2 + z * (p3 + z * (p4 + z * p5))));
                double Q = q0 + z * (q1 + z * (q2 + z * (q3 + z * (q4 + z * q5))));

                double y = (2 / Math.PI) * Math.Log(x / Zero1) * _Bessel.J0(x);
                factor = (x + Zero1) * ((x - Zero1a) - Zero1b);
                value = y + factor * (P / Q);

            } else if (x <= 5.5) {                  // x in (3, 5.5]

                const double Zero2 = 3.9576784193148578684e+00;
                const double Zero2a = 16998099379.0 / 4294967296;
                const double Zero2b = 9.846238317388612698651281418603765563630625507511794841e-12;

                const double p0 = -2.2213976967566192242e+13;
                const double p1 = -5.5107435206722644429e+11;
                const double p2 = 4.3600098638603061642e+10;
                const double p3 = -6.9590439394619619534e+08;
                const double p4 = 4.6905288611678631510e+06;
                const double p5 = -1.4566865832663635920e+04;
                const double p6 = 1.7427031242901594547e+01;

                const double q0 = 4.3386146580707264428e+14;
                const double q1 = 5.4266824419412347550e+12;
                const double q2 = 3.4015103849971240096e+10;
                const double q3 = 1.3960202770986831075e+08;
                const double q4 = 4.0669982352539552018e+05;
                const double q5 = 8.3030857612070288823e+02;
                const double q6 = 1.0;

                double z = x * x;
                double P = p0 + z * (p1 + z * (p2 + z * (p3 + z * (p4 + z * (p5 + z * p6)))));
                double Q = q0 + z * (q1 + z * (q2 + z * (q3 + z * (q4 + z * (q5 + z * q6)))));

                double y = (2 / Math.PI) * Math.Log(x / Zero2) * _Bessel.J0(x);
                factor = (x + Zero2) * ((x - Zero2a) - Zero2b);
                value = y + factor * (P / Q);

            } else if (x <= 8) {                // x in (5.5, 8]

                const double Zero3 = 7.0860510603017726976e+00;
                const double Zero3a = 15217178781.0 / 2147483648;
                const double Zero3b = -5.07017534414388797421475310284896188222355301448323476e-11;

                const double p0 = -8.0728726905150210443e+15;
                const double p1 = 6.7016641869173237784e+14;
                const double p2 = -1.2829912364088687306e+11;
                const double p3 = -1.9363051266772083678e+11;
                const double p4 = 2.1958827170518100757e+09;
                const double p5 = -1.0085539923498211426e+07;
                const double p6 = 2.1363534169313901632e+04;
                const double p7 = -1.7439661319197499338e+01;

                const double q0 = 3.4563724628846457519e+17;
                const double q1 = 3.9272425569640309819e+15;
                const double q2 = 2.2598377924042897629e+13;
                const double q3 = 8.6926121104209825246e+10;
                const double q4 = 2.4727219475672302327e+08;
                const double q5 = 5.3924739209768057030e+05;
                const double q6 = 8.7903362168128450017e+02;
                const double q7 = 1.0;

                double z = x * x;
                double P = p0 + z * (p1 + z * (p2 + z * (p3 + z * (p4 + z * (p5 + z * (p6 + z * p7))))));
                double Q = q0 + z * (q1 + z * (q2 + z * (q3 + z * (q4 + z * (q5 + z * (q6 + z * q7))))));

                double y = (2 / Math.PI) * Math.Log(x / Zero3) * _Bessel.J0(x);
                factor = (x + Zero3) * ((x - Zero3a) - Zero3b);
                value = y + factor * (P / Q);

            } else {                                // x in (8, infinity)

                double y = 8 / x;
                double z = y * y;

                double rc, rs;
                {
                    const double p0 = 2.2779090197304684302e+04;
                    const double p1 = 4.1345386639580765797e+04;
                    const double p2 = 2.1170523380864944322e+04;
                    const double p3 = 3.4806486443249270347e+03;
                    const double p4 = 1.5376201909008354296e+02;
                    const double p5 = 8.8961548424210455236e-01;

                    const double q0 = 2.2779090197304684318e+04;
                    const double q1 = 4.1370412495510416640e+04;
                    const double q2 = 2.1215350561880115730e+04;
                    const double q3 = 3.5028735138235608207e+03;
                    const double q4 = 1.5711159858080893649e+02;
                    const double q5 = 1.0;


                    double P = p0 + z * (p1 + z * (p2 + z * (p3 + z * (p4 + z * p5))));
                    double Q = q0 + z * (q1 + z * (q2 + z * (q3 + z * (q4 + z * q5))));

                    rc = P / Q;

                }

                {
                    const double p0 = -8.9226600200800094098e+01;
                    const double p1 = -1.8591953644342993800e+02;
                    const double p2 = -1.1183429920482737611e+02;
                    const double p3 = -2.2300261666214198472e+01;
                    const double p4 = -1.2441026745835638459e+00;
                    const double p5 = -8.8033303048680751817e-03;

                    const double q0 = 5.7105024128512061905e+03;
                    const double q1 = 1.1951131543434613647e+04;
                    const double q2 = 7.2642780169211018836e+03;
                    const double q3 = 1.4887231232283756582e+03;
                    const double q4 = 9.0593769594993125859e+01;
                    const double q5 = 1.0;

                    double P = p0 + z * (p1 + z * (p2 + z * (p3 + z * (p4 + z * p5))));
                    double Q = q0 + z * (q1 + z * (q2 + z * (q3 + z * (q4 + z * q5))));

                    rs = P / Q;

                }


                // The following code is:
                //double z = x - 0.25f * Math.PI;
                //value = (Constants.RecipSqrtHalfPI / Math.Sqrt(x)) * (rc * Math.Sin(z) + y * rs * Math.Cos(z));

                // Which can be written as:
                // (double, double) sincos = Math2.SinCos(x,-0.25);
                // value = (Constants.RecipSqrtHalfPI / Math.Sqrt(x)) * (rc * sincos.Item1 + y * rs * sincos.Item2);

                // Or, using trig addition rules, simplified to:
                value = (Constants.RecipSqrtPI / Math.Sqrt(x)) * ((rc + y * rs) * Math.Sin(x) + (-rc + y * rs) * Math.Cos(x));


            }

            return value;
        }


        /// <summary>
        /// Returns Y{1}(x)
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        public static double Y1(double x)
        {
            if (x <= 0) {
                Policies.ReportDomainError("BesselY(v: 1, x: {0}): Requires x >= 0; complex number result not supported", x);
                return double.NaN;
            }

            if (double.IsInfinity(x))
                return 0;


            // Bessel function of the second kind of order one
            // x <= 8, minimax rational approximations on root-bracketing intervals
            // x > 8, Hankel asymptotic expansion in Hart, Computer Approximations, 1968

            double value, factor;

            if (x < DoubleLimits.MachineEpsilon)
                return (x / Math.PI) * Math.Log(x / 2) - 2 / (Math.PI * x) - (x / (2 * Math.PI)) * (1 - 2 * Constants.EulerMascheroni);

            if (x <= 4) {                      // x in (0, 4]
                const double Zero1 = 2.1971413260310170351e+00;
                const double Zero1a = 2359162535.0 / 1073741824;
                const double Zero1b = -1.56191001087854667603372694698166849966580723533404956e-12;

                const double p0 = 4.0535726612579544093e+13;
                const double p1 = 5.4708611716525426053e+12;
                const double p2 = -3.7595974497819597599e+11;
                const double p3 = 7.2144548214502560419e+09;
                const double p4 = -5.9157479997408395984e+07;
                const double p5 = 2.2157953222280260820e+05;
                const double p6 = -3.1714424660046133456e+02;

                const double q0 = 3.0737873921079286084e+14;
                const double q1 = 4.1272286200406461981e+12;
                const double q2 = 2.7800352738690585613e+10;
                const double q3 = 1.2250435122182963220e+08;
                const double q4 = 3.8136470753052572164e+05;
                const double q5 = 8.2079908168393867438e+02;
                const double q6 = 1.0;


                double z = x * x;
                double P = p0 + z * (p1 + z * (p2 + z * (p3 + z * (p4 + z * (p5 + z * p6)))));
                double Q = q0 + z * (q1 + z * (q2 + z * (q3 + z * (q4 + z * (q5 + z * q6)))));

                double y = (2 / Math.PI) * Math.Log(x / Zero1) * _Bessel.J1(x);
                factor = (x + Zero1) * ((x - Zero1a) - Zero1b) / x;
                value = y + factor * (P/ Q);

            } else if (x <= 8) {                 // x in (4, 8]

                const double Zero2 = 5.4296810407941351328e+00;
                const double Zero2a = 11660151249.0 / 2147483648;
                const double Zero2b = -1.81486215767496919599158034162425239708541634963014629e-11;

                const double p0 = 1.1514276357909013326e+19;
                const double p1 = -5.6808094574724204577e+18;
                const double p2 = -2.3638408497043134724e+16;
                const double p3 = 4.0686275289804744814e+15;
                const double p4 = -5.9530713129741981618e+13;
                const double p5 = 3.7453673962438488783e+11;
                const double p6 = -1.1957961912070617006e+09;
                const double p7 = 1.9153806858264202986e+06;
                const double p8 = -1.2337180442012953128e+03;

                const double q0 = 5.3321844313316185697e+20;
                const double q1 = 5.6968198822857178911e+18;
                const double q2 = 3.0837179548112881950e+16;
                const double q3 = 1.1187010065856971027e+14;
                const double q4 = 3.0221766852960403645e+11;
                const double q5 = 6.3550318087088919566e+08;
                const double q6 = 1.0453748201934079734e+06;
                const double q7 = 1.2855164849321609336e+03;
                const double q8 = 1.0;


                double z = x * x;
                double P = p0 + z * (p1 + z * (p2 + z * (p3 + z * (p4 + z * (p5 + z * (p6 + z * (p7 + z * p8)))))));
                double Q = q0 + z * (q1 + z * (q2 + z * (q3 + z * (q4 + z * (q5 + z * (q6 + z * (q7 + z * q8)))))));

                double y = (2 / Math.PI) * Math.Log(x / Zero2) * _Bessel.J1(x);
                factor = (x + Zero2) * ((x - Zero2a) - Zero2b) / x;
                value = y + factor * (P/ Q);
            } else {                              // x in (8, infinity)

                double y = 8 / x;
                double z = y * y;
                double rc, rs;
                {
                    const double p0 = -4.4357578167941278571e+06;
                    const double p1 = -9.9422465050776411957e+06;
                    const double p2 = -6.6033732483649391093e+06;
                    const double p3 = -1.5235293511811373833e+06;
                    const double p4 = -1.0982405543459346727e+05;
                    const double p5 = -1.6116166443246101165e+03;
                    const double p6 = 0.0;

                    const double q0 = -4.4357578167941278568e+06;
                    const double q1 = -9.9341243899345856590e+06;
                    const double q2 = -6.5853394797230870728e+06;
                    const double q3 = -1.5118095066341608816e+06;
                    const double q4 = -1.0726385991103820119e+05;
                    const double q5 = -1.4550094401904961825e+03;
                    const double q6 = 1.0;

                    double P = p0 + z * (p1 + z * (p2 + z * (p3 + z * (p4 + z * (p5 + z * p6)))));
                    double Q = q0 + z * (q1 + z * (q2 + z * (q3 + z * (q4 + z * (q5 + z * q6)))));
                    rc = P / Q;
                }
                {

                    const double p0 = 3.3220913409857223519e+04;
                    const double p1 = 8.5145160675335701966e+04;
                    const double p2 = 6.6178836581270835179e+04;
                    const double p3 = 1.8494262873223866797e+04;
                    const double p4 = 1.7063754290207680021e+03;
                    const double p5 = 3.5265133846636032186e+01;
                    const double p6 = 0.0;

                    const double q0 = 7.0871281941028743574e+05;
                    const double q1 = 1.8194580422439972989e+06;
                    const double q2 = 1.4194606696037208929e+06;
                    const double q3 = 4.0029443582266975117e+05;
                    const double q4 = 3.7890229745772202641e+04;
                    const double q5 = 8.6383677696049909675e+02;
                    const double q6 = 1.0;

                    double P = p0 + z * (p1 + z * (p2 + z * (p3 + z * (p4 + z * (p5 + z * p6)))));
                    double Q = q0 + z * (q1 + z * (q2 + z * (q3 + z * (q4 + z * (q5 + z * q6)))));
                    rs = P / Q;

                }

                // The following code is:
                // double z = x - 0.75f * Math.PI;
                // value = (Constants.RecipSqrtHalfPI / Math.Sqrt(x)) * (rc * Math.Sin(z) + y * rs * Math.Cos(z));

                // Which can be written as:
                //(double, double) sincos = Math2.SinCos(x,-0.75);
                //value = (Constants.RecipSqrtHalfPI / Math.Sqrt(x)) * (rc * sincos.Item1 + y * rs * sincos.Item2);

                // Or, using trig addition rules, simplified to:
                value = (Constants.RecipSqrtPI / Math.Sqrt(x)) * ((-rc + y * rs) * Math.Sin(x) - (rc + y * rs) * Math.Cos(x));


            }

            return value;
        }


        /// <summary>
        /// Returns Y{n}(z) for positive integer order and z &lt;= sqrt(eps)
        /// </summary>
        /// <param name="n"></param>
        /// <param name="z"></param>
        /// <returns></returns>
        static DoubleX YN_SmallArg(int n, double z)
        {
            Debug.Assert(n >= 0);
            Debug.Assert(z <= DoubleLimits.RootMachineEpsilon._2);


            //
            // See http://functions.wolfram.com/Bessel-TypeFunctions/BesselY/06/01/04/01/02/
            //
            // Note that when called we assume that x < epsilon and n is a positive integer.
            //

            if (n == 0)
                return new DoubleX((2 / Math.PI) * (Math.Log(z / 2) + Constants.EulerMascheroni));

            if (n == 1) {
                //return (z / Math.PI) * Math.Log(z / 2) - 2 / (Math.PI * z) - (z / (2 * Math.PI)) * (1 - 2 * Constants.EulerMascheroni);
                return new DoubleX((z/2 * (2*Math.Log(z / 2)  - (1 - 2 * Constants.EulerMascheroni)) - 2/z) / Math.PI);
            }

            if (n == 2) {
                double z2 = z * z;
                return new DoubleX((z2/4 * (Math.Log(z / 2) - (3 - 4 * Constants.EulerMascheroni)) - (4/z2) ) / Math.PI);
            }

            // To double check the maximum x value on Mathematica
            // BesselYDiff[n_, x_] := Block[{ a = BesselY[n, x], b = -(Factorial[n - 1]/Pi)*(x/2)^(-n)}, (b - a)/a]
            // Table[Block[{eps = 2^-52}, N[BesselYDiff[n, Sqrt[eps]], 20]], {n, 3, 100}]

            double num = -((Math2.Factorial(n - 1) / Math.PI));
            double den = Math.Pow(z / 2, n);
            if (double.IsInfinity(num) || den == 0) 
                return DoubleX.NegativeInfinity;

            DoubleX result = new DoubleX(num) / new DoubleX(den);

            return result;

        }



        /// <summary>
        /// Computes Y{v}(x), where v is an integer
        /// </summary>
        /// <param name="v">Integer order</param>
        /// <param name="x">Argument. Requires x &gt; 0</param>
        /// <returns></returns>
        public static DoubleX YN(double v, double x)
        {
            if ((x == 0) && (v == 0))
                return double.NegativeInfinity;
            if (x <= 0) {
                Policies.ReportDomainError("BesselY(v: {0}, x: {1}): Complex number result not supported. Requires x >= 0", v, x);
                return double.NaN;
            }


            //
            // Reflection comes first:
            //

            double sign = 1;
            if (v < 0) {
                // Y_{-n}(z) = (-1)^n Y_n(z)
                if (Math2.IsOdd(v))
                    sign = -1;
                v = -v;
            }

            Debug.Assert(v >= 0 && x >= 0);

            if (v > int.MaxValue) {
                Policies.ReportNotImplementedError("BesselY(v: {0}, x: {1}): Large integer values not yet implemented", v, x);
                return double.NaN;
            }

            int n = (int)v;

            if (n == 0)
                return Y0(x);

            if (n == 1)
                return sign * Y1(x);

            if (x < DoubleLimits.RootMachineEpsilon._2) {
                DoubleX smallArgValue = YN_SmallArg(n, x);
                return smallArgValue * sign;
            }

            // if v > MinAsymptoticV, use the asymptotics
            const int MinAsymptoticV = 7; // arbitrary constant to reduce iterations
            if (v > MinAsymptoticV) {
                double Jv;
                DoubleX Yv;
                if (JY_TryAsymptotics(v, x, out Jv, out Yv, false, true))
                    return sign * Yv;
            }


            // forward recurrence always OK (though unstable)
            double prev = Y0(x);
            double current = Y1(x);

            var (Yvpn, Yvpnm1, YScale) = Recurrence.ForwardJY_B(1.0, x, n - 1, current, prev);

            return DoubleX.Ldexp(sign * Yvpn, YScale);
        }
    }


} // namespaces


