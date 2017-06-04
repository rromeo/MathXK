//  Copyright (c) 2013 Rocco Romeo
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

//  Ported to .NET based on the original:
//      Copyright (c) 2006 Xiaogang Zhang, Boost Software License v1.0
//  History:
//      XZ wrote original code.
//      RR improved results in for large arguments in J0, J1  


using System;
using System.Diagnostics;
using MathXK.Numerics;

namespace MathXK
{

    static partial class _Bessel
    {

        /// <summary>
        /// Returns J{0}(x)
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        public static double J0(double x)
        {
            if (double.IsInfinity(x))
                return 0;

            // even function
            if (x < 0)
                x = -x;

            // Bessel function of the first kind of order zero
            // x <= 8, minimax rational approximations on root-bracketing intervals
            // x > 8, Hankel asymptotic expansion in Hart, Computer Approximations, 1968

            if (x == 0)
                return 1;

            if (x <= 4) {                      // x in (0, 4]

                const double Zero1 = 2.4048255576957727686e+00;
                const double Zero1a = 10328647123.0 / 4294967296.0;
                const double Zero1b = -2.60059817979848198922953568757550908540328642930009094e-11;

                const double p0 = -4.1298668500990866786e+11;
                const double p1 = 2.7282507878605942706e+10;
                const double p2 = -6.2140700423540120665e+08;
                const double p3 = 6.6302997904833794242e+06;
                const double p4 = -3.6629814655107086448e+04;
                const double p5 = 1.0344222815443188943e+02;
                const double p6 = -1.2117036164593528341e-01;

                const double q0 = 2.3883787996332290397e+12;
                const double q1 = 2.6328198300859648632e+10;
                const double q2 = 1.3985097372263433271e+08;
                const double q3 = 4.5612696224219938200e+05;
                const double q4 = 9.3614022392337710626e+02;
                const double q5 = 1.0;


                double z = x * x;
                double P = p0 + z * (p1 + z * (p2 + z * (p3 + z * (p4 + z * (p5 + z * p6)))));
                double Q = q0 + z * (q1 + z * (q2 + z * (q3 + z * (q4 + z * q5))));

                double factor = (x + Zero1) * ((x - Zero1a) - Zero1b); // (x^2 - Zero1^2)
                return factor * (P / Q);
            }

            if (x <= 8.0) {                 // x in (4, 8]

                const double Zero2 = 5.5200781102863106496e+00;
                const double Zero2a = 23708554955.0 / 4294967296.0; // divided by 2^32
                const double Zero2b = 1.052055426731944484427742522186547878290985375755203814e-11;

                const double p0 = -1.8319397969392084011e+03;
                const double p1 = -1.2254078161378989535e+04;
                const double p2 = -7.2879702464464618998e+03;
                const double p3 = 1.0341910641583726701e+04;
                const double p4 = 1.1725046279757103576e+04;
                const double p5 = 4.4176707025325087628e+03;
                const double p6 = 7.4321196680624245801e+02;
                const double p7 = 4.8591703355916499363e+01;

                const double q0 = -3.5783478026152301072e+05;
                const double q1 = 2.4599102262586308984e+05;
                const double q2 = -8.4055062591169562211e+04;
                const double q3 = 1.8680990008359188352e+04;
                const double q4 = -2.9458766545509337327e+03;
                const double q5 = 3.3307310774649071172e+02;
                const double q6 = -2.5258076240801555057e+01;
                const double q7 = 1.0;


                double z = 1 - (x * x) / 64;
                double P = p0 + z * (p1 + z * (p2 + z * (p3 + z * (p4 + z * (p5 + z * (p6 + z * p7))))));
                double Q = q0 + z * (q1 + z * (q2 + z * (q3 + z * (q4 + z * (q5 + z * (q6 + z * q7))))));

                double factor = (x + Zero2) * ((x - Zero2a) - Zero2b); // (x^2 - Zero2^2)
                return factor * (P / Q);


            } else {

                // x in (8, infinity)

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
                //double z = x - 0.25 * Math.PI;
                //double value = (rc * Math.Cos(z) - y * rs * Math.Sin(z))/ Math.Sqrt(x * (0.5 * Math.PI));
                // using trig addition rules

                double value = (Constants.RecipSqrtPI / Math.Sqrt(x)) * ((rc - y * rs) * Math.Sin(x) + (rc + y * rs) * Math.Cos(x));
                return value;
            }
        }

        /// <summary>
        /// Returns J{1}(x)
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        public static double J1(double x)
        {
            if (double.IsInfinity(x))
                return 0;


            double value;

            // Bessel function of the first kind of order one
            // x <= 8, minimax rational approximations on root-bracketing intervals
            // x > 8, Hankel asymptotic expansion in Hart, Computer Approximations, 1968

            double w = Math.Abs(x);
            if (x == 0)
                return 0;

            if (w <= 4) {                       // w in (0, 4]
                const double Zero1 = 3.8317059702075123156e+00;
                const double Zero1a = 8228525915.0 / 2147483648;
                const double Zero1b = -1.64807473585865746293392334354547257121980712377010100e-11;

                const double p0 = -1.4258509801366645672e+11;
                const double p1 = 6.6781041261492395835e+09;
                const double p2 = -1.1548696764841276794e+08;
                const double p3 = 9.8062904098958257677e+05;
                const double p4 = -4.4615792982775076130e+03;
                const double p5 = 1.0650724020080236441e+01;
                const double p6 = -1.0767857011487300348e-02;

                const double q0 = 4.1868604460820175290e+12;
                const double q1 = 4.2091902282580133541e+10;
                const double q2 = 2.0228375140097033958e+08;
                const double q3 = 5.9117614494174794095e+05;
                const double q4 = 1.0742272239517380498e+03;
                const double q5 = 1.0;


                double z = x * x;
                double P = p0 + z * (p1 + z * (p2 + z * (p3 + z * (p4 + z * (p5 + z * p6)))));
                double Q = q0 + z * (q1 + z * (q2 + z * (q3 + z * (q4 + z * q5))));

                double factor = w * (w + Zero1) * ((w - Zero1a) - Zero1b);
                value = factor * (P/Q);

            } else if (w <= 8) {                  // w in (4, 8]
                const double Zero2 = 7.0155866698156187535e+00;
                const double Zero2a = 30131715309.0 / 4294967296;
                const double Zero2b = 2.599155913244793069527474327631150291131389605537782698e-11;

                const double p0 = -1.7527881995806511112e+16;
                const double p1 = 1.6608531731299018674e+15;
                const double p2 = -3.6658018905416665164e+13;
                const double p3 = 3.5580665670910619166e+11;
                const double p4 = -1.8113931269860667829e+09;
                const double p5 = 5.0793266148011179143e+06;
                const double p6 = -7.5023342220781607561e+03;
                const double p7 = 4.6179191852758252278e+00;

                const double q0 = 1.7253905888447681194e+18;
                const double q1 = 1.7128800897135812012e+16;
                const double q2 = 8.4899346165481429307e+13;
                const double q3 = 2.7622777286244082666e+11;
                const double q4 = 6.4872502899596389593e+08;
                const double q5 = 1.1267125065029138050e+06;
                const double q6 = 1.3886978985861357615e+03;
                const double q7 = 1.0;

                double z = x * x;
                double P = p0 + z * (p1 + z * (p2 + z * (p3 + z * (p4 + z * (p5 + z * (p6 + z * p7))))));
                double Q = q0 + z * (q1 + z * (q2 + z * (q3 + z * (q4 + z * (q5 + z * (q6 + z * q7))))));

                double factor = w * (w + Zero2) * ((w - Zero2a) - Zero2b);
                value = factor * (P / Q);

            } else {                               // w in (8, infinity)

                double y = 8 / w;
                double z = y * y;
                double rc, rs;
                {
                    const double p0 = -4.4357578167941278571e+06;
                    const double p1 = -9.9422465050776411957e+06;
                    const double p2 = -6.6033732483649391093e+06;
                    const double p3 = -1.5235293511811373833e+06;
                    const double p4 = -1.0982405543459346727e+05;
                    const double p5 = -1.6116166443246101165e+03;

                    const double q0 = -4.4357578167941278568e+06;
                    const double q1 = -9.9341243899345856590e+06;
                    const double q2 = -6.5853394797230870728e+06;
                    const double q3 = -1.5118095066341608816e+06;
                    const double q4 = -1.0726385991103820119e+05;
                    const double q5 = -1.4550094401904961825e+03;
                    const double q6 = 1.0;

                    double P = p0 + z * (p1 + z * (p2 + z * (p3 + z * (p4 + z * p5))));
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

                    const double q0 = 7.0871281941028743574e+05;
                    const double q1 = 1.8194580422439972989e+06;
                    const double q2 = 1.4194606696037208929e+06;
                    const double q3 = 4.0029443582266975117e+05;
                    const double q4 = 3.7890229745772202641e+04;
                    const double q5 = 8.6383677696049909675e+02;
                    const double q6 = 1.0;

                    double P = p0 + z * (p1 + z * (p2 + z * (p3 + z * (p4 + z * p5))));
                    double Q = q0 + z * (q1 + z * (q2 + z * (q3 + z * (q4 + z * (q5 + z * q6)))));

                    rs = P / Q;
                }

                // The following code is:
                // double z = w - 0.75 * Math.PI;
                // value = (rc * Math.Cos(z) - y * rs * Math.Sin(z))/ Math.Sqrt((w * 0.5 * Math.PI));
                // using trig addition rules

                value = (Constants.RecipSqrtPI / Math.Sqrt(w)) * ((rc + y * rs) * Math.Sin(w) + (-rc + y * rs) * Math.Cos(w));


            }

            // odd function
            return (x < 0) ? -value : value;
        }


        /// <summary>
        /// Returns J{v}(x) for integer v
        /// </summary>
        /// <param name="v"></param>
        /// <param name="x"></param>
        /// <returns></returns>
        public static double JN(double v, double x)
        {
            Debug.Assert(Math2.IsInteger(v));


            // For integer orders only, we can use two identities:
            // #1: J{-n}(x) = (-1)^n * J{n}(x) 
            // #2: J{n}(-x) = (-1)^n * J{n}(x)

            double sign = 1;
            if (v < 0) {
                v = -v;
                if (Math2.IsOdd(v))
                    sign = -sign;
            }

            if (x < 0) {
                x = -x;
                if (Math2.IsOdd(v))
                    sign = -sign;
            }

            Debug.Assert(v >= 0 && x >= 0);

            if (v > int.MaxValue) {
                Policies.ReportNotImplementedError("BesselJ(v: {0}): Large integer values not yet implemented", v);
                return double.NaN;
            }

            int n = (int)v;

            //
            // Special cases:
            //
            if (n == 0)
                return sign * J0(x);

            if (n == 1)
                return sign * J1(x);

            // n >= 2
            Debug.Assert(n >= 2);

            if (x == 0)
                return 0;

            if (x < 5 || (v > x * x / 4))
                return sign * J_SmallArg(v, x);


            // if v > MinAsymptoticV, use the asymptotics
            const int MinAsymptoticV = 7; // arbitrary constant to reduce iterations
            if (v > MinAsymptoticV) {
                double Jv;
                DoubleX Yv;
                if (JY_TryAsymptotics(v, x, out Jv, out Yv, true, false))
                    return sign * Jv;
            }


            // v < abs(x), forward recurrence stable and usable
            // v >= abs(x), forward recurrence unstable, use Miller's algorithm

            if (v < x) {
                // use forward recurrence
                double prev = J0(x);
                double current = J1(x);
                return sign * Recurrence.ForwardJY(1.0, x, n - 1, current, prev).JYvpn;
            } else {

                // Steed is somewhat more accurate when n gets large
                if (v >= 200) {
                    var (J, Y) = JY_Steed(n, x, true, false);
                    return J;
                }

                // J{n+1}(x) / J{n}(x)
                var (fv, s) = J_CF1(n, x);
                var (Jvmn, Jvmnp1, scale) = Recurrence.BackwardJY_B(n, x, n, s, fv * s);

                // scale the result
                double Jv = (J0(x) / Jvmn);
                Jv = Math2.Ldexp(Jv, -scale);

                return (s * sign) * Jv;      // normalization

            }

        }
    }

} // namespaces


