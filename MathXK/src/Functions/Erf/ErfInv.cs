//  Copyright (c) 2013 Rocco Romeo
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

//  Ported to .NET based on the original:
//      Copyright (C) John Maddock 2006, Boost Software License, Version 1.0


using System;

namespace MathXK
{

    public static partial class Math2
    {
        /// <summary>
        /// Common implementation for the inverse erf and erfc functions.
        /// This version is for 80-bit long double's and smaller
        /// </summary>
        /// <param name="p"></param>
        /// <param name="q"></param>
        /// <returns></returns>
        private static double ErfInv_Imp(double p, double q)
        {


            double result = 0;

            if (p <= 0.5) {
                //
                // Evaluate inverse erf using the rational approximation:
                //
                // x = p(p+10)(Y+R(p))
                //
                // Where Y is a constant, and R(p) is optimised for a low
                // absolute error compared to |Y|.
                //
                // double: Max error found: 2.001849e-18
                // long double: Max error found: 1.017064e-20
                // Maximum Deviation Found (actual error term at infinite precision) 8.030e-21
                //
                const double Y = 0.0891314744949340820313;
    
                const double p0 = -0.000508781949658280665617;
                const double p1 = -0.00836874819741736770379;
                const double p2 = 0.0334806625409744615033;
                const double p3 = -0.0126926147662974029034;
                const double p4 = -0.0365637971411762664006;
                const double p5 = 0.0219878681111168899165;
                const double p6 = 0.00822687874676915743155;
                const double p7 = -0.00538772965071242932965;
                
                const double q0 = 1;
                const double q1 = -0.970005043303290640362;
                const double q2 = -1.56574558234175846809;
                const double q3 = 1.56221558398423026363;
                const double q4 = 0.662328840472002992063;
                const double q5 = -0.71228902341542847553;
                const double q6 = -0.0527396382340099713954;
                const double q7 = 0.0795283687341571680018;
                const double q8 = -0.00233393759374190016776;
                const double q9 = 0.000886216390456424707504;

                double z = p;
                double P = p0 + z * (p1 + z * (p2 + z * (p3 + z * (p4 + z * (p5 + z * (p6 + z * p7))))));
                double Q = q0 + z * (q1 + z * (q2 + z * (q3 + z * (q4 + z * (q5 + z * (q6 + z * (q7 + z * (q8 + z * q9))))))));

                double g = p * (p + 10);
                result = g * (Y + P/Q);

            } else if (q >= 0.25) {
                //
                // Rational approximation for 0.5 > q >= 0.25
                //
                // x = Math.Sqrt(-2*Math.Log(q)) / (Y + R(q))
                //
                // Where Y is a constant, and R(q) is optimised for a low
                // absolute error compared to Y.
                //
                // double : Max error found: 7.403372e-17
                // long double : Max error found: 6.084616e-20
                // Maximum Deviation Found (error term) 4.811e-20
                //
                const double Y = 2.249481201171875;
                
                const double p0 = -0.202433508355938759655;
                const double p1 = 0.105264680699391713268;
                const double p2 = 8.37050328343119927838;
                const double p3 = 17.6447298408374015486;
                const double p4 = -18.8510648058714251895;
                const double p5 = -44.6382324441786960818;
                const double p6 = 17.445385985570866523;
                const double p7 = 21.1294655448340526258;
                const double p8 = -3.67192254707729348546;
                
                    
                const double q0 = 1;
                const double q1 = 6.24264124854247537712;
                const double q2 = 3.9713437953343869095;
                const double q3 = -28.6608180499800029974;
                const double q4 = -20.1432634680485188801;
                const double q5 = 48.5609213108739935468;
                const double q6 = 10.8268667355460159008;
                const double q7 = -22.6436933413139721736;
                const double q8 = 1.72114765761200282724;

                double z = q - 0.25;
                double P = p0 + z * (p1 + z * (p2 + z * (p3 + z * (p4 + z * (p5 + z * (p6 + z * (p7 + z * p8)))))));
                double Q = q0 + z * (q1 + z * (q2 + z * (q3 + z * (q4 + z * (q5 + z * (q6 + z * (q7 + z * q8)))))));

                double g = Math.Sqrt(-2 * Math.Log(q));
                result = g / (Y + P/Q);

            } else {
                //
                // For q < 0.25 we have a series of rational approximations all
                // of the general form:
                //
                // let: x = Math.Sqrt(-Math.Log(q))
                //
                // Then the result is given by:
                //
                // x(Y+R(x-B))
                //
                // where Y is a constant, B is the lowest value of x for which 
                // the approximation is valid, and R(x-B) is optimised for a low
                // absolute error compared to Y.
                //
                // Note that almost all code will really go through the first
                // or maybe second approximation.  After than we're dealing with very
                // small input values indeed: 80 and 128 bit long double's go all the
                // way down to ~ 1e-5000 so the "tail" is rather long...
                //
                double x = Math.Sqrt(-Math.Log(q));
                if (x < 3) {

                    // Max error found: 1.089051e-20
                    const double Y = 0.807220458984375;
                       
                    const double p0 = -0.131102781679951906451;
                    const double p1 = -0.163794047193317060787;
                    const double p2 = 0.117030156341995252019;
                    const double p3 = 0.387079738972604337464;
                    const double p4 = 0.337785538912035898924;
                    const double p5 = 0.142869534408157156766;
                    const double p6 = 0.0290157910005329060432;
                    const double p7 = 0.00214558995388805277169;
                    const double p8 = -0.679465575181126350155e-6;
                    const double p9 = 0.285225331782217055858e-7;
                    const double p10 = -0.681149956853776992068e-9;
                    
                    const double q0 = 1;
                    const double q1 = 3.46625407242567245975;
                    const double q2 = 5.38168345707006855425;
                    const double q3 = 4.77846592945843778382;
                    const double q4 = 2.59301921623620271374;
                    const double q5 = 0.848854343457902036425;
                    const double q6 = 0.152264338295331783612;
                    const double q7 = 0.01105924229346489121;
                    
                    double z = x - 1.125;
                    double P = p0 + z * (p1 + z * (p2 + z * (p3 + z * (p4 + z * (p5 + z * (p6 + z * (p7 + z * (p8 + z * (p9 + z * p10)))))))));
                    double Q = q0 + z * (q1 + z * (q2 + z * (q3 + z * (q4 + z * (q5 + z * (q6 + z * q7))))));

                    result = x * (Y + P/Q);

                } else if (x < 6) {

                    // Max error found: 8.389174e-21
                    const double Y = 0.93995571136474609375;
                    
                    const double p0 = -0.0350353787183177984712;
                    const double p1 = -0.00222426529213447927281;
                    const double p2 = 0.0185573306514231072324;
                    const double p3 = 0.00950804701325919603619;
                    const double p4 = 0.00187123492819559223345;
                    const double p5 = 0.000157544617424960554631;
                    const double p6 = 0.460469890584317994083e-5;
                    const double p7 = -0.230404776911882601748e-9;
                    const double p8 = 0.266339227425782031962e-11;
                    
                        
                    const double q0 = 1;
                    const double q1 = 1.3653349817554063097;
                    const double q2 = 0.762059164553623404043;
                    const double q3 = 0.220091105764131249824;
                    const double q4 = 0.0341589143670947727934;
                    const double q5 = 0.00263861676657015992959;
                    const double q6 = 0.764675292302794483503e-4;
                    

                    double z = x - 3;
                    double P = p0 + z * (p1 + z * (p2 + z * (p3 + z * (p4 + z * (p5 + z * (p6 + z * (p7 + z * p8)))))));
                    double Q = q0 + z * (q1 + z * (q2 + z * (q3 + z * (q4 + z * (q5 + z * q6)))));

                    result = x * (Y + P/Q);

                } else if (x < 18) {

                    // Max error found: 1.481312e-19
                    const double Y = 0.98362827301025390625;
                    
                    const double p0 = -0.0167431005076633737133;
                    const double p1 = -0.00112951438745580278863;
                    const double p2 = 0.00105628862152492910091;
                    const double p3 = 0.000209386317487588078668;
                    const double p4 = 0.149624783758342370182e-4;
                    const double p5 = 0.449696789927706453732e-6;
                    const double p6 = 0.462596163522878599135e-8;
                    const double p7 = -0.281128735628831791805e-13;
                    const double p8 = 0.99055709973310326855e-16;
                    
                        
                    const double q0 = 1;
                    const double q1 = 0.591429344886417493481;
                    const double q2 = 0.138151865749083321638;
                    const double q3 = 0.0160746087093676504695;
                    const double q4 = 0.000964011807005165528527;
                    const double q5 = 0.275335474764726041141e-4;
                    const double q6 = 0.282243172016108031869e-6;
                    
                    double z = x - 6;
                    double P = p0 + z * (p1 + z * (p2 + z * (p3 + z * (p4 + z * (p5 + z * (p6 + z * (p7 + z * p8)))))));
                    double Q = q0 + z * (q1 + z * (q2 + z * (q3 + z * (q4 + z * (q5 + z * q6)))));

                    result = x * (Y + P/Q);

                } else if (x < 44) {

                    // Max error found: 5.697761e-20
                    const double Y = 0.99714565277099609375;
                        
                    const double p0 = -0.0024978212791898131227;
                    const double p1 = -0.779190719229053954292e-5;
                    const double p2 = 0.254723037413027451751e-4;
                    const double p3 = 0.162397777342510920873e-5;
                    const double p4 = 0.396341011304801168516e-7;
                    const double p5 = 0.411632831190944208473e-9;
                    const double p6 = 0.145596286718675035587e-11;
                    const double p7 = -0.116765012397184275695e-17;
                    
                    const double q0 = 1;
                    const double q1 = 0.207123112214422517181;
                    const double q2 = 0.0169410838120975906478;
                    const double q3 = 0.000690538265622684595676;
                    const double q4 = 0.145007359818232637924e-4;
                    const double q5 = 0.144437756628144157666e-6;
                    const double q6 = 0.509761276599778486139e-9;

                    double z = x - 18;
                    double P = p0 + z * (p1 + z * (p2 + z * (p3 + z * (p4 + z * (p5 + z * (p6 + z * p7))))));
                    double Q = q0 + z * (q1 + z * (q2 + z * (q3 + z * (q4 + z * (q5 + z * q6)))));

                    result = x * (Y + P/Q);

                } else {
                    // Max error found: 1.279746e-20
                    const double Y = 0.99941349029541015625;
                        
                    const double p0 = -0.000539042911019078575891;
                    const double p1 = -0.28398759004727721098e-6;
                    const double p2 = 0.899465114892291446442e-6;
                    const double p3 = 0.229345859265920864296e-7;
                    const double p4 = 0.225561444863500149219e-9;
                    const double p5 = 0.947846627503022684216e-12;
                    const double p6 = 0.135880130108924861008e-14;
                    const double p7 = -0.348890393399948882918e-21;
                    
                        
                    const double q0 = 1;
                    const double q1 = 0.0845746234001899436914;
                    const double q2 = 0.00282092984726264681981;
                    const double q3 = 0.468292921940894236786e-4;
                    const double q4 = 0.399968812193862100054e-6;
                    const double q5 = 0.161809290887904476097e-8;
                    const double q6 = 0.231558608310259605225e-11;

                    double z = x - 44;
                    double P = p0 + z * (p1 + z * (p2 + z * (p3 + z * (p4 + z * (p5 + z * (p6 + z * p7))))));
                    double Q = q0 + z * (q1 + z * (q2 + z * (q3 + z * (q4 + z * (q5 + z * q6)))));

                    result = x * (Y + P/Q);
                }
            }
            return result;
        }

        /// <summary>
        /// Returns the value of the Inverse of the Complemetary Error Function at <paramref name="x"/>
        /// </summary>
        /// <param name="x">The argument. Requires x in [0,2]</param>
        public static double ErfcInv(double x)
        {

            if (!(x >= 0 && x <= 2)) {
                Policies.ReportDomainError("ErfcInv(x: {0}) Requires x in [0,2]", x);
                return double.NaN;
            }

            if (x == 0) 
                return double.PositiveInfinity;
            if (x == 2) 
                return double.NegativeInfinity;
            //
            // Normalise the input, so it's in the range [0,1], we will
            // negate the result if x is outside that range.  This is a simple
            // application of the erfc reflection formula: erfc(-x) = 2 - erfc(x)
            //
            double p, q, s;
            if (x > 1) {
                q = 2 - x;
                p = 1 - q;
                s = -1;
            } else {
                p = 1 - x;
                q = x;
                s = 1;
            }

            //
            // And get the result, negating where required:
            //
            return s * ErfInv_Imp(p, q);
        }

        /// <summary>
        /// Returns the value of the Inverse Error Function at <paramref name="x"/>
        /// </summary>
        /// <param name="x">The argument. Requires x in [-1,1]</param>
        public static double ErfInv(double x)
        {

            if (!(x >= -1 && x <= 1)) {
                Policies.ReportDomainError("ErfInv(x: {0}): Requires x in [-1, 1]", x);
                return double.NaN;
            }
            if (x == 1 || x == -1) 
                return x * double.PositiveInfinity;

            if (x == 0)
                return 0;
            //
            // Normalise the input, so it's in the range [0,1], we will
            // negate the result if x is outside that range.  This is a simple
            // application of the erf reflection formula: erf(-x) = -erf(x)
            //
            double p, q, s;
            if (x < 0) {
                p = -x;
                q = 1 - p;
                s = -1;
            } else {
                p = x;
                q = 1 - x;
                s = 1;
            }
            //
            // And get the result, negating where required:
            //
            return s * ErfInv_Imp(p, q);
        }

    }


} // namespaces


