//  Copyright (c) 2013 Rocco Romeo
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

//  Ported to .NET based on the original:
//      Copyright (C) John Maddock 2007, Boost Software License, Version 1.0

using System;
using System.Diagnostics;

namespace MathXK
{
    static class _Expint
    {

#if false
        // not used but kept for future reference

        /// <summary>
        /// Ei(x) = γ + log(x) + Σ (x^k/(k*k!)) k={1,∞}
        /// </summary>
        /// <param name="z"></param>
        /// <returns></returns>
        public static double Ei_Series(double z)
        {
            const double DefaultTolerance = 2 * DoubleLimits.MachineEpsilon;

            double term = 1.0;
            double sum = Math.Log(z) + Constants.EulerMascheroni;

            for (int n = 1; n < Policies.MaxSeriesIterations; n++) {

                double prevSum = sum;
                term *= (z/n);
                sum += term/n;

                if (Math.Abs(term) <= Math.Abs(prevSum) * DefaultTolerance)
                    return sum;
            }

            Policies.ReportConvergenceError("Series did not converge after {0} iterations", Policies.MaxSeriesIterations);
            return double.NaN;
        }
#endif

        /// <summary>
        /// Return En(x) using continued fractions:
        /// <para>En(x) = e^-x/(x + n + ContinuedFraction(-k * (k + n - 1), x + n + 2*k) k={1, ∞}</para>
        /// </summary>
        /// <param name="n"></param>
        /// <param name="x">Requires x &gt; 0</param>
        /// <returns></returns>
        /// <seealso href="http://functions.wolfram.com/GammaBetaErf/ExpIntegralE/10/0004/"/>
        public static double En_Fraction(int n, double x)
        {
            double b = n + x;
            double k = 0;
            Func<(double, double)> fracF = () => {
                ++k;
                return (-k * (k + n - 1), b + 2 * k);
            };

            double frac = ContinuedFraction.Eval(b, fracF);

            return Math.Exp(-x) / frac;
        }

        /// <summary>
        /// Implements the Exponential Series for integer n
        /// <para>Sum = ((-z)^(n-1)/(n-1)!) * (ψ(n)-Log(z)) - Σ( (k==n-1) ? 0 : ((-1)^k * z^k)/((k-n+1) * k!)) k={0,∞}</para>
        /// </summary>
        /// <param name="n"></param>
        /// <param name="z"></param>
        /// <returns></returns>
        /// <see href="http://functions.wolfram.com/GammaBetaErf/ExpIntegralE/06/01/04/01/02/0005/"/>
        public static double En_Series(int n, double z)
        {
            Debug.Assert(n > 0);

            const double DefaultTolerance = 2 * DoubleLimits.MachineEpsilon;

            // the sign of the summation is negative
            // so negate these
            double sum = (n > 1) ? -1.0 / (1.0 - n) : 0;
            double term = -1.0;

            int k = 1;
            for (; k < n - 1; k++) {
                term *= -z / k;
                if (term == 0)
                    break;
                sum += term / (k - n + 1);
            }

            sum += Math.Pow(-z, n - 1) * (Math2.Digamma(n) - Math.Log(z)) / Math2.Factorial(n - 1);

            if (term == 0)
                return sum;

            // skip n-1 term
            term *= -z / k;

            for (int i = 0; i < Policies.MaxSeriesIterations; i++) {
                double prevSum = sum;

                k++;
                term *= -z / k;
                double delta = term / (k - n + 1);
                sum += delta;

                if (Math.Abs(delta) <= Math.Abs(prevSum) * DefaultTolerance)
                    return sum;

            }


            Policies.ReportConvergenceError("Series did not converge in {0} iterations", Policies.MaxSeriesIterations);
            return double.NaN;

        }

    }



    public static partial class Math2
    {
        /// <summary>
        /// Returns E1(x) using a rational approximation
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        private static double E1(double x)
        {
            Debug.Assert(x >= 0);

            double result;
            if (x <= 1) {

                // Maximum Deviation Found:                     2.006e-18
                // Expected Error Term:                         2.006e-18
                // Max error found at double precision:         2.760e-17
                const double Y = 0.66373538970947265625;

                const double p0 = 0.0865197248079397976498;
                const double p1 = 0.0320913665303559189999;
                const double p2 = -0.245088216639761496153;
                const double p3 = -0.0368031736257943745142;
                const double p4 = -0.00399167106081113256961;
                const double p5 = -0.000111507792921197858394;

                const double q0 = 1;
                const double q1 = 0.37091387659397013215;
                const double q2 = 0.056770677104207528384;
                const double q3 = 0.00427347600017103698101;
                const double q4 = 0.000131049900798434683324;
                const double q5 = -0.528611029520217142048e-6;

                double z = x;
                double P = p0 + z * (p1 + z * (p2 + z * (p3 + z * (p4 + z * p5))));
                double Q = q0 + z * (q1 + z * (q2 + z * (q3 + z * (q4 + z * q5))));

                result = P / Q + (x - Math.Log(x) - Y);

            } else if (x <= -DoubleLimits.MinLogValue) {

                // Maximum Deviation Found (interpolated):      1.444e-17
                // Max error found at double precision:         3.119e-17

                const double p0 = -0.121013190657725568138e-18;
                const double p1 = -0.999999999999998811143;
                const double p2 = -43.3058660811817946037;
                const double p3 = -724.581482791462469795;
                const double p4 = -6046.8250112711035463;
                const double p5 = -27182.6254466733970467;
                const double p6 = -66598.2652345418633509;
                const double p7 = -86273.1567711649528784;
                const double p8 = -54844.4587226402067411;
                const double p9 = -14751.4895786128450662;
                const double p10 = -1185.45720315201027667;

                const double q0 = 1;
                const double q1 = 45.3058660811801465927;
                const double q2 = 809.193214954550328455;
                const double q3 = 7417.37624454689546708;
                const double q4 = 38129.5594484818471461;
                const double q5 = 113057.05869159631492;
                const double q6 = 192104.047790227984431;
                const double q7 = 180329.498380501819718;
                const double q8 = 86722.3403467334749201;
                const double q9 = 18455.4124737722049515;
                const double q10 = 1229.20784182403048905;
                const double q11 = -0.776491285282330997549;

                double z = 1 / x;
                double P = p0 + z * (p1 + z * (p2 + z * (p3 + z * (p4 + z * (p5 + z * (p6 + z * (p7 + z * (p8 + z * (p9 + z * p10)))))))));
                double Q = q0 + z * (q1 + z * (q2 + z * (q3 + z * (q4 + z * (q5 + z * (q6 + z * (q7 + z * (q8 + z * (q9 + z * (q10 + z * q11))))))))));

                result = Math.Exp(-x) * ((1 + P / Q) / x);

            } else {
                // the asymptotic series from http://dlmf.nist.gov/6.12#E1
                result = Math.Exp(-x) / x;
                if ( result > 0 ) {
                    double z = 1 / x;
                    double P = (1 + z * (-1 + z * (2 + z * (-6 + z * (24 + z * (-120 + z * (720 - 5040 * z)))))));
                    result *= P;
                }

            }

            return result;
        }

        /// <summary>
        /// Returns the value of the Exponential Integral
        /// <para>Ei(x) = ∫ e^t/t dt, t={-∞,x}</para>
        /// </summary>
        public static double Expint(double x)
        {
            if (x < 0)
                return -E1(-x);

            if (double.IsNaN(x)) {
                Policies.ReportDomainError("Expint(x: {0}) NaN not allowed", x);
                return double.NaN;
            }

            if (x == 0)
                return double.NegativeInfinity;

            double result;

            if (x <= 6) {
                // Maximum Deviation Found:                     2.852e-18
                // Expected Error Term:                         2.852e-18
                // Max Error found at double precision =        Poly: 2.636335e-16   Cheb: 4.187027e-16

                const double R = 0.372507410781366634461991866580119133535689497771654051555657435242200120636201854384926049951548942392;
                const double R1 = 1677624236387711.0 / 4503599627370496.0;
                const double R2 = 0.131401834143860282009280387409357165515556574352422001206362e-16;

                const double p0 = 2.98677224343598593013;
                const double p1 = 0.356343618769377415068;
                const double p2 = 0.780836076283730801839;
                const double p3 = 0.114670926327032002811;
                const double p4 = 0.0499434773576515260534;
                const double p5 = 0.00726224593341228159561;
                const double p6 = 0.00115478237227804306827;
                const double p7 = 0.000116419523609765200999;
                const double p8 = 0.798296365679269702435e-5;
                const double p9 = 0.2777056254402008721e-6;

                const double q0 = 1;
                const double q1 = -1.17090412365413911947;
                const double q2 = 0.62215109846016746276;
                const double q3 = -0.195114782069495403315;
                const double q4 = 0.0391523431392967238166;
                const double q5 = -0.00504800158663705747345;
                const double q6 = 0.000389034007436065401822;
                const double q7 = -0.138972589601781706598e-4;

                double z = (x / 3) - 1;
                double P = p0 + z * (p1 + z * (p2 + z * (p3 + z * (p4 + z * (p5 + z * (p6 + z * (p7 + z * (p8 + z * p9))))))));
                double Q = q0 + z * (q1 + z * (q2 + z * (q3 + z * (q4 + z * (q5 + z * (q6 + z * q7))))));

                double xmr = (x - R1) - R2;
                result = (P / Q) * xmr;
                result += (Math.Abs(xmr) < 0.1) ? Math2.Log1p(xmr / R) : Math.Log(x / R);

            } else if (x <= 10) {
                // Maximum Deviation Found:                     6.546e-17
                // Expected Error Term:                         6.546e-17
                // Max Error found at double precision =        Poly: 6.890169e-17   Cheb: 6.772128e-17

                const double Y = 1.158985137939453125;

                const double p0 = 0.00139324086199402804173;
                const double p1 = -0.0349921221823888744966;
                const double p2 = -0.0264095520754134848538;
                const double p3 = -0.00761224003005476438412;
                const double p4 = -0.00247496209592143627977;
                const double p5 = -0.000374885917942100256775;
                const double p6 = -0.554086272024881826253e-4;
                const double p7 = -0.396487648924804510056e-5;

                const double q0 = 1;
                const double q1 = 0.744625566823272107711;
                const double q2 = 0.329061095011767059236;
                const double q3 = 0.100128624977313872323;
                const double q4 = 0.0223851099128506347278;
                const double q5 = 0.00365334190742316650106;
                const double q6 = 0.000402453408512476836472;
                const double q7 = 0.263649630720255691787e-4;

                double z = x / 2 - 4;
                double P = p0 + z * (p1 + z * (p2 + z * (p3 + z * (p4 + z * (p5 + z * (p6 + z * p7))))));
                double Q = q0 + z * (q1 + z * (q2 + z * (q3 + z * (q4 + z * (q5 + z * (q6 + z * q7))))));

                result = Math.Exp(x) / x * (Y + P / Q) + x;

            } else if (x <= 20) {
                // Maximum Deviation Found:                     1.843e-17
                // Expected Error Term:                         -1.842e-17
                // Max Error found at double precision =        Poly: 4.375868e-17   Cheb: 5.860967e-17

                const double Y = 1.0869731903076171875;

                const double p0 = -0.00893891094356945667451;
                const double p1 = -0.0484607730127134045806;
                const double p2 = -0.0652810444222236895772;
                const double p3 = -0.0478447572647309671455;
                const double p4 = -0.0226059218923777094596;
                const double p5 = -0.00720603636917482065907;
                const double p6 = -0.00155941947035972031334;
                const double p7 = -0.000209750022660200888349;
                const double p8 = -0.138652200349182596186e-4;

                const double q0 = 1;
                const double q1 = 1.97017214039061194971;
                const double q2 = 1.86232465043073157508;
                const double q3 = 1.09601437090337519977;
                const double q4 = 0.438873285773088870812;
                const double q5 = 0.122537731979686102756;
                const double q6 = 0.0233458478275769288159;
                const double q7 = 0.00278170769163303669021;
                const double q8 = 0.000159150281166108755531;


                double z = x / 5 - 3;
                double P = p0 + z * (p1 + z * (p2 + z * (p3 + z * (p4 + z * (p5 + z * (p6 + z * (p7 + z * p8)))))));
                double Q = q0 + z * (q1 + z * (q2 + z * (q3 + z * (q4 + z * (q5 + z * (q6 + z * (q7 + z * q8)))))));

                result = Math.Exp(x) / x * (Y + P / Q) + x;

            } else if (x <= 40) {
                // Maximum Deviation Found:                     5.102e-18
                // Expected Error Term:                         5.101e-18
                // Max Error found at double precision =        Poly: 1.441088e-16   Cheb: 1.864792e-16

                const double Y = 1.03937530517578125;

                const double p0 = -0.00356165148914447597995;
                const double p1 = -0.0229930320357982333406;
                const double p2 = -0.0449814350482277917716;
                const double p3 = -0.0453759383048193402336;
                const double p4 = -0.0272050837209380717069;
                const double p5 = -0.00994403059883350813295;
                const double p6 = -0.00207592267812291726961;
                const double p7 = -0.000192178045857733706044;
                const double p8 = -0.113161784705911400295e-9;

                const double q0 = 1;
                const double q1 = 2.84354408840148561131;
                const double q2 = 3.6599610090072393012;
                const double q3 = 2.75088464344293083595;
                const double q4 = 1.2985244073998398643;
                const double q5 = 0.383213198510794507409;
                const double q6 = 0.0651165455496281337831;
                const double q7 = 0.00488071077519227853585;

                double z = x / 10 - 3;
                double P = p0 + z * (p1 + z * (p2 + z * (p3 + z * (p4 + z * (p5 + z * (p6 + z * (p7 + z * p8)))))));
                double Q = q0 + z * (q1 + z * (q2 + z * (q3 + z * (q4 + z * (q5 + z * (q6 + z * q7))))));

                result = Math.Exp(x) / x * (Y + P / Q) + x;

            } else {
                // Max Error found at double precision =        3.381886e-17

                const double Y = 1.013065338134765625;

                const double p0 = -0.0130653381347656243849;
                const double p1 = 0.19029710559486576682;
                const double p2 = 94.7365094537197236011;
                const double p3 = -2516.35323679844256203;
                const double p4 = 18932.0850014925993025;
                const double p5 = -38703.1431362056714134;

                const double q0 = 1;
                const double q1 = 61.9733592849439884145;
                const double q2 = -2354.56211323420194283;
                const double q3 = 22329.1459489893079041;
                const double q4 = -70126.245140396567133;
                const double q5 = 54738.2833147775537106;
                const double q6 = 8297.16296356518409347;


                double z = 1 / x;
                double P = p0 + z * (p1 + z * (p2 + z * (p3 + z * (p4 + z * p5))));
                double Q = q0 + z * (q1 + z * (q2 + z * (q3 + z * (q4 + z * (q5 + z * q6)))));

                result = Y + P / Q;
                if (x < 41)
                    result *= Math.Exp(x) / x;
                else {
                    // Avoid premature overflow if we can:
                    const double exp40 = 2.35385266837019985407899910749034804508871617254555467236651e17;
                    result *= (Math.Exp(x - 40) / x);
                    result *= exp40;
                }
                result += x;
            }

            return result;
        }


        /// <summary>
        /// Returns the Exponential Integral
        /// <para>E<sub>n</sub>(x) = ∫ e^(-x*t)/t^n dt, t={1,∞}</para>
        /// </summary>
        /// <param name="n">Requires n ≥ 0</param>
        /// <param name="x">Requires x ≥ 0</param>
        public static double Expint(int n, double x)
        {
            if (!(n >= 0) || !(x >= 0)) {
                Policies.ReportDomainError("Expint(n: {0}, x: {1}): Requires n,x >= 0", n, x);
                return double.NaN;
            }

            // E_{n}(∞) = 0 
            // see: http://functions.wolfram.com/06.34.03.0014.01 
            if (double.IsInfinity(x))
                return 0;

            if (x == 0) {
                if (n <= 1) {
                    Policies.ReportPoleError("Expint(n: {0}, x: {1}): Requires x > 0 when n <= 1", n, x);
                    return double.NaN;
                }
                return (1.0 / (n - 1));
            }

            if (n == 0)
                return Math.Exp(-x) / x;

            if (n == 1)
                return E1(x);

            if (n < 3) {
                if (x < 0.5)
                    return _Expint.En_Series(n, x);

                return _Expint.En_Fraction(n, x);
            }

            if (x < (n - 2.0) / (n - 1.0) )
                return _Expint.En_Series(n, x);

            return _Expint.En_Fraction(n, x);

        }
    }

} // namespaces



