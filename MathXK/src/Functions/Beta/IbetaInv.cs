//  Copyright (c) 2013 Rocco Romeo
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

//  Ported to .NET based on the original:
//      Copyright (c) John Maddock 2006, Boost Software License v1.0
//      Copyright (c) Paul A. Bristow 2007, Boost Software License v1.0

#if DEBUG
//#define EXTRA_DEBUG
#endif

using System;
using System.Diagnostics;
using MathXK;
using MathXK.Roots;

namespace MathXK
{

    internal static class _IbetaInv
    {

        private static readonly double[] _Co1 = { -1, -5, 5 };
        private static readonly double[] _Co2 = { 1, 21, -69, 46 };
        private static readonly double[] _Co3 = { 7, -2, 33, -62, 31 };
        private static readonly double[] _Co4 = { 25, -52, -17, 88, -115, 46 };
        private static readonly double[] _Co5 = { 7, 12, -78, 52 };
        private static readonly double[] _Co6 = { -7, 2, 183, -370, 185 };
        private static readonly double[] _Co7 = { -533, 776, -1835, 10240, -13525, 5410 };
        private static readonly double[] _Co8 = { -1579, 3747, -3372, -15821, 45588, -45213, 15071 };
        private static readonly double[] _Co9 = { 449, -1259, -769, 6686, -9260, 3704 };
        private static readonly double[] _Co10 = { 63149, -151557, 140052, -727469, 2239932, -2251437, 750479 };
        private static readonly double[] _Co11 = { 29233, -78755, 105222, 146879, -1602610, 3195183, -2554139, 729754 };
        private static readonly double[] _Co12 = { 1, -13, 13 };
        private static readonly double[] _Co13 = { 1, 21, -69, 46 };



        //
        // See:
        // "Asymptotic Inversion of the Incomplete Beta Function"
        // N.M. Temme
        // Journal of Computation and Applied Mathematics 41 (1992) 145-157.
        // Section 2.
        // Link: http://oai.cwi.nl/oai/asset/2294/2294A.pdf

        static double TemmeMethod1(double a, double b, double z)
        {
            const double r2 = Constants.Sqrt2;
            //
            // get the first approximation for eta from the inverse
            // error function (Eq: 2.9 and 2.10).
            //
            double eta0 = -Math2.ErfcInv(2 * z) / Math.Sqrt(a / 2);

            //
            // calculate powers:
            //
            double B = b - a;
            double B_2 = B * B;
            double B_3 = B_2 * B;
            //
            // Calculate correction terms:
            //

            // See eq following 2.15:
            double w0_0 = -B * r2 / 2;
            double w0_1 = (1 - 2 * B) / 8;
            double w0_2 = -(B * r2 / 48);
            double w0_3 = -1.0 / 192;
            double w0_4 = -B * r2 / 3840;

            // Eq Following 2.17:
            double w1_0 = B * r2 * (3 * B - 2) / 12;
            double w1_1 = (20 * B_2 - 12 * B + 1) / 128;
            double w1_2 = B * r2 * (20 * B - 1) / 960;
            double w1_3 = (16 * B_2 + 30 * B - 15) / 4608;
            double w1_4 = B * r2 * (21 * B + 32) / 53760;
            double w1_5 = (-32 * B_2 + 63) / 368640;
            double w1_6 = -B * r2 * (120 * B + 17) / 25804480;

            // Eq Following 2.17:
            double w2_0 = B * r2 * (-75 * B_2 + 80 * B - 16) / 480;
            double w2_1 = (-1080 * B_3 + 868 * B_2 - 90 * B - 45) / 9216;
            double w2_2 = B * r2 * (-1190 * B_2 + 84 * B + 373) / 53760;
            double w2_3 = (-2240 * B_3 - 2508 * B_2 + 2100 * B - 165) / 368640;

            // Using the coefficents above
            // Compute the polynomials w0 + w1*eta + w2*eta^2 + v3*eta^3 + ...
            double t0 = eta0;
            double t1 = w0_0 + t0 * (w0_1 + t0 * (w0_2 + t0 * (w0_3 + t0 * w0_4)));
            double t2 = w1_0 + t0 * (w1_1 + t0 * (w1_2 + t0 * (w1_3 + t0 * (w1_4 + t0 * (w1_5 + t0 * w1_6)))));
            double t3 = w2_0 + t0 * (w2_1 + t0 * (w2_2 + t0 * w2_3));


            //
            // Bring them together to get a final estimate for eta:
            //
            double ra = 1 / a;
            double eta = t0 + ra * (t1 + ra * (t2 + ra * t3));


            //
            // now we need to convert eta to x, by solving the appropriate
            // quadratic equation:
            //
            double eta_2 = eta * eta;
            double c = -Math.Exp(-eta_2 / 2);
            double x;
            if (eta < 0)
                x = -Math2.Sqrt1pm1(c) / 2; // (1 - Math.Sqrt(1 + c)) / 2;
            else if (eta == 0)
                x = 0.5;
            else
                x = (1 + Math.Sqrt(1 + c)) / 2;


            Debug.Assert(x >= 0 && x <= 1);
            Debug.Assert(eta * (x - 0.5) >= 0);

#if EXTRA_DEBUG
            Debug.WriteLine("Estimating x with Temme method 1: " + x);
#endif


            return x;
        }

        //
        // See: "Asymptotic Inversion of the Incomplete Beta Function"
        // N.M. Temme
        // Journal of Computation and Applied Mathematics 41 (1992) 145-157. Section 3.
        //

        static double TemmeMethod2(double a, double b, double p, double q)
        {

            double z = p;
            double r = a + b;


            //
            // Get first estimate for eta, see Eq 3.9 and 3.10,
            // but note there is a typo in Eq 3.10:
            //
            double eta0 = -Math2.ErfcInv(2 * z) / Math.Sqrt(r / 2);

            double theta = Math.Asin(Math.Sqrt(a / r));
            double s = Math.Sin(theta);
            double c = Math.Cos(theta);
            //
            // Now we need to purturb eta0 to get eta, which we do by
            // evaluating the polynomial in 1/r at the bottom of page 151,
            // to do this we first need the error terms e1, e2 e3
            // which we'll fill into the array "terms".  Since these
            // terms are themselves polynomials, we'll need another
            // array "workspace" to calculate those...
            //
            double[] terms = new double[4];
            terms[0] = eta0;

            double[] workspace = new double[6];
            //
            // some powers of Sin(theta)Cos(theta) that we'll need later:
            //
            double sc = s * c;
            double sc_2 = sc * sc;
            double sc_3 = sc_2 * sc;
            double sc_4 = sc_2 * sc_2;
            double sc_5 = sc_4 * sc;
            double sc_6 = sc_4 * sc_2;
            double sc_7 = sc_6 * sc;


            //
            // Calculate e1 and put it in terms[1], see the middle of page 151:
            //
            workspace[0] = (2 * s * s - 1) / (3 * s * c);
            workspace[1] = -Polynomial.EvalEven(_Co1, s) / (36 * sc_2);
            workspace[2] = Polynomial.EvalEven(_Co2, s) / (1620 * sc_3);
            workspace[3] = -Polynomial.EvalEven(_Co3, s) / (6480 * sc_4);
            workspace[4] = Polynomial.EvalEven(_Co4, s) / (90720 * sc_5);
            terms[1] = Polynomial.Eval(workspace, eta0, 5);
            //
            // Now evaluate e2 and put it in terms[2]:
            //

            workspace[0] = -Polynomial.EvalEven(_Co5, s) / (405 * sc_3);
            workspace[1] = Polynomial.EvalEven(_Co6, s) / (2592 * sc_4);
            workspace[2] = -Polynomial.EvalEven(_Co7, s) / (204120 * sc_5);
            workspace[3] = -Polynomial.EvalEven(_Co8, s) / (2099520 * sc_6);
            terms[2] = Polynomial.Eval(workspace, eta0, 4);
            //
            // And e3, and put it in terms[3]:
            //

            workspace[0] = Polynomial.EvalEven(_Co9, s) / (102060 * sc_5);
            workspace[1] = -Polynomial.EvalEven(_Co10, s) / (20995200 * sc_6);
            workspace[2] = Polynomial.EvalEven(_Co11, s) / (36741600 * sc_7);
            terms[3] = Polynomial.Eval(workspace, eta0, 3);
            //
            // Bring the correction terms together to evaluate eta,
            // this is the last equation on page 151:
            //
            double eta = Polynomial.Eval(terms, 1 / r, 4);
            //
            // Now that we have eta we need to back solve for x,
            // we seek the value of x that gives eta in Eq 3.2.
            // The two methods used are described in section 5.
            //
            // Begin by defining a few variables we'll need later:
            //
            double x;
            double s_2 = s * s;
            double c_2 = c * c;
            double alpha = c / s;
            alpha *= alpha;
            double lu = (-(eta * eta) / (2 * s_2) + Math.Log(s_2) + c_2 * Math.Log(c_2) / s_2);
            //
            // Temme doesn't specify what value to switch on here,
            // but this seems to work pretty well:
            //
            if (Math.Abs(eta) < 0.7) {
                //
                // Small eta use the expansion Temme gives in the second equation
                // of section 5, it's a polynomial in eta:
                //
                workspace[0] = s_2;
                workspace[1] = s * c;
                workspace[2] = (1 - 2 * s_2) / 3;
                workspace[3] = Polynomial.Eval(_Co12, s_2) / (36 * s * c);
                workspace[4] = Polynomial.Eval(_Co13, s_2) / (270 * s_2 * c_2);
                x = Polynomial.Eval(workspace, eta, 5);

#if EXTRA_DEBUG
                Debug.WriteLine("Estimating x with Temme method 2 (small eta): "  + x);
#endif
            } else {
                //
                // If eta is large we need to solve Eq 3.2 more directly,
                // begin by getting an initial approximation for x from
                // the last equation on page 155, this is a polynomial in u:
                //
                double u = Math.Exp(lu);
                workspace[0] = u;
                workspace[1] = alpha;
                workspace[2] = 0;
                workspace[3] = 3 * alpha * (3 * alpha + 1) / 6;
                workspace[4] = 4 * alpha * (4 * alpha + 1) * (4 * alpha + 2) / 24;
                workspace[5] = 5 * alpha * (5 * alpha + 1) * (5 * alpha + 2) * (5 * alpha + 3) / 120;
                x = Polynomial.Eval(workspace, u, 6);
                //
                // At this point we may or may not have the right answer, Eq-3.2 has
                // two solutions for x for any given eta, however the mapping in 3.2
                // is 1:1 with the sign of eta and x-sin^2(theta) being the same.
                // So we can check if we have the right root of 3.2, and if not
                // switch x for 1-x.  This transformation is motivated by the fact
                // that the distribution is *almost* symetric so 1-x will be in the right
                // ball park for the solution:
                //
                if ((x - s_2) * eta < 0)
                    x = 1 - x;

#if EXTRA_DEBUG
                Debug.WriteLine("Estimating x with Temme method 2 (large eta): " + x);
#endif

            }
            //
            // The final step is a few Newton-Raphson iterations to
            // clean up our approximation for x, this is pretty cheap
            // in general, and very cheap compared to an incomplete beta
            // evaluation.  The limits set on x come from the observation
            // that the sign of eta and x-sin^2(theta) are the same.
            //
            double lower, upper;
            if (eta < 0) {
                lower = 0;
                upper = s_2;
            } else {
                lower = s_2;
                upper = 1;
            }
            //
            // If our initial approximation is out of bounds then bisect:
            //
            if ((x < lower) || (x > upper))
                x = (lower + upper) / 2;
            //
            // And iterate:
            //

            return SolveTemmeEqn(-lu, alpha, x, lower, upper);

        }


        //
        // See:
        // "Asymptotic Inversion of the Incomplete Beta Function"
        // N.M. Temme
        // Journal of Computation and Applied Mathematics 41 (1992) 145-157.
        // Section 4.
        //
        static double TemmeMethod3(double a, double b, double p, double q)
        {
            Debug.Assert(a > 0 && b > 0);
            Debug.Assert(p >= 0 && p <= 1);
            Debug.Assert(q >= 0 && q <= 1);

            //
            // Begin by getting an initial approximation for the quantity
            // eta from the dominant part of the incomplete beta:
            //
            double eta0;
            if (p < q)
                eta0 = Math2.GammaQInv(b, p);
            else
                eta0 = Math2.GammaPInv(b, q);
            eta0 /= a;
            //
            // Define the variables and powers we'll need later on:
            //
            double mu = b / a;
            double w = Math.Sqrt(1 + mu);
            double w_2 = 1 + mu; //w * w;
            double w_3 = w_2 * w;
            double w_4 = w_2 * w_2;
            double w_5 = w_4 * w;
            double w_6 = w_4 * w_2;
            double w_7 = w_6 * w;
            double w_8 = w_6 * w_2;
            double w_9 = w_8 * w;
            double w_10 = w_8 * w_2;
            double d = eta0 - mu;
            double d_2 = d * d;
            double d_3 = d_2 * d;
            double d_4 = d_2 * d_2;
            double w1 = w + 1;
            double w1_2 = w1 * w1;
            double w1_3 = w1 * w1_2;
            double w1_4 = w1_2 * w1_2;
            //
            // Now we need to compute the purturbation error terms that
            // convert eta0 to eta, these are all polynomials of polynomials.
            // Probably these should be re-written to use tabulated data
            // (see examples above), but it's less of a win in this case as we
            // need to calculate the individual powers for the denominator terms
            // anyway, so we might as well use them for the numerator-polynomials
            // as well....
            //
            // Refer to p154-p155 for the details of these expansions:
            //
            double e1 = (w + 2) * (w - 1) / (3 * w);
            e1 += (w_3 + 9 * w_2 + 21 * w + 5) * d / (36 * w_2 * w1);
            e1 -= (w_4 - 13 * w_3 + 69 * w_2 + 167 * w + 46) * d_2 / (1620 * w1_2 * w_3);
            e1 -= (7 * w_5 + 21 * w_4 + 70 * w_3 + 26 * w_2 - 93 * w - 31) * d_3 / (6480 * w1_3 * w_4);
            e1 -= (75 * w_6 + 202 * w_5 + 188 * w_4 - 888 * w_3 - 1345 * w_2 + 118 * w + 138) * d_4 / (272160 * w1_4 * w_5);

            double e2 = (28 * w_4 + 131 * w_3 + 402 * w_2 + 581 * w + 208) * (w - 1) / (1620 * w1 * w_3);
            e2 -= (35 * w_6 - 154 * w_5 - 623 * w_4 - 1636 * w_3 - 3983 * w_2 - 3514 * w - 925) * d / (12960 * w1_2 * w_4);
            e2 -= (2132 * w_7 + 7915 * w_6 + 16821 * w_5 + 35066 * w_4 + 87490 * w_3 + 141183 * w_2 + 95993 * w + 21640) * d_2 / (816480 * w_5 * w1_3);
            e2 -= (11053 * w_8 + 53308 * w_7 + 117010 * w_6 + 163924 * w_5 + 116188 * w_4 - 258428 * w_3 - 677042 * w_2 - 481940 * w - 105497) * d_3 / (14696640 * w1_4 * w_6);

            double e3 = -((3592 * w_7 + 8375 * w_6 - 1323 * w_5 - 29198 * w_4 - 89578 * w_3 - 154413 * w_2 - 116063 * w - 29632) * (w - 1)) / (816480 * w_5 * w1_2);
            e3 -= (442043 * w_9 + 2054169 * w_8 + 3803094 * w_7 + 3470754 * w_6 + 2141568 * w_5 - 2393568 * w_4 - 19904934 * w_3 - 34714674 * w_2 - 23128299 * w - 5253353) * d / (146966400 * w_6 * w1_3);
            e3 -= (116932 * w_10 + 819281 * w_9 + 2378172 * w_8 + 4341330 * w_7 + 6806004 * w_6 + 10622748 * w_5 + 18739500 * w_4 + 30651894 * w_3 + 30869976 * w_2 + 15431867 * w + 2919016) * d_2 / (146966400 * w1_4 * w_7);
            //
            // Combine eta0 and the error terms to compute eta (Second eqaution p155):
            //

            // The following is:
            // double eta = eta0 + e1 / a + e2 / (a * a) + e3 / (a * a * a);
            double ra = 1 / a;
            double eta = eta0 + ra * (e1 + ra * (e2 + ra * e3));

            //
            // Now we need to solve Eq 4.2 to obtain x.  For any given value of
            // eta there are two solutions to this equation, and since the distribtion
            // may be very skewed, these are not related by x ~ 1-x we used when
            // implementing section 3 above.  However we know that:
            //
            //  cross < x <= 1       ; iff eta < mu
            //          x == cross   ; iff eta == mu
            //     0 <= x < cross    ; iff eta > mu
            //
            // Where cross == 1 / (1 + mu)
            // Many thanks to Prof Temme for clarifying this point.
            //
            // Therefore we'll just jump straight into Newton iterations
            // to solve Eq 4.2 using these bounds, and simple bisection
            // as the first guess, in practice this converges pretty quickly
            // and we only need a few digits correct anyway:
            //
            if (eta <= 0)
                eta = DoubleLimits.MinNormalValue;
            //double u = eta - mu * Math.Log(eta) + (1 + mu) * Math2.Log1p(mu) - mu;
            //double u = eta - mu * Math.Log(eta/(1+mu)) + Math2.Log1p(mu) - mu;
            double u = eta - mu * Math.Log(eta / (1 + mu)) + Math2.Log1pmx(mu);

            double cross = 1 / (1 + mu); // a/(a+b);

            double lower, upper;
            if (eta < mu) {
                // find the root closest to 1
                lower = cross;
                upper = 1;
            } else {
                // find the root closest to 0
                lower = 0;
                upper = cross;

                // we don't have the resolution to have cross (the maxima) different from 1
                if (upper == 1)
                    upper -= DoubleLimits.MachineEpsilon;
            }

            double x = (lower + upper) / 2;
            if (x == 1)
                return x;

            return SolveTemmeEqn(u, mu, x, lower, upper);

        }

        /// <summary>
        /// Uses the root finder to solve Log(x) + a*Log(1-x) == -t, for x.
        /// Or equivalently, x*(1-x)^a == e^-t, for x
        /// </summary>
        /// <param name="t"></param>
        /// <param name="a"></param>
        /// <param name="guess"></param>
        /// <param name="lower"></param>
        /// <param name="upper"></param>
        /// <returns></returns>
        public static double SolveTemmeEqn(double t, double a, double guess, double lower, double upper)
        {
            // check if our resolution is too low 
            if (lower == 0 && -t < DoubleLimits.MinLogValue)
                return 0;
            if (upper == 1 && -t < a * DoubleLimits.MinLogValue)
                return 1;

            Func<double, ValueTuple<double, double, double>> halleyF = (double x) => {
                const double inf = double.MaxValue / 4;

                if (x == 0)
                    return (-inf, inf, -inf);

                if (x == 1)
                    return (-inf, -inf, -inf);

                // log(x*(1-x)^a)+t
                double y = 1 - x;
                double f = Math.Log(x) + a * Math2.Log1p(-x) + t;
                double f1 = (1 / x) - (a / y);
                double f2 = -(1 / (x * x) + a / (y * y));
                return (f, f1, f2);
            };


            RootResults rr = RootFinder.Halley(halleyF, guess, lower, upper);
            if (rr == null) {
                Policies.ReportRootNotFoundError("Invalid parameter in root solver");
                return double.NaN;
            }
            if (!rr.Success) {
                Policies.ReportRootNotFoundError("Root not found after {0} iterations", rr.Iterations);
                return double.NaN;
            }

            return rr.SolutionX;

        }


        /// <summary>
        /// Returns Pow(a*p*Beta(a,b),1/a) 
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="p"></param>
        /// <param name="q"></param>
        /// <returns></returns>
        /// <remarks>When "a" gets tiny, the function can be very sensitive to variations in the base</remarks>
        static double IbetaInvPrefix(double a, double b, double p, double q)
        {
            Debug.Assert(a > 0 && b > 0);
            Debug.Assert(p >= 0 && p <= 1);
            Debug.Assert(q >= 0 && q <= 1);

            double c = a + b;

            // When x < TinyX, Γ(x) ~= 1/x 
            const double TinyX = DoubleLimits.MachineEpsilon / Constants.EulerMascheroni;
            if (c < TinyX) {
                // c < TinyX => a, b < TinyX 
                // Pow(a*p*Beta(a,b),1/a) = Pow(p*(a+b)/b,1/a)

                // When the inner term of Pow ~== 1, the direct computation is error prone 
                // innerM1 = p*(a+b)/b - 1 = (pa - qb)/b;
                double innerM1 = (p * a - q * b) / b;
                if (Math.Abs(innerM1) < 0.5)
                    return Math.Exp(Math2.Log1p(innerM1) / a);

                // otherwise, compute it normally
                return Math.Pow(p * (c / b), 1 / a);
            }

            // When x < SqrtMachineEpsilon, Γ(x) = 1/x  - γ + O(x)
            if (c < DoubleLimits.RootMachineEpsilon._2) {
                // c < SqrtMachineEpsilon => a, b < SqrtMachineEpsilon

                double gammaC = 1 / c - Constants.EulerMascheroni;
                double innerM1 = (p * a - q * b) / b + (p * a * (Constants.EulerMascheroni * Constants.EulerMascheroni) / gammaC);
                if (Math.Abs(innerM1) < 0.5) {
                    return Math.Exp(Math2.Log1p(innerM1) / a);
                }

                // otherwise, compute it normally
                return Math.Pow(p * (c / b) * (1 - a * Constants.EulerMascheroni) * (1 - b * Constants.EulerMascheroni) / (1 - c * Constants.EulerMascheroni), 1 / a);

            }

            if (a < TinyX) {


                // With a very small a, the direct computation of Pow(a*p*Beta(a,b),1/a)
                // can be very error prone, particulary when the base ~= 1. Small errors in the base may be magnified. 
                // Because of cancellation errors, using Γ(a) = 1/a is inadequate; use Γ(a) = 1/a - γ instead
                // So, using Γ(z)/Γ(z + ϵ) = Γ(z)*(1 + ψ(z)ϵ + O(ϵ^2)), from
                // http://functions.wolfram.com/GammaBetaErf/Gamma/06/01/04/01/
                // p * a * (1/a - γ) * Γ(b)/Γ(b + a) - 1 = p(1 - γa)/(1 + ψ(b)a) - 1 = -(q + γap + ψ(b)*a)/(1 + ψ(b)*a)
                // when b is small ψ(b) ~= 1/b

                Debug.Assert(a < b);
                Debug.Assert(b >= DoubleLimits.RootMachineEpsilon._2);

                double dgb = Math2.Digamma(b);
                double dgba = a * dgb;
                double innerM1 = -(q + Constants.EulerMascheroni * p * a + dgba) / (1 + dgba);
                if (Math.Abs(innerM1) < 0.5) {
                    return Math.Exp(Math2.Log1p(innerM1) / a);
                }

                return Math.Pow(p / (1 + dgba), 1 / a);
            }

            // Lanczos calculation. Compute:
            // Pow(a * p * (agh/cgh)^a * (bgh/cgh)^b * Sqrt((cgh * e)/(agh * bgh) * (SumGScaled(a)*SumGScaled(b)/SumGScaled(c))), 1/a)
            // = (a* p * Sqrt((cgh * e)/(agh * bgh) * (SumGScaled(a)*SumGScaled(b)/SumGScaled(c)))^(1/a) * (agh/cgh) * (bgh/cgh)^(b/a)

            double agh = a + (Lanczos.G - 0.5);
            double bgh = b + (Lanczos.G - 0.5);
            double cgh = c + (Lanczos.G - 0.5);


            double inner = (Lanczos.SeriesExpGScaled(a) / Lanczos.SeriesExpGScaled(c)) * Lanczos.SeriesExpGScaled(b);
            inner *= Math.Sqrt((cgh / agh) * (Math.E / bgh));
            inner *= a;

            double factor;
            const double cutoff = 0.90625; // estimate; TODO: verify
            if (p <= cutoff) {
                // watch for underflows or denorms when computing Pow(p*inner,1/a)
                double pinner = p * inner;
                if ((a > 1) && (pinner == 0 || pinner < DoubleLimits.MinNormalValue))
                    factor = Math.Pow(p, 1 / a) * Math.Pow(inner, 1 / a);
                else
                    factor = Math.Pow(pinner, 1 / a);
            } else {
                // use p = 1-q instead
                factor = Math.Exp((Math2.Log1p(-q) + Math.Log(inner)) / a);
            }

            // because the two power terms are always less than 1,
            // try to maximize result to prevent unneccessary underflows.
            double result = 1;
            if (factor > 1) {
                result = factor;
                factor = 1;
            }

            result *= (agh / cgh);

            // calculate (bgh/cgh)^(b/a) or the equivalent (1 - a/cgh)^(b/a) 
            // choose the latter method if bgh/cgh is sufficiently close to 1

            double t2 = a / cgh;
            if (t2 < 0.5)
                result *= Math.Exp((b / a) * Math2.Log1p(-t2));
            else
                result *= Math.Pow(bgh / cgh, b / a);

            result *= factor;

            return result;

        }


        static double IbetaInvPowerSeriesEstimate(double a, double b, double w)
        {
            Debug.Assert(a > 0 && b > 0);

            //double w = Math.Pow(a*p*Math2.Beta(a,b), 1/a);

            // compute a^n
            double a_2 = a * a;
            double a_3 = a_2 * a;
            double a_4 = a_3 * a;

            // term 1 : sum = w
            double mult = w;
            double term = mult;
            double sum = term;


            // term 2 : sum += w^2 * (b-1)/(a+1)
            double lastTerm = Math.Abs(term);
            mult *= w * ((b - 1) / (a + 1));
            term = mult;
            if (term == 0 || Math.Abs(term) >= lastTerm)
                return sum;
            sum += term;

            // term 3: sum += w^3*(b-1)*(a^2 + 3*a*b - a + 5*b - 4)/(2*(a+1)^2*(a+2))
            lastTerm = Math.Abs(term);
            mult *= w / ((a + 1) * (a + 2));
            term = mult * (a_2 + 3 * a * b - a + 5 * b - 4) / 2;
            if (term == 0 || Math.Abs(term) >= lastTerm)
                return sum;
            sum += term;

            // term 4: w^4*(b-1)*(a^4 + (6*b-1) * a^3 + (b+2)*(8*b-5)*a^2 + (33*b*b - 30 * b + 4)*a + b*(31*b-47)+18)/(3 * (a+1)^3 * (a+2)*(a+3))
            lastTerm = Math.Abs(term);
            mult *= w / (a + 1);
            term = mult * (a_4 + (6 * b - 1) * a_3 + (b + 2) * (8 * b - 5) * a_2 + (33 * b * b - 30 * b + 4) * a + b * (31 * b - 47) + 18) / (3 * (a + 3));
            if (term == 0 || Math.Abs(term) >= lastTerm)
                return sum;
            sum += term;


            return sum;
        }


        // Returns a power series estimate to IbetaInv using the first five terms of:
        // http://functions.wolfram.com/06.23.06.0004.01 
        // Note: this is an asymptotic series
        static double IbetaInvPowerSeriesEstimate(double a, double b, double p, double q)
        {
            Debug.Assert(a > 0 && b > 0);
            Debug.Assert(p >= 0 && p <= 1);
            Debug.Assert(q >= 0 && q <= 1);

            //double w = Math.Pow(a*p*Math2.Beta(a,b), 1/a);
            double w = IbetaInvPrefix(a, b, p, q);
            if (w == 0 || double.IsInfinity(w))
                return w;

            return IbetaInvPowerSeriesEstimate(a, b, w);
        }

#if false
        // When, q is tiny, this code suffers from rounding errors.
        // It is left here for future reference


        double find_ibeta_inv_from_t_dist(double a, double p, double q,out double py )
        {
            Debug.Assert( a > 0 );
            Debug.Assert( p >= 0 && p <= 1 );
            Debug.Assert( q >= 0 && q <= 1 );

            //double u = (p > q) ? double(0.5f - q) / double(2) : double(p / 2);
            //double v = 1 - u; // u < 0.5 so no cancellation error

            double u, v;
            if ( p > q ) {
                u = 0.5-q/2; //(1-q)/2;
                v = 0.5+q/2; //1-u;
            } else {
                u = p/2;
                v = 1-u;
            }

            double df = a * 2;
            bool exact;
            double t = Math2.StudentsTInv(df, u, v, exact);
            double den = (df + t * t);
            double x = df / den;
            py = t * t / den;
            return x;
        }
#endif

        /// <summary>
        /// Use the root finder to solve I{x}(a,b) == p or I{1-x}(a,b) == q for x
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="p"></param>
        /// <param name="q"></param>
        /// <param name="guess"></param>
        /// <param name="min"></param>
        /// <param name="max"></param>
        /// <returns></returns>
        static double FindRoot(double a, double b, double p, double q, double guess, double min, double max)
        {
            // Define our root finding equation for Halley iterations
            Func<double, ValueTuple<double, double, double>> halleyF;
            if (p < q) {
                halleyF = (double x) => {
                    double y = 1 - x;

                    double f = Math2.Ibeta(a, b, x) - p;
                    double f1 = Math2.IbetaDerivative(a, b, x);

                    // f2 = -(1-x)^(b-2) * x^(a-2) * (1 + (x - 1)*a + (b - 2)* x) / Beta(a,b)
                    // f2 = -(1-x)^(b-2) * x^(a-2) * ( (1-a) + (a + b - 2) * x) / Beta(a,b)

                    // adjust x and y on second derivative to avoid overflows
                    if (x == 0) {
                        x = DoubleLimits.MinNormalValue * 64;
                    } else if (y == 0) {
                        y = DoubleLimits.MinNormalValue * 64;
                    }


                    //double f2 = -f1 * ( (1-_a) + (_a + _b -2)*x )/(y * x);
                    double f2 = -f1 * (-y * a + (b - 2) * x + 1) / (y * x);
                    return (f, f1, f2);
                };
            } else {
                halleyF = (double x) => {
                    double y = 1 - x;

                    double f = Math2.Ibetac(a, b, x) - q;
                    double f1 = -Math2.IbetaDerivative(a, b, x);

                    // adjust x and y on second derivative to avoid overflows
                    if (x == 0) {
                        x = DoubleLimits.MinNormalValue * 64;
                    } else if (y == 0) {
                        y = DoubleLimits.MinNormalValue * 64;
                    }

                    // f2 = -(1-x)^(b-2) * x^(a-2) * (1 + (x - 1)*a + (b - 2)* x) /Beta(a,b)
                    // =  -(1-x)^(b-2) * x^(a-2) * ( (1-a) + (a + b -2)*x )

                    //double f2 = -f1 * ( (1-_a) + (_a + _b -2)*x )/(y * x);
                    double f2 = -f1 * (-y * a + (b - 2) * x + 1) / (y * x);

                    return (f, f1, f2);
                };

            }


            //Debug.WriteLine("BetaInv(a={0}, b={1}, p={2}) guess = {3}, range = [{4}, {5}]", a, b, p, guess, min, max);

            RootResults rr = RootFinder.Halley(halleyF, guess, min, max);
            if (rr == null) {
                Policies.ReportRootNotFoundError("Invalid parameter in root solver");
                return double.NaN;
            }

#if EXTRA_DEBUG
            if (!(guess == 0 && rr.SolutionX == 0) && rr.Iterations >= 10)
            {
                Debug.WriteLine("BetaInvPQ(a={0}, b={1}, p={2}, q={3}) guess: {4}, result: {5}, n: {6}", a, b, p, q, guess, rr.SolutionX, rr.Iterations);
            }
#endif

            if (!rr.Success) {
                Policies.ReportRootNotFoundError("Root not found after {0} iterations", rr.Iterations);
                return double.NaN;
            }

            return rr.SolutionX;
        }



        // returns a pair (result, 1-result) to minimize cancellation errors
        public static (double, double) Imp(double a, double b, double p, double q)
        {
            Debug.Assert(a > 0 && b > 0);
            Debug.Assert(p >= 0 && p <= 1);
            Debug.Assert(q >= 0 && q <= 1);

            //
            // Handle trivial cases first:
            //
            if (q == 0)
                return (1.0, 0.0);

            if (p == 0)
                return (0.0, 1.0);

            if ((a == 1) && (b == 1))
                return (p, 1.0 - p);

            const double PowerSeriesLimit = DoubleLimits.RootMachineEpsilon._5; // ~=7.4e-4

            // For reference, the plots with a > b look (very roughly) as follows:
            //
            // I{x}(a,b)        I{x}(b, a)
            //    |       /        |  /-------
            //    |      /         | /
            //    |_____/____      |/__________
            //
            //    I{1-x}(a,b)    I{1-x}(b, a}
            //    |\               |-------\       
            //    | \              |        \
            //    |__\_______      |_________\__
            //


            //
            // Depending upon which approximation method we use, we may end up
            // calculating either x or y initially (where y = 1-x):
            //
            double x, y;

#if false
            // When, p is tiny in the a==0.5 case, this code suffers from rounding errors.
            // It is left here for future reference

            //
            // Handle Student's T case with b = 0.5 gets handled as a special case
            // swap if a == 0.5:
            // 
            if((b == 0.5) && (a >= 0.5)) {
                // We have a Student's T distribution:
                x = find_ibeta_inv_from_t_dist(a, p, q, y);
                return (x, y); 
            }

            if((a == 0.5) && (b > 0.5)) {
                // swap arguments and invert the result
                x = find_ibeta_inv_from_t_dist(b, q, p, y);
                return (y, x); 
            }
#endif

            //
            // The flag invert is set to true if we swap a for b and p for q,
            // in which case the result has to be subtracted from 1:
            //
            bool invert = false;

            // For some of the methods we can put tighter bounds
            // on the result than simply [0,1]:
            //
            double lower = 0;
            double upper = 1;

            //
            // Select calculation method for the initial estimate:
            //

            if (a + b > 5) {
                //
                // When a+b is large then we can use one of Prof Temme's
                // asymptotic expansions, begin by swapping things around
                // so that p < 0.5, we do this to avoid cancellations errors
                // when p is large.
                //
                if (p > 0.5) {
                    Utility.Swap(ref a, ref b);
                    Utility.Swap(ref p, ref q);
                    invert = !invert;
                }
                double minv = Math.Min(a, b);
                double maxv = Math.Max(a, b);
                if ((minv > 5) && (Math.Sqrt(minv) > (maxv - minv))) {
                    //
                    // When a and b differ by a small amount
                    // the curve is quite symmetrical and we can use an error
                    // function to approximate the inverse. This is the cheapest
                    // of the three Temme expantions, and the calculated value
                    // for x will never be much larger than p, so we don't have
                    // to worry about cancellation as long as p is small.
                    //


                    double w = IbetaInvPrefix(a, b, p, q);
                    if (w < PowerSeriesLimit) {
                        x = IbetaInvPowerSeriesEstimate(a, b, w);
                        y = 1 - x;
                    } else {
                        x = TemmeMethod1(a, b, p);
                        y = 1 - x;
                    }
                } else {
                    double r = a + b;
                    double lambda = minv / r;
                    if ((r >= 10) && ((lambda >= 0.2) && (lambda <= 0.8))) {
                        //
                        // The second error function case is the next cheapest
                        // to use, it breaks down when the result is likely to be
                        // very small, if a+b is also small, but we can use a
                        // cheaper expansion there in any case.  As before x won't
                        // be much larger than p, so as long as p is small we should
                        // be free of cancellation error.
                        //

                        double ppa = Math.Pow(p, 1 / a);
                        if ((ppa < 0.0025) && (a + b < 200)) {
                            x = ppa * Math.Pow(a * Math2.Beta(a, b), 1 / a);
                            //x = IbetaInvPrefix(a, b, p, q);
                            //x = IbetaInvPowerSeriesEstimate(a, b, p, q);
                        } else {
                            x = TemmeMethod2(a, b, p, q);
                        }
                        y = 1 - x;
                    } else {
                        //
                        // If we get here then a and b are very different in magnitude
                        // and we need to use the third of Temme's methods which
                        // involves inverting the incomplete gamma.  This is much more
                        // expensive than the other methods.  We also can only use this
                        // method when a > b, which can lead to cancellation errors
                        // if we really want y (as we will when x is close to 1), so
                        // a different expansion is used in that case.
                        //
                        if (a < b) {
                            Utility.Swap(ref a, ref b);
                            Utility.Swap(ref p, ref q);
                            invert = !invert;
                        }

                        // at this point a > b

                        // our implementation of temme method 3 looses accuracy when p is small
                        // so use the power series estimate
                        double w = IbetaInvPrefix(a, b, p, q);
                        if (w < PowerSeriesLimit) {
                            x = IbetaInvPowerSeriesEstimate(a, b, w);
                            y = 1 - x;
                        } else {

                            // Try to compute the easy way first:
                            y = IbetaInvPowerSeriesEstimate(b, a, q, p);
                            x = 1 - y;

                            if (y > 1e-5) {
                                x = TemmeMethod3(a, b, p, q);
                                y = 1 - x;
                            }
                        }

                    }
                }
            } else if ((a < 1) && (b < 1)) {
                //
                // Both a and b less than 1,
                // there is a point of inflection at xs:
                // xsc = 1-xs
                double c = a + b;
                double xs = (1 - a) / (2 - c);
                double xsc = (1 - b) / (2 - c);
                //
                // Now we need to ensure that we start our iteration from the
                // right side of the inflection point:
                //
                double fs = Math2.Ibeta(a, b, xs) - p;
                if (fs == 0) {
                    // The result is at the point of inflection, best just return it:
                    return (invert) ? (xsc, xs) : (xs, xsc);
                }

                if (fs < 0) {
                    Utility.Swap(ref a, ref b);
                    Utility.Swap(ref p, ref q);
                    Utility.Swap(ref xs, ref xsc);
                    invert = !invert;
                }


                //double xg = Math.Pow(a * p * Math2.Beta(a, b), 1/a);
                double w = IbetaInvPrefix(a, b, p, q);
                if (double.IsInfinity(w)) {
                    x = 1;
                    y = 0;
                } else if (w < PowerSeriesLimit) {
                    x = IbetaInvPowerSeriesEstimate(a, b, w);
                    y = 1 - x;
                } else {
                    x = w / (1 + w);
                    y = 1 / (1 + w);
                }

                //
                // And finally we know that our result is below the inflection
                // point, so set an upper limit on our search:
                //

                if (x > xs) {
                    x = xs;
                    y = xsc;
                }
                upper = xs;

            } else if ((a > 1) && (b > 1)) {
                //
                // Small a and b, both greater than 1,
                // there is a point of inflection at xs,
                // and it's complement is xsc, we must always
                // start our iteration from the right side of the
                // point of inflection.
                //
                double c = a + b;
                double xs = (a - 1) / (c - 2);
                double xsc = (b - 1) / (c - 2);
                double ps = Math2.Ibeta(a, b, xs) - p;
                if (ps == 0) {
                    // The result is at the point of inflection, best just return it:
                    return (invert) ? (xsc, xs) : (xs, xsc);
                }

                if (ps < 0) {
                    Utility.Swap(ref a, ref b);
                    Utility.Swap(ref p, ref q);
                    Utility.Swap(ref xs, ref xsc);
                    invert = !invert;
                }
                //
                // Estimate x and y, using expm1 to get a good estimate
                // for y when it's very small:
                //

                x = IbetaInvPowerSeriesEstimate(a, b, p, q);
                const double cutoff = 0.90625;
                if (x < cutoff)
                    y = 1 - x;
                else {
                    y = -Math2.Expm1(Math.Log(p * a * Math2.Beta(a, b)) / a);
                    x = 1 - y;
                }

                //
                // And finally we know that our result is below the inflection
                // point, so set an upper limit on our search:
                //
                if (x > xs) {
                    x = xs;
                    y = xsc;
                }
                upper = xs;
            } else { //if((a <= 1) != (b <= 1))
        
                //
                // If all else fails we get here, only one of a and b
                // is above 1, and a+b is small.  Start by swapping
                // things around so that we have a concave curve with b > a
                // and no points of inflection in [0,1].  As long as we expect
                // x to be small then we can use the simple (and cheap) power
                // term to estimate x, but when we expect x to be large then
                // this greatly underestimates x and leaves us trying to
                // iterate "round the corner" which may take almost forever...
                //
                // We could use Temme's inverse gamma function case in that case,
                // this works really rather well (albeit expensively) even though
                // strictly speaking we're outside it's defined range.
                //
                // However it's expensive to compute, and an alternative approach
                // which models the curve as a distorted quarter circle is much
                // cheaper to compute, and still keeps the number of iterations
                // required down to a reasonable level.  With thanks to Prof Temme
                // for this suggestion.
                //
                if (b < a) {
                    Utility.Swap(ref a, ref b);
                    Utility.Swap(ref p, ref q);
                    invert = !invert;
                }

                // at this point a <= 1 && b >= 1; also b < 5

                double w = IbetaInvPrefix(a, b, p, q);
                if (w < PowerSeriesLimit) {
                    x = IbetaInvPowerSeriesEstimate(a, b, w);
                    y = 1 - x;
                } else if (Math.Pow(p, 1 / a) < 0.5) {
                    x = IbetaInvPowerSeriesEstimate(a, b, w);
                    y = 1 - x;
                } else { //if(pow(q, 1/b) < 0.1)

                    // model a distorted quarter circle:
                    // y = Math.Pow(1 - Math.Pow(p, b * Math2.Beta(a, b)), 1/b);
                    double base_ = b * Math2.Beta(a, b);
                    base_ *= (p > q) ? Math2.Log1p(-q) : Math.Log(p);
                    y = Math.Pow(-Math2.Expm1(base_), 1 / b);
                    x = 1 - y;
                }

            }

            //
            // Now we have a guess for x (and for y) we can set things up for
            // iteration.  If x > 0.5 it pays to swap things round:
            //

            if (x > 0.5) {
                Utility.Swap(ref a, ref b);
                Utility.Swap(ref p, ref q);
                Utility.Swap(ref x, ref y);
                invert = !invert;

                double l = 1 - upper;
                double u = 1 - lower;
                lower = l;
                upper = u;

                // adjust for any minor floating point errors in the complements
                if (x < lower)
                    x = lower;
                else if (x > upper)
                    x = upper;

            }

            x = FindRoot(a, b, p, q, x, lower, upper);

            return (invert) ? (1 - x, x) : (x, 1 - x);
        }

    };


    public static partial class Math2
    {

        /// <summary>
        /// Returns the value "x" that is the solution to: I<sub>x</sub>(a, b) == p. Sets <paramref name="xc"/> to 1-x.
        /// </summary>
        /// <param name="a">Requires a &gt; 0</param>
        /// <param name="b">Requires b &gt; 0</param>
        /// <param name="p">Requires 0 ≤ p ≤ 1</param>
        /// <param name="xc">Returns 1-x without cancellation errors</param>
        public static double IbetaInv(double a, double b, double p, out double xc)
        {
            xc = double.NaN;

            if ((!(a > 0) || double.IsInfinity(a))
            || (!(b > 0) || double.IsInfinity(b))
            || (!(p >= 0 && p <= 1))) {
                Policies.ReportDomainError("IbetaInv(a: {0}, b: {1}, p: {2}, py): Requires finite a,b > 0; p in [0,1]", a, b, p);
                return double.NaN;
            }

            var result = _IbetaInv.Imp(a, b, p, 1 - p);

            xc = result.Item2;
            return result.Item1;
        }

        /// <summary>
        /// Returns the value "x" that is the solution to: 1 - I<sub>x</sub>(a, b) == q. Sets <paramref name="xc"/> to 1-x. 
        /// </summary>
        /// <param name="a">Requires a &gt; 0</param>
        /// <param name="b">Requires b &gt; 0</param>
        /// <param name="q">Requires 0 ≤ q ≤ 1</param>
        /// <param name="xc">Returns 1-x without cancellation errors</param>
        public static double IbetacInv(double a, double b, double q, out double xc)
        {
            xc = double.NaN;

            if ((!(a > 0) || double.IsInfinity(a))
            || (!(b > 0) || double.IsInfinity(b))
            || (!(q >= 0 && q <= 1))) {
                Policies.ReportDomainError("IbetacInv(a: {0}, b: {1}, q: {2}, py): Requires finite a,b > 0; q in [0,1]", a, b, q);
                return double.NaN;
            }

            var result = _IbetaInv.Imp(a, b, 1 - q, q);

            xc = result.Item2;
            return result.Item1;
        }

        /// <summary>
        /// Returns the value "x" that is the solution to: I<sub>x</sub>(a, b) == p
        /// </summary>
        /// <param name="a">Requires a &gt; 0</param>
        /// <param name="b">Requires b &gt; 0</param>
        /// <param name="p">Requires 0 ≤ p ≤ 1</param>
        public static double IbetaInv(double a, double b, double p)
        {
            if ((!(a > 0) || double.IsInfinity(a))
            || (!(b > 0) || double.IsInfinity(b))
            || (!(p >= 0 && p <= 1))) {
                Policies.ReportDomainError("IbetaInv(a: {0}, b: {1}, p: {2}): Requires finite a,b > 0; p in [0,1]", a, b, p);
                return double.NaN;
            }


            return _IbetaInv.Imp(a, b, p, 1 - p).Item1;
        }

        /// <summary>
        /// Returns the value "x" that is the solution to: 1 - I<sub>x</sub>(a, b) == q 
        /// </summary>
        /// <param name="a">Requires a &gt; 0</param>
        /// <param name="b">Requires b &gt; 0</param>
        /// <param name="q">Requires 0 ≤ q ≤ 1</param>
        public static double IbetacInv(double a, double b, double q)
        {
            if ((!(a > 0) || double.IsInfinity(a))
            || (!(b > 0) || double.IsInfinity(b))
            || (!(q >= 0 && q <= 1))) {
                Policies.ReportDomainError("IbetacInv(a: {0}, b: {1}, q: {2}): Requires finite a,b > 0; q in [0,1]", a, b, q);
                return double.NaN;
            }

            return _IbetaInv.Imp(a, b, 1 - q, q).Item1;
        }
    }

} // namespaces






