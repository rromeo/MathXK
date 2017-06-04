//  Copyright (c) 2013 Rocco Romeo
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

//  Ported to .NET based on the original:
//      Copyright (c) John Maddock 2006, Boost Software License, Version 1.0
//      Copyright (c) John Maddock 2007, Boost Software License, Version 1.0


using System;
using System.Diagnostics;
using System.Runtime.InteropServices;

namespace MathXK
{

    internal static class _TDistInvSeries
    {

        private static readonly double[][] _Polynomial = {
            new double[] { 0 },
            new double[] { 1 },   
            new double[] { 0.16666666666666666667, 0.16666666666666666667 },
            new double[] { 0.058333333333333333333, 0.066666666666666666667, 0.0083333333333333333333 },
            new double[] { 0.025198412698412698413, 0.026785714285714285714, 0.0017857142857142857143, 0.00019841269841269841270 },
            new double[] { 0.012039792768959435626, 0.010559964726631393298, -0.0011078042328042328042, 0.00037477954144620811287, 2.7557319223985890653e-6 },
            new double[] { 0.0038370059724226390893, 0.0061039211560044893378, -0.0016095979637646304313, 0.00059458674042007375341, -0.000062705427288760622094, 2.5052108385441718775e-8 },
            new double[] { 0.0032177478835464946576, 0.0010898206731540064873, -0.0012579159844784844785, 0.00069084207973096861986, -0.00016376804137220803887, 0.000015401265401265401265, 1.6059043836821614599e-10 },
            new double[] { 0.0017438262298340009980, 0.000033530976880017885309, -0.00076245135440323932387, 0.00064513046951456342991, -0.00024947258047043099953, 0.000049255746366361445727, -3.9851014346715404916e-6, 7.6471637318198164759e-13 },
            new double[] { 0.00096472747321388644237, -0.00031101086326318780412, -0.00036307660358786885787, 0.00051406605788341121363, -0.00029133414466938067350, 0.000090867107935219902229, -0.000015303004486655377567, 1.0914179173496789432e-6, 2.8114572543455207632e-15 },
            new double[] { 0.00054229262813129686486, -0.00036942667800009661203, -0.00010230378073700412687, 0.00035764655430568632777, -0.00028690924218514613987, 0.00012645437628698076975, -0.000033202652391372058698, 4.8903045291975346210e-6, -3.1239569599829868045e-7, 8.2206352466243297170e-18 }
        };


        public static double[] EvaluatePolynomials(double x)
        {

            double[] c = new double[_Polynomial.Length];
            for (int i = 0; i < _Polynomial.Length; i++)
                c[i] = Polynomial.Eval(_Polynomial[i], x);

            return c;
        }

    }

    internal static class _StudentsTInv
    {

        //
        // The main method used is due to Hill:
        //
        // G. W. Hill, Algorithm 396, Student's t-Quantiles,
        // Communications of the ACM, 13(10): 619-620, Oct., 1970.
        //
        public static double Hill(double ndf, double u)
        {
            Debug.Assert(u <= 0.5);

            if (ndf > 1e20)
                return -Math2.ErfcInv(2 * u) * Constants.Sqrt2;

            double a = 1 / (ndf - 0.5);
            double b = 48 / (a * a);
            double c = ((20700 * a / b - 98) * a - 16) * a + 96.36;
            double d = ((94.5 / (b + c) - 3) / b + 1) * Math.Sqrt(a * Math.PI / 2) * ndf;
            double y = Math.Pow(d * 2 * u, 2 / ndf);

            if (y > (0.05 + a)) {
                //
                // Asymptotic inverse expansion about normal:
                //
                double x = -Math2.ErfcInv(2 * u) * Constants.Sqrt2;
                y = x * x;

                if (ndf < 5)
                    c += 0.3 * (ndf - 4.5) * (x + 0.6);
                c += (((0.05 * d * x - 5) * x - 7) * x - 2) * x + b;
                y = (((((0.4 * y + 6.3) * y + 36) * y + 94.5) / c - y - 3) / b + 1) * x;
                y = Math2.Expm1(a * y * y);
            } else {
                y = ((1 / (((ndf + 6) / (ndf * y) - 0.089 * d - 0.822)
                        * (ndf + 2) * 3) + 0.5 / (ndf + 4)) * y - 1)
                        * (ndf + 1) / (ndf + 2) + 1 / y;
            }
            double q = Math.Sqrt(ndf * y);

            return -q;
        }
        //
        // Tail and body series are due to Shaw:
        //
        // www.mth.kcl.ac.uk/~shaww/web_page/papers/Tdistribution06.pdf
        //
        // Shaw, W.double., 2006, "Sampling Student's double distribution - use of
        // the inverse cumulative distribution function."
        // Journal of Computational Finance, Vol 9 Issue 4, pp 37-73, Summer 2006
        //
        public static double TailSeries(double df, double v)
        {

            // Tail series expansion, see section 6 of Shaw's paper.
            // w is calculated using Eq 60:
            double w = Math2.TgammaDeltaRatio(df / 2, 0.5) * Math.Sqrt(df * Math.PI) * v;

            // define some variables:
            double np2 = df + 2;
            double np4 = df + 4;
            double np6 = df + 6;
            //
            // Calculate the coefficients d(k), these depend only on the
            // number of degrees of freedom df, so at least in theory
            // we could tabulate these for fixed df, see p15 of Shaw:
            //
            double[] d = new double[7];
            d[0] = 1;
            d[1] = -(df + 1) / (2 * np2);
            np2 *= (df + 2);
            d[2] = -df * (df + 1) * (df + 3) / (8 * np2 * np4);
            np2 *= df + 2;
            d[3] = -df * (df + 1) * (df + 5) * (((3 * df) + 7) * df - 2) / (48 * np2 * np4 * np6);
            np2 *= (df + 2);
            np4 *= (df + 4);
            d[4] = -df * (df + 1) * (df + 7) *
               ((((((15 * df) + 154) * df + 465) * df + 286) * df - 336) * df + 64)
               / (384 * np2 * np4 * np6 * (df + 8));
            np2 *= (df + 2);
            d[5] = -df * (df + 1) * (df + 3) * (df + 9)
                     * (((((((35 * df + 452) * df + 1573) * df + 600) * df - 2020) * df) + 928) * df - 128)
                     / (1280 * np2 * np4 * np6 * (df + 8) * (df + 10));
            np2 *= (df + 2);
            np4 *= (df + 4);
            np6 *= (df + 6);
            d[6] = -df * (df + 1) * (df + 11)
                     * ((((((((((((945 * df) + 31506) * df + 425858) * df + 2980236) * df + 11266745) * df + 20675018) * df + 7747124) * df - 22574632) * df - 8565600) * df + 18108416) * df - 7099392) * df + 884736)
                     / (46080 * np2 * np4 * np6 * (df + 8) * (df + 10) * (df + 12));
            //
            // Now bring everthing together to provide the result,
            // this is Eq 62 of Shaw:
            //
            double rn = Math.Sqrt(df);
            double div = Math.Pow(rn * w, 1 / df);
            double power = div * div;
            double result = Polynomial.Eval(d, power);
            result *= rn;
            result /= div;
            return -result;
        }



        public static double BodySeries(double df, double u)
        {

            //
            // Body series for small N:
            //
            // Start with Eq 56 of Shaw:
            //
            double v = Math2.TgammaDeltaRatio(df / 2, 0.5) * Math.Sqrt(df * Math.PI) * (u - 0.5);
            //
            // Workspace for the polynomial coefficients:
            //
            //
            // Figure out what the coefficients are, note these depend
            // only on the degrees of freedom (Eq 57 of Shaw):
            //

            double[] c = _TDistInvSeries.EvaluatePolynomials(1 / df);


            //
            // The result is then a polynomial in v (see Eq 56 of Shaw):
            //
            return Polynomial.EvalOdd(c, v);
        }


        //
        // The following functions below are StudentsTInv_{df}(u, v)
        // They are exact solutions to the Inverse student T for the particualr df
        // u = probability
        // v = 1 -u
        //

        /// <summary>
        /// Student T Inverse with df = 1 
        /// </summary>
        /// <param name="u"></param>
        /// <param name="v"></param>
        /// <returns></returns>
        public static double Df1(double u, double v)
        {
            Debug.Assert(u <= v);

            // df = 1 is the same as the Cauchy distribution, see Shaw Eq 35:
            //if(u == 0.5)
            //  return 0;

            // check endpoints - watch the reverse in sign
            if (u == 0)
                return double.NegativeInfinity;
            if (u == 1) 
                return double.PositiveInfinity;

            return -Math2.CotPI(u);
        }

        /// <summary>
        /// Student T Inverse with df = 2
        /// </summary>
        /// <param name="u"></param>
        /// <param name="v"></param>
        /// <returns></returns>
        public static double Df2(double u, double v)
        {
            Debug.Assert(u <= v);

            // df = 2 has an exact result, see Shaw Eq 36:
            return (2 * u - 1) / Math.Sqrt(2 * u * v);
        }

        /// <summary>
        /// Student T Inverse with df = 4
        /// </summary>
        /// <param name="u"></param>
        /// <param name="v"></param>
        /// <returns></returns>
        public static double Df4(double u, double v)
        {
            // df = 4 has an exact result, see Shaw Eq 38 & 39:
            double alpha = 4 * u * v;
            double root_alpha = Math.Sqrt(alpha);
            double r = 4 * Math.Cos(Math.Acos(root_alpha) / 3) / root_alpha;
            double x = Math.Sqrt(r - 4);
            return (u - 0.5 < 0) ? -x : x;
        }

        /// <summary>
        /// Student T Inverse with df = 6
        /// </summary>
        /// <param name="u"></param>
        /// <param name="v"></param>
        /// <returns></returns>
        public static double Df6(double u, double v)
        {
            const double tolerance = 10 * DoubleLimits.MachineEpsilon;
            const double df = 6;

            //
            // Newton-Raphson iteration of a polynomial case,
            // choice of seed value is taken from Shaw's online
            // supplement:
            //
            double a = 4 * u * v;
            double b = Math2.Cbrt(a);
            const double c = 0.85498797333834849467655443627193;
            double p = 6 * (1 + c * (1 / b - 1));
            double p0;
            do {
                double p2 = p * p;
                double p4 = p2 * p2;
                double p5 = p * p4;
                p0 = p;
                // next term is given by Eq 41:
                p = 2 * (8 * a * p5 - 270 * p2 + 2187) / (5 * (4 * a * p4 - 216 * p - 243));
            } while (!Math2.AreNearRel(p0, p, tolerance));
            //
            // Use Eq 45 to extract the result:
            //
            p = Math.Sqrt(p - df);
            return ((u - 0.5f) < 0) ? -p : p;
        }

#if false
                                                                                                                                                                                                                                            // unrefenced function

        // 
        /// <summary>
        /// Student T Inverse with df = 8
        /// </summary>
        /// <param name="u"></param>
        /// <param name="v"></param>
        /// <returns></returns>
        public static double Df8(double u, double v) 
        {
            const double tolerance = 10 * DoubleLimits.MachineEpsilon;
            const double df = 8;

            //
            // Newton-Raphson iteration of a polynomial case,
            // choice of seed value is taken from Shaw's online
            // supplement:
            //
            const double c8 = 0.85994765706259820318168359251872;
            double a = 4 * (u - u * u); //1 - 4 * (u - 0.5f) * (u - 0.5f);
            double b = Math.Pow(a, 1.0 / 4);
            double p = 8 * (1 + c8 * (1 / b - 1));
            double p0 = p;
            do{
                double p5 = p * p;
                p5 *= p5 * p;
                p0 = p;
                // Next term is given by Eq 42:
                p = 2 * (3 * p + (640 * (160 + p * (24 + p * (p + 4)))) / (-5120 + p * (-2048 - 960 * p + a * p5))) / 7;
            }while( !Math2.AreNearRel(p0, p, tolerance) );
            //
            // Use Eq 45 to extract the result:
            //
            p = Math.Sqrt(p - df);
            return (u - 0.5f < 0) ? -p : p;

        }

        /// <summary>
        /// Student T Inverse with df = 10
        /// </summary>
        /// <param name="u"></param>
        /// <param name="v"></param>
        /// <returns></returns>
        public static double Df10(double u, double v) 
        {
            const double tolerance = 10 * DoubleLimits.MachineEpsilon;
            const double df = 8;

            //
            // Newton-Raphson iteration of a polynomial case,
            // choice of seed value is taken from Shaw's online
            // supplement:
            //
            const double c10 = 0.86781292867813396759105692122285;
            double a = 4 * (u - u * u); //1 - 4 * (u - 0.5f) * (u - 0.5f);
            double b = Math.Pow(a, 1.0 / 5);
            double p = 10 * (1 + c10 * (1 / b - 1));
            double p0;
            do{
                double p6 = p * p;
                p6 *= p6 * p6;
                p0 = p;
                // Next term given by Eq 43:
                p = (8 * p) / 9 + (218750 * (21875 + 4 * p * (625 + p * (75 + 2 * p * (5 + p))))) /
                    (9 * (-68359375 + 8 * p * (-2343750 + p * (-546875 - 175000 * p + 8 * a * p6))));
            } while( !Math2.AreNearRel(p0, p, tolerance) );
            //
            // Use Eq 45 to extract the result:
            //
            p = Math.Sqrt(p - df);
            return (u - 0.5f) < 0 ? -p : p;
        }
#endif

        /// <summary>
        /// Returns the Students T quantile
        /// </summary>
        /// <param name="df">Number of degrees of freedom</param>
        /// <param name="p">Probablity</param>
        /// <param name="q">Probability complement (1-p)</param>
        /// <param name="pexact"></param>
        /// <returns></returns>
        static double StudentsTInv(double df, double p, double q, out bool pexact)
        {

            // if u > v, function is symmetric, invert it:
            if (p > q)
                return -StudentsTInv(df, q, p, out pexact);

            pexact = false;
            double result;

            if (df < 20 && Math2.IsInteger(df)) {
                //
                // we have integer degrees of freedom, 
                // try for the special cases first: 
                //

                switch ((int)df) {
                case 1:
                    pexact = true;
                    return _StudentsTInv.Df1(p, q);
                case 2:
                    pexact = true;
                    return _StudentsTInv.Df2(p, q);
                case 4:
                    pexact = true;
                    return _StudentsTInv.Df4(p, q);
                case 6:
                    // We get numeric overflow in this area:
                    if (p < 1e-150)
                        return _StudentsTInv.Hill(df, p);
                    pexact = true;
                    return _StudentsTInv.Df6(p, q);
#if false
                //
                // These are Shaw's "exact" but iterative solutions
                // for even df, the numerical accuracy of these is
                // rather less than Hill's method, so these are disabled
                // for now, which is a shame because they are reasonably
                // quick to evaluate...
                //
                case 8:
                    return _StudentsTInv.Df8(u,v);
            
                case 10:
                    return _StudentsTInv.Df10(u,v);
            
#endif
                default:
                    break;
                }
            }

            // calculate real

            if (df < 3) {

                // Use a roughly linear scheme to choose between Shaw's tail series and body series:

                double crossover = 0.2742 - df * 0.0242143;
                if (p > crossover)
                    result = _StudentsTInv.BodySeries(df, p);
                else
                    result = _StudentsTInv.TailSeries(df, p);


            } else {

                // Use Hill's method except in the exteme tails
                // where we use Shaw's tail series.
                // The crossover point is roughly exponential in -df:
                double crossover = Math2.Ldexp(1, (int)-Math.Floor(df / 0.654 + 0.5));
                if (p > crossover)
                    result = _StudentsTInv.Hill(df, p);
                else
                    result = _StudentsTInv.TailSeries(df, p);
            }

            return result;
        }



        static double FastQuantileImp2(double df, double p)
        {

            //
            // Need to use inverse incomplete beta to get
            // required precision so not so fast:
            //
            double probability = (p > 0.5) ? 1 - p : p;
            double t, x, y = 0;
            x = Math2.IbetaInv(df / 2, 0.5, 2 * probability, out y);
            t = df * y / x;
            if ( double.IsInfinity(t) )
                t = Math.Sqrt(df)*(Math.Sqrt(y)/Math.Sqrt(x)); // let's do our best to preserve the value 
            else
                t = Math.Sqrt(t);
            //
            // Figure out sign based on the size of p:
            //
            if (p < 0.5)
                t = -t;

            return t;
        }


        public static double FastQuantileImp(double df, double p)
        {

            bool invert = false;
            if ((df < 2) && (!Math2.IsInteger(df)))
                return FastQuantileImp2(df, p);
            if (p > 0.5) {
                p = 1 - p;
                invert = true;
            }
            //
            // Get an estimate of the result:
            //
            bool exact;
            double t = StudentsTInv(df, p, 1 - p, out exact);
            if ((t == 0) || exact)
                return invert ? -t : t; // can't do better!
            //
            // Change variables to inverse incomplete beta:
            //
            double t2 = t * t;
            double xb = df / (df + t2);
            double y = t2 / (df + t2);
            double a = df / 2;
            //
            // t can be so large that x underflows,
            // just return our estimate in that case:
            //
            if (xb == 0)
                return t;
            //
            // Get incomplete beta and it's derivative:
            //
            double f1;
            double f0;
            if (xb < y) {
                f0 = Math2.Ibeta(a, 0.5, xb);
                f1 = Math2.IbetaDerivative(a, 0.5, xb);
            } else {
                f0 = Math2.Ibetac(0.5, a, y);
                f1 = -Math2.IbetaDerivative(0.5, a, y);
            }

            // Get cdf from incomplete beta result:
            double p0 = f0 / 2 - p;
            // Get pdf from derivative:
            double p1 = f1 * Math.Sqrt(y * xb * xb * xb / df);
            //
            // Second derivative divided by p1:
            //
            // yacas gives:
            //
            // In> PrettyForm(Simplify(D(t) (1 + t^2/v) ^ (-(v+1)/2)))
            //
            //  |                        | v + 1     |     |
            //  |                       -| ----- + 1 |     |
            //  |                        |   2       |     |
            // -|             |  2     |                   |
            //  |             | t      |                   |
            //  |             | -- + 1 |                   |
            //  | ( v + 1 ) * | v      |               * t |
            // ---------------------------------------------
            //                       v
            //
            // Which after some manipulation is:
            //
            // -p1 * t * (df + 1) / (t^2 + df)
            //
            double p2 = t * (df + 1) / (t * t + df);
            // Halley step:
            t = Math.Abs(t);
            t += p0 / (p1 + p0 * p2 / 2);
            return !invert ? -t : t;
        }



    }

    public static partial class Math2
    {
        /// <summary>
        /// Returns the Students T quantile
        /// </summary>
        /// <param name="df">Number of degrees of freedom</param>
        /// <param name="p">Probablity</param>
        /// <returns></returns>
        internal static double FastStudentsTQuantile(double df, double p)
        {

            return _StudentsTInv.FastQuantileImp(df, p);
        }

    }

} // namespaces





