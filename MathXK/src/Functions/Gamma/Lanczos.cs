//  Copyright (c) 2013 Rocco Romeo
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

//  Ported to .NET based on the original:
//      Copyright (c) John Maddock 2006, Boost Software License, Version 1.0

using System;
using System.Diagnostics;

namespace MathXK
{

    /// <summary>
    /// Lanczos approximation of the gamma function
    /// <para>Γ(z) = Sqrt(2*π) * (z+g-0.5)^(z-0.5) * e^-(z+g-0.5) * A_g(z)</para>
    /// </summary>
    internal static class Lanczos
    {


        // Optimal values for G for each N are taken from
        // http://web.mala.bc.ca/pughg/phdThesis/phdThesis.pdf,
        // as are the theoretical error bounds.
        //
        // Constants calculated using the method described by Godfrey
        // http://my.fit.edu/~gabdo/gamma.txt and elaborated by Toth at
        // http://www.rskey.org/gamma.htm using NTL.RR at 1000 bit precision.
        //


        // Lanczos Coefficients for N=13 G=6.024680040776729583740234375
        // Max experimental error (with arbitary precision arithmetic) 1.196214e-17
        // Generated with compiler: Microsoft Visual C++ version 8.0 on Win32 at Mar 23 2006


        // Sum


        private static readonly double[] _N_series = {
            23531376880.41075968857200767445163675473,
            42919803642.64909876895789904700198885093,
            35711959237.35566804944018545154716670596,
            17921034426.03720969991975575445893111267,
            6039542586.35202800506429164430729792107,
            1439720407.311721673663223072794912393972,
            248874557.8620541565114603864132294232163,
            31426415.58540019438061423162831820536287,
            2876370.628935372441225409051620849613599,
            186056.2653952234950402949897160456992822,
            8071.672002365816210638002902272250613822,
            210.8242777515793458725097339207133627117,
            2.506628274631000270164908177133837338626
        };
        private static readonly double[] _D_series = {
            0,
            39916800,
            120543840,
            150917976,
            105258076,
            45995730,
            13339535,
            2637558,
            357423,
            32670,
            1925,
            66,
            1
        };


        // Sum_expG_scaled

        private static readonly double[] _N_SumExpGScaled = {
            56906521.91347156388090791033559122686859,
            103794043.1163445451906271053616070238554,
            86363131.28813859145546927288977868422342,
            43338889.32467613834773723740590533316085,
            14605578.08768506808414169982791359218571,
            3481712.15498064590882071018964774556468,
            601859.6171681098786670226533699352302507,
            75999.29304014542649875303443598909137092,
            6955.999602515376140356310115515198987526,
            449.9445569063168119446858607650988409623,
            19.51992788247617482847860966235652136208,
            0.5098416655656676188125178644804694509993,
            0.006061842346248906525783753964555936883222
        };
        private static readonly double[] _D_SumExpGScaled = {
            0,
            39916800,
            120543840,
            150917976,
            105258076,
            45995730,
            13339535,
            2637558,
            357423,
            32670,
            1925,
            66,
            1
        };

#if false
        // not being used

        // near 1
        private static readonly double[] _N1 = {
            2.208709979316623790862569924861841433016,
            -3.327150580651624233553677113928873034916,
            1.483082862367253753040442933770164111678,
            -0.1993758927614728757314233026257810172008,
            0.004785200610085071473880915854204301886437,
            -0.1515973019871092388943437623825208095123e-5,
            -0.2752907702903126466004207345038327818713e-7,
            0.3075580174791348492737947340039992829546e-7,
            -0.1933117898880828348692541394841204288047e-7,
            0.8690926181038057039526127422002498960172e-8,
            -0.2499505151487868335680273909354071938387e-8,
            0.3394643171893132535170101292240837927725e-9
       };

        // near 2
        private static readonly double[] _N2 = {
             6.565936202082889535528455955485877361223,
             -9.8907772644920670589288081640128194231,
             4.408830289125943377923077727900630927902,
             -0.5926941084905061794445733628891024027949,
             0.01422519127192419234315002746252160965831,
             -0.4506604409707170077136555010018549819192e-5,
             -0.8183698410724358930823737982119474130069e-7,
             0.9142922068165324132060550591210267992072e-7,
             -0.5746670642147041587497159649318454348117e-7,
             0.2583592566524439230844378948704262291927e-7,
             -0.7430396708998719707642735577238449585822e-8,
             0.1009141566987569892221439918230042368112e-8
          };


        private static double Sum_near_1(double dz)
        {
            double result = 0;
            for(int k = 1; k <= _N1.Length; ++k)
            {
                result += (-_N1[k-1]*dz)/(k*dz + k*k);
            }
            return result;
        }

        private static double Sum_near_2(double dz)
        {
            double result = 0;
            double z = dz + 2;
            for(int k = 1; k <= _N2.Length; ++k)
            {
                result += (-_N2[k-1]*dz)/(z + k*z + k*k - 1);
            }
            return result;
        }
#endif


        public const double G = 6.024680040776729583740234375;

        /// <summary>e^(G-0.5)</summary>
        public const double ExpGmHalf = 250.8060774637396673987542292884223801881291602213100597290193;

        /// <summary>e^-(G-0.5)</summary>
        public const double RecipExpGmHalf = 0.003987144211625315081300873894929618620683773340927465353086;

        /// <summary>
        /// Choose x such that (x+G-0.5)^(x-0.5) &lt; double.MaxValue
        /// <para>G-0.5 ~= 5.525</para>
        /// </summary>
        public const double MaxPowerTermX = 142; 

        /// <summary>
        /// Returns the Lanczos polynomial: A_g(z) * Sqrt(2*π)
        /// <para>Note: as z->inf, Sum-> ~2.5;</para> 
        /// <para>As z->0+, Sum->2.35e10 + 1075/z</para>
        /// </summary>
        public static double Series(double z)
        {
            Debug.Assert(z >= 256 * DoubleLimits.MinNormalValue);

            return Polynomial.EvalRational(_N_series, _D_series, z);
        }


        /// <summary>
        /// Returns the scaled Lanczos polynomial: A_g(z) * Sqrt(2*π) * e^-g
        /// <para>As z->inf, Sum-> 0.006...</para>
        /// <para>At z = 1, Sum = 0.217...</para>
        /// <para>As z->0+, Sum->5.7e7 + 2.6/z</para>
        /// </summary>
        public static double SeriesExpGScaled(double z)
        {
            return Polynomial.EvalRational(_N_SumExpGScaled, _D_SumExpGScaled, z);
        }


#if false
        // not being used but kept for future reference

        /// <summary>
        /// Returns ln(|Γ(2+dz)|) using the Lanczos approzimation 
        /// </summary>
        public static double Lgamma2p(double dz) 
        {
            double result = dz * Math.Log( (dz + (Lanczos.G + 1.5))/Math.E );
            result += Math2.Log1p(dz / (Lanczos.G + 1.5)) * 1.5;
            result += Math2.Log1p(Lanczos.Sum_near_2(dz));
            return result;        
        }

        /// <summary>
        /// Returns ln(|Γ(1+dz)|) using the Lanczos approzimation 
        /// </summary>
        public static double Lgamma1p(double dz) 
        {
            double result = dz * Math.Log((dz + (Lanczos.G + 0.5)) / Math.E);
            result += Math2.Log1p(dz / (Lanczos.G + 0.5)) / 2;
            result += Math2.Log1p(Lanczos.Sum_near_1(dz));
            return result;        
        }
#endif

        /// <summary>
        /// Returns the Laczos approximation of Γ(z) 
        /// </summary>
        public static double Tgamma(double x)
        {
            Debug.Assert(x > 0);

            // check for large x to avoid potential Inf/Inf below
            if (x > DoubleLimits.MaxGamma + 1)
                return double.PositiveInfinity;

            // use the Lanczos approximation
            double xgh = x + (Lanczos.G - 0.5);
            double result = Lanczos.Series(x);
            if (x > MaxPowerTermX) {
                // we're going to overflow unless this is done with care:
                double hp = Math.Pow(xgh, (x / 2) - 0.25);
                result *= (hp / Math.Exp(xgh));
                result *= hp;
            } else {
                result *= (Math.Pow(xgh, x - 0.5) / Math.Exp(xgh));
            }

            return result;
        }


        /// <summary>
        /// Returns the Laczos approximation of ln(Γ(z)) 
        /// </summary>
        public static double Lgamma(double z)
        {
            Debug.Assert(z > 0);

            double zgh = z + (Lanczos.G - 0.5);
            double result = ((Math.Log(zgh) - 1) * (z - 0.5)) + Math.Log(SeriesExpGScaled(z));
            return result;
        }


        /// <summary>
        /// Returns Γ(z)/Γ(z+delta) using the Lanczos approximation
        /// </summary>
        public static double TgammaDeltaRatio(double x, double delta)
        {
            Debug.Assert(x > 0 && delta + x > 0);

            // Compute Γ(x)/Γ(x+Δ), which reduces to:
            //
            // ((x+g-0.5)/(x+Δ+g-0.5))^(x-0.5) * (x+Δ+g-0.5)^-Δ * e^Δ * (Sum(x)/Sum(x+Δ))
            // = (1+Δ/(x+g-0.5))^-(x-0.5) * (x+Δ+g-0.5)^-Δ * e^Δ * (Sum(x)/Sum(x+Δ))

            double xgh = x + (Lanczos.G - 0.5);
            double xdgh = x + delta + (Lanczos.G - 0.5);

            double mu = delta / xgh;
            double result = Lanczos.Series(x) / Lanczos.Series(x + delta);
            result *= (Math.Abs(mu) < 0.5) ? Math.Exp((0.5 - x) * Math2.Log1p(mu)) : Math.Pow(xgh / xdgh, x - 0.5);
            result *= Math.Pow(Math.E /xdgh, delta);

            return result;
        }

        /// <summary>
        /// Returns x^a/Γ(a) using the Lanczos approximation  
        /// </summary>
        /// <param name="a"></param>
        /// <param name="x"></param>
        public static double PowDivGamma(double x, double a) 
        {
            // Lanczos calculation. Compute:
            // (x/agh)^a * e^a * Sqrt(agh/e) / SumExpGScaled(a)

            double agh = a + (Lanczos.G - 0.5);
            double factor = Constants.RecipSqrtE * Math.Sqrt(agh) / SeriesExpGScaled(a);

            // when a <= 1, there is no danger of overflow or underflow
            // when a >= 1 && t >= 1, factor > 1, so we can only overflow the whole expression

            double t = x/agh;
            if ( a <= 1 || t >= 1 )  
                return Math.Pow(t, a) * Math.Exp(a) * factor;

            // at this point a > 1 && 0 < t < 1
            double alogt = a * Math.Log(t);
            if (a <= DoubleLimits.MaxLogValue && alogt >= DoubleLimits.MinLogValue ) {
                return Math.Pow(t, a) * Math.Exp(a) * factor;
            }
            if (a <= 2 * DoubleLimits.MaxLogValue && alogt >= 2 * DoubleLimits.MinLogValue) {
                double hp = Math.Pow(t, a / 2) * Math.Exp(a / 2);
                return hp * factor * hp;
            }

            return Math.Exp(alogt + a + Math.Log(factor));
       
        }


        /// <summary>
        /// Returns (x^a)(e^-x)/Γ(a) using the Lanczos approximation
        /// </summary>
        /// <param name="a"></param>
        /// <param name="x"></param>
        /// <returns></returns>
        public static double PowExpDivGamma(double a, double x)
        {
            Debug.Assert(a > 0 && x >= 0);

            // Lanczos
            // Compute: (x/agh)^a * e^(a-x) * Sqrt(agh/e)  / SumExpGScaled(a)

            double agh = a + (Lanczos.G - 0.5);
            double factor = Constants.RecipSqrtE * Math.Sqrt(agh) / SeriesExpGScaled(a);

            // for a ~ x
            // a * (ln(1+(x-agh)/agh) - ((x-agh)/agh)) === a*ln(x/agh) + a - a * x/agh 
            var d = (x - agh) / agh;
            if (Math.Abs(d) < 0.5)
                return factor * Math.Exp(a * Math2.Log1pmx(d) + x / agh * (0.5 - Lanczos.G));

            double t = x / agh;
            double alogt = a * Math.Log(t);
            double amx = a - x;
            if (t > 1) {
                // x > agh

                if (amx >= DoubleLimits.MinLogValue && alogt <= DoubleLimits.MaxLogValue)
                    return Math.Pow(t, a) * Math.Exp(amx) * factor;
                if (amx >= 2 * DoubleLimits.MinLogValue && alogt <= 2 * DoubleLimits.MaxLogValue) {
                    double h = Math.Pow(t, a / 2) * Math.Exp(amx / 2);
                    return h * factor * h;
                }

            } else {
                // x < agh
                if (amx <= DoubleLimits.MaxLogValue && alogt >= DoubleLimits.MinLogValue)
                    return Math.Pow(t, a) * Math.Exp(amx) * factor;
                if (amx <= 2 * DoubleLimits.MaxLogValue && alogt >= 2 * DoubleLimits.MinLogValue) {
                    double h = Math.Pow(t, a / 2) * Math.Exp(amx / 2);
                    return h * factor * h;
                }

            }

            return Math.Exp(alogt + amx + Math.Log(factor));
        }


        /// <summary>
        /// Returns Beta(a, b) using the Lanczos approximation
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <returns></returns>
        public static double Beta(double a, double b)
        {
            // Lanczos calculation. Compute:
            // (agh/cgh)^a * (bgh/cgh)^b * Sqrt((cgh * e)/(agh * bgh))
            //  = (1-b/cgh)^a * (1 - a/cgh)^b * Sqrt((cgh * e)/(agh * bgh)) 

            double c = a + b;
            double agh = a + (Lanczos.G - 0.5);
            double bgh = b + (Lanczos.G - 0.5);
            double cgh = c + (Lanczos.G - 0.5);


            double factor = (SeriesExpGScaled(a) / SeriesExpGScaled(c)) * SeriesExpGScaled(b);
            double result = Constants.SqrtE * Math.Sqrt((cgh/agh)/bgh) * factor;

            // calculate (agh/cgh)^a or the equivalent (1 - b/cgh)^a 
            // calculate (bgh/cgh)^b or the equivalent (1 - a/cgh)^b 

            double t1 = b / cgh;
            result *= (t1 <= 0.5) ? Math.Exp(a * Math2.Log1p(-t1)) : Math.Pow(agh / cgh, a);

            double t2 = a / cgh;
            result *= (t2 <= 0.5) ? Math.Exp(b * Math2.Log1p(-t2)) : Math.Pow(bgh / cgh, b);

            return result;

        }

        /// <summary>
        /// Returns the Log(Beta(a,b)) using the Lanczos approximation
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <returns></returns>
        public static double LogBeta(double a, double b)
        {
            // B(a,b) == B(b, a)
            if (a < b)
                Utility.Swap(ref a, ref b);

            // from this point a >= b

            // Lanczos calculation. Compute:
            // ln((agh/cgh)^a * (bgh/cgh)^b * Sqrt((cgh * e)/(agh * bgh)))
            // = ln((1-b/cgh)^a * (1 - a/cgh)^b * Sqrt((cgh * e)/(agh * bgh))) 

            double c = a + b;
            double agh = a + (Lanczos.G - 0.5);
            double bgh = b + (Lanczos.G - 0.5);
            double cgh = c + (Lanczos.G - 0.5);


            double logValue = 0;
            double factor = (SeriesExpGScaled(a) / SeriesExpGScaled(c)) * SeriesExpGScaled(b);
            factor *= Constants.SqrtE * Math.Sqrt((cgh / agh) / bgh);

            // calculate a*ln(agh/cgh) or the equivalent a*ln(1 - b/cgh) 
            // calculate b*ln(bgh/cgh) or the equivalent b*ln(1 - a/cgh) 

            double t1 = b / cgh;
            logValue += (t1 <= 0.5) ? a * Math2.Log1p(-t1) : a * Math.Log(agh / cgh);

            double t2 = a / cgh;
            logValue += (t2 <= 0.5) ? b * Math2.Log1p(-t2) : b * Math.Log(bgh / cgh);

            logValue += Math.Log(factor);

            return logValue;
        }


    };


} 






