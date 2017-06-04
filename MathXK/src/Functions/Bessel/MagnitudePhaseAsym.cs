//  Copyright (c) 2013 Rocco Romeo
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

//  Ported to .NET based on the original:
//      Copyright (c) 2007 John Maddock, Boost Software License v1.0



using System;

namespace MathXK
{

    internal partial class _Bessel
    {
        /// <summary>
        /// Asymptotic magnitude phase approach
        /// </summary>
        private static class MagnitudePhase
        {
            static double Amplitude(double v, double x)
            {
                // see: http://dlmf.nist.gov/10.18.iii


                double k = 1;
                double mu = 4 * v * v;
                double txq = 4 * x * x;
                double sq = 1;
                double term = 1;


                // The following is the generalized form of the following
                // as per DLMF:
                // double sum = 1;
                // sum += (mu - 1) / (2 * txq);
                // sum += 3 * (mu - 1) * (mu - 9) / (txq * txq * 8);
                // sum += 15 * (mu - 1) * (mu - 9) * (mu - 25) / (txq * txq * txq * 8 * 6);



                double sum = 1, lastSum = sum;

                for(int n = 0; n < Policies.MaxSeriesIterations ; n++) {
                    // check to see when the series starts to diverge        
                    double mult = (k/(k + 1.0)) * ((mu - sq * sq) / txq);
                    if (Math.Abs(mult) >= 1.0) {
                        Policies.ReportConvergenceError("Series diverges after {0} iterations", k);
                        return double.NaN;
                    }


                    lastSum = sum;
                    term *= mult;
                    sum += term;
                    k += 2;
                    sq += 2;

                    if (Math.Abs(term) <= Math.Abs(lastSum) * Policies.SeriesTolerance) 
                        return Constants.RecipSqrtHalfPI * Math.Sqrt(sum / x);
                    

                }

                Policies.ReportConvergenceError("Series did not converge after {0} iterations", Policies.MaxSeriesIterations);
                return double.NaN;

            }

            static double PhaseAdjustment(double v, double x)
            {
                // see: http://dlmf.nist.gov/10.18.iii

                //
                // Calculate the phase of J(v, x) and Y(v, x) for large x.
                // See A&S 9.2.29.
                // Note must add (x - (1/2v+1/4)PI) to get the actual phase.
                //
                double mu = 4 * v * v;
                double denom = 4 * x;

                double m = (mu-1)/denom;
                double rd = 1/(denom * denom);

                double p1 = (mu - 25)/6;
                double p2 = (1073 + mu*(-114 + mu))/5;
                double p3 = (-375733 + mu*(54703 + mu*(-1535 + mu*5)))/14;

                return m * (0.5 + rd*(p1 + rd*(p2 + rd*p3)));
            }

            // Usability Limits:
            // Using the magnitude/phase equations: http://dlmf.nist.gov/10.18#iii
            // Mathematica min phase computation
            // Solving ((375733 - 430436 mu + 56238 mu^2 - 1540 mu^3 + 5 mu^4)/(14*(4x)^7) / x ) < epsilon, for x
            // F[v_] := Block[{epsilon = 2^-52, mu = 4 v^2}, (Abs[(375733 - 430436 mu + 56238 mu^2 - 1540 mu^3 + 5 mu^4)/(229376 epsilon)])^(1/8)]
            // Limit[F[v]/v, v -> Infinity] = 32 * (5/7)^(1/8) * 2^(5/8)

            // These are the approximate limits for the phase portion
            public const double SmallestX = 250;
            public const double V_Multiple = 47.318146131550666944983629673470458177564741046220; //32 * (5/7)^(1/8) * 2^(5/8)


            ///<summary>
            /// Returns the minimum x to use Asymptotic Magnitude/Phase
            /// <para>Currently set to x = Max(250, 47.32v)</para>
            /// </summary>
            public static double MinX(double v)
            {
                return Math.Max(SmallestX, V_Multiple * v);
            }

            /// <summary>
            /// Compute Bessel J and Y using the asymptotic Amplitude Phase approach
            /// </summary>
            /// <param name="v"></param>
            /// <param name="x"></param>
            /// <returns></returns>
            public static (double J, double Y) BesselJY(double v, double x)
            {
                // http://dlmf.nist.gov/10.18#iii
                // See A&S 9.2.19.

                // compute the sin and cos of the base angle
                var (sin, cos) = Math2.SinCos(x, -(Math2.Mod(v / 2, 2) + 0.25));

                // Get the phase and amplitude:
                double ampl = Amplitude(v, x);
                double phaseAdj = PhaseAdjustment(v, x);

                // Calculate the sin and cos of the phase, using:
                // sin(x+p) = sin(x)cos(p) + cos(x)sin(p)
                // cos(x+p) = cos(x)cos(p) - sin(x)sin(p)

                double sin_phase = Math.Sin(phaseAdj) * cos + Math.Cos(phaseAdj) * sin;
                double cos_phase = Math.Cos(phaseAdj) * cos - Math.Sin(phaseAdj) * sin;

                double J = cos_phase * ampl;
                double Y = sin_phase * ampl;

                return (J, Y);
            }

        };

    }

} // namespaces


