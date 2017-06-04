//  Copyright (c) 2015 Rocco Romeo
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

//  Ported to .NET based on the original:
//      Copyright (c) 2015 John Maddock, Boost Software License v1.0


using System;
using System.Diagnostics;

namespace MathXK
{

    public partial class Math2
    {
        /// <summary>
        /// Carlson's symmetric form of elliptical integrals 
        /// <para>RG = 1/4 ∫ sqrt((t+x)(t+y)(t+z))*(x/(t+x)+y/(t+y)+z/(t+z))*t*dt t={0,∞}</para>
        /// <para>RG = 1/(4π) ∫∫ sqrt(x*sin^2(θ)*cos^2(φ)+y*sin^2(θ)*sin^2(φ)+z*cos^2(θ))*sin(θ) dθdφ θ={0,π} φ={0,2π}</para>
        /// </summary>
        /// <param name="x">Argument. Requires x >= 0</param>
        /// <param name="y">Argument. Requires y >= 0</param>
        /// <param name="z">Argument. Requires z >= 0</param>
        /// <returns></returns>
        /// <remarks>
        /// Numerical Computation of Real or Complex Elliptic Integrals, B.C. Carlson, 6 Sep 1994
        /// </remarks>
        /// <see href="http://arxiv.org/abs/math/9409227v1"/>
        public static double EllintRG(double x, double y, double z)
        {
            if (!(x >= 0) && !(y >= 0) && !(z >= 0)) {
                Policies.ReportDomainError("EllintRG(x: {0}, y: {1}, z: {2}) requires finite x >= 0, y >= 0, z >= 0", x, y, z);
                return double.NaN;
            }

            // reorder x,y,z such that x <= y <= z 
            if (x > y)
                Utility.Swap(ref x, ref y);
            if (y > z)
                Utility.Swap(ref y, ref z);
            if (x > y)
                Utility.Swap(ref x, ref y);

            if (x == 0) {
                if ( y == 0 )
                    return Math.Sqrt(z) / 2; // RG(0, 0, z)
                if ( y == z )
                    return (Math.PI / 4) * Math.Sqrt(y); // RG(0, y, y)


                double xn = Math.Sqrt(z);
                double yn = Math.Sqrt(y);
                double a = (xn + yn) / 2;
                double sum = 0;
                double sum_pow = 0.25;

                double c = xn - yn;
                while (Math.Abs(c) >= 2.7 * DoubleLimits.RootMachineEpsilon._2 * xn) {
                    double t = Math.Sqrt(xn * yn);
                    xn = (xn + yn) / 2;
                    yn = t;

                    c = xn - yn;
                    sum_pow *= 2;
                    sum += sum_pow * (c*c);
                }
                return (a * a - sum) * ((Math.PI/2) /(xn + yn));
            }

            if (x == y) {
                if (x == z)
                    return Math.Sqrt(x); // RG(x, x, x)
                return (x * EllintRC(z, x) + Math.Sqrt(z)) / 2; //RG(x, x, z)
            }

            if (y == z) 
                return (y * EllintRC(x, y) + Math.Sqrt(x)) / 2; //RG(x, y, y)

            // Swap y and z to avoid cancellation error in the equation:
            // We want the multiplier (x - z)(y - z) < 0 so that the term is additive. 
            // So let x < z and let y > z 

            Utility.Swap(ref y, ref z);
            return (z * EllintRF(x, y, z) - (x - z) * (y - z) * EllintRD(x, y, z) / 3 + Math.Sqrt(x * y / z)) / 2;
        }


    }
}