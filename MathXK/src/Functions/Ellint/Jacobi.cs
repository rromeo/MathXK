//  Copyright (c) 2015 Rocco Romeo
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

//  Ported to .NET based on the original:
//      Copyright (c) John Maddock 2012, Boost Software License v1.0


using System;
using System.Diagnostics;

namespace MathXK
{
    public partial class Math2
    {
        /// <summary>
        /// Results of Jacobi Elliptical Sn, Cn, Dn
        /// </summary>
        struct JacobiSnCnDn
        {
            public double Sn { get; internal set; }
            public double Cn { get; internal set; }
            public double Dn { get; internal set; }

            public double Ns { get { return 1 / Sn; } }
            public double Nc { get { return 1 / Cn; } }
            public double Nd { get { return 1 / Dn; } }

            public double Cd { get { return Cn / Dn; } }
            public double Dc { get { return Dn / Cn; } }

            public double Sd { get { return Sn / Dn; } }
            public double Ds { get { return Dn / Sn; } }

            public double Sc { get { return Sn / Cn; } }
            public double Cs { get { return Cn / Sn; } }
        };

        /// <summary>
        /// Computes the Jacobi elliptic functions using the arithmetic-geometric mean
        /// </summary>
        /// <param name="k">The modulus</param>
        /// <param name="u">The argument</param>
        /// <param name="anm1"></param>
        /// <param name="bnm1"></param>
        /// <param name="N"></param>
        /// <param name="pTn"></param>
        /// <returns></returns>
        /// <see href="http://dlmf.nist.gov/22.20#ii"/>
        private static double Jacobi_Recurse(double k, double u, double anm1, double bnm1, int N, out double pTn)
        {

            ++N;
            double Tn;
            double cn = (anm1 - bnm1) / 2;
            double an = (anm1 + bnm1) / 2;

            if (cn < DoubleLimits.MachineEpsilon) {
                Tn = Math2.Ldexp(1, N) * u * an;
            } else {
                Tn = Jacobi_Recurse(k, u, an, Math.Sqrt(anm1 * bnm1), N, out double unused);
            }

            pTn = Tn;
            return (Tn + Math.Asin((cn / an) * Math.Sin(Tn))) / 2;
        }


        /// <summary>
        /// Returns the Jacobi elliptic functions: sn, cn, dn
        /// </summary>
        /// <param name="k">The modulus</param>
        /// <param name="u">The argument</param>
        /// <returns></returns>
        static JacobiSnCnDn JacobiElliptic(double k, double u)
        {
                if ( double.IsNaN(k) || double.IsInfinity(k) || double.IsNaN(u) || double.IsInfinity(u)) {
                    Policies.ReportDomainError("Jacobi(k: {0}, u: {1}): Requires finite u, k", k, u);
                    return new JacobiSnCnDn { Sn = double.NaN, Cn = double.NaN, Dn = double.NaN };
                }

                // For k < 0 or k > 1, see: http://dlmf.nist.gov/22.17#i

                if (k < 0)
                    k = -k;

                if (k > 1) {
                    var results = JacobiElliptic(1 / k, u * k);
                    return new JacobiSnCnDn { Sn = results.Sn / k, Cn = results.Dn, Dn = results.Cn };
                }

                if (u == 0) 
                    return new JacobiSnCnDn { Sn = 0, Cn = 1, Dn = 1 };

                if ( k == 0 ) 
                    return new JacobiSnCnDn { Sn = Math.Sin(u), Cn = Math.Cos(u), Dn = 1 };

                if (k == 1) {
                    double sech = 1 / Math.Cosh(u);
                    return new JacobiSnCnDn { Sn = Math.Tanh(u), Cn = sech, Dn = sech };
                }

                double sn, cn, dn;
                if (k < DoubleLimits.RootMachineEpsilon._4) {

                    // Asymptotic form from A&S 16.13:
                    // and http://dlmf.nist.gov/22.10#ii

                    double tu = Math.Tanh(u);
                    double su = Math.Sin(u);
                    double cu = Math.Cos(u);
                    double m = k * k;
                    sn = su - m * (u - cu*su) * cu / 4;
                    cn = cu + m * (u - cu*su) * su / 4;          
                    dn = 1 - m * su * su / 2;

                } else if (k*k > 1 - DoubleLimits.RootMachineEpsilon._2) {
                    // Asymptotic form from A&S 16.15:
                    // and http://dlmf.nist.gov/22.10#ii
                    // O(k'^4) = O((1-k^2)^2)

                    double tu = Math.Tanh(u);
                    double cu = Math.Cosh(u);
                    double su = Math.Sinh(u);
                    double sech = 1 / Math.Cosh(u);

                    double kc = 1 - k;
                    double kp2 = kc * (2 - kc); // k'^2 = 1-k^2

                    dn = (sech + kp2 * (u * sech + su) * tu) / 4;
                    cn = (sech + kp2 * (u * sech - su) * tu) / 4;
                    sn = (tu + kp2 * (tu - u/(cu*cu))) / 4;
                    //sn -= (72 * u * cu + 4 * (8 * u * u - 5) * su - 19 * Math.Sinh(3 * u) + Math.Sinh(5 * u)) * sec * sec * sec * kp2 * kp2 / 512;

                } else {

                    double kc = 1 - k;
                    double k_prime = k < 0.5 ? Math.Sqrt((1 - k) * (1 + k)) : Math.Sqrt(kc * (2 - kc));

                    double T1;
                    double T0 = Jacobi_Recurse(k, u, 1, k_prime, 0, out T1);

                    sn = Math.Sin(T0);
                    cn = Math.Cos(T0);
                    dn = Math.Cos(T0) / Math.Cos(T1 - T0);
                }


                return new JacobiSnCnDn { Sn = sn, Cn = cn, Dn = dn };

        }

        /// <summary>
        /// Returns the Jacobi Elliptic sn function
        /// <para>If u = F(φ, k), then sn(u, k) = sin(φ)</para>
        /// </summary>
        /// <param name="k">The modulus</param>
        /// <param name="u">The argument</param>
        /// <returns></returns>
        public static double JacobiSN(double k, double u)
        {
            return JacobiElliptic(k, u).Sn;
        }

        /// <summary>
        /// Returns the Jacobi Elliptic cn function
        /// <para>If u = F(φ, k), then cn(u, k) = cos(φ)</para>
        /// </summary>
        /// <param name="k">The modulus</param>
        /// <param name="u"></param>
        /// <returns></returns>
        public static double JacobiCN(double k, double u)
        {
            return JacobiElliptic(k, u).Cn;
        }

        /// <summary>
        /// Returns the Jacobi Elliptic dn function
        /// <para>If u = F(φ, k), then dn(u, k) = Sqrt(1-k^2*sin^2(φ))</para>
        /// </summary>
        /// <param name="k">The modulus</param>
        /// <param name="u">The argument</param>
        /// <returns></returns>
        public static double JacobiDN(double k, double u)
        {
            return JacobiElliptic(k, u).Dn;
        }

        /// <summary>
        /// Returns the Jacobi Elliptic cd function
        /// </summary>
        /// <param name="k">The modulus</param>
        /// <param name="u">The argument</param>
        /// <returns></returns>
        public static double JacobiCD(double k, double u)
        {
            return JacobiElliptic(k, u).Cd;
        }

        /// <summary>
        /// Returns the Jacobi Elliptic dc function
        /// </summary>
        /// <param name="k">The modulus</param>
        /// <param name="u">The argument</param>
        /// <returns></returns>
        public static double JacobiDC(double k, double u)
        {
            return JacobiElliptic(k, u).Dc;
        }

        /// <summary>
        /// Returns the Jacobi Elliptic ns function
        /// <para>ns = 1/sn</para>
        /// </summary>
        /// <param name="k">The modulus</param>
        /// <param name="u">The argument</param>
        /// <returns></returns>
        public static double JacobiNS(double k, double u)
        {
            return JacobiElliptic(k, u).Ns;
        }

        /// <summary>
        /// Returns the Jacobi Elliptic sd function
        /// </summary>
        /// <param name="k">The modulus</param>
        /// <param name="u">The argument</param>
        /// <returns></returns>
        public static double JacobiSD(double k, double u)
        {
            return JacobiElliptic(k, u).Sd;
        }

        /// <summary>
        /// Returns the Jacobi Elliptic ds function
        /// </summary>
        /// <param name="k">The modulus</param>
        /// <param name="u">The argument</param>
        /// <returns></returns>
        public static double JacobiDS(double k, double u)
        {
            return JacobiElliptic(k, u).Ds;
        }

        /// <summary>
        /// Returns the Jacobi Elliptic nc function
        /// <para>nc = 1/cn</para>
        /// </summary>
        /// <param name="k">The modulus</param>
        /// <param name="u">The argument</param>
        /// <returns></returns>
        public static double JacobiNC(double k, double u)
        {
            return JacobiElliptic(k, u).Nc;
        }

        /// <summary>
        /// Returns the Jacobi Elliptic nd function
        /// <para>nd = 1/dn</para>
        /// </summary>
        /// <param name="k">The modulus</param>
        /// <param name="u">The argument</param>
        /// <returns></returns>
        public static double JacobiND(double k, double u)
        {
            return JacobiElliptic(k, u).Nd;
        }

        /// <summary>
        /// Returns the Jacobi Elliptic sc function
        /// </summary>
        /// <param name="k">The modulus</param>
        /// <param name="u">The argument</param>
        /// <returns></returns>
        public static double JacobiSC(double k, double u)
        {
            return JacobiElliptic(k, u).Sc;
        }

        /// <summary>
        /// Returns the Jacobi Elliptic cs function
        /// </summary>
        /// <param name="k">The modulus</param>
        /// <param name="u">The argument</param>
        /// <returns></returns>
        public static double JacobiCS(double k, double u)
        {
            return JacobiElliptic(k, u).Cs;
        }
    }

} // namespaces

