//  Copyright (c) 2013 Rocco Romeo
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

using System;
using System.Diagnostics;

namespace MathXK
{

    static class Trig
    {

        ///<summary>Computed three stage Cody-Waite range reduction constants for π</summary>
        static class Pi
        {
            ///<summary>Double precision π. Note that this value &gt; true π</summary>
            public const double Value = Math.PI;

            ///<summary>1/π</summary>
            public const double Reciprocal = 0.318309886183790671537767526745028724068919291480912897495334;

            ///<summary>Floor[Pi*2^32]/2^32 = 1686629713/536870912</summary>
            public const double A = 3.14159265346825122833251953125;

            ///<summary>Floor[(Pi - 1686629713/536870912)*(2^67)]/2^67 = 2242054355/18446744073709551616</summary>
            public const double B = 1.215420101260793195319109827323700301349163055419921875e-10;

            ///<summary>Pi - (1686629713/536870912 + 2242054355/18446744073709551616)</summary>
            public const double C = 4.044532497591901464799369240189515432952404807816406286e-21;

            // If we go higher, than MaxLimit, k*A and k*B are no longer exact, so the result would be inexact
            public const int MaxMultiple = 1048576; // 2^20
            public const double MaxLimit = MaxMultiple * Value; // 2^20 * Pi ~= 3.29e6


            static Pi()
            {
                // Check that constants are exact
                Debug.Assert(A == BitConverter.Int64BitsToDouble(0x400921fb54400000L));
                Debug.Assert(B == BitConverter.Int64BitsToDouble(0x3de0b4611a600000L));
                Debug.Assert((A + (B + C)) == Value);
            }

        }

        ///<summary>Computed three stage Cody-Waite range reduction constants for π/2</summary>
        static class HalfPi
        {
            ///<summary>Double precision π/2. Note that this value &gt; true π/2</summary>
            public const double Value = Math.PI / 2;

            // Constants A,B and reduction range (2^20) for Cody-Waite are based on Sun's Netlib approach.
            // See: http://www.netlib.org/fdlibm/e_rem_pio2.c
            // ====================================================
            // Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
            //
            // Developed at SunSoft, a Sun Microsystems, Inc. business.
            // Permission to use, copy, modify, and distribute this
            // software is freely granted, provided that this notice 
            // is preserved.
            // ====================================================

            ///<summary>1/(π/2)</summary>
            public const double Reciprocal = 0.636619772367581343075535053490057448137838582961825794990669;

            ///<summary>Floor[(Pi/2)*2^33]/2^33 = 1686629713/1073741824</summary>
            public const double A = 1.570796326734125614166259765625;

            ///<summary>Floor[((Pi/2) - 1686629713/1073741824)*(2^68)]/2^68 = 2242054355/36893488147419103232</summary>
            public const double B = 6.077100506303965976595549136618501506745815277099609375e-11;

            ///<summary>Pi/2 - (1686629713/1073741824 + 2242054355/36893488147419103232)</summary>
            public const double C = 2.022266248795950732399684620094757716476202403908203143e-21;

            // If we go higher, than MaxLimit, k*PI_A and k*PI_B are no longer exact, so the result would be inexact
            public const int MaxMultiple = 1048576; // 2^20
            public const double MaxLimit = MaxMultiple * Value; //2^20 * (Pi/2) ~= 1.6e6

            static HalfPi()
            {
                // Check that PI/2 constants are exact
                Debug.Assert(A == BitConverter.Int64BitsToDouble(0x3ff921fb54400000L));
                Debug.Assert(B == BitConverter.Int64BitsToDouble(0x3dd0b4611a600000L));
                Debug.Assert((A + (B + C)) == Value);
            }
        }

        ///<summary>The maximum angle allowed for Pi reductions</summary>
        public const double PiReductionLimit = Pi.MaxLimit;

        ///<summary>The maximum angle allowed for Half Pi reductions</summary>
        public const double HalfPiReductionLimit = HalfPi.MaxLimit;



        /// <summary>
        /// Returns x/π and its remainder.
        /// Note that both the multiple and remainder will have the same sign as the dividend: <paramref name="x"/>
        /// <para>Multiple = Sign(x) * Floor(|x|/π)</para>
        /// <para>Remainder = Sign(x) * (|x|- π*Floor(|x|/π))</para>
        /// </summary>
        public static (double Multiple, double Remainder) DivRemPI(double x)
        {

            if (double.IsNaN(x) || double.IsInfinity(x)) {
                Policies.ReportDomainError("DivRemPI(x: {0}) Requires finite argument", x);
                return (double.NaN, double.NaN);
            }

            double absX = Math.Abs(x);
            if (absX > Pi.MaxLimit) {
                Policies.ReportNotImplementedError("DivRemPI(x: {0}) |x| > {1} is not implemented", x, Pi.MaxLimit);
                return (double.NaN, double.NaN);
            }

            if (absX < Pi.Value)
                return (0, x);

            double k = Math.Floor(absX * Pi.Reciprocal);
            double rem;
            for (; ; ) {
                rem = ((absX - (k * Pi.A)) - k * Pi.B) - k * Pi.C;

                if (rem >= 0)
                    break;

                k--;
            }


            if (x < 0) {
                k = -k;
                rem = -rem;
            }

            Debug.Assert(Math.Abs(rem) < Pi.Value);

            return (k, rem);
        }

        /// <summary>
        /// Returns a pair representing x/π ensuring that |remainder| ≤ π/2
        /// <para>Multiple = Sign(x) * Floor(|x|/π + 1/2)</para>
        /// <para>Remainder = Sign(x) * (|x|- π*Floor(|x|/π + 1/2))</para>
        /// </summary>
        public static (double Multiple, double Remainder) RangeReducePI(double x)
        {
            if (double.IsNaN(x) || double.IsInfinity(x)) {
                Policies.ReportDomainError("RangeReducePI(x: {0}) Requires finite argument", x);
                return (double.NaN, double.NaN);
            }

            double absX = Math.Abs(x);
            if (absX > Pi.MaxLimit) {
                Policies.ReportNotImplementedError("RangeReducePI(x: {0}) |x| > {1} is not implemented", x, Pi.MaxLimit);
                return (double.NaN, double.NaN);
            }

            if (absX < Pi.Value / 2)
                return (0, x);

            double k = Math.Floor(absX * Pi.Reciprocal + 0.5);
            double rem;
            for (; ; ) {
                rem = ((absX - (k * Pi.A)) - k * Pi.B) - k * Pi.C;

                if (rem >= -Pi.Value / 2)
                    break;

                k--;
            }

            if (x < 0) {
                k = -k;
                rem = -rem;
            }

            Debug.Assert(Math.Abs(rem) <= Pi.Value / 2);


            return (k, rem);

        }

        /// <summary>
        /// Returns x/(π/2) and its remainder.
        /// Note that both the multiple and remainder will have the same sign as the dividend: <paramref name="x"/>
        /// <para>Multiple = Sign(x) * Floor(|x|/(π/2))</para>
        /// <para>Remainder = Sign(x) * (|x|- (π/2)*Floor(|x|/(π/2)))</para>
        /// </summary>
        public static (double Multiple, double Remainder) DivRemHalfPI(double x)
        {

            if (double.IsNaN(x) || double.IsInfinity(x)) {
                Policies.ReportDomainError("DivRemHalfPI(x: {0}) Requires finite argument", x);
                return (double.NaN, double.NaN);
            }

            double absX = Math.Abs(x);
            if (absX > HalfPi.MaxLimit) {
                Policies.ReportNotImplementedError("DivRemHalfPI(x: {0}) |x| > {1} is not implemented", x, HalfPi.MaxLimit);
                return (double.NaN, double.NaN);
            }

            if (absX < HalfPi.Value)
                return (0, x);

            double k = Math.Floor(absX * HalfPi.Reciprocal);
            double rem;
            for (; ; ) {
                rem = ((absX - (k * HalfPi.A)) - k * HalfPi.B) - k * HalfPi.C;

                if (rem >= 0)
                    break;

                k--;
            }


            if (x < 0) {
                k = -k;
                rem = -rem;
            }

            Debug.Assert(Math.Abs(rem) < HalfPi.Value);

            return (k, rem);
        }

        /// <summary>
        /// Returns a pair representing x/(π/2) ensuring that |remainder| ≤ π/4
        /// <para>Multiple = Sign(x) * Floor(|x|/(π/2) + 1/2)</para>
        /// <para>Remainder = Sign(x) * (|x|-(π/2)*Floor(|x|/(π/2) + 1/2))</para>
        /// </summary>
        public static (double Multiple, double Remainder) RangeReduceHalfPI(double x)
        {
            if (double.IsNaN(x) || double.IsInfinity(x)) {
                Policies.ReportDomainError("RangeReduceHalfPI(x: {0}) Requires finite argument", x);
                return (double.NaN, double.NaN);
            }

            double absX = Math.Abs(x);
            if (absX > HalfPi.MaxLimit) {
                Policies.ReportNotImplementedError("RangeReduceHalfPI(x: {0}) |x| > {1} is not implemented", x, HalfPi.MaxLimit);
                return (double.NaN, double.NaN);
            }

            if (absX < HalfPi.Value / 2)
                return (0, x);

            double k = Math.Floor(absX * HalfPi.Reciprocal + 0.5);
            double rem;
            for (; ; ) {
                rem = ((absX - (k * HalfPi.A)) - k * HalfPi.B) - k * HalfPi.C;

                if (rem >= -HalfPi.Value / 2)
                    break;

                k--;
            }

            if (x < 0) {
                k = -k;
                rem = -rem;
            }

            Debug.Assert(Math.Abs(rem) <= HalfPi.Value / 2);


            return (k, rem);

        }


    };


}