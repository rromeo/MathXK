//  (C) Copyright Rocco Romeo 2013.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

namespace MathXK
{

    /// <summary>
    /// Limits for double precision values
    /// </summary>
    public static class DoubleLimits
    {

        /// <summary>
        /// Number of decimal digits of precision
        /// <para>Value = 15</para>
        /// </summary>
        public const int Digits = 15;

        /// <summary>
        /// Smallest value such that 1.0+x != 1.0
        /// <para>Value = 2^(-52) = 2.2204460492503131e-016</para>
        /// </summary>
        public const double MachineEpsilon = 2.2204460492503131e-016;

        /// <summary>
        /// Number of bits in mantissa (including the implied bit)
        /// <para>Value = 53</para>
        /// </summary>
        public const int MantissaBits = 53;

        /// <summary>
        /// Maximum value 
        /// <para>Value = 1.7976931348623158e+308</para>
        /// </summary>
        public const double MaxValue = double.MaxValue; // 1.7976931348623158e+308

        /// <summary>
        /// Minimim (i.e. negative) value 
        /// <para>Value = -1.79769e+308</para>
        /// </summary>
        public const double MinValue = double.MinValue;

        /// <summary>
        /// Minimum normalized positive value 
        /// <para>Value = 2^-1022 = 2.2250738585072014e-308</para>
        /// </summary>
        public const double MinNormalValue = 2.2250738585072014e-308; // 2^-1022

        /// <summary>
        /// Minimum denormalized positive value 
        /// <para>Value = 2^-1074</para>
        /// </summary>
        public const double MinDenormValue = double.Epsilon; // 2^-1074

        /// <summary>
        /// Maximum decimal exponent (ie. exp in 1.0e+exp )
        /// <para>Value = 308</para>
        /// </summary>
        public const int MaxDecimalExponent = 308;

        /// <summary>
        /// Minimum decimal exponent (ie. exp in 1.0e+exp )
        /// <para>Value = -307</para>
        /// </summary>
        public const int MinDecimalExponent = -307;


        /// <summary>
        /// The maximum binary exponent such that 2^(exp-1) is normal
        /// <para>Value = 1024</para>
        /// </summary>
        public const int MaxExponent = 1024;

        /// <summary>
        /// The minimum binary exponent such that 2^(exp-1) is normal
        /// <para>Value = -1021</para>
        /// </summary>
        public const int MinExponent = -1021;

        /// <summary>
        /// Log( MaxValue )
        /// <para>Value = 709.5</para>
        /// </summary>
        public const double MaxLogValue = 709.5;

        /// <summary>
        /// Log( MinNormalizedValue )
        /// <para>Value = -708.0</para>
        /// </summary>
        public const double MinLogValue = -708.0;

        /// <summary>
        /// Maximum x for x^x
        /// <para>Value = 143.0</para>
        /// </summary>
        public const double MaxArgXPowX = 143.0;

        /// <summary>
        /// Maximum x for Γ(x)
        /// <para>Value = 171.0</para>
        /// </summary>
        public const double MaxGamma = 171.0;

        /// <summary>
        /// The exponent radix
        /// <para> Value = 2</para>
        /// </summary>
        public const int Radix = 2;


        //      #define _DBL_ROUNDS     1                       /* addition rounding: near */

        // Unfortunately the only way to make them truly constant :-(

        /// <summary>
        /// Each element = Power(MachineEpsilon, 1/n)
        /// </summary>
        public static class RootMachineEpsilon
        {
            // Not ideal, but this is the only way these can behave as true constants

            /// <summary>MachineEpsilon</summary>
            public const double _1 = 2.220446049250313080847263336181640625000e-16;

            /// <summary>MachineEpsilon^(1/2)</summary>
            public const double _2 = 1.490116119384765625000000000000000000000e-8;

            /// <summary>MachineEpsilon^(1/3)</summary>
            public const double _3 = 6.055454452393339060789892727936966935698e-6;

            /// <summary>MachineEpsilon^(1/4)</summary>
            public const double _4 = 0.0001220703125000000000000000000000000000000;

            /// <summary>MachineEpsilon^(1/5)</summary>
            public const double _5 = 0.0007400959797414053136461229498566627141998;

            /// <summary>MachineEpsilon^(1/6)</summary>
            public const double _6 = 0.002460783300575924149935958217340289747208;

            /// <summary>MachineEpsilon^(1/7)</summary>
            public const double _7 = 0.005804665191941204784374426901228877368452;

            /// <summary>MachineEpsilon^(1/8)</summary>
            public const double _8 = 0.01104854345603980506876319315788826623883;

            /// <summary>MachineEpsilon^(1/9)</summary>
            public const double _9 = 0.01822701624337682157646394741657293845965;

            /// <summary>MachineEpsilon^(1/10)</summary>
            public const double _10 = 0.02720470510300387934800843804624206559310;

            /// <summary>MachineEpsilon^(1/11)</summary>
            public const double _11 = 0.03775279513763897854028881907738273640637;

            /// <summary>MachineEpsilon^(1/12)</summary>
            public const double _12 = 0.04960628287400623358599080122725963313723;

            /// <summary>MachineEpsilon^(1/13)</summary>
            public const double _13 = 0.06250000000000000000000000000000000000000;

            /// <summary>MachineEpsilon^(1/14)</summary>
            public const double _14 = 0.076188353387779715056980687662035053579838245991248;

            /// <summary>MachineEpsilon^(1/15)</summary>
            public const double _15 = 0.090454327340023629372158355263346229039798221062208;

            /// <summary>MachineEpsilon^(1/16)</summary>
            public const double _16 = 0.10511205190671431787889068452915186188000428279460;

            /// <summary>MachineEpsilon^(1/17)</summary>
            public const double _17 = 0.12000583585684915290754539848978514770432884000696;

            /// <summary>MachineEpsilon^(1/18)</summary>
            public const double _18 = 0.13500746736153827123411635391107461409218334530821;

            /// <summary>MachineEpsilon^(1/19)</summary>
            public const double _19 = 0.15001283994726287950319688532791978492409255065997;

            /// <summary>MachineEpsilon^(1/20)</summary>
            public const double _20 = 0.16493848884661178242175024640370501662918362664918;

            /// <summary>MachineEpsilon^(1/21)</summary>
            public const double _21 = 0.17971833718062722098312581088618386639490548832824;

            /// <summary>MachineEpsilon^(1/22)</summary>
            public const double _22 = 0.19430078522136491494691180354193947182568105836439;

            /// <summary>MachineEpsilon^(1/23)</summary>
            public const double _23 = 0.20864618324226299057132000202231164317279676701859;

        };

    }


} // namespaces


