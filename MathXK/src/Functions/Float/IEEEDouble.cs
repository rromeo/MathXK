//  (C) Copyright Rocco Romeo 2013.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

using System;
using System.Diagnostics;

namespace MathXK
{

    /*
        Thanks to Thomas Finley's website for this explanation: 
        http://www.tfinley.net/notes/cps104/floating.html

        Single Precision        Double Precision        Object Represented 
        Exponent    Mantissa    Exponent    Mantissa 
            0        0              0            0          zero 
            0        nonzero        0           nonzero     +/- denormalized number 
            1-254    anything       1-2046      anything    +/- normalized number 
            255      0              2047        0           +/- infinity 
            255      nonzero        2047        nonzero     NaN (Not a Number) 


        see also: http://en.wikipedia.org/wiki/Double_precision_floating-point_format
        and: http://en.wikipedia.org/wiki/IEEE_754-1985
    */



    /// <summary>
    /// The purpose of this class is to facilitate the work with IEEE 754 doubles
    /// </summary>
    internal struct IEEEDouble
    {

        #region Constants
        /// <summary>
        /// The sign mask
        /// </summary>
        public const long SignMask = unchecked((long)0x8000000000000000UL);

        /// <summary>
        /// The exponent mask
        /// </summary>
        public const long ExponentMask = 0x7ff0000000000000L;

        /// <summary>
        /// The significand mask
        /// </summary>
        public const long SignificandMask = 0x000fffffffffffffL;

        /// <summary>
        /// The infinity exponent bits
        /// </summary>
        public const long InfinityExponentBits = 0x7ff0000000000000L; // (2047-1023) << 52

        /// <summary>
        /// Add the bias to create the exponent
        /// </summary>
        public const int ExponentBias = 1023; // 0x3ff

        /// <summary>
        /// The maximum value of the floating point exponent (mantissa * 2^exponent);
        /// </summary>
        public const int MaxBinaryExponent = 1023;


        /// <summary>
        /// The minimum physical value of the floating point exponent (mantissa * 2^exponent);
        /// </summary>
        public const int MinBinaryExponent = -1022;

        /// <summary>
        /// The minimum implied value of the floating point exponent (mantissa * 2^exponent),
        /// when mantissa is denormalized.
        /// </summary>
        public const int MinBinaryDenormExponent = -1074;

        /// <summary>
        /// The number of bits in the mantissa
        /// </summary>
        public const int MantissaBits = 52;
        #endregion



        #region Member

        private long _bits;

        public IEEEDouble(double x)
        {
            _bits = BitConverter.DoubleToInt64Bits(x);
        }

        public IEEEDouble(long x)
        {
            _bits = x;
        }

        public static implicit operator long(IEEEDouble x) { return x._bits; }
        public static implicit operator double(IEEEDouble x) { return BitConverter.Int64BitsToDouble(x._bits); }
        public static implicit operator IEEEDouble(double x) { return new IEEEDouble(x); }
        public static implicit operator IEEEDouble(long x) { return new IEEEDouble(x); }

        /// <summary>
        /// Returns the signbit of the specified argument
        /// </summary>
        public long Sign
        {
            get { return (_bits & SignMask); }
        }

        /// <summary>
        /// Returns the exponent bits of the specified argument
        /// </summary>
        public long Exponent
        {
            get { return (_bits & ExponentMask); }
        }

        /// <summary>
        /// Returns the significand bits  of the specified argument
        /// </summary>
        public long Significand
        {
            get { return (_bits & SignificandMask); }
        }


        /// <summary>
        /// The binary exponent value (that is, the exponent in  m*2^exponent)
        /// </summary>
        public int ExponentValue
        {
            get { return (int) ((_bits & ExponentMask) >> MantissaBits) - ExponentBias; }
            set {
                Debug.Assert(value + ExponentBias >= 0);
                _bits &= ~ExponentMask;
                _bits |= ((long)(value + ExponentBias) << MantissaBits);
            }
        }

        /// <summary>
        /// Returns 1 if the sign bit is set, otherwise 0;
        /// </summary>
        public int Signbit
        {
            get { return ((_bits & SignMask) == 0) ? 0 : 1; }
        }

        /// <summary>
        /// Returns true if the number represents NaN
        /// </summary>
        public bool IsNaN
        {
            // ExponentBits = 0x7ff, Significand !=0 
            get { return (_bits & (ExponentMask|SignificandMask)) > ExponentMask; }
        }

        /// <summary>
        /// Returns true if the number represents positive or negative infinity
        /// </summary>
        public bool IsInfinity
        {
            // ExponentBits = 0x7ff, Significand == 0 
            get { return (_bits & (ExponentMask | SignificandMask)) == ExponentMask; }
        }

        /// <summary>
        /// Returns true if |number| &lt; Infinity || !IsNaN(number)
        /// </summary>
        public bool IsFinite
        {
            // ExponentBits != 0x7ff 
            get { return ((_bits & ExponentMask) != ExponentMask); }
        }

        /// <summary>
        /// Returns true if the number is subnormal
        /// </summary>
        /// <returns><c>true</c> if the specified bits is subnormal; otherwise, <c>false</c>.</returns>
        public bool IsSubnormal
        {
            // ExponentBits = 0, Significand !=0 
            get { return ((_bits & ExponentMask) == 0) && ((_bits & SignificandMask) != 0); }
        }

        /// <summary>
        /// Flips the sign bit
        /// </summary>
        public void ChangeSign()
        {
            _bits ^= SignMask;
        }

        /// <summary>
        /// Sets this object's signbit to that of x
        /// </summary>
        /// <param name="x">The argument that contains the desired signbit</param>
        public void CopySign(double x)
        {
            long xbits = BitConverter.DoubleToInt64Bits(x) & SignMask;
            if (xbits == 0)
                _bits &= ~SignMask;
            else
                _bits |= SignMask;
        }

        public override string ToString()
        {
            return _bits.ToString("x8");
        }

        #endregion

        #region Statics


        /// <summary>
        /// Returns the value of x with the signbit set to the signbit of y
        /// </summary>
        /// <param name="x">The argument</param>
        /// <param name="y">The argument that contains the desired signbit</param>
        /// <returns>System.Double.</returns>
        public static double CopySign(double x, double y)
        {
            long xbits = BitConverter.DoubleToInt64Bits(x);
            long ybits = BitConverter.DoubleToInt64Bits(y);

            if ((xbits & SignMask) == (ybits & SignMask))
                return x;

            return BitConverter.Int64BitsToDouble(xbits ^ SignMask);
        }

        #endregion
    }



} // namespace