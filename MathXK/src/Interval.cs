//  (C) Copyright Rocco Romeo 2013.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

using System;

namespace MathXK
{

    // Would have preferred this class to be a generic/template
    // but needs double arithmetic

    /// <summary>
    /// Represents an interval with [lower, upper] limits
    /// </summary>
    public struct Interval : IEquatable<Interval>
    {
        private readonly double _upper;
        private readonly double _lower;


        /// <summary>
        /// Constructs an interval with [lower, upper] limits
        /// </summary>
        /// <param name="lower">The lower limit, or NaN if unknown</param>
        /// <param name="upper">The upper limit, or NaN if unknown</param>
        /// <exception cref="MathXK.Exceptions.DomainException">if lower &gt; upper</exception>
        public Interval(double lower, double upper)
        {
            if (lower > upper)
                throw new Exceptions.DomainException("Interval(lower: {0}, upper: {1}): lower and upper are interchanged", lower, upper);

            _lower = lower;
            _upper = upper;
        }


        /// <summary>
        /// Constructs an interval with [lower, upper] limits
        /// </summary>
        /// <param name="lower">The lower limit, or null if unknown</param>
        /// <param name="upper">The upper limit, or null if unknown</param>
        /// <exception cref="MathXK.Exceptions.DomainException">if lower &gt; upper</exception>
        public Interval(double? lower, double? upper)
        {
            if (lower.HasValue && upper.HasValue && lower.Value > upper.Value)
                throw new Exceptions.DomainException("Interval(lower: {0}, upper: {1}): lower and upper are interchanged", lower.Value, upper.Value);

            _lower = (lower.HasValue) ? lower.Value : double.NaN;
            _upper = (upper.HasValue) ? upper.Value : double.NaN;
        }


        /// <summary>
        /// The upper bound of the range
        /// </summary>
        public double Upper
        {
            get { return _upper; }
        }

        /// <summary>
        /// The lower bound of the range
        /// </summary>
        public double Lower
        {
            get { return _lower; }
        }

        /// <summary>
        /// The width of the range ( upper - lower )
        /// </summary>
        public double Width
        {
            get { return _upper - _lower; }
        }


        /// <summary>
        /// Returns an undefined (NaN) interval
        /// </summary>
        public static Interval NaN
        {
            get { return new Interval(double.NaN, double.NaN); }
        }


        /// <summary>
        /// Returns true if the argument is in [lower,upper] (upper &gt;= x &gt;= lower)
        /// </summary>
        /// <param name="x">Argument</param>
        /// <returns><c>true</c> if x is in [lower,upper]; otherwise, <c>false</c>.</returns>
        public bool IsInRange(double x)
        {
            return (x >= _lower && x <= _upper);
        }

        /// <summary>
        /// Returns x clipped to the boundary [lower,upper]
        /// </summary>
        /// <param name="x">Argument</param>
        /// <returns>Min(Max(lower,x),upper)</returns>
        public double Clip(double x)
        {
            if (x < _lower)
                return _lower;
            if (x > _upper)
                return _upper;
            return x;
        }



        /// <summary>
        /// Returns true if both ranges are equal
        /// </summary>
        /// <param name="range1">The first interval</param>
        /// <param name="range2">The second interval</param>
        /// <returns><c>true</c> if range1 == range2; if either range contains a NaN, the result will be <c>false</c></returns>
        public static bool operator ==(Interval range1, Interval range2)
        {
            // we treat infinities as equal, but NaNs are still not equal
            return ((range1._lower == range2._lower) && (range1._upper == range2._upper));
        }

        /// <summary>
        /// Returns true if both ranges are not equal
        /// </summary>
        /// <param name="range1">The first interval</param>
        /// <param name="range2">The second interval</param>
        /// <returns><c>true</c> if range1 != range2; if either range contains a NaN, the result will be <c>true</c></returns>
        public static bool operator !=(Interval range1, Interval range2)
        {
            return !(range1 == range2);
        }

        /// <summary>
        /// Unary plus. Returns +range
        /// </summary>
        /// <param name="range">The argument</param>
        /// <returns>The argument</returns>
        public static Interval operator +(Interval range) { return range; }


        /// <summary>
        /// Unary minus. Returns -range
        /// </summary>
        /// <param name="range">The argument</param>
        /// <returns>The negated argument</returns>
        public static Interval operator -(Interval range)
        {
            return new Interval(-range._upper, -range._lower);
        }

        /// <summary>
        /// Returns the interval [lower+x, upper+x]
        /// </summary>
        /// <param name="range">The interval</param>
        /// <param name="x">The value to add to the interval</param>
        public static Interval operator +(Interval range, double x)
        {
            return new Interval(range._lower + x, range._upper + x);
        }

        /// <summary>
        /// Returns the interval [lower+x, upper+x]
        /// </summary>
        /// <param name="x">The value to add to the interval</param>
        /// <param name="range">The interval</param>
        public static Interval operator +(double x, Interval range)
        {
            return new Interval(range._lower + x, range._upper + x);
        }

        /// <summary>
        /// Returns the interval [lower-x, upper-x]
        /// </summary>
        /// <param name="range">The interval</param>
        /// <param name="x">The value to subtract from the interval</param>
        public static Interval operator -(Interval range, double x)
        {
            return new Interval(range._lower - x, range._upper - x);
        }


        /// <summary>
        /// Returns the scaled interval
        /// </summary>
        /// <param name="range">The interval</param>
        /// <param name="x">The value to multiply to the interval</param>
        public static Interval operator *(Interval range, double x)
        {
            if (x < 0)
                return new Interval(x * range._upper, x * range._lower);
            return new Interval(range._lower * x, range._upper * x);
        }

        /// <summary>
        /// Returns the scaled interval
        /// </summary>
        /// <param name="x">The value to multiply to the interval</param>
        /// <param name="range">The interval</param>
        public static Interval operator *(double x, Interval range) { return range * x; }


        /// <summary>
        /// Returns the scaled interval
        /// </summary>
        /// <param name="range">The interval</param>
        /// <param name="x">The value to divide the interval by</param>
        public static Interval operator /(Interval range, double x)
        {
            if (x == 0) {
                Policies.ReportDomainError("Requires x != 0");
                return new Interval(double.NaN, double.NaN);
            }

            if (x < 0)
                return new Interval(range._upper / x, range._lower / x);
            return new Interval(range._lower / x, range._upper / x);
        }

        /// <summary>
        /// Returns a string representation of the interval "[lower, upper]"
        /// </summary>
        /// <returns>A <see cref="System.String" /> that represents this instance.</returns>
        public override string ToString()
        {
            return $"[{_lower}, {_upper}]";
        }

        /// <summary>
        /// Gets the hash code
        /// </summary>
        /// <returns>A hash code for this instance, suitable for use in hashing algorithms and data structures like a hash table.</returns>
        /// <see cref="System.Object.GetHashCode" />
        public override int GetHashCode()
        {
            return _lower.GetHashCode() ^ _upper.GetHashCode();
        }

        /// <summary>
        /// Returns true if this object equals the specified object
        /// </summary>
        /// <param name="obj">Another object to compare to.</param>
        /// <returns><c>true</c> if the specified <see cref="System.Object" /> is equal to this instance; otherwise, <c>false</c>.</returns>
        public override bool Equals(object obj)
        {
            if (obj == null || GetType() != obj.GetType())
                return false;

            return this.Equals((Interval)obj);
        }

        /// <summary>
        /// Returns true if this object equals the specified object. Unlike operator==, NaNs are treated as equal.
        /// </summary>
        /// <param name="other">Another object to compare to.</param>
        /// <returns><c>true</c> if the specified <see cref="Interval" /> is equal to this instance; otherwise, <c>false</c>.</returns>
        public bool Equals(Interval other)
        {
            return this._lower.Equals(other._lower) && this._upper.Equals(other._upper);
        }
    };

} // namespaces