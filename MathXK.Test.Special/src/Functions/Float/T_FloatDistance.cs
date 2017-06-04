//  Copyright (c) 2013 Rocco Romeo
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

//  Test Data:
//      Copyright (c) 2006-7 John Maddock, Boost Software License v1.0

using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace MathXK.Test.Functions
{
    [TestClass]
    public class T_FloatDistance
    {

        private const double MachineEpsilon = 1.0 / (1L << 51); // 2^-52 ~= 1.11e-16

        struct Data
        {
            public readonly double X;
            public readonly double Result;
            public readonly int Distance;

            public Data(double x, int distance, double result)
            {
                X = x;
                Distance = distance;
                Result = result;
            }

            public Data(ulong x, int distance, ulong result)
            {
                X = BitConverter.Int64BitsToDouble((long)x);
                Distance = distance;
                Result = BitConverter.Int64BitsToDouble((long)result);
            }
        }

        // note: routines depend on the reverse test being true

        private static readonly Data[] _FloatData = {

            // cannot use NaN because equality test does not work
            // (i.e NaN == NaN fails)
            // new Data(double.NaN, 1, double.NaN),


            // denorms

            new Data(double.Epsilon,        1, double.Epsilon*2.0),
            new Data(0.0,                    1, double.Epsilon),
            new Data(-double.Epsilon,        1, 0.0),            
            new Data(-double.Epsilon*2.0,    1, -double.Epsilon),

            new Data(double.Epsilon,        2, double.Epsilon*3.0),
            new Data(0.0,                    2, double.Epsilon*2.0),
            new Data(-double.Epsilon,        2, double.Epsilon),            
            new Data(-double.Epsilon*2.0,    2, 0.0),


            new Data(0x000fffffffffffffUL,    1, 0x0010000000000000UL), // max denorm, min norm


            new Data(0x8010000000000000UL,    1, 0x800fffffffffffffUL), // -max denorm, -max denorm+1eps

            // new Data(double.MaxValue, 1, double.PositiveInfinity),

            new Data(0x800fffffffffffffUL,    1, 0x800ffffffffffffeUL),


            new Data(0x3ff0000000000000UL, 1, 0x3ff0000000000001UL), // 1.0, 1.0+
            new Data(0x3ff0000000000001UL, 1, 0x3ff0000000000002UL), // 1.0

            new Data(0x7feffffffffffffeUL, 1, 0x7fefffffffffffffUL) // MaxValue-1eps, MaxValue

        
        };


        static void TestValue(double val)
        {
            double upper = double.MaxValue;
            double lower = double.MinValue;

            Assert.IsTrue(Math2.FloatNext(val) > val);
            Assert.AreEqual(1.0, Math2.FloatDistance(Math2.FloatNext(val), val));

            Assert.IsTrue(Math2.FloatPrior(val) < val);
            Assert.AreEqual(-1.0, Math2.FloatDistance(Math2.FloatPrior(val), val));

            Assert.IsTrue(Math2.Nextafter(val, upper) > val);
            Assert.AreEqual(1.0, Math2.FloatDistance(Math2.Nextafter(val, upper), val));

            Assert.IsTrue(Math2.Nextafter(val, lower) < val);
            Assert.AreEqual(-1.0, Math2.FloatDistance(Math2.Nextafter(val, lower), val));

            Assert.AreEqual(2.0, Math2.FloatDistance(Math2.FloatNext(Math2.FloatNext(val)), val));
            Assert.AreEqual(-2.0, Math2.FloatDistance(Math2.FloatPrior(Math2.FloatPrior(val)), val));
            Assert.AreEqual(-4.0, Math2.FloatDistance(Math2.FloatPrior(Math2.FloatPrior(val)), Math2.FloatNext(Math2.FloatNext(val))));
            Assert.AreEqual(0.0, Math2.FloatDistance(Math2.FloatPrior(Math2.FloatNext(val)), val));
            Assert.AreEqual(0.0, Math2.FloatDistance(Math2.FloatNext(Math2.FloatPrior(val)), val));
            Assert.AreEqual(Math2.FloatPrior(Math2.FloatNext(val)), val);
            Assert.AreEqual(Math2.FloatNext(Math2.FloatPrior(val)), val);

            Assert.AreEqual(4.0, Math2.FloatDistance(Math2.FloatAdvance(val, 4), val));
            Assert.AreEqual(-4.0, Math2.FloatDistance(Math2.FloatAdvance(val, -4), val));
            Assert.AreEqual(4.0, Math2.FloatDistance(Math2.FloatAdvance(Math2.FloatNext(Math2.FloatNext(val)), 4), Math2.FloatNext(Math2.FloatNext(val))));
            Assert.AreEqual(-4.0, Math2.FloatDistance(Math2.FloatAdvance(Math2.FloatNext(Math2.FloatNext(val)), -4), Math2.FloatNext(Math2.FloatNext(val))));
        }

        static readonly int[] _Primes = {
            11,     13,     17,     19,     23,     29, 
            31,     37,     41,     43,     47,     53,     59,     61,     67,     71, 
            73,     79,     83,     89,     97,    101,    103,    107,    109,    113, 
            127,    131,    137,    139,    149,    151,    157,    163,    167,    173, 
            179,    181,    191,    193,    197,    199,    211,    223,    227,    229, 
            233,    239,    241,    251,    257,    263,    269,    271,    277,    281, 
            283,    293,    307,    311,    313,    317,    331,    337,    347,    349, 
            353,    359,    367,    373,    379,    383,    389,    397,    401,    409, 
            419,    421,    431,    433,    439,    443,    449,    457,    461,    463, 
        };

        void test_values(double val)
        {
            const double a = 1.3456724e22;
            const double b = 1.3456724e-22;
            const double z = 0;
            const double one = 1;
            const double two = 2;


            TestValue(a);
            TestValue(-a);
            TestValue(b);
            TestValue(-b);
            TestValue(MachineEpsilon);
            TestValue(-MachineEpsilon);
            // TestValue(double.MinValue);
            // TestValue(-double.MinValue);
            TestValue(z);
            TestValue(-z);
            TestValue(one);
            TestValue(-one);
            TestValue(two);
            TestValue(-two);
            TestValue(double.Epsilon);
            TestValue(-double.Epsilon);
            TestValue(2 * double.Epsilon);
            TestValue(-2 * double.Epsilon);



            for (int i = 0; i < _Primes.Length; ++i) {
                double v1 = val;
                double v2 = val;
                for (int j = 0; j < _Primes[i]; ++j) {
                    v1 = Math2.FloatNext(v1);
                    v2 = Math2.FloatPrior(v2);
                }
                Assert.AreEqual(_Primes[i], Math2.FloatDistance(v1, val));
                Assert.AreEqual(-_Primes[i], Math2.FloatDistance(v2, val));
                Assert.AreEqual(v1, Math2.FloatAdvance(val, _Primes[i]));
                Assert.AreEqual(v2, Math2.FloatAdvance(val, -_Primes[i]));
            }
        }

        [TestMethod]
        public void FloatTestAll()
        {
            test_values(1.0);
        }

        [TestMethod]
        public void FloatNextTest()
        {
            // edges

            Assert2.AreEqual<Exceptions.DomainException>(double.NaN, () => Math2.FloatNext(double.NaN));
            Assert2.AreEqual<Exceptions.DomainException>(double.NaN, () => Math2.FloatNext(double.PositiveInfinity));
            Assert2.AreEqual<Exceptions.DomainException>(double.NaN, () => Math2.FloatNext(double.NegativeInfinity));
            Assert.AreEqual(double.PositiveInfinity, Math2.FloatNext(double.MaxValue));
            Assert.AreEqual(double.Epsilon, Math2.FloatNext(-0.0), 1); // -0.0, double.Epsilon


            foreach (var item in _FloatData) {

                if (item.Distance != 1)
                    continue;

                double x = item.X;
                double expected = item.Result;
                double actual = Math2.FloatNext(x);

                Assert2.AreNear(expected, actual, 0,
                        () => { return string.Format("FloatNext({0:E16}) = {1:E16}; expected = {2:E16}; ulps = {3}", x, actual, expected, 0); });

            }

        }

        [TestMethod]
        public void FloatPriorTest()
        {
            // edges


            Assert2.AreEqual<Exceptions.DomainException>(double.NaN, () => Math2.FloatPrior(double.NaN));
            Assert2.AreEqual<Exceptions.DomainException>(double.NaN, () => Math2.FloatPrior(double.PositiveInfinity));
            Assert2.AreEqual<Exceptions.DomainException>(double.NaN, () => Math2.FloatPrior(double.NegativeInfinity));
            Assert.AreEqual(double.NegativeInfinity, Math2.FloatPrior(double.MinValue));
            Assert.AreEqual(-double.Epsilon, Math2.FloatPrior(-0.0), 1); // -0.0, double.Epsilon

            foreach (var item in _FloatData) {
                // swap result and x

                if (item.Distance != 1)
                    continue;

                double x = item.Result;
                double expected = item.X;
                double actual = Math2.FloatPrior(x);

                Assert2.AreNear(expected, actual, 0,
                        () => { return string.Format("FloatPrior({0:E16}) = {1:E16}; got {2:E16}; ulps = {3}", x, expected, actual, 0); });

            }
        }

        [TestMethod]
        public void FloatAdvanceTest()
        {

            // test non-symettrical behaviors - NaN, -0.0, Infinity

            for (int i = 1; i < 255; i++) {

                Assert2.AreEqual<Exceptions.DomainException>(double.NaN, () => Math2.FloatAdvance(double.NaN, i));
                Assert2.AreEqual<Exceptions.DomainException>(double.NaN, () => Math2.FloatAdvance(double.NaN, -i));
                Assert2.AreEqual<Exceptions.DomainException>(double.NaN, () => Math2.FloatAdvance(double.PositiveInfinity, i));
                Assert2.AreEqual<Exceptions.DomainException>(double.NaN, () => Math2.FloatAdvance(double.PositiveInfinity, -i));
                Assert2.AreEqual<Exceptions.DomainException>(double.NaN, () => Math2.FloatAdvance(double.NegativeInfinity, i));
                Assert2.AreEqual<Exceptions.DomainException>(double.NaN, () => Math2.FloatAdvance(double.NegativeInfinity, -i));

                Assert.AreEqual(double.PositiveInfinity, Math2.FloatAdvance(double.MaxValue, i));
                Assert.AreEqual(double.NegativeInfinity, Math2.FloatAdvance(double.MinValue, -i));

                Assert.AreEqual(double.Epsilon * (double)i, Math2.FloatAdvance(-0.0, i)); // -0.0, double.Epsilon
                Assert.AreEqual(-double.Epsilon * (double)i, Math2.FloatAdvance(-0.0, -i)); // -0.0, -double.Epsilon
            }



            foreach (var item in _FloatData) {
                double x = item.X;
                int distance = item.Distance;
                double expected = item.Result;
                double actual = Math2.FloatAdvance(x, distance);

                Assert2.AreNear(expected, actual, 0,
                        () => { return string.Format("FloatAdvance({0:E16},{1}) = {2:E16}; got = {3:E16}; ulps = {4}", x, distance, expected, actual, 0); });

                // check to see that the reverse is true

                x = item.Result;
                distance = -item.Distance;
                expected = item.X;
                actual = Math2.FloatAdvance(x, distance);

                Assert2.AreNear(expected, actual, 0,
                        () => { return string.Format("FloatAdvance({0:E16},{1}) = {2:E16}; got = {3:E16}; ulps = {4}", x, distance, expected, actual, 0); });



            }


        }

        [TestMethod]
        public void FloatDistanceTest()
        {

            // test non-symettrical behaviors - NaN, -0.0, Infinity

            Assert2.AreEqual<Exceptions.DomainException>(double.NaN, () => Math2.FloatDistance(double.NaN, double.NaN));
            Assert2.AreEqual<Exceptions.DomainException>(double.NaN, () => Math2.FloatDistance(double.PositiveInfinity, double.NaN));
            Assert2.AreEqual<Exceptions.DomainException>(double.NaN, () => Math2.FloatDistance(double.NegativeInfinity, double.NaN));
            Assert2.AreEqual<Exceptions.DomainException>(double.NaN, () => Math2.FloatDistance(double.MaxValue, double.NaN));
            Assert2.AreEqual<Exceptions.DomainException>(double.NaN, () => Math2.FloatDistance(double.MinValue, double.NaN));
            Assert2.AreEqual<Exceptions.DomainException>(double.NaN, () => Math2.FloatDistance(0.0, double.NaN));
            Assert2.AreEqual<Exceptions.DomainException>(double.NaN, () => Math2.FloatDistance(-0.0, double.NaN));

            Assert2.AreEqual<Exceptions.DomainException>(double.NaN, () => Math2.FloatDistance(double.NaN, double.NaN));
            Assert2.AreEqual<Exceptions.DomainException>(double.NaN, () => Math2.FloatDistance(double.NaN, double.PositiveInfinity));
            Assert2.AreEqual<Exceptions.DomainException>(double.NaN, () => Math2.FloatDistance(double.NaN, double.NegativeInfinity));
            Assert2.AreEqual<Exceptions.DomainException>(double.NaN, () => Math2.FloatDistance(double.NaN, double.MaxValue));
            Assert2.AreEqual<Exceptions.DomainException>(double.NaN, () => Math2.FloatDistance(double.NaN, double.MinValue));
            Assert2.AreEqual<Exceptions.DomainException>(double.NaN, () => Math2.FloatDistance(double.NaN, 0.0));
            Assert2.AreEqual<Exceptions.DomainException>(double.NaN, () => Math2.FloatDistance(double.NaN, -0.0));

            Assert2.AreEqual<Exceptions.DomainException>(double.NaN, () => Math2.FloatDistance(double.PositiveInfinity, double.MaxValue));
            Assert2.AreEqual<Exceptions.DomainException>(double.NaN, () => Math2.FloatDistance(double.MaxValue, double.PositiveInfinity));
            Assert2.AreEqual<Exceptions.DomainException>(double.NaN, () => Math2.FloatDistance(double.NegativeInfinity, double.MinValue));
            Assert2.AreEqual<Exceptions.DomainException>(double.NaN, () => Math2.FloatDistance(double.MinValue, double.NegativeInfinity));

            Assert.AreEqual(1.0, Math2.FloatDistance(double.Epsilon, -0.0));
            Assert.AreEqual(-1.0, Math2.FloatDistance(-double.Epsilon, -0.0));
            Assert.AreEqual(-1.0, Math2.FloatDistance(-0.0, double.Epsilon));
            Assert.AreEqual(1.0, Math2.FloatDistance(-0.0, -double.Epsilon));

            // large distance

            double ZeroToMax = (double)BitConverter.DoubleToInt64Bits(double.MaxValue);
            double ZeroToMin = ZeroToMax + 1.0;
            double MaxToMin = ZeroToMax + ZeroToMin;

            Assert.AreEqual(ZeroToMax, Math2.FloatDistance(double.MaxValue, 0));
            Assert.AreEqual(-ZeroToMax, Math2.FloatDistance(0, double.MaxValue));
            Assert.AreEqual(-ZeroToMin, Math2.FloatDistance(double.MinValue, 0));
            Assert.AreEqual(ZeroToMin, Math2.FloatDistance(0, double.MinValue));
            Assert.AreEqual(MaxToMin, Math2.FloatDistance(double.MaxValue, double.MinValue));
            Assert.AreEqual(-MaxToMin, Math2.FloatDistance(double.MinValue, double.MaxValue));


            // now test symmetrical behavior

            foreach (var item in _FloatData) {
                double x = item.X;
                double y = item.Result;

                double expected = item.Distance;
                double actual = Math2.FloatDistance(y, x);

                Assert2.AreNear(expected, actual, 0,
                        () => { return string.Format("FloatDistance({0:E16},{1:E16}) = {2:E16}; got = {3:E16}; ulps = {4}", y, x, expected, actual, 0); });

                // Reverse parameters should yield negative answer

                expected = -item.Distance;
                actual = Math2.FloatDistance(x, y);

                Assert2.AreNear(expected, actual, 0,
                        () => { return string.Format("FloatDistance({0:E16},{1:E16}) = {2:E16}; got = {3:E16}; ulps = {4}", x, y, expected, actual, 0); });
            }

        }
    }
}
