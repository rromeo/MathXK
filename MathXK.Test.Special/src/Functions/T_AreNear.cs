//  Copyright (c) 2013 Rocco Romeo
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

using Microsoft.VisualStudio.TestTools.UnitTesting;


namespace MathXK.Test.Functions
{

    [TestClass()]
    public class T_AreNear
    {

        private TestContext testContextInstance;

        /// <summary>
        ///Gets or sets the test context which provides
        ///information about and functionality for the current test run.
        ///</summary>
        public TestContext TestContext
        {
            get
            {
                return testContextInstance;
            }
            set
            {
                testContextInstance = value;
            }
        }

        #region Additional test attributes
        // 
        //You can use the following additional attributes as you write your tests:
        //
        //Use ClassInitialize to run code before running the first test in the class
        //[ClassInitialize()]
        //public static void MyClassInitialize(TestContext testContext)
        //{
        //}
        //
        //Use ClassCleanup to run code after all tests in a class have run
        //[ClassCleanup()]
        //public static void MyClassCleanup()
        //{
        //}
        //
        //Use TestInitialize to run code before running each test
        //[TestInitialize()]
        //public void MyTestInitialize()
        //{
        //}
        //
        //Use TestCleanup to run code after each test has run
        //[TestCleanup()]
        //public void MyTestCleanup()
        //{
        //}
        //
        #endregion


        struct Data
        {
            public readonly double X;
            public readonly double Y;
            public readonly int Ulps;
            public readonly bool Result;

            public Data(double x, double y, int ulps, bool result)
            {
                X = x;
                Y = y;
                Ulps = ulps;
                Result = result;
            }
        }

        private static readonly Data[] _Cases = {

            // border cases
            new Data(double.NaN , double.NaN, 10, false ),
            new Data(double.PositiveInfinity , double.PositiveInfinity, 1, true ),
            new Data(double.NegativeInfinity , double.NegativeInfinity, 1, true ),
            new Data(double.MaxValue , double.PositiveInfinity, 100, false ),
            new Data(double.MinValue , double.NegativeInfinity, 100, false ),
            new Data(double.PositiveInfinity , double.NegativeInfinity, 4*1024*1024, false ),

            // note that double.Epsilon is subnormal
            new Data(0.0 , 0.0, 1, true ),
            new Data(0.0 , double.Epsilon, 1, true ),
            new Data(0.0 , -double.Epsilon, 1, false ),
            new Data(-0.0 , double.Epsilon, 1, false ),
            new Data(-0.0 , -double.Epsilon, 1, true ),
            new Data(double.Epsilon , -double.Epsilon, 2, false ),

            new Data(0.1 , 0.0999999999999999, 10, true ),
            new Data(0.1 , 0.0999999999999950, 10, false ),
            new Data(-0.1 , -0.0999999999999999, 10, true ),
            new Data(-0.1 , -0.0999999999999950, 10, false ),

            new Data(1.0 , 0.9999999999999999, 10, true ),
            new Data(1.0 , 0.9999999999999950, 10, false ),
            new Data(-1.0 , -0.9999999999999999, 10, true ),
            new Data(-1.0 , -0.9999999999999950, 10, false ),


            new Data(2.0 , 1.9999999999999999, 10, true ),
            new Data(2.0 , 1.9999999999999950, 10, false ),
            new Data(-2.0 , -1.9999999999999999, 10, true ),
            new Data(-2.0 , -1.9999999999999950, 10, false ),

            new Data(0.01 , 0.009999999999999999, 10, true ),
            new Data(0.01 , 0.0099999999999950, 10, false ),
            new Data(-0.01 , -0.009999999999999999, 10, true ),
            new Data(-0.01 , -0.0099999999999950, 10, false )

            // TODO: large ulps
        };



        /// <summary>
        ///A test for AreNearUlps
        ///</summary>
        [TestMethod()]
        public void AreNearUlpsTest()
        {

            // Uses Bruce Dawson's technique
            // See http://www.cygnus-software.com/papers/comparingfloats/comparingfloats.htm
            // This only works for IEEE
            // Adapted from his tests


            foreach (var item in _Cases) {
                double x = item.X;
                double y = item.Y;
                int maxUlps = item.Ulps;
                bool result = Math2.AreNearUlps(x, y, maxUlps);
                Assert.AreEqual(item.Result, result, $"Wrong comparison - Value({x}) compared to({y})");

                // reverse
                result = Math2.AreNearUlps(y, x, maxUlps);
                Assert.AreEqual(item.Result, result, $"Wrong comparison - Value({y}) compared to({x})");

            }



        }
    }
}
