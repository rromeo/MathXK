//  Copyright (c) 2013 Rocco Romeo
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

//  Test Data:
//      Copyright (c) 2006-7 John Maddock, Boost Software License v1.0

using System;
using System.Numerics;
using System.Collections.Generic;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace MathXK.Test.Functions
{

    using TestCase1 = TestCase<int>;
    using TestCase2 = TestCase<double, int>;


    [TestClass]
    public class T_Factorial
    {


        private static readonly TestCase2[] _RisingFactorialData = {


            // added: n = 0 
           new TestCase2( 0, 0, 1),
           new TestCase2( 100, 0, 1),
            
           // x = 0
           new TestCase2( 0, 1, 0),
           new TestCase2( 0, -1, -1),

           // x = 1
           new TestCase2( 1, 1, 1),
           new TestCase2( -1, 2, 0),


           new TestCase2(3, 4, 360),
           new TestCase2(7, -4, 0.00277777777777777777777777777777777777777777777777777777777778),
           new TestCase2(120.5, 8, 5.58187566784927180664062500e16),
           new TestCase2(120.5, -4, 5.15881498170104646868208445266116850161120996179812063177241e-9),
           new TestCase2(5000.25, 8, 3.92974581976666067544013393509103775024414062500000e29),
           new TestCase2(5000.25, -7, 1.28674092710208810281923019294164707555099052561945725535047e-26),
           new TestCase2(30.25, 21, 3.93286957998925490693364184100209193343633629069699964020401e33), 
           new TestCase2(30.25, -21, 3.35010902064291983728782493133164809108646650368560147505884e-27),
           new TestCase2(-30.25, 21,-9.76168312768123676601980433377916854311706629232503473758698e26), 
           new TestCase2(-30.25, -21,-1.50079704000923674318934280259377728203516775215430875839823e-34), 
           new TestCase2(-30.25, 5, -1.78799177197265625000000e7),
           new TestCase2(-30.25, -5,-2.47177487004482195012362027432181137141899692171397467859150e-8),
           new TestCase2(-30.25, 6, 4.5146792242309570312500000e8),
           new TestCase2(-30.25, -6, 6.81868929667537089689274558433603136943171564610751635473516e-10),
           new TestCase2(-3, 6, 0),
           new TestCase2(-3.25, 6, 2.99926757812500),
           new TestCase2(-5.25, 6, 50.987548828125000000000000),
           new TestCase2(-5.25, 13, 127230.91046623885631561279296875000),
           new TestCase2(-3.25, -6, 0.0000129609865918182348202632178291407500332449622510474437452125),
           new TestCase2(-5.25, -6, 2.50789821857946332294524052303699065683926911849535903362649e-6),
           new TestCase2(-5.25, -13,-1.38984989447269128946284683518361786049649013886981662962096e-14),

           new TestCase2(-3.5, 3, -13.125),

           // Added x < -n
           new TestCase2(0.5, -1, -2),
           new TestCase2(1.5, -2, -4),

           new TestCase2(40.5, -41, -2.75643016796662963097096639854040835565778207128865739e-47),
           new TestCase2(-Math.Pow(2, -52), 1, -2.220446049250313080847263336181640625e-16),
           new TestCase2(-Math.Pow(2, -52), -1, -0.9999999999999997779553950749687412190802426950628280887828),

           // small cases
           new TestCase2(-Math.Pow(2, -53), 1, -1.1102230246251565404236316680908203125e-16),
           new TestCase2(-Math.Pow(2, -53), -1, -0.9999999999999998889776975374843582835884772692260598527273),

           new TestCase2(-1-Math.Pow(2, -51), 89, 9.360322731930654460042192362644943510156327796891048366e116),
           new TestCase2(-1-Math.Pow(2, -26), 89, 3.140802936573047366577676836664604784483364254405888260e124),
           new TestCase2(-1-Math.Pow(2, -13), 89, 2.5716745174369217373494582134530875678121794192478503659e128),
           new TestCase2(-1-Math.Pow(2, -12), 89, 5.1408072827703783910055602044257354416556253465987239465e128),
           new TestCase2(-1.25, 89, 1.7580834079618036144965982333807443862834744883311477160e131),

           

                                                                   };


        // results from fuctions.wolfram.com

        private static readonly TestCase2[] _FallingFactorialData = {
           new TestCase2( 0, 0, 1),
           new TestCase2( 100, 0, 1),

           new TestCase2( 1, 1, 1),
           new TestCase2( 1, 2, 0),
           new TestCase2( 1, 3, 0),

           new TestCase2( 30.25, 0, 1),
           new TestCase2( 30.25, 1, 30.25),
           new TestCase2( 30.25, 2, 884.8125),
           new TestCase2( 30.25, 5, 1.78799177197265625e7),
           new TestCase2( 30.25, 22, 9.0295568931051440085683190087457309023832863204006571e27),
           new TestCase2( 100.5, 6, 8.85035464185234375e11),
           new TestCase2( 30.75, 30, 3.7803017610368795281377414370846070478617117090045689e33),
           new TestCase2(-12.0, 6,8910720),
           new TestCase2( -12, 5, -524160),
           new TestCase2(-3.0, 6, 20160),
           new TestCase2( -3.0, 5, -2520),
           new TestCase2( 3.0, 6, 0),
           new TestCase2( 3.0, 5, 0),
           new TestCase2( 3.25, 4, 2.28515625),
           new TestCase2( 3.25, 5, -1.7138671875),
           new TestCase2( 3.25, 6, 2.999267578125),
           new TestCase2( 3.25, 7, -8.24798583984375),
           new TestCase2( 8.25, 12, -68796.780131757259368896484375),


           // small cases
           new TestCase2(Math.Pow(2, -53), 1, 1.1102230246251565404236316680908203125e-16),
           new TestCase2(-Math.Pow(2, -53), 1, -1.1102230246251565404236316680908203125e-16),

           new TestCase2(Math.Pow(2, -52), 173, 4.7396555142284711816788142871761393464147030759327225634e295),
           new TestCase2(-Math.Pow(2, -52), 173, -4.73965551422848323733249323035900320877988995463670310e295),



        };

        /// <summary>
        /// Computes RisingFactorial(x:[-170...170], n:[-170...170]) using big integers
        /// </summary>
        /// <returns></returns>
        public TestCase2[] GetRisingFactorialIntData()
        {

            double DoubleMinNormal = Math.Pow(2, -1022);

            int Len = Math2.FactorialTable.Length - 1;
            var data = new List<TestCase2>(4 * Len * Len);

            for (int x = -Len; x < Len; x++) {
                for (int n = -Len; n < Len; n++) {

                    // divide by zero
                    if (x > 0 && x + n <= 0)
                        continue;

                    // our big integer routines don't work well with denorms
                    double result = BigInt.FactorialRising((BigInteger)x, n);
                    if (Math.Abs(result) < DoubleMinNormal)
                        result = 0;

                    data.Add(new TestCase2(x, n, result));
                }
            }

            return data.ToArray();

        }

        public TestCase2[] GetFallingFactorialIntData()
        {

            int Len = Math2.FactorialTable.Length - 1;
            TestCase2[] data = new TestCase2[2 * Len * Len];

            int index = 0;
            for (int x = -Len; x < Len; x++) {
                for (int n = 0; n < Len; n++) {
                    double result = BigInt.FactorialFalling((BigInteger)x, n);
                    data[index++] = new TestCase2(x, n, result);
                }
            }

            return data;

        }

        /// <summary>
        ///A test for FactorialTable
        ///</summary>
        [TestMethod()]
        public void FactorialTableTest()
        {

            TestCase1[] factorial_data = new TestCase1[Math2.FactorialTable.Length];

            // step 1: compute factorials using big integer to avoid trunctaion errors


            factorial_data[0] = new TestCase1(0, 1);
            for (int i = 1; i < Math2.FactorialTable.Length; i++) {
                factorial_data[i] = new TestCase1(i, BigInt.Factorial(i));
            }

            // Run the test case

            TestCaseSet<int>[] testCases = {
                new TestCaseSet<int>(factorial_data,    "Factorial: Complete table",    1),
            };

            NumericFunctionTest.RunSet((int i) => Math2.FactorialTable[i], "Factorial Table", testCases);


        }


        /// <summary>
        ///A test for Factorial
        ///</summary>
        [TestMethod()]
        public void FactorialTest()
        {

            TestCase1[] factorial_data = new TestCase1[Math2.FactorialTable.Length];

            // step 1: compute factorials using big integer to avoid trunctaion errors


            factorial_data[0] = new TestCase1(0, 1);
            for (int i = 1; i < Math2.FactorialTable.Length; i++) {
                factorial_data[i] = new TestCase1(i, BigInt.Factorial(i));
            }

            // Run the test case

            TestCaseSet<int>[] testCases = {
                new TestCaseSet<int>(factorial_data,    "Factorial: Complete table",    1),
            };

            NumericFunctionTest.RunSet(Math2.Factorial, "Factorial", testCases);

            Assert.AreEqual(double.PositiveInfinity, Math2.Factorial(Math2.FactorialTable.Length));
            Assert.AreEqual(double.PositiveInfinity, Math2.Factorial(int.MaxValue));

            // TODO: negative
        }

        /// <summary>
        ///A test for Factorial2
        ///</summary>
        [TestMethod()]
        public void Factorial2Test()
        {

            TestCase1[] factorial2_data = new TestCase1[250];

            // step 1: compute factorials using big integer to avoid trunctaion errors


            factorial2_data[0] = new TestCase1(-1, 1);
            factorial2_data[1] = new TestCase1(0, 1);
            for (int i = 2; i < factorial2_data.Length; i++) {
                factorial2_data[i] = new TestCase1(i - 1, BigInt.Factorial2(i - 1));
            }

            // Run the test case

            TestCaseSet<int>[] testCases = {
                new TestCaseSet<int>(factorial2_data,    "Factorial: Complete table",    10),
            };

            NumericFunctionTest.RunSet(Math2.Factorial2, "Factorial2", testCases);

        }


        private static double RisingFactorialNoDenorm(double x, int n)
        {
            double DoubleMinNormal = Math.Pow(2, -1022);
            double result = Math2.FactorialRising(x, n);
            if (Math.Abs(result) < DoubleMinNormal)
                result = 0;
            return result;
        }


        /// <summary>
        ///A test for Rising Factorial
        ///</summary>
        [TestMethod()]
        public void RisingFactorialTest()
        {

            Assert2.AreEqual<Exceptions.PoleException>(double.NaN, () => Math2.FactorialRising(1, -1));


            // Our BigInteger routines don't return denorms


            TestCase2[] intdata = GetRisingFactorialIntData();

            TestCaseSet<double, int>[] testCasesND = {
                new TestCaseSet<double,int>(intdata, "Integer Data", 280),
                
            };

            NumericFunctionTest.RunSet(RisingFactorialNoDenorm, "FactorialRising", testCasesND);

            TestCaseSet<double, int>[] testCases = {
                new TestCaseSet<double,int>(_RisingFactorialData, "Spots", 60),
                
            };


            NumericFunctionTest.RunSet(Math2.FactorialRising, "FactorialRising", testCases);

        }



        /// <summary>
        ///A test for Falling Factorial
        ///</summary>
        [TestMethod()]
        public void FallingFactorialTest()
        {

            TestCase2[] intdata = GetFallingFactorialIntData();

            TestCaseSet<double, int>[] testCases = {
                new TestCaseSet<double,int>(_FallingFactorialData, "Spots", 100),
                new TestCaseSet<double,int>(intdata, "Integer Data", 250),
            };

            NumericFunctionTest.RunSet(Math2.FactorialFalling, "FactorialFalling", testCases);


        }


    }
}
