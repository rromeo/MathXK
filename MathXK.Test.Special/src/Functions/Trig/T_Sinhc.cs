//  Copyright (c) 2013 Rocco Romeo
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

using System;
using System.Collections.Generic;
using System.Linq;
using Microsoft.VisualStudio.TestTools.UnitTesting;


namespace MathXK.Test.Functions
{

    using TestCase1 = TestCase<double>;

    [TestClass]
    public class T_Sinhc
    {


        #region Data

        // added to boost tests
        // function values calculated on http://functions.wolfram.com/
        private static readonly TestCase1[] _SinhcData = { 
            // range1 >= 2^-13
            new TestCase1(2, 1.8134302039235093838341069914006308524431710061605678606547422374),

            // note this case was from http://keisan.casio.com - wolfram incorrect
            new TestCase1(1.5, 1.4195196367298783312229249964517544325545260796695), 

            new TestCase1(1, 1.175201193643801456882381850595600815155717981334095870229565413013),
            new TestCase1(1.0/32, 1.0001627683641374254502179525624697694363611639241748061844803212),
            new TestCase1(1.0/64, 1.0000406906008749270894797563306041343633691296582229589660292671),
            new TestCase1(1.0/256, 1.000002543133450672735352570661278414388537404837808446669754768),
            new TestCase1(1.0/1024, 1.00000015894572698016435353665386298708456878244813511254845238),
            new TestCase1(1.0/4096, 1.00000000993410749217105153201969762337571516585781459646172047),

            // 2^-13 <= Range2 <= 2^-26
            new TestCase1(Math.Pow(2,-16), 1.0000000000388051072760966890562464645011704817235098284832832), // 2^-16
            new TestCase1(Math.Pow(2,-18), 1.000000000002425319204729573286413551237112723890170099295167), // 2^-18
            new TestCase1(Math.Pow(2,-20), 1.00000000000015158245029549493282427778579655971898637707276), // 2^-20
            new TestCase1(Math.Pow(2,-22), 1.00000000000000947390314346802940473388844494773404420724708), // 2^-22
            new TestCase1(Math.Pow(2,-24), 1.0000000000000005921189464667502600740575857833537853635822), // 2^-24
            new TestCase1(Math.Pow(2,-26), 1.0000000000000000370074341541718850916527770723043279074269), // 2^-26

            // Range3 < 2^-26
            new TestCase1(Math.Pow(2,-28), 1.000000000000000002312964634635742794154174262178572195785),
            new TestCase1(Math.Pow(2,-30), 1.00000000000000000014456028966473392454059634332037776157),
            new TestCase1(Math.Pow(2,-40), 1.00000000000000000000000013786343542550461247856812104689),
            new TestCase1(Math.Pow(2,-52), 1.000000000000000000000000000000008217301096052206306372172555),
            new TestCase1(Math.Pow(2,-100),1.000000000000000000000000000000000000000000000000000000000000103716921297685695119),
            new TestCase1(0.0 , 1),
        };

        #endregion

        static IEnumerable<TestCase<double>> GetNegative(TestCase1[] data)
        {
            return data.Select(d => TestCase.Create(-d.Param1, d.Result));
        }

        [TestMethod]
        public void SinhcTest()
        {
            TestCaseSet<double>[] testCases = {
                new TestCaseSet<double>(_SinhcData, "Sinhc spot data", 5),
                new TestCaseSet<double>(GetNegative(_SinhcData), "Sinhc spot data", 5),
            };

            NumericFunctionTest.RunSet(Math2.Sinhc, "Sinhc", testCases);

        }
    }
}
