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
    public class T_Sinc
    {

        #region Test Data

        // added to boost tests
        // function values calculated on http://functions.wolfram.com/
        private static readonly TestCase1[] _SincData = { 
            // range1 >= 2^-13
            new TestCase1(2, 0.4546487134128408476980099329558724213511274857239451),
            new TestCase1(1.5, 0.6649966577360362872944822474276582151377676172814105),
            new TestCase1(1, 0.8414709848078965066525023216302989996225630607983710),
            new TestCase1(0.03125, 0.9998372475304345196665958802388628680344998015657932),
            new TestCase1(0.015625, 0.9999593103925358191866113133277558808093941468711926),
            new TestCase1(0.00390625, 0.9999974568704298379922122219308794575027612890345612),
            new TestCase1(9.765625e-4, 0.9999998410542881780806760121546612396056760793172566),
            new TestCase1(2.44140625e-4, 0.9999999900658925670408431146553179354490941888458829),

            // 2^-13 <= Range2 <= 2^-26
            new TestCase1(6.103515625e-5, 0.999999999379118283705329221304862776873118511660449), // 2^-16
            new TestCase1(1.52587890625e-5, 0.999999999961194892724806812754158122527169078473952), // 2^-18
            new TestCase1(3.814697265625e-6, 0.999999999997574680795273956017533341680966727453830), // 2^-20
            new TestCase1(9.5367431640625e-7, 0.999999999999848417549704518853519264764664688137825), // 2^-22
            new TestCase1(2.384185791015625e-7, 0.999999999999990526096856532024448170574642791515396), // 2^-24
            new TestCase1(1.490116119384765625e-8, 0.9999999999999999629925658458281157300773325329163027), // 2^-26

            // Range3 < 2^-26
            new TestCase1(Math.Pow(2,-28), 0.9999999999999999976870353653642572090557089784668209),
            new TestCase1(Math.Pow(2,-30), 0.9999999999999999998554397103352660754719422630883933),
            new TestCase1(Math.Pow(2,-40), 0.9999999999999999999999998621365645744953875214318904),
            new TestCase1(Math.Pow(2,-52), 0.9999999999999999999999999999999917826989039477936936),
            new TestCase1(Math.Pow(2,-100),1.0),
            new TestCase1(0.0 ,1.0),
        };

        #endregion

        static IEnumerable<TestCase<double>> GetNegative(TestCase1[] data)
        {
            return data.Select(d => TestCase.Create(-d.Param1, d.Result));
        }

        [TestMethod]
        public void SincTest()
        {
            TestCaseSet<double>[] testCases = {
                new TestCaseSet<double>(_SincData, "Sinc spot data", 5),
                new TestCaseSet<double>(GetNegative(_SincData), "Sinc spot data", 5),
            };

            NumericFunctionTest.RunSet(Math2.Sinc, "Sinc", testCases);
        }
    }
}
