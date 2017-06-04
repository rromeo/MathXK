//  Copyright (c) 2013 Rocco Romeo
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)


using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace MathXK.Test.Functions
{
    [TestClass]
    public class T_Ldexp
    {

        struct FrexpData
        {
            public readonly double X;
            public readonly double Fraction;
            public readonly int Exponent;

            public FrexpData(double x, double fraction, int exponent)
            {
                X = x;
                Fraction = fraction;
                Exponent = exponent;
            }

        }

        struct LdexpData
        {
            public readonly double X;
            public readonly int N;
            public readonly double Result;

            public LdexpData(double x, int n, double result)
            {
                X = x;
                N = n;
                Result = result;
            }

        }

        // computed using MS Visual C++ 2010 frexp
        static readonly FrexpData[] _FrexpData = {
            new FrexpData(double.MaxValue, 9.9999999999999989E-001, 1024),
            new FrexpData(double.MinValue, -9.9999999999999989E-001, 1024),

            new FrexpData(2.0499389086188863E-283, 9.5259181826961781E-001, -939),
            new FrexpData(6.7027964222990170E-283, 7.7868528234711165E-001, -937),
            new FrexpData(1.6305665152541377E-270, 8.6141937266170021E-001, -896),
            new FrexpData(2.4542247330227735E-253, 8.9966472093931615E-001, -839),
            new FrexpData(7.3249864038251399E-250, 6.5556119319776129E-001, -827),
            new FrexpData(2.2943213475729025E-224, 5.3077569628636156E-001, -742),
            new FrexpData(6.2943669415884398E-210, 5.1733164885888416E-001, -694),
            new FrexpData(6.5833774382011824E-202, 8.0627994649404655E-001, -668),
            new FrexpData(4.6041617255076366E-199, 5.5066530851212581E-001, -658),
            new FrexpData(1.1249819372582208E-166, 8.2922647652645898E-001, -551),
            new FrexpData(1.0167358245407663E-159, 8.9339977078763755E-001, -528),
            new FrexpData(3.8505598738252232E-154, 6.4534459013740741E-001, -509),
            new FrexpData(1.3163952896517198E-148, 8.4161640184075404E-001, -491),
            new FrexpData(6.3549969286379542E-136, 9.2381157210274212E-001, -449),
            new FrexpData(9.2885314589490588E-116, 9.1496693525228967E-001, -382),
            new FrexpData(3.3104889444870933E-093, 8.6317889013475690E-001, -307),
            new FrexpData(1.0560301943378588E-078, 9.7823954000055768E-001, -259),
            new FrexpData(1.8451425843247698E-055, 5.6553455654789442E-001, -181),
            new FrexpData(4.8142124295707140E-040, 6.5527664015780973E-001, -130),
            new FrexpData(7.0200971222798428E-034, 9.1126070586557528E-001, -110),
            new FrexpData(7.1932686572502240E-018, 5.1832963270988763E-001, -56),
            new FrexpData(5.6353182901909533E-007, 5.9090595114552691E-001, -20),
            new FrexpData(1.2612818445873301E+010, 7.3416265926052005E-001, 34),
            new FrexpData(2.1000075322937882E+043, 9.4167594562362689E-001, 144),
            new FrexpData(8.7104881397019384E+054, 7.1048233970556518E-001, 183),
            new FrexpData(1.0850128584832956E+063, 6.5938003857591188E-001, 210),
            new FrexpData(2.1525640682932284E+080, 9.0771025973731201E-001, 267),
            new FrexpData(3.0401416230826577E+081, 8.0124436914978747E-001, 271),
            new FrexpData(1.0002008415788584E+088, 6.2849016516864775E-001, 293),
            new FrexpData(1.5355734211978261E+109, 8.1730124485553302E-001, 363),
            new FrexpData(1.3876349671583673E+118, 6.8783922618620064E-001, 393),
            new FrexpData(3.5745193225673791E+122, 5.4072869640808963E-001, 408),
            new FrexpData(6.7122261111968920E+146, 8.3990215176711891E-001, 488),
            new FrexpData(4.5682138894646269E+148, 8.9315857453884495E-001, 494),
            new FrexpData(3.8494746483456993E+156, 5.6075573366170550E-001, 521),
            new FrexpData(3.1186371786513974E+157, 5.6786764392995637E-001, 524),
            new FrexpData(1.1766197550315303E+165, 6.3851119514485410E-001, 549),
            new FrexpData(3.4843878103394662E+184, 5.1251802989864625E-001, 614),
            new FrexpData(2.2238167941757738E+252, 6.0664223581861809E-001, 839),
            new FrexpData(3.6319636147663703E+254, 7.7404312452955315E-001, 846),
            new FrexpData(1.4050570125936581E+255, 7.4861344683385145E-001, 848),
            new FrexpData(4.9885830448476960E+263, 9.9015002930078189E-001, 876),
            new FrexpData(7.9108090258027163E+268, 5.9896959438872044E-001, 894),
            new FrexpData(9.1606671321958320E+271, 6.7734670484314985E-001, 904),
            new FrexpData(4.0376916254603294E+275, 7.2888182137575086E-001, 916),
            new FrexpData(3.9961336849962276E+283, 5.3746983806484838E-001, 943),
            new FrexpData(2.2169343761368105E+286, 5.8236727890668771E-001, 952),
            new FrexpData(2.0227575024349430E+291, 8.1078931727017711E-001, 968),
            new FrexpData(8.3966578664809865E+292, 5.2588487580692622E-001, 974),
            new FrexpData(4.0425944074561342E+299, 6.0364900557494205E-001, 996),
                                                 };
        // computed using MS Visual C++ 2010 ldexp
        static readonly LdexpData[] _LdexpData = {

            new LdexpData(4.8636085050476759E+003, -1010, 4.4326453832068602E-301),
            new LdexpData(2.7529941698317390E+002, -1004, 1.6057933129096461E-300),
            new LdexpData(7.3055576567098290E+007, -1014, 4.1613837731231149E-298),
            new LdexpData(4.2801973010786796E+001, -924, 3.0181984748043184E-277),
            new LdexpData(9.6452305231453991E+008, -912, 2.7858425430882980E-266),
            new LdexpData(4.6767922838576101E+006, -901, 2.7664446678328244E-265),
            new LdexpData(1.0469126613097789E+004, -842, 3.5698839896969531E-250),
            new LdexpData(1.6441768368911821E+007, -813, 3.0099691940511503E-238),
            new LdexpData(4.0927807121038342E+004, -794, 3.9282753241059916E-235),
            new LdexpData(6.0253840526684106E+005, -792, 2.3132798121687962E-233),
            new LdexpData(2.5885049684850984E+009, -800, 3.8819726539901126E-232),
            new LdexpData(3.8908226061103973E+006, -787, 4.7800764545571703E-231),
            new LdexpData(4.1300937599176978E+004, -719, 1.4975903348899463E-212),
            new LdexpData(4.4917026695291053E+004, -677, 7.1631481871299143E-200),
            new LdexpData(2.4933161132472175E+000, -654, 3.3354941205130282E-197),
            new LdexpData(9.2748247921815255E+002, -646, 3.1763512015497479E-192),
            new LdexpData(5.5562303520535259E+009, -658, 4.6456137202740359E-189),
            new LdexpData(7.7303186746921018E+000, -612, 4.5482052094033925E-184),
            new LdexpData(3.1903132392048434E+001, -559, 1.6906968751050579E-167),
            new LdexpData(7.5194197276231926E+003, -536, 3.3427711981250299E-158),
            new LdexpData(2.8637652531563140E+000, -514, 5.3397342580529021E-155),
            new LdexpData(5.2296142118189550E+000, -485, 5.2350611039060131E-146),
            new LdexpData(7.3759574198052096E+008, -502, 5.6332701343469135E-143),
            new LdexpData(5.7349615384521712E+002, -450, 1.9725701681655826E-133),
            new LdexpData(9.0344022582445308E+004, -448, 1.2429720591730910E-130),
            new LdexpData(4.3336268208611886E+000, -426, 2.5007677359046331E-128),
            new LdexpData(8.1862985427427557E+004, -431, 1.4762484258584842E-125),
            new LdexpData(1.6622910144097598E+009, -445, 1.8296119442680684E-125),
            new LdexpData(3.9432482020665184E+006, -419, 2.9126341067441768E-120),
            new LdexpData(2.6189615636220954E+007, -409, 1.9808924563640733E-116),
            new LdexpData(7.7371705359486729E+007, -394, 1.9176258599248126E-111),
            new LdexpData(3.4685122284425944E+008, -393, 1.7193129487901162E-110),
            new LdexpData(1.0411650506039918E+006, -346, 7.2634187322029691E-099),
            new LdexpData(5.5703066734458813E+006, -304, 1.7090722556448433E-085),
            new LdexpData(2.1460126257250145E+008, -291, 5.3939079973852593E-080),
            new LdexpData(2.3066625801411750E+001, -245, 4.0797648572064224E-073),
            new LdexpData(2.9994927453806123E+005, -250, 1.6578639954489569E-070),
            new LdexpData(4.8218237677690685E+008, -257, 2.0821041400707124E-069),
            new LdexpData(6.8314304839966030E+009, -241, 1.9332263160120962E-063),
            new LdexpData(4.7295838243931456E+001, -175, 9.8758318275475655E-052),
            new LdexpData(1.0264552203130234E+006, -189, 1.3081900068962855E-051),
            new LdexpData(9.7753270574730475E+000, -169, 1.3063574594411673E-050),
            new LdexpData(4.8594714630672112E+007, -176, 5.0735249360102080E-046),
            new LdexpData(5.2660256974706543E+004, -150, 3.6896368614803258E-041),
            new LdexpData(8.6922672012319162E+007, -121, 3.2696675170835115E-029),
            new LdexpData(6.8816088125489693E+009, -116, 8.2834352985294678E-026),
            new LdexpData(8.8767943712308968E+005, -101, 3.5012780215749938E-025),
            new LdexpData(1.3969066140227515E+000, -78, 4.6219762746668616E-024),
            new LdexpData(2.2116576285155758E+005, -86, 2.8585004873642039E-021),
            new LdexpData(6.9693228495164413E+005, -58, 2.4179695917435494E-012),
            new LdexpData(7.6508202113447813E+004, -23, 9.1204884187516943E-003),
            new LdexpData(8.7454025767489355E+006, -18, 3.3361063296313993E+001),
            new LdexpData(5.1913331426639729E+002, 1, 1.0382666285327946E+003),
            new LdexpData(3.7973696197370859E+000, 37, 5.2190650648303168E+011),
            new LdexpData(4.7288703893911419E+005, 60, 5.4520163644276017E+023),
            new LdexpData(9.6274369720497343E+003, 96, 7.6276414101739332E+032),
            new LdexpData(6.1318248538914258E+003, 102, 3.1092045826119364E+034),
            new LdexpData(1.6212772958079714E+005, 138, 5.6493268546838017E+046),
            new LdexpData(4.1985775270492663E+009, 126, 3.5717547465134622E+047),
            new LdexpData(3.3120295513942975E+008, 153, 3.7816692283211732E+054),
            new LdexpData(5.3427227006027110E+001, 182, 3.2750794859026456E+056),
            new LdexpData(9.5618045048610316E+000, 210, 1.5733992888963266E+064),
            new LdexpData(8.4981670549596508E+005, 212, 5.5935090471252150E+069),
            new LdexpData(4.8273788927250447E+004, 247, 1.0917427490795997E+079),
            new LdexpData(3.7236138636821945E+000, 280, 7.2337488196358557E+084),
            new LdexpData(2.3821273596874100E+008, 256, 2.7583150380757699E+085),
            new LdexpData(5.2058509156135151E+003, 273, 7.9009723682243787E+085),
            new LdexpData(4.4560000740252435E+004, 270, 8.4536452417620963E+085),
            new LdexpData(7.9170872494192272E+006, 272, 6.0079215298764432E+088),
            new LdexpData(6.2637871020770569E+002, 304, 2.0415295480049438E+094),
            new LdexpData(3.7453894334590942E+001, 318, 2.0000258185859468E+097),
            new LdexpData(7.5917618896773842E+005, 306, 9.8974029395442645E+097),
            new LdexpData(1.1638551164250147E+005, 373, 2.2391712162638598E+117),
            new LdexpData(1.7563684013468998E+009, 379, 2.1626387072819378E+123),
            new LdexpData(2.8332544410756111E+002, 495, 2.8982338991240432E+151),
            new LdexpData(8.2307754588456703E+000, 503, 2.1554034466149980E+152),
            new LdexpData(8.1499477420700472E+006, 524, 4.4758193751959365E+164),
            new LdexpData(3.1049695252929673E+007, 545, 3.5760614986180533E+171),
            new LdexpData(2.1811057683923774E+000, 575, 2.6972688267589542E+173),
            new LdexpData(4.0672931373432712E+003, 568, 3.9295517365934147E+174),
            new LdexpData(4.1161784502287295E+000, 618, 4.4774579407464826E+186),
            new LdexpData(7.4898143380367267E+008, 594, 4.8561095630799196E+187),
            new LdexpData(4.5092963168906499E+002, 621, 3.9240640016700872E+189),
            new LdexpData(7.7258451784615609E+003, 651, 7.2189357768035579E+199),
            new LdexpData(3.8502315403196174E+006, 707, 2.5923508677771665E+219),
            new LdexpData(2.8067516408892121E+009, 717, 1.9351332550044749E+225),
            new LdexpData(8.9073347714298097E+002, 795, 1.8560699004495797E+242),
            new LdexpData(2.7364946169482982E+008, 799, 9.1234928006546404E+248),
            new LdexpData(5.0080924923569393E+002, 845, 1.1749480833306311E+257),
            new LdexpData(3.1585249366045581E+001, 920, 2.7994990268824671E+278),
            new LdexpData(8.3470604141927600E+007, 902, 2.8222120754413030E+279),
            new LdexpData(6.7883159531225991E+006, 923, 4.8133566867746270E+284),
            new LdexpData(8.0827721525369078E+007, 927, 9.1699362626615514E+286),
            new LdexpData(1.9563946369832358E+001, 995, 6.5509177893143440E+300),
            new LdexpData(6.1337641608592889E+002, 991, 1.2836681821897227E+301),
            new LdexpData(1.7945476816010418E+002, 1007, 2.4612778071573471E+305),
            new LdexpData(8.3544960424092121E+009, 984, 1.3659537381202816E+306),
            new LdexpData(5.3895575019482321E+003, 1011, 1.1827112452634314E+308),
            new LdexpData(7.3332627570877135E+007, 1013, double.PositiveInfinity),
            new LdexpData(3.3561006855946500E+006, 1007, double.PositiveInfinity),
                                                 };

        [TestMethod]
        public void FrexpTest()
        {

            const int MaxNormalExponent = 1023;
            const int MinNormalExponent = -1022;
            const int MinDenormExponent = -1074;


            int exponent;
            double mantissa;

            // Zero
            mantissa = Math2.Frexp(0.0, out exponent);
            Assert.AreEqual(0, mantissa);
            Assert.AreEqual(0, exponent);

            // NaN
            Assert2.AreEqual<Exceptions.DomainException>(double.NaN, () => Math2.Frexp(double.NaN, out exponent));

            mantissa = Math2.Frexp(double.PositiveInfinity, out exponent);
            Assert.AreEqual(double.PositiveInfinity, mantissa);

            mantissa = Math2.Frexp(double.NegativeInfinity, out exponent);
            Assert.AreEqual(double.NegativeInfinity, mantissa);

            mantissa = Math2.Frexp(double.Epsilon, out exponent);
            Assert.AreEqual(0.5, mantissa);
            Assert.AreEqual(MinDenormExponent + 1, exponent);

            mantissa = Math2.Frexp(double.MaxValue, out exponent);
            // Assert.AreEqual(0.5, mantissa);
            Assert.AreEqual(1024, exponent);

            // check normal numbers -- powers of 2
            for (int i = MinNormalExponent; i < MaxNormalExponent; i++) {
                double x = Math.Pow(2, i);
                mantissa = Math2.Frexp(x, out exponent);
                Assert2.AreNear(0.5, mantissa, 1);
                Assert.AreEqual(i + 1, exponent);

                x = -x;
                mantissa = Math2.Frexp(x, out exponent);
                Assert2.AreNear(-0.5, mantissa, 1);
                Assert.AreEqual(i + 1, exponent);

            }

            // denorm
            for (int i = 0; i < 52; i++) {
                double x = double.Epsilon * (double)(1UL << i);
                Assert.IsTrue(x > 0);

                mantissa = Math2.Frexp(x, out exponent);
                Assert2.AreNear(0.5, mantissa, 0);
                Assert.AreEqual(MinDenormExponent + 1 + i, exponent);

                x = -x;
                mantissa = Math2.Frexp(x, out exponent);
                Assert2.AreNear(-0.5, mantissa, 0);
                Assert.AreEqual(MinDenormExponent + 1 + i, exponent);

            }

            foreach (var item in _FrexpData) {
                double x = item.X;
                double fraction = item.Fraction;
                int exp = item.Exponent;

                int computed_exponent;
                double computed_fraction = Math2.Frexp(x, out computed_exponent);
                Assert2.AreNear(fraction, computed_fraction, 1);
                Assert.AreEqual(exp, computed_exponent);

            }

        }


        [TestMethod]
        public void LdexpTest()
        {
            Assert.AreEqual(0.75 * Math.Pow(2, -1073), Math2.Ldexp(0.75, -1073));

            const int MaxNormalExponent = 1023;
            const int MinNormalExponent = -1022;
            const int MinDenormExponent = -1074;

            // NaN
            int[] nValues = { -1, 0, 1, 2 };
            foreach (int n in nValues)
                Assert2.AreEqual<Exceptions.DomainException>(double.NaN, () => Math2.Ldexp(double.NaN, n));

            // edge values
            double[] xValues = { double.PositiveInfinity, double.NegativeInfinity, 0 };
            nValues = new int[] { -1, 0, 1, 2 };
            foreach (double x in xValues)
                foreach (int n in nValues)
                    Assert.AreEqual(Math2.Ldexp(x, n), x);

            // overflow and underflow
            Assert.AreEqual(0.0, Math2.Ldexp(1, MinDenormExponent - 1));
            Assert.AreEqual(-0.0, Math2.Ldexp(-1, MinDenormExponent - 1));
            Assert.AreEqual(double.PositiveInfinity, Math2.Ldexp(1, MaxNormalExponent + 1));
            Assert.AreEqual(double.NegativeInfinity, Math2.Ldexp(-1, MaxNormalExponent + 1));


            // normal powers of 2
            for (int i = MinNormalExponent; i < MaxNormalExponent; i++) {
                double expected = Math.Pow(2, i);
                double calculated = Math2.Ldexp(1, i);
                Assert2.AreNear(expected, calculated, 2,
                    () => string.Format("Ldexp(1, {0}) = {1:E16}; got {2:E16}", i, expected, calculated));

                expected = -expected;
                calculated = Math2.Ldexp(-1, i);
                Assert2.AreNear(expected, calculated, 2,
                    () => string.Format("Ldexp(-1, {0}) = {1:E16}; got {2:E16}", i, expected, calculated));

            }

            // denorm
            for (int i = 0; i < 52; i++) {
                double expected = double.Epsilon * (1UL << i);
                Assert.IsTrue(expected > 0);

                int n = MinDenormExponent + i;

                double calculated = Math2.Ldexp(1, n);
                Assert2.AreNear(expected, calculated, 0,
                    () => string.Format("Ldexp(1, {0}) = {1:E16}; got {2:E16}", n, expected, calculated));


                expected = -expected;
                calculated = Math2.Ldexp(-1, n);
                Assert2.AreNear(expected, calculated, 0,
                    () => string.Format("Ldexp(-1, {0}) = {1:E16}; got {2:E16}", n, expected, calculated));

            }


            // denorm
            for (int i = 0; i < 52; i++) {
                double expected = double.Epsilon * (1UL << i);
                Assert.IsTrue(expected > 0);

                double calculated = Math2.Ldexp(double.Epsilon, i);
                Assert2.AreNear(expected, calculated, 0,
                    () => string.Format("Ldexp(double.Epsilon, {0}) = {1:E16}; got {2:E16}", i, expected, calculated));

                expected = -expected;
                calculated = Math2.Ldexp(-double.Epsilon, i);
                Assert2.AreNear(expected, calculated, 0,
                    () => string.Format("Ldexp(-double.Epsilon, {0}) = {1:E16}; got {2:E16}", i, expected, calculated));




            }


            foreach (var item in _LdexpData) {
                double x = item.X;
                int n = item.N;
                double expected = item.Result;

                double computed = Math2.Ldexp(x, n);
                Assert2.AreNear(expected, computed, 1);

            }


        }

        [TestMethod]
        public void LogbTest()
        {

            const int MaxNormalExponent = 1023;
            const int MinNormalExponent = -1022;
            const int MinDenormExponent = -1074;


            // NaN, 0, +/- Infinity
            Assert2.AreEqual<Exceptions.DomainException>(double.NaN, () => Math2.Logb(double.NaN));
            Assert2.AreEqual<Exceptions.PoleException>(double.NegativeInfinity, () => Math2.Logb(0));
            Assert.AreEqual(double.PositiveInfinity, Math2.Logb(double.PositiveInfinity));
            Assert.AreEqual(double.PositiveInfinity, Math2.Logb(double.NegativeInfinity));


            Assert.AreEqual(MinDenormExponent, Math2.Logb(double.Epsilon));
            Assert.AreEqual(MaxNormalExponent, Math2.Logb(double.MaxValue));
            Assert.AreEqual(MaxNormalExponent, Math2.Logb(double.MaxValue));


            // check normal numbers -- powers of 2
            for (int i = MinNormalExponent; i < MaxNormalExponent; i++) {
                double x = Math.Pow(2, i);
                Assert.AreEqual(i, Math2.Logb(x));
                Assert.AreEqual(i, Math2.Logb(-x));

            }

            // denorm
            for (int i = 0; i < 52; i++) {
                double x = double.Epsilon * (1UL << i);
                Assert.IsTrue(x > 0);

                Assert.AreEqual(MinDenormExponent+i, Math2.Logb(x));
                Assert.AreEqual(MinDenormExponent+i, Math2.Logb(-x));
            }

            foreach (var item in _FrexpData) {
                Assert.AreEqual(item.Exponent - 1, Math2.Logb(item.X));

            }

        }

    }
}
