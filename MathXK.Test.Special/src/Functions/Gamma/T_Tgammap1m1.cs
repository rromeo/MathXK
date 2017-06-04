﻿//  Copyright (c) 2013 Rocco Romeo
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

//  Test cases ported from Boost:
//      Copyright (c) John Maddock 2006, Boost Software License v1.0


using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace MathXK.Test.Functions
{

    using TestCase1 = TestCase<double>;

    [TestClass]
    public class T_Tgammap1m1
    {


        #region Test Case Data

        static readonly TestCase1[] _Tgamma1pm1Data = { 
            new TestCase1(-0.4952165186405181884765625, 0.7559827693907095754807809442951050489732), 
            new TestCase1(-0.4642883241176605224609375, 0.6574328869566978138139138799311066062094), 
            new TestCase1(-0.4024595916271209716796875, 0.4948624198600628575485791492257182331098), 
            new TestCase1(-0.3901382386684417724609375, 0.466994695902624837582482771934607569504), 
            new TestCase1(-0.3875354826450347900390625, 0.4612760821854033810460931227828015343156), 
            new TestCase1(-0.373013198375701904296875, 0.4303940312894863907637421435297095229576), 
            new TestCase1(-0.364522993564605712890625, 0.4131121158463471061555369932726425070959), 
            new TestCase1(-0.35811364650726318359375, 0.4004260060282522370607617707987066040664), 
            new TestCase1(-0.342386901378631591796875, 0.3705491389143732767395411311340356537958), 
            new TestCase1(-0.311618030071258544921875, 0.3168394199636867612664810682360364987851), 
            new TestCase1(-0.2880756855010986328125, 0.2795634828731819295069040561405863199118), 
            new TestCase1(-0.27896595001220703125, 0.2659494215075372040950689557832739194382), 
            new TestCase1(-0.221501767635345458984375, 0.1892800230490205183297725744929501015544), 
            new TestCase1(-0.202970564365386962890625, 0.1675837744862475402389911769714649852773), 
            new TestCase1(-0.191832959651947021484375, 0.1551782986717939911912971517092154248689), 
            new TestCase1(-0.1387059986591339111328125, 0.1019358076187795495751950704491293377337), 
            new TestCase1(-0.1012614667415618896484375, 0.06964849559574637322009665674759450635903), 
            new TestCase1(-0.0782387256622314453125, 0.05168944861403749320872745190411594799105), 
            new TestCase1(-0.0146243572235107421875, 0.008655823217365420051945171708362679166147), 
            new TestCase1(0.1431564604442703013402649929484277542953e-29, -0.8263215150028947028222714725368343067077e-30), 
            new TestCase1(0.1791466932348087634896446282571611213266e-29, -0.1034062776504410791508686126702507453693e-29), 
            new TestCase1(0.6013619202535540063110633226832922483532e-29, -0.3471155206456177563387123499681270341265e-29), 
            new TestCase1(0.115805324961653822428570241697281798758e-28, -0.6684464764687909153739175517973613571261e-29), 
            new TestCase1(0.1422457400834001098175711728787848259007e-28, -0.8210646944165041873250036406993552111264e-29), 
            new TestCase1(0.4970121018327539153628705477876439795096e-28, -0.2868831708235014101261168518349749234889e-28), 
            new TestCase1(0.9660079415057497591758174164417478444323e-28, -0.5575949162564024099347963251548294827165e-28), 
            new TestCase1(0.1232929313253182131376331095427391968754e-27, -0.7116661133260258147780879126862696392479e-28), 
            new TestCase1(0.3296523285617759312781860549364832953326e-27, -0.1902804880171240659658521971091414147824e-27), 
            new TestCase1(0.528364435768055252017009628713605422886e-27, -0.3049802291021812635024061181080113889506e-27), 
            new TestCase1(0.886586057273120049620324386849842094685e-27, -0.5117513605413324831805933701656815521236e-27), 
            new TestCase1(0.2499669674831043259218157022821422146034e-26, -0.1442848493391799075204112945868205575921e-26), 
            new TestCase1(0.4131050397232622964314362671638736040881e-26, -0.2384507001780369908735157814349257650306e-26), 
            new TestCase1(0.7679738097881433551381658732998641759182e-26, -0.4432865132438264916724504940164877771849e-26), 
            new TestCase1(0.199929739820949207249437007767740538737e-25, -0.1154025777043396680338534232061371984153e-25), 
            new TestCase1(0.5151477415246978459754129800826163591626e-25, -0.2973513461467014566158126478567428007026e-25), 
            new TestCase1(0.101200734533556026342258477595279955025e-24, -0.5841464927231005972586520769367701801272e-25), 
            new TestCase1(0.2064292695896540981798546456623054911033e-24, -0.1191542081013299677372796818685763621433e-24), 
            new TestCase1(0.4063294332896333395257434433879773416284e-24, -0.2345397140053387487331153502526122099635e-24), 
            new TestCase1(0.8138195767936862452966745688936976428456e-24, -0.4697494081288516881709792363579127898501e-24), 
            new TestCase1(0.9575550627132253801929510132578249716542e-24, -0.5527157822038433842858279196451987644551e-24), 
            new TestCase1(0.2855160956298500804375620841706273850616e-23, -0.1648043629790735548426012109835882111878e-23), 
            new TestCase1(0.65201444297915461398563707001320281266e-23, -0.3763529502296153146061378052316774680657e-23), 
            new TestCase1(0.1310988374636350038320977491775043421995e-22, -0.7567230263439006455136418784154172627278e-23), 
            new TestCase1(0.2590288837798696209228010176465529547374e-22, -0.1495155293796993236481538654450094828034e-22), 
            new TestCase1(0.2937779542193655202274099291941187976629e-22, -0.1695732371781431498672748630538267989951e-22), 
            new TestCase1(0.7863513178004503049754083414326234074965e-22, -0.4538942987503834952636621253829968594981e-22), 
            new TestCase1(0.1903818607087388763706780167350761726053e-21, -0.109891392314185724619218859435147914583e-21), 
            new TestCase1(0.3812242142377350870566942975497647799754e-21, -0.2200485882977986685255989911406099731464e-21), 
            new TestCase1(0.5493133580141330277178034419485741501887e-21, -0.3170722751854215599980143372792630827825e-21), 
            new TestCase1(0.9672153634284186955666772243312215295852e-21, -0.5582918591043124472447977064584587087382e-21), 
            new TestCase1(0.1702169477623814384559878647986894129041e-20, -0.9825188868017248805929057833465618943933e-21), 
            new TestCase1(0.4817114569977399785676754474621208412799e-20, -0.2780513989416366360400948471972900830387e-20), 
            new TestCase1(0.7538352992756463183303278501219690799218e-20, -0.4351255434976382024407146526549699565796e-20), 
            new TestCase1(0.2596305715949999708394617609422128090557e-19, -0.1498628330119729391523576917959646915902e-19), 
            new TestCase1(0.4444587480324321591032923385589104015025e-19, -0.2565485517668431887674899070788880246089e-19), 
            new TestCase1(0.9715574921498573937069095571295029856174e-19, -0.5607982038213457280620337061143959968827e-19), 
            new TestCase1(0.2036598542733453787268262970278076551267e-18, -0.1175556581981383412558675525031327289634e-18), 
            new TestCase1(0.4248971931658660264162106698360155121463e-18, -0.2452573158680304024128136998709195441154e-18), 
            new TestCase1(0.6521097487613458963613731825259556273977e-18, -0.3764079622200518159115220576707372289309e-18), 
            new TestCase1(0.1436126164096190058281493628911107407475e-17, -0.8289545186912702312307596583951073368993e-18), 
            new TestCase1(0.3118908901459261162419055180006211003274e-17, -0.1800283075323116855097674712567357424907e-17), 
            new TestCase1(0.3593346613595175715618300349429858897565e-17, -0.2074135954788010816821644850649799519675e-17), 
            new TestCase1(0.9445874854124767215374919304693435151421e-17, -0.5452306934500297136990208691880314661067e-17), 
            new TestCase1(0.2566182432094081539023303073498993853718e-16, -0.1481240698799817909729203530750507585117e-16), 
            new TestCase1(0.3363765695149349330660137891158001366421e-16, -0.1941618252298598450712101459398214835952e-16), 
            new TestCase1(0.1073581901339262605326457800103412409953e-15, -0.6196882910077942026172170471644452240087e-16), 
            new TestCase1(0.186668406231853462907965823802669547149e-15, -0.1077479282192287023344331697438668504475e-15), 
            new TestCase1(0.3727540802657755688795382376099496468669e-15, -0.2151594942853688563255514533064804448857e-15), 
            new TestCase1(0.6211646767866855090717281839829411183018e-15, -0.3585459819247720488263755672803206966864e-15), 
            new TestCase1(0.1561186859754253464932505224282976996619e-14, -0.9011415112885851355703262394946230181441e-15), 
            new TestCase1(0.3092010764722992466335682593125966377556e-14, -0.1784757049442269726369719145429686196432e-14), 
            new TestCase1(0.6192850577371690132255643845837767003104e-14, -0.3574610363653403859138387348420733335946e-14), 
            new TestCase1(0.1047879028014987723427253740737796761096e-13, -0.6048521898920322580919537460083385958307e-14), 
            new TestCase1(0.1978473638988408750405412206418986897916e-13, -0.1142005977018810929368580238123776497134e-13), 
            new TestCase1(0.4041816252346730475863978426787070930004e-13, -0.2332999655507978182935820166341591465938e-13), 
            new TestCase1(0.9410302262901834580155480125540634617209e-13, -0.5431773877604405885692410233977439758357e-13), 
            new TestCase1(0.1334530223958893535574077304772799834609e-12, -0.7703117505534481432167486795079037095016e-13), 
            new TestCase1(0.266297021326439287136622624529991298914e-12, -0.1537108122261681912683838906901126667239e-12), 
            new TestCase1(0.5920415525016708979677559909760020673275e-12, -0.3417356583762410657228951913532598329392e-12), 
            new TestCase1(0.155163989296047688526414276566356420517e-11, -0.8956308525005437046951559037146364454014e-12), 
            new TestCase1(0.326923297461201300961874949280172586441e-11, -0.1887052485148118271327652592003779278951e-11), 
            new TestCase1(0.3753785910581841633870681107509881258011e-11, -0.2166744030260566997416579153839836535144e-11), 
            new TestCase1(0.9579165585749116473834874341264367103577e-11, -0.5529244432689301568675461523208969657637e-11), 
            new TestCase1(0.1858167439361402273334533674642443656921e-10, -0.1072563353975220567027440177997243696045e-10), 
            new TestCase1(0.5449485307451595872407779097557067871094e-10, -0.3145528284818088265306644806429650505592e-10), 
            new TestCase1(0.6089519166696533147842274047434329986572e-10, -0.3514965854368603547185748339618264962e-10), 
            new TestCase1(0.1337744776064297980155970435589551925659e-9, -0.7721672402075083279577276691763777798079e-10), 
            new TestCase1(0.2554458866654840676346793770790100097656e-9, -0.1474473672534405166880467009401481006846e-9), 
            new TestCase1(0.9285605062636648199259070679545402526855e-9, -0.5359796691714968347006887450221011388011e-9), 
            new TestCase1(0.1698227447555211711005540564656257629395e-8, -0.980243482442200345890191122898173278899e-9), 
            new TestCase1(0.339355921141759608872234821319580078125e-8, -0.1958815525210918994679361737180445367869e-8), 
            new TestCase1(0.6313728651008432279922999441623687744141e-8, -0.3644383041872783822741274272616254061398e-8), 
            new TestCase1(0.8383264749056706932606175541877746582031e-8, -0.4838951666662356901644608730186546413827e-8), 
            new TestCase1(0.1962631124285962869180366396903991699219e-7, -0.1132861391263510827635198807950784606364e-7), 
            new TestCase1(0.5256384838503436185419559478759765625e-7, -0.3034067396263077516870264512967026960836e-7), 
            new TestCase1(0.116242290459922514855861663818359375e-6, -0.6709685761311096545819618213723265238395e-7), 
            new TestCase1(0.1776920584006802528165280818939208984375e-6, -0.1025666084085592538025029643343934584008e-6), 
            new TestCase1(0.246631174150024889968335628509521484375e-6, -0.142359317011220184093241753035251343271e-6), 
            new TestCase1(0.7932688959044753573834896087646484375e-6, -0.4578866108069128537168607745624614794722e-6), 
            new TestCase1(0.1372093493046122603118419647216796875e-5, -0.791991995861101927407880386227298092424e-6), 
            new TestCase1(0.214747751670074649155139923095703125e-5, -0.1239553101482841564081262476002391564028e-5), 
            new TestCase1(0.527022712049074470996856689453125e-5, -0.3042030180348038757030022358659056016131e-5), 
            new TestCase1(0.9233162927557714283466339111328125e-5, -0.5329441960781667799420300524501852822265e-5), 
            new TestCase1(0.269396477960981428623199462890625e-4, -0.1554926893050895772896453957379474045684e-4), 
            new TestCase1(0.3208058114978484809398651123046875e-4, -0.1851639610824637542072201121307897455015e-4), 
            new TestCase1(0.00010957030463032424449920654296875, -0.6323382317252073303439803607481865268353e-4), 
            new TestCase1(0.000126518702018074691295623779296875, -0.7301274674392621832822527575061291232133e-4), 
            new TestCase1(0.00028976381872780621051788330078125, -0.0001671731931845284705043154051351730706298), 
            new TestCase1(0.000687857042066752910614013671875, -0.0003965741858362509385654340724912705823236), 
            new TestCase1(0.00145484809763729572296142578125, -0.0008376704829300962209037156359622989708524), 
            new TestCase1(0.00366270542144775390625, -0.002100946766981816464155333982357752904999), 
            new TestCase1(0.046881496906280517578125, -0.02497588947336944943591732868959193811385), 
            new TestCase1(0.04722058773040771484375, -0.02514197077286474061941968460870171431946), 
            new TestCase1(0.1323592662811279296875, -0.06091072639354085529687250076608036167983), 
            new TestCase1(0.139763355255126953125, -0.06350237679012056526225278765134891026831), 
            new TestCase1(0.155740678310394287109375, -0.06883441875617310404657299803607175885806), 
            new TestCase1(0.17873513698577880859375, -0.07590347542832235360535014668688262042636), 
            new TestCase1(0.1813595294952392578125, -0.07666619009216679212771233211297326572625), 
            new TestCase1(0.225838959217071533203125, -0.08828121583826801836069461556215827733624), 
            new TestCase1(0.292207300662994384765625, -0.1013142147509511845410398110728244632327), 
            new TestCase1(0.297928631305694580078125, -0.1022125377703106514269725748156778633457), 
            new TestCase1(0.29810583591461181640625, -0.1022398124529602122123098193784073635017), 
            new TestCase1(0.30028045177459716796875, -0.1025718475377491669113360177547935208456), 
            new TestCase1(0.314723670482635498046875, -0.1046527148786850527567817015912636515928), 
            new TestCase1(0.335008561611175537109375, -0.1072166100764849221291653907100739024018), 
            new TestCase1(0.34912931919097900390625, -0.1087597918806720999316678995289547437333), 
            new TestCase1(0.378430664539337158203125, -0.1113472284261879798379642546535625292991), 
            new TestCase1(0.405791938304901123046875, -0.1130363541503776146472634379575969350382), 
            new TestCase1(0.4133758544921875, -0.1133834162443148976453781329005327033269), 
            new TestCase1(0.415735542774200439453125, -0.1134808289122708542397946052703635909288), 
            new TestCase1(0.43399322032928466796875, -0.1140666251072475325460987089801937766512), 
            new TestCase1(0.457166969776153564453125, -0.1143882508090213200510081514983929178479), 
            new TestCase1(0.457506835460662841796875, -0.1143895043011533114503086722896019260013), 
            new TestCase1(0.4594924449920654296875, -0.1143948425572877686750521180317509748309), 
            new TestCase1(0.464888513088226318359375, -0.1143922664386800836288584840918904877633), 
            new TestCase1(0.467694938182830810546875, -0.1143810844126184461958335384166685988825), 
            new TestCase1(0.468867778778076171875, -0.114374421494214380845931849419009841586), 
            new TestCase1(0.470592796802520751953125, -0.1143624939837521390307936225841266176506), 
            new TestCase1(0.481109678745269775390625, -0.1142351916017016455044428590981136786372), 
            new TestCase1(0.492881298065185546875, -0.1139822218283491394475161618135268104636), 
            new TestCase1(0.496461331844329833984375, -0.1138823102618022293565213085801050842487) 
        };


        /// <summary>
        /// Spot values computed from wolfram
        /// </summary>
        private static readonly TestCase1[] _Tgammap1m1Spot = {

            new TestCase1(0.93708860874176, -0.02498709936699429764758222555630719056199926841436 ),
            new TestCase1(1.00243771076202, 0.0010330744419057840257992588440244305136114079594494 ),
            new TestCase1(1.07191634178162, 0.0325674540584143101790070409280318699126211188593092),
           
           
        };


        #endregion

        [TestMethod]
        public void Tgammap1m1Test()
        {

            TestCaseSet<double>[] testCases = {
                new TestCaseSet<double>(_Tgamma1pm1Data,"Tgamma1pm1: Random Values",30),
 
                // added to test set: Rocco Romeo
                new TestCaseSet<double>(_Tgammap1m1Spot,"Tgamma1pm1: Spot Values",140) 
            };

            NumericFunctionTest.RunSet(Math2.Tgamma1pm1, "Tgamma1pm1", testCases);

        }
    }
}