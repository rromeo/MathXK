﻿//  Copyright (c) 2013 Rocco Romeo
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

//  Test Data:
//      Copyright (c) 2008 John Maddock, Boost Software License v1.0

using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace MathXK.Test.Functions
{

    using TestCase1 = TestCase<double>;


    /// <summary>
    ///This is a test class for Math2Test and is intended
    ///to contain all Math2Test Unit Tests
    ///</summary>
    [TestClass()]
    public class T_Asinh
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



        #region Test Case Data

        private static readonly TestCase1[] _AsinhData = {
            new TestCase1(-497.181640625, -6.902103625349694816896128061929344488774), 
            new TestCase1(-495.216552734375, -6.898143347143858285101511567975488446084), 
            new TestCase1(-488.0980224609375, -6.883664481302669222996695409880306523751), 
            new TestCase1(-486.4609375, -6.880304842490272781295755933349565066089), 
            new TestCase1(-482.2261962890625, -6.87156154650904602462361416438487324277), 
            new TestCase1(-468.167236328125, -6.841973895837549314760586943747302390125), 
            new TestCase1(-465.553955078125, -6.836376331805492593765664636569933371943), 
            new TestCase1(-464.288330078125, -6.833654100575195198545828306587445927099), 
            new TestCase1(-463.558837890625, -6.832081663500904935919737735423649299577), 
            new TestCase1(-453.82861328125, -6.81086801736630853656142245668310545817), 
            new TestCase1(-448.7835693359375, -6.799689165151486910757237507031024212372), 
            new TestCase1(-446.0499267578125, -6.793579326246197169098958913912928750651), 
            new TestCase1(-432.4046630859375, -6.762510387544996144846476710732722287747), 
            new TestCase1(-424.145751953125, -6.74322572098922178423984261142063288429), 
            new TestCase1(-402.8682861328125, -6.69175839599430693315848859453612846314), 
            new TestCase1(-402.4595947265625, -6.690743430063694316478378930442947372479), 
            new TestCase1(-390.1383056640625, -6.65965012921145083591836824643494979965), 
            new TestCase1(-387.5355224609375, -6.652956360641761371501126846246743833405), 
            new TestCase1(-381.0023193359375, -6.635954365364266311169729037072694814964), 
            new TestCase1(-374.8172607421875, -6.619587562578274050554380148151080691039), 
            new TestCase1(-374.1033935546875, -6.617681179427804403332414367961774779791), 
            new TestCase1(-373.01318359375, -6.614762741096184771534770781005436075363), 
            new TestCase1(-370.0938720703125, -6.606905687537060691122483933557635629252), 
            new TestCase1(-364.5230712890625, -6.591738907156093293807291816593779744849), 
            new TestCase1(-361.3756103515625, -6.583066984213974024130269349536598335519), 
            new TestCase1(-358.1136474609375, -6.573999516974134646764605942478671320562), 
            new TestCase1(-350.8861083984375, -6.553610904389895936746972395658913458292), 
            new TestCase1(-350.7060546875, -6.553097634736137347416591397792124912806), 
            new TestCase1(-345.5616455078125, -6.538320325468201861196672683728074373679), 
            new TestCase1(-342.386962890625, -6.529090881007076202584033383216307968614), 
            new TestCase1(-341.9425048828125, -6.527791927233787108981347318463639017897), 
            new TestCase1(-337.3883056640625, -6.51438388615078099259979168949837683455), 
            new TestCase1(-328.8133544921875, -6.488639771044976384845381893407049569445), 
            new TestCase1(-326.1348876953125, -6.480460592697476652284072277108613530788), 
            new TestCase1(-313.12744140625, -6.439759999015992467577771895808599408942), 
            new TestCase1(-311.6180419921875, -6.434927968512049412470577925892037200619), 
            new TestCase1(-303.40478515625, -6.408217734896572191669656774970329190625), 
            new TestCase1(-291.9320068359375, -6.369671035834964490622480119650197877027), 
            new TestCase1(-289.791015625, -6.362310184909174825124230502026784750431), 
            new TestCase1(-288.07568359375, -6.356373428913315483326002974850146028147), 
            new TestCase1(-282.76220703125, -6.337756593913613854161781726501299103316), 
            new TestCase1(-278.9659423828125, -6.3242400970614694374776334482330766038), 
            new TestCase1(-276.1881103515625, -6.314232650754295158800596164826483845909), 
            new TestCase1(-269.843994140625, -6.290994606392703105943420882029653368811), 
            new TestCase1(-256.47509765625, -6.240182555852785458496224304526376585373), 
            new TestCase1(-248.91619873046875, -6.210267503979360726887449195178911415313), 
            new TestCase1(-245.71783447265625, -6.197335184435549180000272192577714732671), 
            new TestCase1(-244.9049072265625, -6.194021350132334791819182281057947070343), 
            new TestCase1(-242.49176025390625, -6.184119163536405964259882205805060501935), 
            new TestCase1(-223.97491455078125, -6.104686221071835053635371385862353306867), 
            new TestCase1(-223.0770263671875, -6.100669325836893543736746195548392269128), 
            new TestCase1(-221.50177001953125, -6.093582856519022569016903564242082589239), 
            new TestCase1(-214.1610107421875, -6.059880750068793807232821057354816859259), 
            new TestCase1(-202.9705810546875, -6.006214296526251845395554109089589227751), 
            new TestCase1(-200.1683349609375, -5.992312107336994557313003250781997964749), 
            new TestCase1(-198.0869140625, -5.981859446096082920408103292664994574137), 
            new TestCase1(-191.8330078125, -5.949779216585290285138095990044779123562), 
            new TestCase1(-183.4495849609375, -5.905094497458789993605346377453537455651), 
            new TestCase1(-182.9005126953125, -5.90209701227578839648368537690039781634), 
            new TestCase1(-167.5517578125, -5.814448391006795100253762785669751772879), 
            new TestCase1(-162.87738037109375, -5.786154254111213380715127555746974140457), 
            new TestCase1(-159.6142578125, -5.765917008989404365321507899157382728927), 
            new TestCase1(-150.01629638671875, -5.703902219845273396365649927451080691209), 
            new TestCase1(-148.34051513671875, -5.69266895044603959347064515800579166058), 
            new TestCase1(-147.23760986328125, -5.685206387751923160853276372644826819905), 
            new TestCase1(-143.65484619140625, -5.660572815631807047475426769799365771506), 
            new TestCase1(-138.70599365234375, -5.625516713960632970798501811037458042606), 
            new TestCase1(-119.554168701171875, -5.476934234171878793825519703714137456694), 
            new TestCase1(-118.441558837890625, -5.467584665632571047898569146378815818543), 
            new TestCase1(-112.7041015625, -5.41793267560343365856590869885928647495), 
            new TestCase1(-111.430206298828125, -5.406565756574078635952822819651240567772), 
            new TestCase1(-107.772979736328125, -5.373195678988387579981120804358907493552), 
            new TestCase1(-107.6795654296875, -5.372328571218373638952712194621628880899), 
            new TestCase1(-105.091796875, -5.348004040102252742967063201190196082691), 
            new TestCase1(-101.261474609375, -5.310877589708960332430067923458969169195), 
            new TestCase1(-95.79150390625, -5.255348419702703704863588329233038165262), 
            new TestCase1(-91.26885986328125, -5.206986845736275651134721240774128212633), 
            new TestCase1(-87.33349609375, -5.162914035396618398952584430311544798966), 
            new TestCase1(-78.238739013671875, -5.052952927749896114494197793427686940384), 
            new TestCase1(-77.912353515625, -5.048772883924985058524898899156261070079), 
            new TestCase1(-76.83489990234375, -5.03484848764480888947313271747141240657), 
            new TestCase1(-61.255645751953125, -4.808269821238498732397629948430226353629), 
            new TestCase1(-54.4138031005859375, -4.689849459883310559788337051110370703783), 
            new TestCase1(-43.967193603515625, -4.476720236388958724591705984989483327501), 
            new TestCase1(-42.0108489990234375, -4.431216695067421800858619793836241888979), 
            new TestCase1(-30.609375, -4.114720236218123790586146467787932912866), 
            new TestCase1(-26.7111663818359375, -3.978579083165602241813343216612781617222), 
            new TestCase1(-25.2413177490234375, -3.922021583095348294972916890287113874009), 
            new TestCase1(-14.624359130859375, -3.377002632462029261128559390722023174755), 
            new TestCase1(-12.431087493896484375, -3.214961448471211148851788664569558466737), 
            new TestCase1(-10.235607147216796875, -3.021397455139020950019259089284989476735), 
            new TestCase1(-9.41094970703125, -2.937831931335705068892682801940661221814), 
            new TestCase1(-1.635939121246337890625, -1.267878515574958901161281414987746394126), 
            new TestCase1(0.165048140085555239409131900174543261528e-11, 0.1650481400855552394091318252402434490969e-11), 
            new TestCase1(0.2065420751096169738048047292977571487427e-11, 0.2065420751096169738048045824476195851805e-11), 
            new TestCase1(0.6933230031758164102484442992135882377625e-11, 0.6933230031758164102484387445779249960378e-11), 
            new TestCase1(0.1335144494962747785393730737268924713135e-10, 0.1335144494962747785393691069885154254178e-10), 
            new TestCase1(0.1639981206391638579589198343455791473389e-10, 0.1639981206391638579589124830249793630233e-10), 
            new TestCase1(0.5730159402528300915946601890027523040771e-10, 0.5730159402528300915943466086387879674642e-10), 
            new TestCase1(0.1113731329382972035091370344161987304688e-9, 0.1113731329382972035089067894949093881495e-9), 
            new TestCase1(0.1421470718909745301061775535345077514648e-9, 0.1421470718909745301056988545527660628072e-9), 
            new TestCase1(0.3800632031314421510614920407533645629883e-9, 0.3800632031314421510523421433949182208556e-9), 
            new TestCase1(0.6091627202664540163823403418064117431641e-9, 0.6091627202664540163446657373156062782485e-9), 
            new TestCase1(0.1022164131114777774200774729251861572266e-8, 0.1022164131114777774022778557990306872369e-8), 
            new TestCase1(0.2881922256392499548383057117462158203125e-8, 0.2881922256392499544393767813667092771012e-8), 
            new TestCase1(0.4762776839584148547146469354629516601563e-8, 0.476277683958414852913996340565203355397e-8), 
            new TestCase1(0.8854133426439148024655878543853759765625e-8, 0.8854133426439147908968245283871083240172e-8), 
            new TestCase1(0.2305032609228874207474291324615478515625e-7, 0.2305032609228874003356918102983865074646e-7), 
            new TestCase1(0.593924909253473742865025997161865234375e-7, 0.5939249092534733936898428443727569533982e-7), 
            new TestCase1(0.116676488914890796877443790435791015625e-6, 0.1166764889148905321500984793498641255822e-6), 
            new TestCase1(0.23799674409019644372165203094482421875e-6, 0.2379967440901941969351979802846104416105e-6), 
            new TestCase1(0.468465941594331525266170501708984375e-6, 0.4684659415943143903171559983024977629855e-6), 
            new TestCase1(0.938269977268646471202373504638671875e-6, 0.9382699772685088034539126253793231155849e-6), 
            new TestCase1(0.11039855962735600769519805908203125e-5, 0.1103985596273335823585612314642277768255e-5), 
            new TestCase1(0.3291776010883040726184844970703125e-5, 0.3291776010877095894302224300936487731914e-5), 
            new TestCase1(0.7517213816754519939422607421875e-5, 0.7517213816683722188794989550887653607132e-5), 
            new TestCase1(0.1511466689407825469970703125e-4, 0.1511466689350275580917503171217159502838e-4), 
            new TestCase1(0.2986399340443313121795654296875e-4, 0.2986399339999405714013629488236926241004e-4), 
            new TestCase1(0.338702811859548091888427734375e-4, 0.3387028117947883430533095948345065139191e-4), 
            new TestCase1(0.906601198948919773101806640625e-4, 0.906601197706988351306957830965994124706e-4), 
            new TestCase1(0.0002194953267462551593780517578125, 0.0002194953249837736286256750985648679367941), 
            new TestCase1(0.000439521507360041141510009765625, 0.0004395214932089767739989257698158711437628), 
            new TestCase1(0.000633315183222293853759765625, 0.0006333151408864353233898099984279240916971), 
            new TestCase1(0.0011151232756674289703369140625, 0.001115123044558274291478926657905120008324), 
            new TestCase1(0.001962467096745967864990234375, 0.001962465837080717523701980763915077400694), 
            new TestCase1(0.00555375404655933380126953125, 0.005553725496786973429128982795141683132844), 
            new TestCase1(0.0086911283433437347412109375, 0.008691018931968293864799414130929206833958), 
            new TestCase1(0.0299333631992340087890625, 0.02992889492062483965221469264905787460865), 
            new TestCase1(0.05124260485172271728515625, 0.05122020579778826952521305025815121247091), 
            new TestCase1(0.1120129525661468505859375, 0.1117800293787827963417928974597546321371), 
            new TestCase1(0.23480379581451416015625, 0.2326980652154337550758180136962670174127), 
            new TestCase1(0.48987305164337158203125, 0.4721357117742938088066477027937692054202), 
            new TestCase1(0.7518312931060791015625, 0.6946115711893359819020679952345318169567), 
            new TestCase1(1.6557407379150390625, 1.278160734826225530236928993772347284054), 
            new TestCase1(3.59585666656494140625, 1.991726234324511503262116200593118895023), 
            new TestCase1(3.66270542144775390625, 2.009484184971721909158848583710336926639), 
            new TestCase1(4.14284515380859375, 2.128787712416204967344704932367445907457), 
            new TestCase1(5.957065582275390625, 2.484696793415547705602607297785951657088), 
            new TestCase1(10.890350341796875, 3.083125584533294091723482211217314707631), 
            new TestCase1(27.3714599609375, 4.002981567623351049359177787856924686562), 
            new TestCase1(29.58606719970703125, 4.080736210902825878904303085045024018186), 
            new TestCase1(30.79753875732421875, 4.12084543001111324730244006549246292804), 
            new TestCase1(38.7815704345703125, 4.351258506393415652318140630603706518075), 
            new TestCase1(46.8814849853515625, 4.540883728536112674069796475595291204506), 
            new TestCase1(47.21551513671875, 4.547981853382592216976253569088895438026), 
            new TestCase1(47.2205810546875, 4.548089117076700023837192332723136228729), 
            new TestCase1(49.7236175537109375, 4.599728302509060806991933759403338520635), 
            new TestCase1(61.557464599609375, 4.813184271185753146544327202950243752583), 
            new TestCase1(67.821624755859375, 4.910082619934557664814376219665476353171), 
            new TestCase1(68.823638916015625, 4.924747230639766605446150280099353765226), 
            new TestCase1(73.754669189453125, 4.993937439635390959095430118863527649672), 
            new TestCase1(80.956695556640625, 5.087099712053553781191118720872632390369), 
            new TestCase1(85.264068603515625, 5.138934697019629394937832600298516485788), 
            new TestCase1(85.2677001953125, 5.138977285472120998185283011836704311053), 
            new TestCase1(92.8238525390625, 5.223879832616765332217842454967156441878), 
            new TestCase1(94.503570556640625, 5.241812789460326774676763952219373084611), 
            new TestCase1(116.044677734375, 5.447141014648796298911597081177174321311), 
            new TestCase1(123.775543212890625, 5.511633288238573314525515498135556594256), 
            new TestCase1(132.3592529296875, 5.578681289305597653175233933020307342597), 
            new TestCase1(139.7633056640625, 5.633110296634630769495301352527335620124), 
            new TestCase1(143.9609375, 5.662701238627724704458477126104574527542), 
            new TestCase1(146.31298828125, 5.67890694100532316423025711195756179793), 
            new TestCase1(155.0980224609375, 5.737214893086865588590011960998256979258), 
            new TestCase1(155.47784423828125, 5.739660763047893353413979379888554852015), 
            new TestCase1(155.74066162109375, 5.741349685869528141606229427222705529931), 
            new TestCase1(163.60546875, 5.790614371552514063824117620171866397141), 
            new TestCase1(178.735107421875, 5.879059869096351492478036425245903640013), 
            new TestCase1(179.70269775390625, 5.884458728291027196042299574564237490922), 
            new TestCase1(179.81976318359375, 5.885109945587401516601219374963261030511), 
            new TestCase1(181.3594970703125, 5.893636014368935823237345043968331464439), 
            new TestCase1(194.82861328125, 5.965274032538233309914029311001854910366), 
            new TestCase1(195.23284912109375, 5.967346683696556361432840609074921098744), 
            new TestCase1(199.07666015625, 5.986843466070591664424697234367005123532), 
            new TestCase1(205.77423095703125, 6.019932686217941432444903969630541836634), 
            new TestCase1(206.04608154296875, 6.021252909681260874009121058600710829746), 
            new TestCase1(209.36480712890625, 6.037231102920488838374618484909263149999), 
            new TestCase1(210.703857421875, 6.043606439928323259236771869386592943524), 
            new TestCase1(215.2139892578125, 6.064785410115010353953388909263073452906), 
            new TestCase1(225.83892822265625, 6.112974120371601210219399663570367666925), 
            new TestCase1(226.95465087890625, 6.117902255760310722524206309196088056026), 
            new TestCase1(232.79864501953125, 6.143325688959409019088830019463644651492), 
            new TestCase1(240.647216796875, 6.176483527820343060796486861138751990657), 
            new TestCase1(243.1324462890625, 6.186757751007361577528719065211514303001), 
            new TestCase1(251.26702880859375, 6.219667373726848075772693817688242362523), 
            new TestCase1(253.72906494140625, 6.22941808808355428724031065923268579187), 
            new TestCase1(254.6866455078125, 6.233184983047428276974209254934385269649), 
            new TestCase1(257.2001953125, 6.243005711460191965779748682542918291165), 
            new TestCase1(257.7401123046875, 6.245102704489326829989715358050375153164), 
            new TestCase1(261.731201171875, 6.260468857392133508067336988271711266899), 
            new TestCase1(263.75, 6.268152459140511369534885468430435001253), 
            new TestCase1(265.5167236328125, 6.27482855458316535698354618677862660051), 
            new TestCase1(273.9171142578125, 6.305976070434008107321296347063621346142), 
            new TestCase1(278.897705078125, 6.323995460699820076921115676535829705602), 
            new TestCase1(279.167236328125, 6.324961403980196616927106861416582346957), 
            new TestCase1(292.207275390625, 6.370613506132747051155136954016356169709), 
            new TestCase1(293.5975341796875, 6.375359978930308922973172601759003753776), 
            new TestCase1(293.9749755859375, 6.37664472001459930030121377025213179339), 
            new TestCase1(295.1998291015625, 6.380802563199264354475951774362027757146), 
            new TestCase1(297.2799072265625, 6.387824152942428878634062028172724987619), 
            new TestCase1(297.9285888671875, 6.390003820200830500292592492268074658657), 
            new TestCase1(298.1058349609375, 6.390598568067900196600610103133530143627), 
            new TestCase1(300.2803955078125, 6.397866642974941387667911791807820745523), 
            new TestCase1(307.531005859375, 6.421725738171607321147138767579512701297), 
            new TestCase1(308.1754150390625, 6.423818963102848059254801023392818581651), 
            new TestCase1(309.7344970703125, 6.42886525591175973950489819292419777646), 
            new TestCase1(314.2847900390625, 6.443449261058926842539512498259875923692), 
            new TestCase1(314.7236328125, 6.444844602076255234209250709648120853169), 
            new TestCase1(320.8406982421875, 6.464094341970106436820739729174428145587), 
            new TestCase1(321.2459716796875, 6.465356699668166045068069215854964871388), 
            new TestCase1(321.9031982421875, 6.467400466944124604717633061439316010097), 
            new TestCase1(323.457763671875, 6.472218114936839319482664927965675017022), 
            new TestCase1(330.82861328125, 6.494749921382326159049677218709809728653), 
            new TestCase1(335.008544921875, 6.507305446835735629461798195367913785079), 
            new TestCase1(340.7171630859375, 6.52420203343567514995053263208100130883), 
            new TestCase1(348.4677734375, 6.546694993078936223411278280975538973745), 
            new TestCase1(349.1292724609375, 6.548591493378012538030220208753185699103), 
            new TestCase1(372.4288330078125, 6.613194950203131899741436452253747959432), 
            new TestCase1(376.7574462890625, 6.624750543633906159543027563737013064495), 
            new TestCase1(378.4306640625, 6.629181796246806383589594972543857520438), 
            new TestCase1(390.9031982421875, 6.661608771130218794653929942801159149921), 
            new TestCase1(405.7918701171875, 6.698989091751706872771950506377356853297), 
            new TestCase1(407.3646240234375, 6.702857353572475083661986001523177171689), 
            new TestCase1(413.3758544921875, 6.717505881986416333938861696792809072989), 
            new TestCase1(415.7354736328125, 6.723197804327891152771611004762978333262), 
            new TestCase1(417.193603515625, 6.726699007993022779613295744938090059081), 
            new TestCase1(420.874755859375, 6.735483889307782773060913517308358177287), 
            new TestCase1(429.2635498046875, 6.755219602793124098109976875531119009337), 
            new TestCase1(429.756103515625, 6.756366380816258251898421006493740831209), 
            new TestCase1(433.9931640625, 6.766177290841292993057892254025351144187), 
            new TestCase1(434.0106201171875, 6.766217511883345263895227495273810179399), 
            new TestCase1(440.073974609375, 6.780091308338912469861271641540800609593), 
            new TestCase1(450.2220458984375, 6.802889310303153472760070274877750654792), 
            new TestCase1(455.017578125, 6.813484439494547291485100306142865313692), 
            new TestCase1(457.1668701171875, 6.818196843455478403993903909497663695391), 
            new TestCase1(457.5068359375, 6.818940201487998386877795327256997263141), 
            new TestCase1(459.2913818359375, 6.822833193143804950038640831090638344206), 
            new TestCase1(459.492431640625, 6.823270835445770195995146570284994476855), 
            new TestCase1(459.743896484375, 6.823817951018000432957797403476399271545), 
            new TestCase1(464.888427734375, 6.834945773756887582002201745232432180165), 
            new TestCase1(464.96630859375, 6.835113285253827054382848168952118996735), 
            new TestCase1(467.6949462890625, 6.840964582694129262617973631406633110533), 
            new TestCase1(468.86767578125, 6.843468905210339769583953244362473799296), 
            new TestCase1(470.5927734375, 6.847141429556456346098564888170408445454), 
            new TestCase1(481.109619140625, 6.869243403190376592226057897975831923528), 
            new TestCase1(487.4595947265625, 6.882355637062963925878474710879054534122), 
            new TestCase1(488.521484375, 6.884531678915821025400961463061609241284), 
            new TestCase1(492.8812255859375, 6.893416432937340181639026008051884214636), 
            new TestCase1(494.0684814453125, 6.895822338701104787981876142689381248646), 
            new TestCase1(496.4613037109375, 6.900653737167637608469868350257964416187), 
            new TestCase1(716.154052734375, 7.2670429692740963323840103680788401489), 
            new TestCase1(1799.92578125, 8.18864796812207220842639194214816612374), 
            new TestCase1(3564.845703125, 8.872023251113288479644702153534943187411), 
            new TestCase1(7139.869140625, 9.566596912986166722124065953497502737374), 
            new TestCase1(12081.22265625, 10.09255486190560839163513867694963651638), 
            new TestCase1(22810.2421875, 10.72811211386442708825132311411659945728), 
            new TestCase1(46598.96875, 11.44248087071561781484413254045580909559), 
            new TestCase1(108493.375, 12.2875915707717694181967755881243802905), 
            new TestCase1(153860.8125, 12.63695083834421849044542638165059455958), 
            new TestCase1(307019.5, 13.32781372303006380727775468144124321804), 
            new TestCase1(682577.25, 14.12677816700977652906247822629132152831), 
            new TestCase1(1788919.0, 15.09026926533497056732451999718664988024), 
            new TestCase1(3769169.0, 15.83551229128394411859348316953904317643), 
            new TestCase1(4327820.0, 15.97372168955474121405849637207319473463), 
            new TestCase1(11044024.0, 16.91054720571544732784968725481341574805), 
            new TestCase1(21423208.0, 17.57313255890322504472542433762806340181), 
            new TestCase1(62828288.0, 18.64906315643796382994112822618240813581), 
            new TestCase1(70207360.0, 18.7601108873651530393019317938988457707), 
            new TestCase1(154231424.0, 19.54711196618087428636344348941647974165), 
            new TestCase1(294509056.0, 20.1939674915675244205666343569672677746), 
            new TestCase1(1070557184.0, 21.48459226315622322979093631189365864251), 
            new TestCase1(1957922816.0, 22.08829714102155481686957732827526425736), 
            new TestCase1(3912507392.0, 22.78059146269991677250675041832419710247), 
            new TestCase1(7279233024.0, 23.40143852031869095098313030835785443555), 
            new TestCase1(9665245184.0, 23.68494949808078517275225570625127768937), 
            new TestCase1(22627590144.0, 24.53558298204260156687347955127694595939), 
            new TestCase1(60601991168.0, 25.52074076759958328225601044752066239618), 
            new TestCase1(134018236416.0, 26.31438890085421876104087786768111623882), 
            new TestCase1(204864946176.0, 26.73876398039978985836947169876252884996), 
            new TestCase1(284346286080.0, 27.06660583008717947194092254982253381033), 
            new TestCase1(914576637952.0, 28.23487428494463574299686107465880501208), 
            new TestCase1(1581915832320.0, 28.7828049610810604762091293247739091717),


            // Added the following test cases to check for overflow
            // Results tested against Wolfram
            new TestCase1(double.MaxValue/2, 709.78271289338399633830450547794019812841279391420417),
            new TestCase1(double.MaxValue, 710.47586007394394164772173759939837469648829404856442)
            
       };

        #endregion

        /// <summary>
        ///A test for Asinh
        ///</summary>
        [TestMethod()]
        public void AsinhTest()
        {

            TestCaseSet<double>[] testCases = {
                new TestCaseSet<double>(_AsinhData, "Asinh: Random Data", 10)
            };

            NumericFunctionTest.RunSet(Math2.Asinh, "Asinh", testCases);

        }


    }
}