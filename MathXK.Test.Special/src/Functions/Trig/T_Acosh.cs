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
    public class T_Acosh
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

        private static readonly TestCase1[] _AcoshData = {
          new TestCase1(1.000001430511474609375, 0.001691455665129294448190238354291760044433), 
          new TestCase1(1.0000019073486328125, 0.001953124689559275021527821917817027620963), 
          new TestCase1(1.000007152557373046875, 0.003782208044661295168504799813496158490314), 
          new TestCase1(1.000013828277587890625, 0.005258943946801101349061072655743616330534), 
          new TestCase1(1.0000171661376953125, 0.005859366618129202398694086527594451883545), 
          new TestCase1(1.00006008148193359375, 0.0109618319921888517811096976159923461784), 
          new TestCase1(1.000116825103759765625, 0.01528547213183042467192017645636643040682), 
          new TestCase1(1.000148773193359375, 0.01724931909352987823311560583970196658141), 
          new TestCase1(1.000398159027099609375, 0.02821817173865537359266716853098519889415), 
          new TestCase1(1.000638484954833984375, 0.03573281468231456624811499142438796295686), 
          new TestCase1(1.001071453094482421875, 0.04628740247287877599360006621134226755174), 
          new TestCase1(1.003021717071533203125, 0.07771996527168971264969648279358369972575), 
          new TestCase1(1.004993915557861328125, 0.09989759308602780434912638996550489375369), 
          new TestCase1(1.00928401947021484375, 0.1361593876803246479600670434737716450022), 
          new TestCase1(1.024169921875, 0.2194227900495835483852561715845842241518), 
          new TestCase1(1.062277317047119140625, 0.3511165938166054588185413287563455693446), 
          new TestCase1(1.12234401702880859375, 0.4897502671128818535428474267470966393752), 
          new TestCase1(1.2495574951171875, 0.6925568837084910405419269283192900693752), 
          new TestCase1(1.491221904754638671875, 0.9545305722214140465734705961617555409538), 
          new TestCase1(1.983847141265869140625, 1.307581416910453029674373062377350125402), 
          new TestCase1(2.15761280059814453125, 1.403518877974133434572205965467405839077), 
          new TestCase1(2.40639781951904296875, 1.525007084542751786819384889715403191957), 
          new TestCase1(3.38695812225341796875, 1.890537201307279875549078665860235683411), 
          new TestCase1(4.4516773223876953125, 2.173567339994825397387431806918552342932), 
          new TestCase1(6.9391326904296875, 2.625091127868242256287879402513352014572), 
          new TestCase1(7.756023406982421875, 2.737434918695162165715461546271032662695), 
          new TestCase1(8.8823699951171875, 2.874031716780194544439172840789924579634), 
          new TestCase1(9.869171142578125, 2.979986393289490221624555712675426245527), 
          new TestCase1(16.848876953125, 3.516549380542481075157493697729821147595), 
          new TestCase1(16.88458251953125, 3.518670034680249794623612003576884164306), 
          new TestCase1(18.18859100341796875, 3.593185165198828891634086870735797352995), 
          new TestCase1(18.82012176513671875, 3.627367214296338506596699897092700261917), 
          new TestCase1(19.18418121337890625, 3.646553244410946142822321573885913155251), 
          new TestCase1(24.039520263671875, 3.872413451393967155852229598937671193827), 
          new TestCase1(26.5569915771484375, 3.97208556893332931613088010770137243142), 
          new TestCase1(27.92110443115234375, 4.022209178119237972634584536383754567227), 
          new TestCase1(32.314666748046875, 4.168428891496629419926716002725343213186), 
          new TestCase1(34.7300872802734375, 4.24054622986100481621284684937772318866), 
          new TestCase1(36.51556396484375, 4.290698214968890417003449585503652329902), 
          new TestCase1(38.851287841796875, 4.352722738491573736218790864551080662126), 
          new TestCase1(49.46875, 4.59438616262944868606926670858880728888), 
          new TestCase1(49.6726531982421875, 4.598500387004538200979463375829093317529), 
          new TestCase1(55.821014404296875, 4.715217340185609026248077357993388963859), 
          new TestCase1(57.119781494140625, 4.738221040019820009132180121068224048563), 
          new TestCase1(60.3798370361328125, 4.793733825338028989056646578701956447074), 
          new TestCase1(63.4661865234375, 4.843592376953016901945130034442741537681), 
          new TestCase1(63.822418212890625, 4.849190310904695081724453083230505499), 
          new TestCase1(64.366424560546875, 4.857678972284480806836897045666147581194), 
          new TestCase1(65.82318115234375, 4.880061548144127001309581646898589070845), 
          new TestCase1(68.60302734375, 4.921430721025434361496543279369773904556), 
          new TestCase1(70.173583984375, 4.944068352080570156852448111271304401145), 
          new TestCase1(71.80126953125, 4.967000841791218009854450984654955748527), 
          new TestCase1(75.407867431640625, 5.016014824864731856945880331601344083118), 
          new TestCase1(75.497711181640625, 5.017205657609766266706283982466292758789), 
          new TestCase1(78.06475830078125, 5.050644871655082237287791451581061020693), 
          new TestCase1(79.64892578125, 5.070736320140527520131044330827542555678), 
          new TestCase1(79.8707275390625, 5.073517411135062274633407326047677797493), 
          new TestCase1(82.14324951171875, 5.101574796209937553992594384708408566472), 
          new TestCase1(86.422149658203125, 5.152357710985635411723269669070507384393), 
          new TestCase1(87.758697509765625, 5.167705692500116617668438212039278403675), 
          new TestCase1(94.249420166015625, 5.239063709802807411774762976165711451289), 
          new TestCase1(95.002593994140625, 5.24702367651990363562135237922593081374), 
          new TestCase1(96.06402587890625, 5.258134994273664228925894855596706922092), 
          new TestCase1(99.10101318359375, 5.289261389093961411474730159621540149161), 
          new TestCase1(104.825958251953125, 5.345425863147170918152692907733404160171), 
          new TestCase1(105.894317626953125, 5.355566478724588805787125327460047194459), 
          new TestCase1(106.750244140625, 5.363617180711894834893456209993366372702), 
          new TestCase1(109.40167236328125, 5.388152468690488330936438928871930644775), 
          new TestCase1(111.295989990234375, 5.40532022596301299871913818031312501789), 
          new TestCase1(112.682159423828125, 5.417698597745428582703028438959417637741), 
          new TestCase1(115.847869873046875, 5.445406415908932972979361284782173288771), 
          new TestCase1(122.518951416015625, 5.501396249028249309478903043850515943644), 
          new TestCase1(126.29083251953125, 5.531718947357248141901016725943097821463), 
          new TestCase1(127.88677978515625, 5.544277233951786833046255699209839124282), 
          new TestCase1(128.29241943359375, 5.547444176085567387868629933507681372844), 
          new TestCase1(129.49658203125, 5.556786759298987489018162005929203901677), 
          new TestCase1(138.73651123046875, 5.625710723366437130110566517151979116591), 
          new TestCase1(139.18450927734375, 5.628934733085021889655473009134778907233), 
          new TestCase1(139.9705810546875, 5.63456668505549125187745033695013277031), 
          new TestCase1(143.6336669921875, 5.660401141376928390068911580516739480767), 
          new TestCase1(149.2176513671875, 5.698541939965667319302405718431306124964), 
          new TestCase1(150.61602783203125, 5.707869896181299042045964962097048877529), 
          new TestCase1(151.65460205078125, 5.714741890601692470131657476637875234333), 
          new TestCase1(154.77532958984375, 5.735111323217677235163096422075933590581), 
          new TestCase1(158.9586181640625, 5.76178119164116056169664536557902282984), 
          new TestCase1(159.23260498046875, 5.76350337802895912932271749968885321291), 
          new TestCase1(166.89166259765625, 5.810483079631769347906173026231591180669), 
          new TestCase1(169.22418212890625, 5.824362807770767372309724658583543187687), 
          new TestCase1(170.85247802734375, 5.833939098607024938518744061528880726076), 
          new TestCase1(175.641845703125, 5.861586030831371111277316749706184493497), 
          new TestCase1(176.47808837890625, 5.866335876872543841804521246497450515591), 
          new TestCase1(177.0284423828125, 5.869449614294116474831852726381675383803), 
          new TestCase1(178.81622314453125, 5.879497954012966336157412623485146460234), 
          new TestCase1(181.28570556640625, 5.893213844044450514365554597289542581697), 
          new TestCase1(190.84246826171875, 5.944588630523773589353947366401268565142), 
          new TestCase1(191.39764404296875, 5.94749352592071287963619396611460598163), 
          new TestCase1(194.2606201171875, 5.962341215900494050109187978733428591562), 
          new TestCase1(194.89630126953125, 5.965608227627600046293288118981009401155), 
          new TestCase1(196.72125244140625, 5.974928484931286802822226328854017998272), 
          new TestCase1(196.76788330078125, 5.97516550017620215691682414383193227447), 
          new TestCase1(198.0592041015625, 5.981706804024238659640680772319516752432), 
          new TestCase1(199.97052001953125, 5.991310884439669647644011374615753552043), 
          new TestCase1(202.70001220703125, 6.004868209578553896003834136537443847497), 
          new TestCase1(204.95684814453125, 6.015940689286515853596140406483086930665), 
          new TestCase1(206.92059326171875, 6.025476453825986378650455090165700700781), 
          new TestCase1(211.4588623046875, 6.047172064627678522672151564001319932574), 
          new TestCase1(211.6217041015625, 6.047941864223159215142104276243607189411), 
          new TestCase1(212.15936279296875, 6.050479329955437674042299338123601544698), 
          new TestCase1(219.93341064453125, 6.086466833749718691309844243695794396376), 
          new TestCase1(223.34747314453125, 6.101870903204913623855234973179987714517), 
          new TestCase1(228.56036376953125, 6.12494274439855238424605675602915271015), 
          new TestCase1(229.53656005859375, 6.129204755426344240698049588914342485212), 
          new TestCase1(231.15753173828125, 6.136241935513705416652882857399712029477), 
          new TestCase1(235.22589111328125, 6.153688953514382528837249519861026992294), 
          new TestCase1(237.17108154296875, 6.16192447986332112027066888322817886098), 
          new TestCase1(237.904541015625, 6.165012268502458223847465413629334545746), 
          new TestCase1(243.202392578125, 6.187036941752031847781237018817248226575), 
          new TestCase1(244.296875, 6.191527178125453588013362995547294190556), 
          new TestCase1(245.39239501953125, 6.196001570568187444354689455736723305464), 
          new TestCase1(245.80389404296875, 6.197677082130341298614824238362609944761), 
          new TestCase1(249.68365478515625, 6.213337906126028634019511362668656048871), 
          new TestCase1(252.32763671875, 6.223871642756904989009941455863848621693), 
          new TestCase1(253.4725341796875, 6.228398760115368978419785861340394334949), 
          new TestCase1(264.1583251953125, 6.269692237869834655812588535960729291645), 
          new TestCase1(265.867919921875, 6.276143287577458434544566538806340395227), 
          new TestCase1(273.893798828125, 6.305884283737175999056149819356227140283), 
          new TestCase1(274.060546875, 6.306492908028796208195612980534832526609), 
          new TestCase1(274.06298828125, 6.306501816321711306242744764166827854238), 
          new TestCase1(275.31201171875, 6.311048924823309539122894248298646711515), 
          new TestCase1(281.2171630859375, 6.332271212543191917887016682596127229871), 
          new TestCase1(284.3428955078125, 6.343324976847916523668021506516982451804), 
          new TestCase1(284.8428955078125, 6.345081883725142172010442438360996090607), 
          new TestCase1(287.3035888671875, 6.353683609448095450690653155116609411345), 
          new TestCase1(290.8973388671875, 6.366114643735996187010064226141367713494), 
          new TestCase1(293.0467529296875, 6.373476431987165195009419682926258568514), 
          new TestCase1(293.048583984375, 6.3734826803404046607923256022891007467), 
          new TestCase1(296.819091796875, 6.386267177599691603449621482941540136447), 
          new TestCase1(297.6572265625, 6.389086936901673374370938395043767464615), 
          new TestCase1(308.40625, 6.424562459508494815578189029717236448261), 
          new TestCase1(316.5472412109375, 6.450617177370153067374226847837008577358), 
          new TestCase1(320.2418212890625, 6.462221144761522067224169792659123404149), 
          new TestCase1(322.33642578125, 6.468740575092417615259160797033210213863), 
          new TestCase1(323.5101318359375, 6.472375224718483717268648919441488851024), 
          new TestCase1(327.8939208984375, 6.485834999462653559161623555380871809288), 
          new TestCase1(328.0833740234375, 6.486412623146554512400295875149968123595), 
          new TestCase1(328.214599609375, 6.486812521370483270051173204357170113696), 
          new TestCase1(332.13916015625, 6.498698952535686590971725425771901435577), 
          new TestCase1(339.6888427734375, 6.521175044233962334855405597829204127829), 
          new TestCase1(340.171630859375, 6.522595306993372841176959822169835487203), 
          new TestCase1(340.22998046875, 6.522766822935215301119395768816736464171), 
          new TestCase1(340.9984130859375, 6.525022854134450397488249525584977969848), 
          new TestCase1(347.719482421875, 6.544541182598698837108180066568285561192), 
          new TestCase1(347.921142578125, 6.545120967585682780158986289403626335039), 
          new TestCase1(349.8392333984375, 6.550618853671590396954212923343777795426), 
          new TestCase1(353.1812744140625, 6.560126626713879038902701972848577987576), 
          new TestCase1(353.3170166015625, 6.560510895819138848557847535694659145744), 
          new TestCase1(354.9730224609375, 6.56518699003913475265728787337345634833), 
          new TestCase1(355.6412353515625, 6.567067660815945254763768403011524099801), 
          new TestCase1(363.193603515625, 6.588081320423385926941395650884632048324), 
          new TestCase1(363.7503662109375, 6.589613116365141119507599520545813820205), 
          new TestCase1(366.66650390625, 6.597598047275183820713128823864925712214), 
          new TestCase1(370.5828857421875, 6.608222493065004674432923629831401595644), 
          new TestCase1(371.822998046875, 6.611563301604297330780249554012149222974), 
          new TestCase1(375.8822021484375, 6.622421213257872616605790226875280590718), 
          new TestCase1(377.1107177734375, 6.625684248051367798967594517111734939065), 
          new TestCase1(377.588623046875, 6.626950731244344518396423054010667385835), 
          new TestCase1(378.8428955078125, 6.630267034079058832474904688332205765807), 
          new TestCase1(379.1123046875, 6.630977920761718188663355128303457269644), 
          new TestCase1(381.1038818359375, 6.636217452968849140199299619634845177402), 
          new TestCase1(382.1112060546875, 6.63885714989915892916806173377813464616), 
          new TestCase1(382.9927978515625, 6.641161660644278254856011296094169025726), 
          new TestCase1(387.1845703125, 6.652047018118426624071335253039983067027), 
          new TestCase1(389.669921875, 6.658445560711747733097363184134660374201), 
          new TestCase1(389.804443359375, 6.658790721334144459483338624834811878532), 
          new TestCase1(396.3114013671875, 6.675345858154136306174834824525273599533), 
          new TestCase1(397.005126953125, 6.677094789236718239327386065258295745882), 
          new TestCase1(397.1934814453125, 6.677569116668089765968288575699796076917), 
          new TestCase1(397.8046875, 6.679106750673113062518913783394959687397), 
          new TestCase1(398.8426513671875, 6.681712590609844816466059398426090860379), 
          new TestCase1(399.1663818359375, 6.682523938576487571986006995202709018951), 
          new TestCase1(399.2547607421875, 6.682745323455160373493982153327465350514), 
          new TestCase1(400.33984375, 6.68545941647717813159478393220239518843), 
          new TestCase1(403.9578857421875, 6.6944562778394978005276529277337054901), 
          new TestCase1(404.279541015625, 6.69525222285407608719250824812439781544), 
          new TestCase1(405.0574951171875, 6.697174677114241660699777628750006369949), 
          new TestCase1(407.328125, 6.702764738337774266407206003201536960489), 
          new TestCase1(407.547119140625, 6.703302231179959306927821131625588949855), 
          new TestCase1(410.5994873046875, 6.710763953621196143396177253879091429692), 
          new TestCase1(410.8016357421875, 6.711256159037372934708585809074885390587), 
          new TestCase1(411.129638671875, 6.712054288828398871484168355341928910672), 
          new TestCase1(411.9053955078125, 6.71393940750234669898410033879058990879), 
          new TestCase1(415.5833740234375, 6.722828986708716597641935076657699092718), 
          new TestCase1(417.669189453125, 6.727835453862131466839432640652077316451), 
          new TestCase1(420.517822265625, 6.734632628835640647510302384260878716087), 
          new TestCase1(424.3853759765625, 6.743787740494532100937084679719803238806), 
          new TestCase1(424.7154541015625, 6.744565219553757289338296980369793898726), 
          new TestCase1(436.3419189453125, 6.771572021268065472972105624143419455948), 
          new TestCase1(438.501953125, 6.776510146304200360440233381235270942742), 
          new TestCase1(439.3369140625, 6.778412462065226061129312213528458641603), 
          new TestCase1(445.5606689453125, 6.792479340600349658156225315328718026016), 
          new TestCase1(452.9901123046875, 6.809016260337228831936738732897410136218), 
          new TestCase1(453.77490234375, 6.810747231716348697176150895548602377257), 
          new TestCase1(456.7745361328125, 6.817335895109250694056391088138248559173), 
          new TestCase1(457.9520263671875, 6.819910421197310978529029433663908470507), 
          new TestCase1(458.6795654296875, 6.821497844004013502448409287800537132756), 
          new TestCase1(460.5164794921875, 6.825494642872147128332763285398798677976), 
          new TestCase1(461.8717041015625, 6.828433164406686771701521855077910926631), 
          new TestCase1(464.7025146484375, 6.834543470287694007429787493802694241702), 
          new TestCase1(467.0626220703125, 6.839609377592375029558701051277174984946), 
          new TestCase1(467.0712890625, 6.83962793384421245687931998241985864452), 
          new TestCase1(470.096923828125, 6.846084943645238808995975141700921096311), 
          new TestCase1(475.1607666015625, 6.856799276049143563229337105263608232029), 
          new TestCase1(477.5537109375, 6.861822721577315741859023908594798480392), 
          new TestCase1(478.626220703125, 6.864066049482580522441871176709541953581), 
          new TestCase1(478.7958984375, 6.8644204973336804815313774177052514308), 
          new TestCase1(479.6864013671875, 6.86627865397306926299976132455469132911), 
          new TestCase1(479.7867431640625, 6.866487814627139251181098007832885765806), 
          new TestCase1(479.9122314453125, 6.866749331118839213169502276407114725047), 
          new TestCase1(482.4793701171875, 6.872084270243208288559827917865564857217), 
          new TestCase1(482.5181884765625, 6.872164723177874204251489651336382789811), 
          new TestCase1(483.8797607421875, 6.874982560453873383658251901541786496761), 
          new TestCase1(484.4649658203125, 6.876191234145179554349993597841191267853), 
          new TestCase1(485.3258056640625, 6.877966548833207350351265217595172972905), 
          new TestCase1(490.57373046875, 6.888721726428236039221890401406352667404), 
          new TestCase1(493.7423095703125, 6.895159895589969853078001763235883415882), 
          new TestCase1(494.272216796875, 6.896232568812717932403585431953843090099), 
          new TestCase1(496.44775390625, 6.900624415355815565025217006388610220335), 
          new TestCase1(497.0401611328125, 6.901816998553274960978831433918910814375), 
          new TestCase1(498.234130859375, 6.904216282287646841845864948036717548583), 
          new TestCase1(665.0791015625, 7.193052598670792558912351099032198333184), 
          new TestCase1(1170.29150390625, 7.758155143419731595147190189684865359835), 
          new TestCase1(2058.7958984375, 8.323023697145112010072186763207801515396), 
          new TestCase1(5824.533203125, 9.362981311610990187535171280837584820237), 
          new TestCase1(9114.30859375, 9.810748008110925084589570791503788344652), 
          new TestCase1(31388.40625, 11.04734105631420306080680845572519805956), 
          new TestCase1(53732.765625, 11.58492543551253544550738607240822985131), 
          new TestCase1(117455.09375, 12.3669585392073964063765790447351364021), 
          new TestCase1(246210.625, 13.10708982832787479175510220013343945066), 
          new TestCase1(513670.125, 13.84248373881161892648765606597953800527), 
          new TestCase1(788353.25, 14.2708487357510805583139091933979145572), 
          new TestCase1(1736171.0, 15.06033985221540896276800056693696357419), 
          new TestCase1(3770530.0, 15.83587331365755605568701554508826253991), 
          new TestCase1(4344090.0, 15.97747403917326507392988269405762331118), 
          new TestCase1(11419360.0, 16.94396789915014469006155507384872743685), 
          new TestCase1(31023240.0, 17.9433943395609684013300439026124645108), 
          new TestCase1(40665424.0, 18.21403593674543079149440340324547029148), 
          new TestCase1(129788064.0, 19.37456058170921543648305825985098346634), 
          new TestCase1(225668224.0, 19.9277236237785460971743917048171789264), 
          new TestCase1(450631936.0, 20.61930863840059519949232655148142693395), 
          new TestCase1(750941952.0, 21.1299860930266968673735341410555546727), 
          new TestCase1(1887358976.0, 22.05159150215413004649732665310838648331), 
          new TestCase1(3738011648.0, 22.73496684263974142690946105862900518113), 
          new TestCase1(7486695424.0, 23.42954051928097083366085869187241049078), 
          new TestCase1(12668080128.0, 23.95549847139166892348970713379911103008), 
          new TestCase1(23918272512.0, 24.59105572458284785473015955570721968617), 
          new TestCase1(48862560256.0, 25.30542448179939429678484289438720749573), 
          new TestCase1(113763549184.0, 26.15053518194943663075151882301281675861), 
          new TestCase1(161334755328.0, 26.49989444953256419131899378956019192586), 
          new TestCase1(321933279232.0, 27.19075733422631778452130090694771384112), 
          new TestCase1(715734122496.0, 27.98972177820814613504868209029528272152), 
          new TestCase1(1875817529344.0, 28.95321287653379863631300659442341812344),

          // added from functions.wolfram.com
          new TestCase1(double.MaxValue, 710.47586007394394164772173759939837469648829404856442),
          
        };

        #endregion

        /// <summary>
        ///A test for Asinh
        ///</summary>
        [TestMethod()]
        public void AcoshTest()
        {

            TestCaseSet<double>[] testCases = {
                new TestCaseSet<double>(_AcoshData, "Acosh: Random Data", 10)
            };

            NumericFunctionTest.RunSet(Math2.Acosh, "Acosh", testCases);

        }



    }
}
