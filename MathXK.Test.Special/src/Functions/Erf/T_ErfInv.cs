﻿//  Copyright (c) 2013 Rocco Romeo
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

//  Test Data:
//      Copyright (c) 2006-7 John Maddock, Boost Software License v1.0

using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace MathXK.Test.Functions
{

    using TestCase1 = TestCase<double>;

    /// <summary>
    ///This is a test class for Math2Test and is intended
    ///to contain all Math2Test Unit Tests
    ///</summary>
    [TestClass()]
    public class T_ErfInv
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


        private static readonly TestCase1[] _ErfInvData = {
            new TestCase1(-0.990433037281036376953125, -1.832184533179510927322805923563700329767), 
            new TestCase1(-0.936334311962127685546875, -1.311339282092737086640055105484822812599), 
            new TestCase1(-0.931107819080352783203125, -1.286316461685184857889337272829270644576), 
            new TestCase1(-0.928576648235321044921875, -1.274755355308702535979544636706295190547), 
            new TestCase1(-0.92711746692657470703125, -1.268242390461126936128446743938528753699), 
            new TestCase1(-0.907657206058502197265625, -1.190178701802009872651528128114629153367), 
            new TestCase1(-0.89756715297698974609375, -1.154826929034364581020764782706068245251), 
            new TestCase1(-0.80573642253875732421875, -0.917873491512746130598510739795171676962), 
            new TestCase1(-0.804919183254241943359375, -0.9161942446913841771580833184012189161559), 
            new TestCase1(-0.780276477336883544921875, -0.867806431357551134410815525153998291827), 
            new TestCase1(-0.775070965290069580078125, -0.8580919924152284867224297625328485999768), 
            new TestCase1(-0.7496345043182373046875, -0.8127924290810407931222432768905514928867), 
            new TestCase1(-0.74820673465728759765625, -0.8103475955423936417307157186107256388648), 
            new TestCase1(-0.74602639675140380859375, -0.806632688757205819001462030577472890805), 
            new TestCase1(-0.72904598712921142578125, -0.7784313874823551106598826200232666844639), 
            new TestCase1(-0.7162272930145263671875, -0.7579355487144890063429377586715787838945), 
            new TestCase1(-0.701772034168243408203125, -0.7355614373173595299453301005770428516518), 
            new TestCase1(-0.68477380275726318359375, -0.7101588399949052872532446197423489724803), 
            new TestCase1(-0.657626628875732421875, -0.6713881533266128640126408255047881638389), 
            new TestCase1(-0.652269661426544189453125, -0.6639738263692763456669198307149427734317), 
            new TestCase1(-0.6262547969818115234375, -0.6289572925573171740836428308746584439674), 
            new TestCase1(-0.62323606014251708984375, -0.6249936843093662471116097431474787933967), 
            new TestCase1(-0.57958185672760009765625, -0.5697131213589784467617578394703976041604), 
            new TestCase1(-0.576151371002197265625, -0.5655172244109430153330195150336141500139), 
            new TestCase1(-0.5579319000244140625, -0.5435569422159360687847790186563654276103), 
            new TestCase1(-0.446154057979583740234375, -0.4186121208546731033057205459902879301199), 
            new TestCase1(-0.44300353527069091796875, -0.4152898953738801984047941692529271391195), 
            new TestCase1(-0.40594112873077392578125, -0.3768620801611051992528860948080812212023), 
            new TestCase1(-0.396173775196075439453125, -0.3669220210390825311547962776125822899061), 
            new TestCase1(-0.38366591930389404296875, -0.3542977152760563782151668726041057557165), 
            new TestCase1(-0.36689913272857666015625, -0.3375493847053488720470432821496358279516), 
            new TestCase1(-0.365801036357879638671875, -0.3364591774366710656954166264654945559873), 
            new TestCase1(-0.277411997318267822265625, -0.2510244067029671790889794981353227476998), 
            new TestCase1(-0.236883103847503662109375, -0.213115119330839975829499967157244997714), 
            new TestCase1(-0.215545952320098876953125, -0.1934073617841803428235669261097060281642), 
            new TestCase1(-0.202522933483123779296875, -0.1814532246720147926398927046057793150106), 
            new TestCase1(-0.18253767490386962890625, -0.1632073953550647568421821286058243218715), 
            new TestCase1(-0.156477451324462890625, -0.1395756320903277910768376053314442757507), 
            new TestCase1(-0.1558246612548828125, -0.1389857795955756484965030151195660030168), 
            new TestCase1(-0.12251126766204833984375, -0.109002961098867662134935094105847496074), 
            new TestCase1(-0.1088275909423828125, -0.09674694516640724629590870677194632943569), 
            new TestCase1(-0.08402168750762939453125, -0.07460044047654119877070700345343119035515), 
            new TestCase1(-0.05048263072967529296875, -0.04476895818328636312384562686715995170129), 
            new TestCase1(-0.029248714447021484375, -0.0259268064334840921104659134138093242797), 
            new TestCase1(-0.02486217021942138671875, -0.02203709146986755832638577823832075055744), 
            new TestCase1(-0.02047121524810791015625, -0.01814413302702029459097557481591610553903), 
            new TestCase1(-0.018821895122528076171875, -0.01668201759439857888105181293763417899072), 
            new TestCase1(0.0073254108428955078125, 0.006492067534753215749601670217642082465642), 
            new TestCase1(0.09376299381256103515625, 0.08328747254794857150987333986733043734817), 
            new TestCase1(0.0944411754608154296875, 0.08389270963798942778622198997355058545872), 
            new TestCase1(0.264718532562255859375, 0.2390787735821979528028028789569770109829), 
            new TestCase1(0.27952671051025390625, 0.2530214201700340392837551955289041822603), 
            new TestCase1(0.29262602329254150390625, 0.2654374523135675523971788948011709467352), 
            new TestCase1(0.3109557628631591796875, 0.282950508503826367238408926581528085458), 
            new TestCase1(0.31148135662078857421875, 0.2834552014554130441860525970673030809536), 
            new TestCase1(0.32721102237701416015625, 0.2986277427848421570858990348074985028421), 
            new TestCase1(0.3574702739715576171875, 0.3282140305634627092431945088114761850208), 
            new TestCase1(0.362719058990478515625, 0.3334035993712283467959295804559099468454), 
            new TestCase1(0.3896572589874267578125, 0.3603304982893212173104266596348905175268), 
            new TestCase1(0.4120922088623046875, 0.3831602323665075533579267768785894144888), 
            new TestCase1(0.41872966289520263671875, 0.3899906753567599452444107492361433402154), 
            new TestCase1(0.45167791843414306640625, 0.4244594733907945411184647153213164209335), 
            new TestCase1(0.48129451274871826171875, 0.4563258063707025027210352963461819167707), 
            new TestCase1(0.4862649440765380859375, 0.4617640058971775089811390737537561779898), 
            new TestCase1(0.50937330722808837890625, 0.4874174763856674076219106695373814892182), 
            new TestCase1(0.5154802799224853515625, 0.4943041993872143888987628020569772222018), 
            new TestCase1(0.52750003337860107421875, 0.5079978091910991117615000459548117088362), 
            new TestCase1(0.53103363513946533203125, 0.5120597685873370942783226077302881881069), 
            new TestCase1(0.58441460132598876953125, 0.5756584292527058478710392476034273328569), 
            new TestCase1(0.5879499912261962890625, 0.5800336103175463592377907341030447077804), 
            new TestCase1(0.59039986133575439453125, 0.5830784871670823806198622501806646778319), 
            new TestCase1(0.59455978870391845703125, 0.588273673825686998734497652983815773459), 
            new TestCase1(0.59585726261138916015625, 0.5899005483108011364541949539839185473259), 
            new TestCase1(0.5962116718292236328125, 0.5903454775096607218832535637355431851718), 
            new TestCase1(0.6005609035491943359375, 0.5958247243549040349587326482492767206448), 
            new TestCase1(0.6150619983673095703125, 0.6143583249050861028039832921829036722514), 
            new TestCase1(0.62944734096527099609375, 0.6331707263097125575937994856370309207836), 
            new TestCase1(0.64380657672882080078125, 0.6524069265890823819975498133014027247554), 
            new TestCase1(0.6469156742095947265625, 0.656635855345815020063728463464436343698), 
            new TestCase1(0.67001712322235107421875, 0.6888269167957872563013714303376548038671), 
            new TestCase1(0.6982586383819580078125, 0.7302336318927408409119676651737758401138), 
            new TestCase1(0.74485766887664794921875, 0.8046505193013635090578266413458426260098), 
            new TestCase1(0.75686132907867431640625, 0.8253191678260588578995203396384711816647), 
            new TestCase1(0.81158387660980224609375, 0.9300427626888758122211127950646282789481), 
            new TestCase1(0.826751708984375, 0.9629665092443368464606966822833571908852), 
            new TestCase1(0.83147108554840087890625, 0.9736479209913771931387473923084901789046), 
            new TestCase1(0.84174954891204833984375, 0.997713670556719074960678197806799852186), 
            new TestCase1(0.8679864406585693359375, 1.065050516333636716777334376076374184102), 
            new TestCase1(0.90044414997100830078125, 1.164612422633086435501625591693259387477), 
            new TestCase1(0.91433393955230712890625, 1.215315881176612875682446995412738776976), 
            new TestCase1(0.91501367092132568359375, 1.217962731073139868794942852653058932976), 
            new TestCase1(0.918984889984130859375, 1.233778505900771488542027767896521427575), 
            new TestCase1(0.92977702617645263671875, 1.28019542575660930623179572273596558907), 
            new TestCase1(0.93538987636566162109375, 1.306695301483797253764522033930388453334), 
            new TestCase1(0.93773555755615234375, 1.318335478463913327121670503572736587296), 
            new TestCase1(0.94118559360504150390625, 1.33613349872692113073358883961598631154), 
            new TestCase1(0.96221935749053955078125, 1.468821071545234761861756248744372345584), 
            new TestCase1(0.98576259613037109375, 1.733272259459038694476413373595347034928), 
            new TestCase1(0.9881370067596435546875, 1.77921769652839903464038407684397479173), 
            new TestCase1(0.99292266368865966796875, 1.904368122482929779094714951471938518496) 
        };



        private static readonly TestCase1[] erfc_inv_data = new TestCase1[] {
      new TestCase1(0.00956696830689907073974609375, 1.832184391051582711731256541599359331735), 
      new TestCase1(0.063665688037872314453125, 1.311339282092737086640055105484822812599), 
      new TestCase1(0.068892158567905426025390625, 1.286316565305373898738127895195338854718), 
      new TestCase1(0.071423359215259552001953125, 1.274755321776344058215704428086211324587), 
      new TestCase1(0.0728825032711029052734375, 1.268242522387371561738393687518023984868), 
      new TestCase1(0.09234277904033660888671875, 1.190178756246535875259766567441510867604), 
      new TestCase1(0.102432854473590850830078125, 1.154826903977823078146497880706118113588), 
      new TestCase1(0.1942635476589202880859375, 0.9178735528443878579995511780412810469667), 
      new TestCase1(0.19508080184459686279296875, 0.9161942752629032646404767631869277618212), 
      new TestCase1(0.21972350776195526123046875, 0.8678064594007161661713535461829067693456), 
      new TestCase1(0.224929034709930419921875, 0.8580919924152284867224297625328485999768), 
      new TestCase1(0.2503655254840850830078125, 0.8127923779477598070926819995714374417663), 
      new TestCase1(0.25179326534271240234375, 0.8103475955423936417307157186107256388648), 
      new TestCase1(0.2539736330509185791015625, 0.8066326381314558738773191719462921998277), 
      new TestCase1(0.27095401287078857421875, 0.7784313874823551106598826200232666844639), 
      new TestCase1(0.2837726771831512451171875, 0.7579355956263440770864088550908621329508), 
      new TestCase1(0.298227965831756591796875, 0.7355614373173595299453301005770428516518), 
      new TestCase1(0.3152261674404144287109375, 0.7101588837290742217667270852502075077668), 
      new TestCase1(0.342373371124267578125, 0.6713881533266128640126408255047881638389), 
      new TestCase1(0.347730338573455810546875, 0.6639738263692763456669198307149427734317), 
      new TestCase1(0.3737452030181884765625, 0.6289572925573171740836428308746584439674), 
      new TestCase1(0.37676393985748291015625, 0.6249936843093662471116097431474787933967), 
      new TestCase1(0.42041814327239990234375, 0.5697131213589784467617578394703976041604), 
      new TestCase1(0.4238486588001251220703125, 0.5655171880456876504494070613171955224472), 
      new TestCase1(0.4420680999755859375, 0.5435569422159360687847790186563654276103), 
      new TestCase1(0.553845942020416259765625, 0.4186121208546731033057205459902879301199), 
      new TestCase1(0.55699646472930908203125, 0.4152898953738801984047941692529271391195), 
      new TestCase1(0.59405887126922607421875, 0.3768620801611051992528860948080812212023), 
      new TestCase1(0.603826224803924560546875, 0.3669220210390825311547962776125822899061), 
      new TestCase1(0.61633408069610595703125, 0.3542977152760563782151668726041057557165), 
      new TestCase1(0.63310086727142333984375, 0.3375493847053488720470432821496358279516), 
      new TestCase1(0.634198963642120361328125, 0.3364591774366710656954166264654945559873), 
      new TestCase1(0.722588002681732177734375, 0.2510244067029671790889794981353227476998), 
      new TestCase1(0.763116896152496337890625, 0.213115119330839975829499967157244997714), 
      new TestCase1(0.784454047679901123046875, 0.1934073617841803428235669261097060281642), 
      new TestCase1(0.797477066516876220703125, 0.1814532246720147926398927046057793150106), 
      new TestCase1(0.81746232509613037109375, 0.1632073953550647568421821286058243218715), 
      new TestCase1(0.843522548675537109375, 0.1395756320903277910768376053314442757507), 
      new TestCase1(0.8441753387451171875, 0.1389857795955756484965030151195660030168), 
      new TestCase1(0.87748873233795166015625, 0.109002961098867662134935094105847496074), 
      new TestCase1(0.8911724090576171875, 0.09674694516640724629590870677194632943569), 
      new TestCase1(0.91597831249237060546875, 0.07460044047654119877070700345343119035515), 
      new TestCase1(0.94951736927032470703125, 0.04476895818328636312384562686715995170129), 
      new TestCase1(0.970751285552978515625, 0.0259268064334840921104659134138093242797), 
      new TestCase1(0.97513782978057861328125, 0.02203709146986755832638577823832075055744), 
      new TestCase1(0.97952878475189208984375, 0.01814413302702029459097557481591610553903), 
      new TestCase1(0.981178104877471923828125, 0.01668201759439857888105181293763417899072), 
      new TestCase1(1.0073254108428955078125, -0.006492067534753215749601670217642082465642), 
      new TestCase1(1.09376299381256103515625, -0.08328747254794857150987333986733043734817), 
      new TestCase1(1.0944411754608154296875, -0.08389270963798942778622198997355058545872), 
      new TestCase1(1.264718532562255859375, -0.2390787735821979528028028789569770109829), 
      new TestCase1(1.27952671051025390625, -0.2530214201700340392837551955289041822603), 
      new TestCase1(1.29262602329254150390625, -0.2654374523135675523971788948011709467352), 
      new TestCase1(1.3109557628631591796875, -0.282950508503826367238408926581528085458), 
      new TestCase1(1.31148135662078857421875, -0.2834552014554130441860525970673030809536), 
      new TestCase1(1.32721102237701416015625, -0.2986277427848421570858990348074985028421), 
      new TestCase1(1.3574702739715576171875, -0.3282140305634627092431945088114761850208), 
      new TestCase1(1.362719058990478515625, -0.3334035993712283467959295804559099468454), 
      new TestCase1(1.3896572589874267578125, -0.3603304982893212173104266596348905175268), 
      new TestCase1(1.4120922088623046875, -0.3831602323665075533579267768785894144888), 
      new TestCase1(1.41872966289520263671875, -0.3899906753567599452444107492361433402154), 
      new TestCase1(1.45167791843414306640625, -0.4244594733907945411184647153213164209335), 
      new TestCase1(1.48129451274871826171875, -0.4563258063707025027210352963461819167707), 
      new TestCase1(1.4862649440765380859375, -0.4617640058971775089811390737537561779898), 
      new TestCase1(1.50937330722808837890625, -0.4874174763856674076219106695373814892182), 
      new TestCase1(1.5154802799224853515625, -0.4943041993872143888987628020569772222018), 
      new TestCase1(1.52750003337860107421875, -0.5079978091910991117615000459548117088362), 
      new TestCase1(1.53103363513946533203125, -0.5120597685873370942783226077302881881069), 
      new TestCase1(1.58441460132598876953125, -0.5756584292527058478710392476034273328569), 
      new TestCase1(1.5879499912261962890625, -0.5800336103175463592377907341030447077804), 
      new TestCase1(1.59039986133575439453125, -0.5830784871670823806198622501806646778319), 
      new TestCase1(1.59455978870391845703125, -0.588273673825686998734497652983815773459), 
      new TestCase1(1.59585726261138916015625, -0.5899005483108011364541949539839185473259), 
      new TestCase1(1.5962116718292236328125, -0.5903454775096607218832535637355431851718), 
      new TestCase1(1.6005609035491943359375, -0.5958247243549040349587326482492767206448), 
      new TestCase1(1.6150619983673095703125, -0.6143583249050861028039832921829036722514), 
      new TestCase1(1.62944734096527099609375, -0.6331707263097125575937994856370309207836), 
      new TestCase1(1.64380657672882080078125, -0.6524069265890823819975498133014027247554), 
      new TestCase1(1.6469156742095947265625, -0.656635855345815020063728463464436343698), 
      new TestCase1(1.67001712322235107421875, -0.6888269167957872563013714303376548038671), 
      new TestCase1(1.6982586383819580078125, -0.7302336318927408409119676651737758401138), 
      new TestCase1(1.74485766887664794921875, -0.8046505193013635090578266413458426260098), 
      new TestCase1(1.75686132907867431640625, -0.8253191678260588578995203396384711816647), 
      new TestCase1(1.81158387660980224609375, -0.9300427626888758122211127950646282789481), 
      new TestCase1(1.826751708984375, -0.9629665092443368464606966822833571908852), 
      new TestCase1(1.83147108554840087890625, -0.9736479209913771931387473923084901789046), 
      new TestCase1(1.84174954891204833984375, -0.997713670556719074960678197806799852186), 
      new TestCase1(1.8679864406585693359375, -1.065050516333636716777334376076374184102), 
      new TestCase1(1.90044414997100830078125, -1.164612422633086435501625591693259387477), 
      new TestCase1(1.91433393955230712890625, -1.215315881176612875682446995412738776976), 
      new TestCase1(1.91501367092132568359375, -1.217962731073139868794942852653058932976), 
      new TestCase1(1.918984889984130859375, -1.233778505900771488542027767896521427575), 
      new TestCase1(1.92977702617645263671875, -1.28019542575660930623179572273596558907), 
      new TestCase1(1.93538987636566162109375, -1.306695301483797253764522033930388453334), 
      new TestCase1(1.93773555755615234375, -1.318335478463913327121670503572736587296), 
      new TestCase1(1.94118559360504150390625, -1.33613349872692113073358883961598631154), 
      new TestCase1(1.96221935749053955078125, -1.468821071545234761861756248744372345584), 
      new TestCase1(1.98576259613037109375, -1.733272259459038694476413373595347034928), 
      new TestCase1(1.9881370067596435546875, -1.77921769652839903464038407684397479173), 
      new TestCase1(1.99292266368865966796875, -1.904368122482929779094714951471938518496) 
   };


    private static readonly TestCase1[] _ErfcInvBigData = {
        new TestCase1(0.1040364789689237257780791904879173249974e-287, 25.67656975866313564221220720210665888424), 
        new TestCase1(0.2450719147172718682892603854216257656211e-282, 25.43473899116627259523415573506406682191), 
        new TestCase1(0.1540635014300999959290867576988859893076e-280, 25.35326739382267978364899192426764093558), 
        new TestCase1(0.1611739209855438499046222515484436760851e-273, 25.03273169525729329624232590800413859205), 
        new TestCase1(0.2781078901633922872827030068281452871724e-257, 24.27512000889987237167959986465966890736), 
        new TestCase1(0.8837820681846044382839818307987545000847e-230, 22.93495456229819340317355998004446065241), 
        new TestCase1(0.2845511846316072211634522373099423142668e-217, 22.29887672117465208331735160799934601323), 
        new TestCase1(0.662011186385949636511655220013184768187e-212, 22.02033498263166080103335807537366542004), 
        new TestCase1(0.1113036584088835972703336205837326051999e-204, 21.63966022176108818726656867477709317659), 
        new TestCase1(0.891278753475701910699424577762349195273e-195, 21.10677883438545215904960900943812596651), 
        new TestCase1(0.1243375606125932867518669565282002251187e-191, 20.93474611896499049089384377937974962862), 
        new TestCase1(0.963481566993839758726177317091390488753e-188, 20.72000418866581389344212197908962674052), 
        new TestCase1(0.1401265332225531798315780528328964001988e-167, 19.56909048260682534539216053987194710384), 
        new TestCase1(0.7560947579095629715718115752684784187336e-154, 18.74494560773490096872968649468667611647), 
        new TestCase1(0.2179858592534951368829782528927445425064e-151, 18.59346814776643084908600446605555382182), 
        new TestCase1(0.5208405693419352209961246987938629746382e-131, 17.28776816936911830963159975455855394785), 
        new TestCase1(0.5376896865772724579278898126622880802843e-127, 17.01882418428649428645068852272615582504), 
        new TestCase1(0.1449951754876351411423239477096839160613e-114, 16.15763261758311680195806763260744183267), 
        new TestCase1(0.3189212836315126918184908007667112228414e-108, 15.70012576642621709878560027430710336585), 
        new TestCase1(0.1297369886602648502231884634897739722895e-107, 15.65546665388935777226018284593032382368), 
        new TestCase1(0.2602821327589283546034644516083144727197e-103, 15.33647743152077165227354683273671003413), 
        new TestCase1(0.1590035594081776210615338926776093483999e-96, 14.81946110818452708436758045284808001241), 
        new TestCase1(0.2720908074515874075065365442027770927174e-96, 14.80136596875698966601900922593357151598), 
        new TestCase1(0.3706750818860597777102228512408751197629e-87, 14.07473190298865115531360325271322276549), 
        new TestCase1(0.3489539664837006993042769801460156282803e-83, 13.74669364199157990777768669458685265472), 
        new TestCase1(0.1947620269151566018960312643171165602404e-79, 13.43009987431500967622734624289673694188), 
        new TestCase1(0.5129832297533069185312436951335270626736e-73, 12.86957639029561313119037589037764922829), 
        new TestCase1(0.9272206234111032518003166213125699095688e-70, 12.57574048509944911276422466059543556873), 
        new TestCase1(0.1460449453933976163637836336160553873706e-68, 12.46599637375040876166599138465603014844), 
        new TestCase1(0.8768204308820258161647995326423963318019e-68, 12.39412855962299302779918858089973570208), 
        new TestCase1(0.3176564512338672251415481009303329209995e-61, 11.77127316747622053532997731227084711173), 
        new TestCase1(0.1172937774279091349820752991146331290881e-54, 11.11297358832626209702547835299261897395), 
        new TestCase1(0.3822700066028246575792241970441049035743e-52, 10.85058828114725350216779162566256390403), 
        new TestCase1(0.4793602326185941396320885530096873382082e-41, 9.607343723542133339487829511083057069273), 
        new TestCase1(0.106201511009181306109371950926034372112e-40, 9.566077632118694585964030456246300276422), 
        new TestCase1(0.1026229395621872573137267432392849554594e-24, 7.413111396559360630306454710337113823294), 
        new TestCase1(0.1433834013674021626308791056277965507134e-19, 6.574536501026477531602186421374903506051), 
        new TestCase1(0.2287484604245551701383736767135709586086e-12, 5.183672646329402197657171322814515766894), 
        new TestCase1(0.4054382118091501052221361988694253843776e-11, 4.903978376537783549083010816442468785255), 
        new TestCase1(0.5517250582467272934968166422809728429772e-9, 4.386634719241854622420456791837035899963), 
        new TestCase1(0.1708920553446964553659243704866375034024e-5, 3.383585107532389209677350898502297644278), 
        new TestCase1(0.2892930668406028955896934416060055101499e-5, 3.308042643013470022433136596357284983766), 
        new TestCase1(0.7029998339557698614154475319387477445787e-5, 3.17687395056466498458108882620439294841), 
        new TestCase1(0.3973417837559780922949941034264323413971e-4, 2.90551585834893682957467771364140403613), 
        new TestCase1(0.04246334897368917726269942070072005435577, 1.434684541881674882915735690265484394337), 
        new TestCase1(0.159920364968560744996532371753339418774, 0.9937250495458864822841753804433979428041), 
        new TestCase1(0.2425390104583346613758947085681683120129, 0.826370276571483139432927467391341431105), 
        new TestCase1(0.7954739362140872788414606986417965117653, 0.1832884885283639734589195991287031281603), 
        new TestCase1(0.9236374866673857278562449757419727802699, 0.06777816075457032123722158048171767182127) 
    };

    #endregion



        /// <summary>
        ///A test for ErfInv
        ///</summary>
        [TestMethod()]
        public void ErfInvTest()
        {

            TestCaseSet<double>[] testCases = {
                new TestCaseSet<double>(_ErfInvData, "ErfInv: Random Data", 20)
            };

            NumericFunctionTest.RunSet(Math2.ErfInv, "ErfInv", testCases);

        }

        /// <summary>
        ///A test for ErfcInv
        ///</summary>
        [TestMethod()]
        public void ErfcInvTest()
        {

            TestCaseSet<double>[] testCases = {
                new TestCaseSet<double>(erfc_inv_data, "ErfcInv: Random Data", 20),
                new TestCaseSet<double>(_ErfcInvBigData, "ErfcInv: Random Data", 20)
            };

            NumericFunctionTest.RunSet(Math2.ErfcInv, "ErfcInv", testCases);
        }
    }
}