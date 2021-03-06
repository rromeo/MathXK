﻿//  Copyright (c) 2013 Rocco Romeo
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

using System.Collections.Generic;
using System.Linq;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace MathXK.Test.Functions
{

    [TestClass]
    public class T_SinCosPI
    {
        struct SinCosPIData
        {
            public readonly double X;
            public readonly double SinPI;
            public readonly double CosPI;

            public SinCosPIData(double x, double sinPI, double cosPI)
            {
                X = x;
                SinPI = sinPI;
                CosPI = cosPI;
            }

        };

        #region Data

        private static readonly SinCosPIData[] _SinCosPIData = {
            new SinCosPIData(0.004394531250000000, 0.013805388528060390723422806806845412778981167540657, 0.99990470108285284456900994449060722831261414989009),
            new SinCosPIData(0.03552246093750000, 0.11136560974333516314402726677449681317545041646980, 0.99377950319298454133085088367286886366556062431413),
            new SinCosPIData(0.05615234375000000, 0.17549425337727142769015884254126488525830802479003, 0.98448045538322092994185044824322327096652825296368),
            new SinCosPIData(0.05822753906250000, 0.18190871836966616979764946668522991204799120205405, 0.98331542151087281612490567299643516206554850405145),
            new SinCosPIData(0.08520507812500000, 0.26449443242780162167711808128837144573353056442639, 0.96438721228285429912390975863694594937938904659940),
            new SinCosPIData(0.1140136718750000, 0.35057454605483757068311792630886273668556154441200, 0.93653482992275550024075034899563495495262308958027),
            new SinCosPIData(0.1166992187500000, 0.35846342063373656486301561775328065594338219200631, 0.93354377297883619205407802547805411171840879627642),
            new SinCosPIData(0.1606445312500000, 0.48352707893291872415908426643412037067824824405302, 0.87532940310411085020444785628534750441072000183455),
            new SinCosPIData(0.1635742187500000, 0.49156291610654994490498767732090361503792028152004, 0.87084206347007890771429282610744800548138820653997),
            new SinCosPIData(0.2637939453125000, 0.73707579299426567506159401172749011067353105162547, 0.67581008825103699779303833296757231585825530510606),
            new SinCosPIData(0.2656250000000000, 0.74095112535495909117561689749516272972895530930909, 0.67155895484701840062537685042742180322875063219979),
            new SinCosPIData(0.2742919921875000, 0.75895953657894241302684322136696409362826467831834, 0.65113779020717033893759995053911962446300117935748),
            new SinCosPIData(0.2868652343750000, 0.78408078850986997193360542171972807117702617870429, 0.62065877669597209764836790572441826374432029151521),
            new SinCosPIData(0.2976074218750000, 0.80457609092630709007006772171622880399477251008901, 0.59384957178543357589572414256595991296525271690820),
            new SinCosPIData(0.3240966796875000, 0.85115395288271536524074373679930597882931516070243, 0.52491613472261295312660286600144415043904324458874),
            new SinCosPIData(0.3762207031250000, 0.92534030782320628581830530355747138051115409171576, 0.37913759338484733784674207434015740440644748280682),
            new SinCosPIData(0.4284667968750000, 0.97485471461870842678472080458577295838211794743741, 0.22284139064742113283539444879401250356219019084154),
            new SinCosPIData(0.4571533203125000, 0.99095417361651848118846134450598112841856173903742, 0.13420069221879202843968745358457088449280628449497),
            new SinCosPIData(0.4969482421875000, 0.99995404142512978778474382612450225933055046547842, 0.0095872330497292246437326384180423519280590486435426),
            new SinCosPIData(0.5366210937500000, 0.99338921114808065075576106525222173191899303047883, -0.11479492660651008639943528225564492447191826382303),
            new SinCosPIData(0.5369873046875000, 0.99325648383484643331057791884051967830564744323691, -0.11593773035572779364479461900463794407881530367420),
            new SinCosPIData(0.5568847656250000, 0.98407404237077648978715833561409116387237524106924, -0.17775904796110716652288466374362040759333995945327),
            new SinCosPIData(0.5842285156250000, 0.96519413117572471231073216952857054913264055912966, -0.26153448939659548653995796488565179814863501406981),
            new SinCosPIData(0.6113281250000000, 0.93945922360218991196266924587042221315839835211967, -0.34266071731199439759278198256128687961155916414389),
            new SinCosPIData(0.6309814453125000, 0.91652573855622815572837358160321765496246259024688, -0.39997571246759535441655989159896709262824979069807),
            new SinCosPIData(0.6590576171875000, 0.87772910922613151685417066189337743931280915239731, -0.47915718800525333383747291668887071598988823160230),
            new SinCosPIData(0.6705322265625000, 0.85988968707660226483289425446642370476603573873443, -0.51047989780137575553796185621778663941264829298073),
            new SinCosPIData(0.6748046875000000, 0.85296060493036365774658808174729556122474293342154, -0.52197529293715434269425831751911061907476660273241),
            new SinCosPIData(0.7031250000000000, 0.80320753148064490980667651296314192387956942717046, -0.59569930449243334346703652882996988951192633843750),
            new SinCosPIData(0.7142333984375000, 0.78193394562693757400738402698116790917498323099189, -0.62336129545897332170073718289646565461168215758923),
            new SinCosPIData(0.7551269531250000, 0.69562632734525488761267352264166204833665445713730, -0.71840379502348976120222439198520745073249928327516),
            new SinCosPIData(0.7583007812500000, 0.68842875278409046079980041353000337555287434538723, -0.72530397237306074237175471494088524758917294674515),
            new SinCosPIData(0.8010253906250000, 0.58517607232673038894627228866703931303544264019835, -0.81090626115245972040705430731849047673018063657583),
            new SinCosPIData(0.8435058593750000, 0.47207302324236866106683942458517966989383616854685, -0.88155944820914378191596769443252628163897896869396),
            new SinCosPIData(0.8828125000000000, 0.35989503653498814877510457232675642020231742112903, -0.93299279883473888771166025554330249829501552051230),
            new SinCosPIData(0.9104003906250000, 0.27778296655185765818331319471173144974417129932018, -0.96064385882263856248939297269826191626003401257321),
            new SinCosPIData(0.9166259765625000, 0.25894251895918054320988156122641853183181946392700, -0.96589273311019087764328463382295642573413327043782),
            new SinCosPIData(0.9182128906250000, 0.25412392304732064745556116586210281709375288669218, -0.96717166611467659033137747683134331219088781639110),
            new SinCosPIData(0.9294433593750000, 0.21984952979877870202795626019310012039214447278137, -0.97553379451829136288289957030202633765542967450611),
            new SinCosPIData(0.9447021484375000, 0.17285081953139408112693747894092544936316451369136, -0.98494801598222707138699250964241055723679153855555),
            new SinCosPIData(0.9781494140625000, 0.068591740687380945866892561978994342531056422788845, -0.99764481310207546695147108225749935382432353595005),
            new SinCosPIData(0.9818115234375000, 0.057109694655158060908995764015099383687442571764999, -0.99836790952854379633866874957690340441344896072513),
            new SinCosPIData(0.9842529296875000, 0.049450703970084666089629188570118582384474331154055, -0.99877656554249562180572858764020193096129521237985),
            new SinCosPIData(1.035156250000000, -0.11022220729388305880789914021567772527447462311499, -0.99390697000235604154692281324779982143559535062268),
            new SinCosPIData(1.051025390625000, -0.15961534723719304545449054877715037851208631001500, -0.98717928509787435188216571376348290095349500376378),
            new SinCosPIData(1.054321289062500, -0.16982822813571985754761592298687481168400081466799, -0.98547367947007183889340818756445635404278245284527),
            new SinCosPIData(1.062988281250000, -0.19659459767008022944893132716464786759977010847402, -0.98048486177346934356204062185012936934077472520764),
            new SinCosPIData(1.081542968750000, -0.25338203699557015475950976522966772394872618830269, -0.96736629222232853104820130629614926272155916329915),
            new SinCosPIData(1.103393554687500, -0.31913861280769589636064305483001488748783782526471, -0.94770804882895215431113586627119286098443423466622),
            new SinCosPIData(1.111694335937500, -0.34374132459779849843926720325234820256759634670687, -0.93906437572924196673598219778240204586222038083134),
            new SinCosPIData(1.115112304687500, -0.35380486100177206446135725272490396429360206025605, -0.93531926117851157452389050726776022395670332454165),
            new SinCosPIData(1.140991210937500, -0.42859483689734440122352140997581841231665665120132, -0.90349679898986844735840821785342106632913960209682),
            new SinCosPIData(1.156005859375000, -0.47072017309907163346193703268586474774361325384110, -0.88228256167600866739108420118332151462109503931178),
            new SinCosPIData(1.241455078125000, -0.68787224916668555581204556764403511943527170856987, -0.72583177722277030564486824926855832976162285896118),
            new SinCosPIData(1.265991210937500, -0.74172325371778417205997052210430754266589846314852, -0.67070605699837210310172965890871658176771354561618),
            new SinCosPIData(1.273071289062500, -0.75645687960083371274415676459120613462083581058214, -0.65404356835349261643773791877980168895302802677008),
            new SinCosPIData(1.288574218750000, -0.78740174702903138992122299741270874876834240255437, -0.61644017453085360756657406356363545047561518849797),
            new SinCosPIData(1.296020507812500, -0.80160550454704610146404887971884084329708867121448, -0.59785333910573390051678947379699891565668999176937),
            new SinCosPIData(1.309448242187500, -0.82610501784466462721075118319831054144946505749197, -0.56351619275036483892287316784816757173786460136956),
            new SinCosPIData(1.312866210937500, -0.83210823743573554216679653590005270307168220094656, -0.55461327174130404434489127134456473149268263209358),
            new SinCosPIData(1.313598632812500, -0.83338218268057974514922802401670466895460221456365, -0.55269714816574981985830970414955916602028124045107),
            new SinCosPIData(1.333618164062500, -0.86647246807174300007182765564459178242288106937812, -0.49922486123355506154608714794439097080837050194777),
            new SinCosPIData(1.425903320312500, -0.97302849033369417389370465529395780374508482889839, -0.23068497350051221660406548641159755104138622432809),
            new SinCosPIData(1.482177734375000, -0.99843295266650845768100912338718169352785980883151, -0.055961049218520569880500359937782035877677393480671),
            new SinCosPIData(1.495117187500000, -0.99988234745421252563304962650591566082925735483662, -0.015339206284988101044151867602462621301745093738657),
            new SinCosPIData(1.508056640625000, -0.99967970176298791385635589130878355635602970346506, 0.025307980620024570359778747908591006267128165919173),
            new SinCosPIData(1.573364257812500, -0.97355671350826554112325660905190331393298655358010, 0.22844545428391647392293810148144747152014356660464),
            new SinCosPIData(1.580810546875000, -0.96794675562898779590194637612057530996528319654642, 0.25115548623774194323605285518811273660146124969746),
            new SinCosPIData(1.582519531250000, -0.96658437447833308519915220431895826946047409625585, 0.25634868248994288975851068915232357159853191738483),
            new SinCosPIData(1.599731445312500, -0.95131689215046557627384843093214777453948525791318, 0.30821448815586110800445981141111355919766072767394),
            new SinCosPIData(1.622558593750000, -0.92678747430458179887576501301638647676905046680464, 0.37558617848921721455761858775327472722706168060315),
            new SinCosPIData(1.623046875000000, -0.92621024213831134197479338843714328234114125061739, 0.37700741021641825672656782319985723230153788364791),
            new SinCosPIData(1.632812500000000, -0.91420975570353065463501482939357740104469111568218, 0.40524131400498987090848130550505246651194775410509),
            new SinCosPIData(1.637817382812500, -0.90772528206767641183018248145385603435819612542142, 0.41956502749294690942605922870540591823044170481709),
            new SinCosPIData(1.648925781250000, -0.89253355540276462147223438162310406169311099093980, 0.45098098904510386908356011414644060392405931383576),
            new SinCosPIData(1.650146484375000, -0.89079750603628151090608689434741957669418194562714, 0.45440048771930362568380998966288866522955338595975),
            new SinCosPIData(1.699462890625000, -0.81000765858164110554058550006063253477607264085195, 0.58641929797636054285359980658890710085639082695044),
            new SinCosPIData(1.769897460937500, -0.66155346876039896600714888553777386395519275201354, 0.74989799837783527008945755655423606303782865304762),
            new SinCosPIData(1.794921875000000, -0.60061647938386892665387589554555953864395455132188, 0.79953726910790503350024623225150776431951023804419),
            new SinCosPIData(1.811035156250000, -0.55939071185913607502653784218704460726181357998191, 0.82890411477186491196279734440419539836682400359610),
            new SinCosPIData(1.829223632812500, -0.51113927471546437880837487317012202363430056395056, 0.85949790101160167823029313490699874713784616354887),
            new SinCosPIData(1.856567382812500, -0.43551189610849202516852206444967423236676505377988, 0.90018297492675679287248205592417732434207584954562),
            new SinCosPIData(1.879150390625000, -0.37060492955905164108271196945182017626890274950545, 0.92879060405805698007605315123046367038234722611183),
            new SinCosPIData(1.884521484375000, -0.35488069794622278666268672277608569021060008287173, 0.93491159487151606617515594113812793861212211904888),
            new SinCosPIData(1.899169921875000, -0.31149607495827589991891854468762662594752412050292, 0.95024743897870525216650109980318294880154610316841),
            new SinCosPIData(1.901367187500000, -0.30492922973540240649072863343652234631925192710263, 0.95237501271976585852989360757100877759109622462244),
            new SinCosPIData(1.919555664062500, -0.25004171147145467401674094848282059189723603946352, 0.96823506573787432318054173897007713798749605591361),
            new SinCosPIData(1.954833984375000, -0.14141756302230303205802218810078196203174007217105, 0.98995003554160901404677939416508732611453798078192),
            new SinCosPIData(1.978149414062500, -0.068591740687380945866892561978994342531056422788845, 0.99764481310207546695147108225749935382432353595005),
            new SinCosPIData(2.010620117187500, 0.033357892543086143778344709122365509295211451368925, 0.99944347064007773117655616329980565767985014719348),
            new SinCosPIData(2.010986328125000, 0.034507715524795753429033897215001134467869302831307, 0.99940443143367128564964301438291766172274306859184),
            new SinCosPIData(2.047119140625000, 0.14748912010315358814312930291520472187908018453416, 0.98906367815788156972653867758731293808134885882172),
            new SinCosPIData(2.047485351562500, 0.14862692575279654766675642836213297701029823167773, 0.98889333951709509042676348734848591543832041080142),
            new SinCosPIData(2.062500000000000, 0.19509032201612826784828486847702224092769161775195, 0.98078528040323044912618223613423903697393373089334),
            new SinCosPIData(2.121704101562500, 0.37309689835856063549092330355735582692659092147755, 0.92779238218214632750946245785451880519834775445407),
            new SinCosPIData(2.193725585937500, 0.57172503454319716795673895919395819390038111375896, 0.82044529669965200250101947575278135359122153344657),
            new SinCosPIData(2.206909179687500, 0.60520579725549652839725500308853227027756452862529, 0.79606905665798795117529708790252252627086146722079),
            new SinCosPIData(2.220947265625000, 0.63971415168764048932067723687628473771512322781137, 0.76861290916205830720889209775471723207877077126050),
            new SinCosPIData(2.235473632812500, 0.67411231056231237722332658309055402720998736759431, 0.73862885994817485778845874464818783496755997477179),
            new SinCosPIData(2.255981445312500, 0.72026858973207713404327722694288187107337794268076, 0.69369529236211827539684608243141343601656519079705),
            new SinCosPIData(2.349853515625000, 0.89079750603628151090608689434741957669418194562714, 0.45440048771930362568380998966288866522955338595975),
            new SinCosPIData(2.356079101562500, 0.89951384919848793705503674606819436585935567876861, 0.43689224655528038161997338360536871933296271804575),
            new SinCosPIData(2.401245117187500, 0.95225800379539958235175825874836019293751308395084, 0.30529443854679166856246436031760123709319375625518),
            new SinCosPIData(2.423583984375000, 0.97132181041978620305447162927297185684476506236132, 0.23776867035593421658391022272144598990293795306227),
            new SinCosPIData(2.426147460937500, 0.97320513727125281871147962787395341347817920703212, 0.22993860221555222312884149397788283626676260381534),
            new SinCosPIData(2.503784179687500, 0.99992933438627605999734223201603863222253311984127, -0.011888071072252092833123258366366978361167906410710),
            new SinCosPIData(2.593750000000000, 0.95694033573220886493579788698026996948284920563004, -0.29028467725446236763619237581739527469147627832415),
            new SinCosPIData(2.602294921875000, 0.94880389496265843806156874468377892831958502073538, -0.31586574506218399658986986510305094045688995329416),
            new SinCosPIData(2.660278320312500, 0.87588511461810375211887565930109212159655034863819, -0.48251970528718436251068191790050966063814311921669),
            new SinCosPIData(2.666625976562500, 0.86608931257458679816530404193518094994802874075851, -0.49988929038746137780345485315745866792596758696835),
            new SinCosPIData(2.668823242187500, 0.86261801283581679391398848969198965608661894228157, -0.50585587268626882372717469280804167875288475781575),
            new SinCosPIData(2.689941406250000, 0.82718402727366911231560708546251750903353089275968, -0.56193112124468941370819180877712253228015407720620),
            new SinCosPIData(2.692504882812500, 0.82263179629451498532588521681730084406675913376412, -0.56857446981486919490413833215420356498902713171518),
            new SinCosPIData(2.723510742187500, 0.76343616653417203738186961986172484252886128266241, -0.64588328638199637137095119515176887214736824646866),
            new SinCosPIData(2.731323242187500, 0.74735546450394024818094695143728970952148723736532, -0.66442441983727517699038765653602383737007157105066),
            new SinCosPIData(2.736450195312500, 0.73655723639791913955873498284333824299842011850576, -0.67637521946761167594849411658463144143086749737566),
            new SinCosPIData(2.739868164062500, 0.72925208705878692814838310933493026784991366094958, -0.68424512677870308009260905850944005393411397853461),
            new SinCosPIData(2.805297851562500, 0.57423943459296786208216093039118311383218667943296, -0.81868740784156964132688152182922399752571520714271),
            new SinCosPIData(2.904907226562500, 0.29431888468262740240431678877743022944059321090326, -0.95570727428390659550766483252857828793717958272412),
            new SinCosPIData(2.916259765625000, 0.26005359301549519442919360781448551744452327050920, -0.96559418430297683158832809149393120827990021830691),
            new SinCosPIData(2.926025390625000, 0.23031180479384545578916264636432756612941585428125, -0.97311688535992510817457181588276250330783779857148),
            new SinCosPIData(2.963500976562500, 0.11441395818328693017185508060764793696178542475818, -0.99343316140182932273974895864445182290578062310842),
            new SinCosPIData(3.002319335937500, -0.0072863442679265220477288039971931799125241142852542, -0.99997345424126597377517955377266345603436152191913),
            new SinCosPIData(3.037109375000000, -0.11631863091190476725254431947051254092385667649227, -0.99321194923479453310460101209278277844699436597554),
            new SinCosPIData(3.068847656250000, -0.21460881099378677183161099535945108175457904340316, -0.97670008612871182956093955174975102655802017149826),
            new SinCosPIData(3.093994140625000, -0.29101855584408506162685298669035799115578046126741, -0.95671740872340310123255484397822798575297217192404),
            new SinCosPIData(3.102172851562500, -0.31550186010755600261264604319888031105927676258547, -0.94892495818619512281647339190312069498723353956200),
            new SinCosPIData(3.111206054687500, -0.34230041402351352995307351288254720555992135402518, -0.93959056325578920850599372540567050067458462731211),
            new SinCosPIData(3.111572265625000, -0.34338117265211504809196022868601901481003169654692, -0.93919612981956987863487235682551164720756088993713),
            new SinCosPIData(3.133666992187500, -0.40769401625328014878109863032848008319187684255554, -0.91311860626715418051901496303712555137402984873811),
            new SinCosPIData(3.171020507812500, -0.51179835093948690328766325679242479443244492312053, -0.85910560932613040715959989484418082499421277319466),
            new SinCosPIData(3.173217773437500, -0.51771644198787112902858124655711666416393411978027, -0.85555226941164691638495304380331657707484677200700),
            new SinCosPIData(3.205566406250000, -0.60184224705858005965004629561674329166957690635464, -0.79861499463476083767048034109114977151346776535820),
            new SinCosPIData(3.231567382812500, -0.66499743881132532197757026500149351906192421352850, -0.74684563758140653301058182784072329699138194531247),
            new SinCosPIData(3.233642578125000, -0.66985227139182102967718549444214745424435237203003, -0.74249440032313923555006954562227411661501874185777),
            new SinCosPIData(3.238891601562500, -0.68200459271844082459138646403815797741149977467372, -0.73134788952382548909378574434824929888351227907222),
            new SinCosPIData(3.305664062500000, -0.81934752007679696082468963724253081583636714099752, -0.57329716669804221282017123894213993360437766249721),
            new SinCosPIData(3.320800781250000, -0.84567324698729905899100660647301780225467275441717, -0.53370100180715295622464514542509334924729484137129),
            new SinCosPIData(3.329467773437500, -0.85988968707660226483289425446642370476603573873443, -0.51047989780137575553796185621778663941264829298073),
            new SinCosPIData(3.370605468750000, -0.91850839432521222554762547029604595682966352016538, -0.39540147894781633694607828110122394363059236212076),
            new SinCosPIData(3.391113281250000, -0.94205973977101735561225214326446936244909534195353, -0.33544514708453163254928496829827344827633460748457),
            new SinCosPIData(3.435546875000000, -0.97956976568544053443932610987989550521323449370167, -0.20110463484209191155844354588206717746552731686488),
            new SinCosPIData(3.442749023437500, -0.98386888192401718922821924206223662941367055341158, -0.17889109307504474568528093663928536750121199985037),
            new SinCosPIData(3.465820312500000, -0.99424044945318794635841344190689881768699897755293, -0.10717242495680884917552914822819675037791937965835),
            new SinCosPIData(3.525634765625000, -0.99675889043081803306413512919055553929101702377102, 0.080446966052950007795265135898733932835758854837756),
            new SinCosPIData(3.534057617187500, -0.99428147645164150519044446362253384600034565322093, 0.10679113064830739485799947892892021505380094438026),
            new SinCosPIData(3.570556640625000, -0.97553379451829136288289957030202633765542967450611, 0.21984952979877870202795626019310012039214447278137),
            new SinCosPIData(3.596557617187500, -0.95434272061471647205852516995004484303744837128140, 0.29871386243310038201085207476042461883731574169521),
            new SinCosPIData(3.624511718750000, -0.92446547432526263363323350297086805977261864765827, 0.38126576922216236948926956486303165475359772232545),
            new SinCosPIData(3.699340820312500, -0.81023248799698232717273262633389614402805378630579, 0.58610862081547639451250761373976604486204442866549),
            new SinCosPIData(3.721923828125000, -0.76664667656531043873299434297041280848103535895455, 0.64206921224379251975056996350966250008260271013597),
            new SinCosPIData(3.756835937500000, -0.69175925836415777490673413208882878382690677052991, 0.72212819392921532124360719766762509982791708677698),
            new SinCosPIData(3.758911132812500, -0.68703673511009563414948166590563200996836293928404, 0.72662268379762290830112158508142224149348337758028),
            new SinCosPIData(3.764282226562500, -0.67467863346558449893811811613382320995171524081368, 0.73811160507406431389100472643631584576020424818282),
            new SinCosPIData(3.807373046875000, -0.56888990334017586801295907775117271793210518920120, 0.82241369022992641192399657096644286231578154030154),
            new SinCosPIData(3.836914062500000, -0.49022648328829115422959844903660872001836233879988, 0.87159508665595103484248143520115068911901646304507),
            new SinCosPIData(3.849365234375000, -0.45576641881943468797102590776166931401930223294796, 0.89009941662519229531471631394569700081986282785949),
            new SinCosPIData(3.877441406250000, -0.37558617848921721455761858775327472722706168060315, 0.92678747430458179887576501301638647676905046680464),
            new SinCosPIData(3.901855468750000, -0.30346794657201130093558009858171686916706922062842, 0.95284164760119868225278490287765409482286613061945),
            new SinCosPIData(3.915649414062500, -0.26190461746922259634146953305223525622142181723721, 0.96509376298279960545886252037035002516044802234760),
            new SinCosPIData(3.924194335937500, -0.23590574814860738121822079333460638042299547003522, 0.97177594021999012708645016220016850774717672506428),
            new SinCosPIData(3.964843750000000, -0.11022220729388305880789914021567772527447462311499, 0.99390697000235604154692281324779982143559535062268),
            new SinCosPIData(3.992309570312500, -0.024157847032299862119380272108097639772880617137513, 0.99970815662710484896318784190108918525305082518986),
            new SinCosPIData(3.998535156250000, -0.0046019261204485707649016992969119673728849251615110, 0.99998941108192837361947235727373283744944631695329),
            new SinCosPIData(4.004516601562500, 0.014188846153786344416273340780549230049368515476147, 0.99989933325551536866332259669356346524594844640949),
            new SinCosPIData(4.069580078125000, 0.21685559964263262680404980471938273544336906226188, 0.97620369232227053287875378090070846953560860476203),
            new SinCosPIData(4.092407226562500, 0.28624530413205715803454626126745498923950115868406, 0.95815636816875882061670537963734817642030876462306),
            new SinCosPIData(4.103637695312500, 0.31986540183563049712809835863203946461291872312236, 0.94746299384647769481893638707053361693832231368042),
            new SinCosPIData(4.104003906250000, 0.32095523242787523926293387394413286987697359283178, 0.94709436635277721540320624972087051206446563018208),
            new SinCosPIData(4.127075195312500, 0.38869841434151920735254523409246641815255184625171, 0.92136504312264238158826074245669411814017304091673),
            new SinCosPIData(4.175415039062500, 0.52360986383422796288920866457916802151273236267829, 0.85195816240910636138627686797698634031118417901889),
            new SinCosPIData(4.196166992187500, 0.57800089298526991148080866263607425164228767672543, 0.81603613137423673737790243147786943380012398142004),
            new SinCosPIData(4.200927734375000, 0.59014068383224889262512658480231576061562551267644, 0.80730042319201446735208789406168340805437402559667),
            new SinCosPIData(4.219360351562500, 0.63587434599469772273199422890666140738335733619646, 0.77179259915201016922233564145979510514466028515153),
            new SinCosPIData(4.235595703125000, 0.67439552160513900371743295491032044856967701656684, 0.73837028680664858621101565815516320046197403766031),
            new SinCosPIData(4.236328125000000, 0.67609270357531596036041922765815255809277314527663, 0.73681656887736987509013252017274694686788445838695),
            new SinCosPIData(4.288452148437500, 0.78716528728765090584871074140433951084893700422872, 0.61674209398203882869758652635018789479637411136229),
            new SinCosPIData(4.313842773437500, 0.83380585091378635382322687224052706968133332330954, 0.55205778953107498402400133241091812740428913669024),
            new SinCosPIData(4.354248046875000, 0.89698578927886398735680590657061438336749537916400, 0.44205937817421474932042156158542146885105827814824),
            new SinCosPIData(4.394409179687500, 0.94548250091445378665100955240537632663345751077647, 0.32567290409941981307625717395156635642470497103022),
            new SinCosPIData(4.416259765625000, 0.96559418430297683158832809149393120827990021830691, 0.26005359301549519442919360781448551744452327050920),
            new SinCosPIData(4.472045898437500, 0.99614626641982466909203063824094414464767378920893, 0.087707558954993673009972928730162808938572840914799),
            new SinCosPIData(4.503173828125000, 0.99995029123649047314916064301127315065152183778988, -0.0099707099074180297611249412847618321489896738216746),
            new SinCosPIData(4.559082031250000, 0.98282355119870527722009930350834495627033988058992, -0.18454773693861961481801819392236258369015286899386),
            new SinCosPIData(4.633422851562500, 0.91343103504855473774216216247846846606267103577356, -0.40699354320446653284605350602283675990108414585621),
            new SinCosPIData(4.654296875000000, 0.88479709843093778010400704058574088975908638194697, -0.46597649576796617790275606488777870342689503673289),
            new SinCosPIData(4.660766601562500, 0.87514390842956032049101534498775533255718164645222, -0.48386272799073224402491355353758267161805180346871),
            new SinCosPIData(4.672607421875000, 0.85654340483771995539032185079708194505428983752125, -0.51607499031536664741422588659539879681603856888577),
            new SinCosPIData(4.680541992187500, 0.84341433969379278981444914384683343631448549649886, -0.53726367046254253880228571009317802068849583260467),
            new SinCosPIData(4.716430664062500, 0.77761234108341997299494785281177046176067347485597, -0.62874402342667474298795315133032904539861092651762),
            new SinCosPIData(4.723632812500000, 0.76318841726338127170483829706586545830731955494282, -0.64617601298331636483280221953658528888217376606774),
            new SinCosPIData(4.800537109375000, 0.58641929797636054285359980658890710085639082695044, -0.81000765858164110554058550006063253477607264085195),
            new SinCosPIData(4.816894531250000, 0.54403852673088388510139720972070268461245984847492, -0.83906023707031269401306993083369298151812048224499),
            new SinCosPIData(4.828002929687500, 0.51443164118322288834043633101998057941491085564124, -0.85753139099949910622628639459821234043621899812435),
            new SinCosPIData(4.892944335937500, 0.33002049619310542025961092291992485856874396752882, -0.94397376663361598462711440944497001173439900956397),
            new SinCosPIData(4.924438476562500, 0.23516033602183473124666530333004541872733353307880, -0.97195659180958167856656353854565326670154546978207),
            new SinCosPIData(4.935180664062500, 0.20223148240144146904027296829768942769259489938194, -0.97933774946425682514090820934857594201470540587103),
            new SinCosPIData(4.949951171875000, 0.15658597269299844050413693763877501398377764347613, -0.98766433222820573147215080395327713233553593859022),
            new SinCosPIData(4.956665039062500, 0.13572063839303992889495640391817444946740282040363, -0.99074714650822270033015722909760980781302793175920),
                                                                  };
        #endregion

        static IEnumerable<TestCase<double>> GetSinPIData(SinCosPIData[] data)
        {
            return data.Select(d => TestCase.Create(d.X, d.SinPI));
        }

        static IEnumerable<TestCase<double>> GetNegativeSinPIData(SinCosPIData[] data)
        {
            return data.Select(d => TestCase.Create(-d.X, -d.SinPI));
        }

        static IEnumerable<TestCase<double>> GetCosPIData(SinCosPIData[] data)
        {
            return data.Select(d => TestCase.Create(d.X, d.CosPI));
        }

        static IEnumerable<TestCase<double>> GetNegativeCosPIData(SinCosPIData[] data)
        {
            return data.Select(d => TestCase.Create(-d.X, d.CosPI));
        }

        /// <summary>
        ///A test for SinPI
        ///</summary>
        [TestMethod]
        public void SinPITest()
        {

            TestCaseSet<double>[] testCases = {
                new TestCaseSet<double>(GetSinPIData(_SinCosPIData),"SinPI: Spot Data",    2,1),
                new TestCaseSet<double>(GetNegativeSinPIData(_SinCosPIData),"SinPI: Spot Data",    2,1),
            };

            NumericFunctionTest.RunSet(Math2.SinPI, "SinPI", testCases);

        }

        /// <summary>
        ///A test for CosPI
        ///</summary>
        [TestMethod]
        public void CosPITest()
        {

            TestCaseSet<double>[] testCases = {
                new TestCaseSet<double>(GetCosPIData(_SinCosPIData),"CosPI: Spot Data", 2,1),
                new TestCaseSet<double>(GetNegativeCosPIData(_SinCosPIData),"CosPI: Spot Data", 2,1),
            };

            NumericFunctionTest.RunSet(Math2.CosPI, "CosPI", testCases);

        }


    }
}
