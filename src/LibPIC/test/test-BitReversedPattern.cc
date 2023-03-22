/*
 * Copyright (c) 2021-2022, Kyungguk Min
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include "catch2/catch.hpp"

#define LIBPIC_INLINE_VERSION 1
#include "../PIC/Random/BitReversedPattern.h"
#include <algorithm>
#include <vector>

TEST_CASE("Test LibPIC::BitReversedPattern::Base2", "[LibPIC::BitReversedPattern::Base2]")
{
    using Gen = BitReversedPattern<2>;
    REQUIRE(Gen::min() == 0);
    REQUIRE(Gen::max() == 9223372036854775808UL);

    const std::vector seq1{
        4611686018427387904UL, 2305843009213693952UL, 6917529027641081856UL, 1152921504606846976UL,
        5764607523034234880UL, 3458764513820540928UL, 8070450532247928832UL, 576460752303423488UL,
        5188146770730811392UL, 2882303761517117440UL, 7493989779944505344UL, 1729382256910270464UL,
        6341068275337658368UL, 4035225266123964416UL, 8646911284551352320UL, 288230376151711744UL,
        4899916394579099648UL, 2594073385365405696UL, 7205759403792793600UL, 1441151880758558720UL,
        6052837899185946624UL, 3746994889972252672UL, 8358680908399640576UL, 864691128455135232UL,
        5476377146882523136UL, 3170534137668829184UL, 7782220156096217088UL, 2017612633061982208UL,
        6629298651489370112UL, 4323455642275676160UL, 8935141660703064064UL, 144115188075855872UL,
        4755801206503243776UL, 2449958197289549824UL, 7061644215716937728UL, 1297036692682702848UL,
        5908722711110090752UL, 3602879701896396800UL, 8214565720323784704UL, 720575940379279360UL,
        5332261958806667264UL, 3026418949592973312UL, 7638104968020361216UL, 1873497444986126336UL,
        6485183463413514240UL, 4179340454199820288UL, 8791026472627208192UL, 432345564227567616UL,
        5044031582654955520UL, 2738188573441261568UL, 7349874591868649472UL, 1585267068834414592UL,
        6196953087261802496UL, 3891110078048108544UL, 8502796096475496448UL, 1008806316530991104UL,
        5620492334958379008UL, 3314649325744685056UL, 7926335344172072960UL, 2161727821137838080UL,
        6773413839565225984UL, 4467570830351532032UL, 9079256848778919936UL, 72057594037927936UL,
        4683743612465315840UL, 2377900603251621888UL, 6989586621679009792UL, 1224979098644774912UL,
        5836665117072162816UL, 3530822107858468864UL, 8142508126285856768UL, 648518346341351424UL,
        5260204364768739328UL, 2954361355555045376UL, 7566047373982433280UL, 1801439850948198400UL,
        6413125869375586304UL, 4107282860161892352UL, 8718968878589280256UL, 360287970189639680UL,
        4971973988617027584UL, 2666130979403333632UL, 7277816997830721536UL, 1513209474796486656UL,
        6124895493223874560UL, 3819052484010180608UL, 8430738502437568512UL, 936748722493063168UL,
        5548434740920451072UL, 3242591731706757120UL, 7854277750134145024UL, 2089670227099910144UL,
        6701356245527298048UL, 4395513236313604096UL, 9007199254740992000UL, 216172782113783808UL,
        4827858800541171712UL, 2522015791327477760UL, 7133701809754865664UL, 1369094286720630784UL
    };
    std::vector<Gen::result_type> seq2(seq1.size());
    std::generate(begin(seq2), end(seq2), [gen = Gen{}]() mutable {
        return gen();
    });
    CHECK(std::equal(begin(seq1), end(seq1), begin(seq2)));
}
TEST_CASE("Test LibPIC::BitReversedPattern::Base3", "[LibPIC::BitReversedPattern::Base3]")
{

    using Gen = BitReversedPattern<3>;
    REQUIRE(Gen::min() == 0);
    REQUIRE(Gen::max() == 12157665459056928801UL);

    const std::vector seq1{
        4052555153018976267UL, 8105110306037952534UL, 1350851717672992089UL,
        5403406870691968356UL, 9455962023710944623UL, 2701703435345984178UL,
        6754258588364960445UL, 10806813741383936712UL, 450283905890997363UL,
        4502839058909973630UL, 8555394211928949897UL, 1801135623563989452UL,
        5853690776582965719UL, 9906245929601941986UL, 3151987341236981541UL,
        7204542494255957808UL, 11257097647274934075UL, 900567811781994726UL,
        4953122964800970993UL, 9005678117819947260UL, 2251419529454986815UL,
        6303974682473963082UL, 10356529835492939349UL, 3602271247127978904UL,
        7654826400146955171UL, 11707381553165931438UL, 150094635296999121UL,
        4202649788315975388UL, 8255204941334951655UL, 1500946352969991210UL,
        5553501505988967477UL, 9606056659007943744UL, 2851798070642983299UL,
        6904353223661959566UL, 10956908376680935833UL, 600378541187996484UL,
        4652933694206972751UL, 8705488847225949018UL, 1951230258860988573UL,
        6003785411879964840UL, 10056340564898941107UL, 3302081976533980662UL,
        7354637129552956929UL, 11407192282571933196UL, 1050662447078993847UL,
        5103217600097970114UL, 9155772753116946381UL, 2401514164751985936UL,
        6454069317770962203UL, 10506624470789938470UL, 3752365882424978025UL,
        7804921035443954292UL, 11857476188462930559UL, 300189270593998242UL,
        4352744423612974509UL, 8405299576631950776UL, 1651040988266990331UL,
        5703596141285966598UL, 9756151294304942865UL, 3001892705939982420UL,
        7054447858958958687UL, 11107003011977934954UL, 750473176484995605UL,
        4803028329503971872UL, 8855583482522948139UL, 2101324894157987694UL,
        6153880047176963961UL, 10206435200195940228UL, 3452176611830979783UL,
        7504731764849956050UL, 11557286917868932317UL, 1200757082375992968UL,
        5253312235394969235UL, 9305867388413945502UL, 2551608800048985057UL,
        6604163953067961324UL, 10656719106086937591UL, 3902460517721977146UL,
        7955015670740953413UL, 12007570823759929680UL, 50031545098999707UL,
        4102586698117975974UL, 8155141851136952241UL, 1400883262771991796UL,
        5453438415790968063UL, 9505993568809944330UL, 2751734980444983885UL,
        6804290133463960152UL, 10856845286482936419UL, 500315450989997070UL,
        4552870604008973337UL, 8605425757027949604UL, 1851167168662989159UL,
        5903722321681965426UL, 9956277474700941693UL, 3202018886335981248UL,
        7254574039354957515UL, 11307129192373933782UL, 950599356880994433UL,
        5003154509899970700UL
    };
    std::vector<Gen::result_type> seq2(seq1.size());
    std::generate(begin(seq2), end(seq2), [gen = Gen{}]() mutable {
        return gen();
    });
    CHECK(std::equal(begin(seq1), end(seq1), begin(seq2)));
}
TEST_CASE("Test LibPIC::BitReversedPattern::Base5", "[LibPIC::BitReversedPattern::Base5]")
{

    using Gen = BitReversedPattern<5>;
    REQUIRE(Gen::min() == 0);
    REQUIRE(Gen::max() == 7450580596923828125UL);

    const std::vector seq1{
        1490116119384765625UL, 2980232238769531250UL, 4470348358154296875UL, 5960464477539062500UL,
        298023223876953125UL, 1788139343261718750UL, 3278255462646484375UL, 4768371582031250000UL,
        6258487701416015625UL, 596046447753906250UL, 2086162567138671875UL, 3576278686523437500UL,
        5066394805908203125UL, 6556510925292968750UL, 894069671630859375UL, 2384185791015625000UL,
        3874301910400390625UL, 5364418029785156250UL, 6854534149169921875UL, 1192092895507812500UL,
        2682209014892578125UL, 4172325134277343750UL, 5662441253662109375UL, 7152557373046875000UL,
        59604644775390625UL, 1549720764160156250UL, 3039836883544921875UL, 4529953002929687500UL,
        6020069122314453125UL, 357627868652343750UL, 1847743988037109375UL, 3337860107421875000UL,
        4827976226806640625UL, 6318092346191406250UL, 655651092529296875UL, 2145767211914062500UL,
        3635883331298828125UL, 5125999450683593750UL, 6616115570068359375UL, 953674316406250000UL,
        2443790435791015625UL, 3933906555175781250UL, 5424022674560546875UL, 6914138793945312500UL,
        1251697540283203125UL, 2741813659667968750UL, 4231929779052734375UL, 5722045898437500000UL,
        7212162017822265625UL, 119209289550781250UL, 1609325408935546875UL, 3099441528320312500UL,
        4589557647705078125UL, 6079673767089843750UL, 417232513427734375UL, 1907348632812500000UL,
        3397464752197265625UL, 4887580871582031250UL, 6377696990966796875UL, 715255737304687500UL,
        2205371856689453125UL, 3695487976074218750UL, 5185604095458984375UL, 6675720214843750000UL,
        1013278961181640625UL, 2503395080566406250UL, 3993511199951171875UL, 5483627319335937500UL,
        6973743438720703125UL, 1311302185058593750UL, 2801418304443359375UL, 4291534423828125000UL,
        5781650543212890625UL, 7271766662597656250UL, 178813934326171875UL, 1668930053710937500UL,
        3159046173095703125UL, 4649162292480468750UL, 6139278411865234375UL, 476837158203125000UL,
        1966953277587890625UL, 3457069396972656250UL, 4947185516357421875UL, 6437301635742187500UL,
        774860382080078125UL, 2264976501464843750UL, 3755092620849609375UL, 5245208740234375000UL,
        6735324859619140625UL, 1072883605957031250UL, 2562999725341796875UL, 4053115844726562500UL,
        5543231964111328125UL, 7033348083496093750UL, 1370906829833984375UL, 2861022949218750000UL,
        4351139068603515625UL, 5841255187988281250UL, 7331371307373046875UL, 238418579101562500UL
    };
    std::vector<Gen::result_type> seq2(seq1.size());
    std::generate(begin(seq2), end(seq2), [gen = Gen{}]() mutable {
        return gen();
    });
    CHECK(std::equal(begin(seq1), end(seq1), begin(seq2)));
}
TEST_CASE("Test LibPIC::BitReversedPattern::Base7", "[LibPIC::BitReversedPattern::Base7]")
{

    using Gen = BitReversedPattern<7>;
    REQUIRE(Gen::min() == 0);
    REQUIRE(Gen::max() == 3909821048582988049UL);

    const std::vector seq1{
        558545864083284007UL, 1117091728166568014UL, 1675637592249852021UL, 2234183456333136028UL,
        2792729320416420035UL, 3351275184499704042UL, 79792266297612001UL, 638338130380896008UL,
        1196883994464180015UL, 1755429858547464022UL, 2313975722630748029UL, 2872521586714032036UL,
        3431067450797316043UL, 159584532595224002UL, 718130396678508009UL, 1276676260761792016UL,
        1835222124845076023UL, 2393767988928360030UL, 2952313853011644037UL, 3510859717094928044UL,
        239376798892836003UL, 797922662976120010UL, 1356468527059404017UL, 1915014391142688024UL,
        2473560255225972031UL, 3032106119309256038UL, 3590651983392540045UL, 319169065190448004UL,
        877714929273732011UL, 1436260793357016018UL, 1994806657440300025UL, 2553352521523584032UL,
        3111898385606868039UL, 3670444249690152046UL, 398961331488060005UL, 957507195571344012UL,
        1516053059654628019UL, 2074598923737912026UL, 2633144787821196033UL, 3191690651904480040UL,
        3750236515987764047UL, 478753597785672006UL, 1037299461868956013UL, 1595845325952240020UL,
        2154391190035524027UL, 2712937054118808034UL, 3271482918202092041UL, 3830028782285376048UL,
        11398895185373143UL, 569944759268657150UL, 1128490623351941157UL, 1687036487435225164UL,
        2245582351518509171UL, 2804128215601793178UL, 3362674079685077185UL, 91191161482985144UL,
        649737025566269151UL, 1208282889649553158UL, 1766828753732837165UL, 2325374617816121172UL,
        2883920481899405179UL, 3442466345982689186UL, 170983427780597145UL, 729529291863881152UL,
        1288075155947165159UL, 1846621020030449166UL, 2405166884113733173UL, 2963712748197017180UL,
        3522258612280301187UL, 250775694078209146UL, 809321558161493153UL, 1367867422244777160UL,
        1926413286328061167UL, 2484959150411345174UL, 3043505014494629181UL, 3602050878577913188UL,
        330567960375821147UL, 889113824459105154UL, 1447659688542389161UL, 2006205552625673168UL,
        2564751416708957175UL, 3123297280792241182UL, 3681843144875525189UL, 410360226673433148UL,
        968906090756717155UL, 1527451954840001162UL, 2085997818923285169UL, 2644543683006569176UL,
        3203089547089853183UL, 3761635411173137190UL, 490152492971045149UL, 1048698357054329156UL,
        1607244221137613163UL, 2165790085220897170UL, 2724335949304181177UL, 3282881813387465184UL,
        3841427677470749191UL, 22797790370746286UL, 581343654454030293UL, 1139889518537314300UL
    };
    std::vector<Gen::result_type> seq2(seq1.size());
    std::generate(begin(seq2), end(seq2), [gen = Gen{}]() mutable {
        return gen();
    });
    CHECK(std::equal(begin(seq1), end(seq1), begin(seq2)));
}
TEST_CASE("Test LibPIC::BitReversedPattern::Base11", "[LibPIC::BitReversedPattern::Base11]")
{

    using Gen = BitReversedPattern<11>;
    REQUIRE(Gen::min() == 0);
    REQUIRE(Gen::max() == 5559917313492231481UL);

    const std::vector seq1{
        505447028499293771UL, 1010894056998587542UL, 1516341085497881313UL, 2021788113997175084UL,
        2527235142496468855UL, 3032682170995762626UL, 3538129199495056397UL, 4043576227994350168UL,
        4549023256493643939UL, 5054470284992937710UL, 45949729863572161UL, 551396758362865932UL,
        1056843786862159703UL, 1562290815361453474UL, 2067737843860747245UL, 2573184872360041016UL,
        3078631900859334787UL, 3584078929358628558UL, 4089525957857922329UL, 4594972986357216100UL,
        5100420014856509871UL, 91899459727144322UL, 597346488226438093UL, 1102793516725731864UL,
        1608240545225025635UL, 2113687573724319406UL, 2619134602223613177UL, 3124581630722906948UL,
        3630028659222200719UL, 4135475687721494490UL, 4640922716220788261UL, 5146369744720082032UL,
        137849189590716483UL, 643296218090010254UL, 1148743246589304025UL, 1654190275088597796UL,
        2159637303587891567UL, 2665084332087185338UL, 3170531360586479109UL, 3675978389085772880UL,
        4181425417585066651UL, 4686872446084360422UL, 5192319474583654193UL, 183798919454288644UL,
        689245947953582415UL, 1194692976452876186UL, 1700140004952169957UL, 2205587033451463728UL,
        2711034061950757499UL, 3216481090450051270UL, 3721928118949345041UL, 4227375147448638812UL,
        4732822175947932583UL, 5238269204447226354UL, 229748649317860805UL, 735195677817154576UL,
        1240642706316448347UL, 1746089734815742118UL, 2251536763315035889UL, 2756983791814329660UL,
        3262430820313623431UL, 3767877848812917202UL, 4273324877312210973UL, 4778771905811504744UL,
        5284218934310798515UL, 275698379181432966UL, 781145407680726737UL, 1286592436180020508UL,
        1792039464679314279UL, 2297486493178608050UL, 2802933521677901821UL, 3308380550177195592UL,
        3813827578676489363UL, 4319274607175783134UL, 4824721635675076905UL, 5330168664174370676UL,
        321648109045005127UL, 827095137544298898UL, 1332542166043592669UL, 1837989194542886440UL,
        2343436223042180211UL, 2848883251541473982UL, 3354330280040767753UL, 3859777308540061524UL,
        4365224337039355295UL, 4870671365538649066UL, 5376118394037942837UL, 367597838908577288UL,
        873044867407871059UL, 1378491895907164830UL, 1883938924406458601UL, 2389385952905752372UL,
        2894832981405046143UL, 3400280009904339914UL, 3905727038403633685UL, 4411174066902927456UL,
        4916621095402221227UL, 5422068123901514998UL, 413547568772149449UL, 918994597271443220UL
    };
    std::vector<Gen::result_type> seq2(seq1.size());
    std::generate(begin(seq2), end(seq2), [gen = Gen{}]() mutable {
        return gen();
    });
    CHECK(std::equal(begin(seq1), end(seq1), begin(seq2)));
}
TEST_CASE("Test LibPIC::BitReversedPattern::Base13", "[LibPIC::BitReversedPattern::Base13]")
{

    using Gen = BitReversedPattern<13>;
    REQUIRE(Gen::min() == 0);
    REQUIRE(Gen::max() == 8650415919381337933UL);

    const std::vector seq1{
        665416609183179841UL, 1330833218366359682UL, 1996249827549539523UL, 2661666436732719364UL,
        3327083045915899205UL, 3992499655099079046UL, 4657916264282258887UL, 5323332873465438728UL,
        5988749482648618569UL, 6654166091831798410UL, 7319582701014978251UL, 7984999310198158092UL,
        51185893014090757UL, 716602502197270598UL, 1382019111380450439UL, 2047435720563630280UL,
        2712852329746810121UL, 3378268938929989962UL, 4043685548113169803UL, 4709102157296349644UL,
        5374518766479529485UL, 6039935375662709326UL, 6705351984845889167UL, 7370768594029069008UL,
        8036185203212248849UL, 102371786028181514UL, 767788395211361355UL, 1433205004394541196UL,
        2098621613577721037UL, 2764038222760900878UL, 3429454831944080719UL, 4094871441127260560UL,
        4760288050310440401UL, 5425704659493620242UL, 6091121268676800083UL, 6756537877859979924UL,
        7421954487043159765UL, 8087371096226339606UL, 153557679042272271UL, 818974288225452112UL,
        1484390897408631953UL, 2149807506591811794UL, 2815224115774991635UL, 3480640724958171476UL,
        4146057334141351317UL, 4811473943324531158UL, 5476890552507710999UL, 6142307161690890840UL,
        6807723770874070681UL, 7473140380057250522UL, 8138556989240430363UL, 204743572056363028UL,
        870160181239542869UL, 1535576790422722710UL, 2200993399605902551UL, 2866410008789082392UL,
        3531826617972262233UL, 4197243227155442074UL, 4862659836338621915UL, 5528076445521801756UL,
        6193493054704981597UL, 6858909663888161438UL, 7524326273071341279UL, 8189742882254521120UL,
        255929465070453785UL, 921346074253633626UL, 1586762683436813467UL, 2252179292619993308UL,
        2917595901803173149UL, 3583012510986352990UL, 4248429120169532831UL, 4913845729352712672UL,
        5579262338535892513UL, 6244678947719072354UL, 6910095556902252195UL, 7575512166085432036UL,
        8240928775268611877UL, 307115358084544542UL, 972531967267724383UL, 1637948576450904224UL,
        2303365185634084065UL, 2968781794817263906UL, 3634198404000443747UL, 4299615013183623588UL,
        4965031622366803429UL, 5630448231549983270UL, 6295864840733163111UL, 6961281449916342952UL,
        7626698059099522793UL, 8292114668282702634UL, 358301251098635299UL, 1023717860281815140UL,
        1689134469464994981UL, 2354551078648174822UL, 3019967687831354663UL, 3685384297014534504UL,
        4350800906197714345UL, 5016217515380894186UL, 5681634124564074027UL, 6347050733747253868UL
    };
    std::vector<Gen::result_type> seq2(seq1.size());
    std::generate(begin(seq2), end(seq2), [gen = Gen{}]() mutable {
        return gen();
    });
    CHECK(std::equal(begin(seq1), end(seq1), begin(seq2)));
}
TEST_CASE("Test LibPIC::BitReversedPattern::Base17", "[LibPIC::BitReversedPattern::Base17]")
{
    using Gen = BitReversedPattern<17>;
    REQUIRE(Gen::min() == 0);
    REQUIRE(Gen::max() == 2862423051509815793UL);

    const std::vector seq1{
        168377826559400929UL, 336755653118801858UL, 505133479678202787UL, 673511306237603716UL,
        841889132797004645UL, 1010266959356405574UL, 1178644785915806503UL, 1347022612475207432UL,
        1515400439034608361UL, 1683778265594009290UL, 1852156092153410219UL, 2020533918712811148UL,
        2188911745272212077UL, 2357289571831613006UL, 2525667398391013935UL, 2694045224950414864UL,
        9904578032905937UL, 178282404592306866UL, 346660231151707795UL, 515038057711108724UL,
        683415884270509653UL, 851793710829910582UL, 1020171537389311511UL, 1188549363948712440UL,
        1356927190508113369UL, 1525305017067514298UL, 1693682843626915227UL, 1862060670186316156UL,
        2030438496745717085UL, 2198816323305118014UL, 2367194149864518943UL, 2535571976423919872UL,
        2703949802983320801UL, 19809156065811874UL, 188186982625212803UL, 356564809184613732UL,
        524942635744014661UL, 693320462303415590UL, 861698288862816519UL, 1030076115422217448UL,
        1198453941981618377UL, 1366831768541019306UL, 1535209595100420235UL, 1703587421659821164UL,
        1871965248219222093UL, 2040343074778623022UL, 2208720901338023951UL, 2377098727897424880UL,
        2545476554456825809UL, 2713854381016226738UL, 29713734098717811UL, 198091560658118740UL,
        366469387217519669UL, 534847213776920598UL, 703225040336321527UL, 871602866895722456UL,
        1039980693455123385UL, 1208358520014524314UL, 1376736346573925243UL, 1545114173133326172UL,
        1713491999692727101UL, 1881869826252128030UL, 2050247652811528959UL, 2218625479370929888UL,
        2387003305930330817UL, 2555381132489731746UL, 2723758959049132675UL, 39618312131623748UL,
        207996138691024677UL, 376373965250425606UL, 544751791809826535UL, 713129618369227464UL,
        881507444928628393UL, 1049885271488029322UL, 1218263098047430251UL, 1386640924606831180UL,
        1555018751166232109UL, 1723396577725633038UL, 1891774404285033967UL, 2060152230844434896UL,
        2228530057403835825UL, 2396907883963236754UL, 2565285710522637683UL, 2733663537082038612UL,
        49522890164529685UL, 217900716723930614UL, 386278543283331543UL, 554656369842732472UL,
        723034196402133401UL, 891412022961534330UL, 1059789849520935259UL, 1228167676080336188UL,
        1396545502639737117UL, 1564923329199138046UL, 1733301155758538975UL, 1901678982317939904UL,
        2070056808877340833UL, 2238434635436741762UL, 2406812461996142691UL, 2575190288555543620UL
    };
    std::vector<Gen::result_type> seq2(seq1.size());
    std::generate(begin(seq2), end(seq2), [gen = Gen{}]() mutable {
        return gen();
    });
    CHECK(std::equal(begin(seq1), end(seq1), begin(seq2)));
}
TEST_CASE("Test LibPIC::BitReversedPattern::Base19", "[LibPIC::BitReversedPattern::Base19]")
{
    using Gen = BitReversedPattern<19>;
    REQUIRE(Gen::min() == 0);
    REQUIRE(Gen::max() == 15181127029874798299UL);

    const std::vector seq1{
        799006685782884121UL, 1598013371565768242UL, 2397020057348652363UL,
        3196026743131536484UL, 3995033428914420605UL, 4794040114697304726UL,
        5593046800480188847UL, 6392053486263072968UL, 7191060172045957089UL,
        7990066857828841210UL, 8789073543611725331UL, 9588080229394609452UL,
        10387086915177493573UL, 11186093600960377694UL, 11985100286743261815UL,
        12784106972526145936UL, 13583113658309030057UL, 14382120344091914178UL,
        42052983462257059UL, 841059669245141180UL, 1640066355028025301UL,
        2439073040810909422UL, 3238079726593793543UL, 4037086412376677664UL,
        4836093098159561785UL, 5635099783942445906UL, 6434106469725330027UL,
        7233113155508214148UL, 8032119841291098269UL, 8831126527073982390UL,
        9630133212856866511UL, 10429139898639750632UL, 11228146584422634753UL,
        12027153270205518874UL, 12826159955988402995UL, 13625166641771287116UL,
        14424173327554171237UL, 84105966924514118UL, 883112652707398239UL,
        1682119338490282360UL, 2481126024273166481UL, 3280132710056050602UL,
        4079139395838934723UL, 4878146081621818844UL, 5677152767404702965UL,
        6476159453187587086UL, 7275166138970471207UL, 8074172824753355328UL,
        8873179510536239449UL, 9672186196319123570UL, 10471192882102007691UL,
        11270199567884891812UL, 12069206253667775933UL, 12868212939450660054UL,
        13667219625233544175UL, 14466226311016428296UL, 126158950386771177UL,
        925165636169655298UL, 1724172321952539419UL, 2523179007735423540UL,
        3322185693518307661UL, 4121192379301191782UL, 4920199065084075903UL,
        5719205750866960024UL, 6518212436649844145UL, 7317219122432728266UL,
        8116225808215612387UL, 8915232493998496508UL, 9714239179781380629UL,
        10513245865564264750UL, 11312252551347148871UL, 12111259237130032992UL,
        12910265922912917113UL, 13709272608695801234UL, 14508279294478685355UL,
        168211933849028236UL, 967218619631912357UL, 1766225305414796478UL,
        2565231991197680599UL, 3364238676980564720UL, 4163245362763448841UL,
        4962252048546332962UL, 5761258734329217083UL, 6560265420112101204UL,
        7359272105894985325UL, 8158278791677869446UL, 8957285477460753567UL,
        9756292163243637688UL, 10555298849026521809UL, 11354305534809405930UL,
        12153312220592290051UL, 12952318906375174172UL, 13751325592158058293UL,
        14550332277940942414UL, 210264917311285295UL, 1009271603094169416UL,
        1808278288877053537UL, 2607284974659937658UL, 3406291660442821779UL,
        4205298346225705900UL
    };
    std::vector<Gen::result_type> seq2(seq1.size());
    std::generate(begin(seq2), end(seq2), [gen = Gen{}]() mutable {
        return gen();
    });
    CHECK(std::equal(begin(seq1), end(seq1), begin(seq2)));
}
