function r8_ei ( x )

!*****************************************************************************80
!
!! R8_EI evaluates the exponential integral Ei for an R8 argument.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) R8_EI, the exponential integral Ei.
!
  implicit none

  real ( kind = 8 ) r8_e1
  real ( kind = 8 ) r8_ei
  real ( kind = 8 ) x

  r8_ei = - r8_e1 ( - x )

  return
end

function r8_e1 ( x )

!*****************************************************************************80
!
!! R8_E1 evaluates the exponential integral E1 for an R8 argument.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Wayne Fullerton.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wayne Fullerton,
!    Portable Special Function Routines,
!    in Portability of Numerical Software,
!    edited by Wayne Cowell,
!    Lecture Notes in Computer Science, Volume 57,
!    Springer 1977,
!    ISBN: 978-3-540-08446-4,
!    LC: QA297.W65.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) R8_E1, the exponential integral E1.
!
  implicit none

  real ( kind = 8 ) ae10cs(50)
  real ( kind = 8 ) ae11cs(60)
  real ( kind = 8 ) ae12cs(41)
  real ( kind = 8 ) ae13cs(50)
  real ( kind = 8 ) ae14cs(64)
  real ( kind = 8 ) e11cs(29)
  real ( kind = 8 ) e12cs(25)
  real ( kind = 8 ) eta
  integer ( kind = 4 ) ntae10
  integer ( kind = 4 ) ntae11
  integer ( kind = 4 ) ntae12
  integer ( kind = 4 ) ntae13
  integer ( kind = 4 ) ntae14
  integer ( kind = 4 ) nte11
  integer ( kind = 4 ) nte12
  real ( kind = 8 ) r8_csevl
  real ( kind = 8 ) r8_e1
  integer ( kind = 4 ) r8_inits
  real ( kind = 8 ) r8_mach
  real ( kind = 8 ) x
  real ( kind = 8 ) xmax

  save ae10cs
  save ae11cs
  save ae12cs
  save ae13cs
  save ae14cs
  save e11cs
  save e12cs
  save ntae10
  save ntae11
  save ntae12
  save ntae13
  save ntae14
  save nte11
  save nte12
  save xmax

  data ae10cs( 1) / +0.3284394579616699087873844201881D-01 /
  data ae10cs( 2) / -0.1669920452031362851476184343387D-01 /
  data ae10cs( 3) / +0.2845284724361346807424899853252D-03 /
  data ae10cs( 4) / -0.7563944358516206489487866938533D-05 /
  data ae10cs( 5) / +0.2798971289450859157504843180879D-06 /
  data ae10cs( 6) / -0.1357901828534531069525563926255D-07 /
  data ae10cs( 7) / +0.8343596202040469255856102904906D-09 /
  data ae10cs( 8) / -0.6370971727640248438275242988532D-10 /
  data ae10cs( 9) / +0.6007247608811861235760831561584D-11 /
  data ae10cs(10) / -0.7022876174679773590750626150088D-12 /
  data ae10cs(11) / +0.1018302673703687693096652346883D-12 /
  data ae10cs(12) / -0.1761812903430880040406309966422D-13 /
  data ae10cs(13) / +0.3250828614235360694244030353877D-14 /
  data ae10cs(14) / -0.5071770025505818678824872259044D-15 /
  data ae10cs(15) / +0.1665177387043294298172486084156D-16 /
  data ae10cs(16) / +0.3166753890797514400677003536555D-16 /
  data ae10cs(17) / -0.1588403763664141515133118343538D-16 /
  data ae10cs(18) / +0.4175513256138018833003034618484D-17 /
  data ae10cs(19) / -0.2892347749707141906710714478852D-18 /
  data ae10cs(20) / -0.2800625903396608103506340589669D-18 /
  data ae10cs(21) / +0.1322938639539270903707580023781D-18 /
  data ae10cs(22) / -0.1804447444177301627283887833557D-19 /
  data ae10cs(23) / -0.7905384086522616076291644817604D-20 /
  data ae10cs(24) / +0.4435711366369570103946235838027D-20 /
  data ae10cs(25) / -0.4264103994978120868865309206555D-21 /
  data ae10cs(26) / -0.3920101766937117541553713162048D-21 /
  data ae10cs(27) / +0.1527378051343994266343752326971D-21 /
  data ae10cs(28) / +0.1024849527049372339310308783117D-22 /
  data ae10cs(29) / -0.2134907874771433576262711405882D-22 /
  data ae10cs(30) / +0.3239139475160028267061694700366D-23 /
  data ae10cs(31) / +0.2142183762299889954762643168296D-23 /
  data ae10cs(32) / -0.8234609419601018414700348082312D-24 /
  data ae10cs(33) / -0.1524652829645809479613694401140D-24 /
  data ae10cs(34) / +0.1378208282460639134668480364325D-24 /
  data ae10cs(35) / +0.2131311202833947879523224999253D-26 /
  data ae10cs(36) / -0.2012649651526484121817466763127D-25 /
  data ae10cs(37) / +0.1995535662263358016106311782673D-26 /
  data ae10cs(38) / +0.2798995808984003464948686520319D-26 /
  data ae10cs(39) / -0.5534511845389626637640819277823D-27 /
  data ae10cs(40) / -0.3884995396159968861682544026146D-27 /
  data ae10cs(41) / +0.1121304434507359382850680354679D-27 /
  data ae10cs(42) / +0.5566568152423740948256563833514D-28 /
  data ae10cs(43) / -0.2045482929810499700448533938176D-28 /
  data ae10cs(44) / -0.8453813992712336233411457493674D-29 /
  data ae10cs(45) / +0.3565758433431291562816111116287D-29 /
  data ae10cs(46) / +0.1383653872125634705539949098871D-29 /
  data ae10cs(47) / -0.6062167864451372436584533764778D-30 /
  data ae10cs(48) / -0.2447198043989313267437655119189D-30 /
  data ae10cs(49) / +0.1006850640933998348011548180480D-30 /
  data ae10cs(50) / +0.4623685555014869015664341461674D-31 /

  data ae11cs( 1) / +0.20263150647078889499401236517381D+00 /
  data ae11cs( 2) / -0.73655140991203130439536898728034D-01 /
  data ae11cs( 3) / +0.63909349118361915862753283840020D-02 /
  data ae11cs( 4) / -0.60797252705247911780653153363999D-03 /
  data ae11cs( 5) / -0.73706498620176629330681411493484D-04 /
  data ae11cs( 6) / +0.48732857449450183453464992488076D-04 /
  data ae11cs( 7) / -0.23837064840448290766588489460235D-05 /
  data ae11cs( 8) / -0.30518612628561521027027332246121D-05 /
  data ae11cs( 9) / +0.17050331572564559009688032992907D-06 /
  data ae11cs(10) / +0.23834204527487747258601598136403D-06 /
  data ae11cs(11) / +0.10781772556163166562596872364020D-07 /
  data ae11cs(12) / -0.17955692847399102653642691446599D-07 /
  data ae11cs(13) / -0.41284072341950457727912394640436D-08 /
  data ae11cs(14) / +0.68622148588631968618346844526664D-09 /
  data ae11cs(15) / +0.53130183120506356147602009675961D-09 /
  data ae11cs(16) / +0.78796880261490694831305022893515D-10 /
  data ae11cs(17) / -0.26261762329356522290341675271232D-10 /
  data ae11cs(18) / -0.15483687636308261963125756294100D-10 /
  data ae11cs(19) / -0.25818962377261390492802405122591D-11 /
  data ae11cs(20) / +0.59542879191591072658903529959352D-12 /
  data ae11cs(21) / +0.46451400387681525833784919321405D-12 /
  data ae11cs(22) / +0.11557855023255861496288006203731D-12 /
  data ae11cs(23) / -0.10475236870835799012317547189670D-14 /
  data ae11cs(24) / -0.11896653502709004368104489260929D-13 /
  data ae11cs(25) / -0.47749077490261778752643019349950D-14 /
  data ae11cs(26) / -0.81077649615772777976249734754135D-15 /
  data ae11cs(27) / +0.13435569250031554199376987998178D-15 /
  data ae11cs(28) / +0.14134530022913106260248873881287D-15 /
  data ae11cs(29) / +0.49451592573953173115520663232883D-16 /
  data ae11cs(30) / +0.79884048480080665648858587399367D-17 /
  data ae11cs(31) / -0.14008632188089809829248711935393D-17 /
  data ae11cs(32) / -0.14814246958417372107722804001680D-17 /
  data ae11cs(33) / -0.55826173646025601904010693937113D-18 /
  data ae11cs(34) / -0.11442074542191647264783072544598D-18 /
  data ae11cs(35) / +0.25371823879566853500524018479923D-20 /
  data ae11cs(36) / +0.13205328154805359813278863389097D-19 /
  data ae11cs(37) / +0.62930261081586809166287426789485D-20 /
  data ae11cs(38) / +0.17688270424882713734999261332548D-20 /
  data ae11cs(39) / +0.23266187985146045209674296887432D-21 /
  data ae11cs(40) / -0.67803060811125233043773831844113D-22 /
  data ae11cs(41) / -0.59440876959676373802874150531891D-22 /
  data ae11cs(42) / -0.23618214531184415968532592503466D-22 /
  data ae11cs(43) / -0.60214499724601478214168478744576D-23 /
  data ae11cs(44) / -0.65517906474348299071370444144639D-24 /
  data ae11cs(45) / +0.29388755297497724587042038699349D-24 /
  data ae11cs(46) / +0.22601606200642115173215728758510D-24 /
  data ae11cs(47) / +0.89534369245958628745091206873087D-25 /
  data ae11cs(48) / +0.24015923471098457555772067457706D-25 /
  data ae11cs(49) / +0.34118376888907172955666423043413D-26 /
  data ae11cs(50) / -0.71617071694630342052355013345279D-27 /
  data ae11cs(51) / -0.75620390659281725157928651980799D-27 /
  data ae11cs(52) / -0.33774612157467324637952920780800D-27 /
  data ae11cs(53) / -0.10479325703300941711526430332245D-27 /
  data ae11cs(54) / -0.21654550252170342240854880201386D-28 /
  data ae11cs(55) / -0.75297125745288269994689298432000D-30 /
  data ae11cs(56) / +0.19103179392798935768638084000426D-29 /
  data ae11cs(57) / +0.11492104966530338547790728833706D-29 /
  data ae11cs(58) / +0.43896970582661751514410359193600D-30 /
  data ae11cs(59) / +0.12320883239205686471647157725866D-30 /
  data ae11cs(60) / +0.22220174457553175317538581162666D-31 /

  data ae12cs( 1) / +0.63629589796747038767129887806803D+00 /
  data ae12cs( 2) / -0.13081168675067634385812671121135D+00 /
  data ae12cs( 3) / -0.84367410213053930014487662129752D-02 /
  data ae12cs( 4) / +0.26568491531006685413029428068906D-02 /
  data ae12cs( 5) / +0.32822721781658133778792170142517D-03 /
  data ae12cs( 6) / -0.23783447771430248269579807851050D-04 /
  data ae12cs( 7) / -0.11439804308100055514447076797047D-04 /
  data ae12cs( 8) / -0.14405943433238338455239717699323D-05 /
  data ae12cs( 9) / +0.52415956651148829963772818061664D-08 /
  data ae12cs(10) / +0.38407306407844323480979203059716D-07 /
  data ae12cs(11) / +0.85880244860267195879660515759344D-08 /
  data ae12cs(12) / +0.10219226625855003286339969553911D-08 /
  data ae12cs(13) / +0.21749132323289724542821339805992D-10 /
  data ae12cs(14) / -0.22090238142623144809523503811741D-10 /
  data ae12cs(15) / -0.63457533544928753294383622208801D-11 /
  data ae12cs(16) / -0.10837746566857661115340539732919D-11 /
  data ae12cs(17) / -0.11909822872222586730262200440277D-12 /
  data ae12cs(18) / -0.28438682389265590299508766008661D-14 /
  data ae12cs(19) / +0.25080327026686769668587195487546D-14 /
  data ae12cs(20) / +0.78729641528559842431597726421265D-15 /
  data ae12cs(21) / +0.15475066347785217148484334637329D-15 /
  data ae12cs(22) / +0.22575322831665075055272608197290D-16 /
  data ae12cs(23) / +0.22233352867266608760281380836693D-17 /
  data ae12cs(24) / +0.16967819563544153513464194662399D-19 /
  data ae12cs(25) / -0.57608316255947682105310087304533D-19 /
  data ae12cs(26) / -0.17591235774646878055625369408853D-19 /
  data ae12cs(27) / -0.36286056375103174394755328682666D-20 /
  data ae12cs(28) / -0.59235569797328991652558143488000D-21 /
  data ae12cs(29) / -0.76030380926310191114429136895999D-22 /
  data ae12cs(30) / -0.62547843521711763842641428479999D-23 /
  data ae12cs(31) / +0.25483360759307648606037606400000D-24 /
  data ae12cs(32) / +0.25598615731739857020168874666666D-24 /
  data ae12cs(33) / +0.71376239357899318800207052800000D-25 /
  data ae12cs(34) / +0.14703759939567568181578956800000D-25 /
  data ae12cs(35) / +0.25105524765386733555198634666666D-26 /
  data ae12cs(36) / +0.35886666387790890886583637333333D-27 /
  data ae12cs(37) / +0.39886035156771301763317759999999D-28 /
  data ae12cs(38) / +0.21763676947356220478805333333333D-29 /
  data ae12cs(39) / -0.46146998487618942367607466666666D-30 /
  data ae12cs(40) / -0.20713517877481987707153066666666D-30 /
  data ae12cs(41) / -0.51890378563534371596970666666666D-31 /

  data e11cs( 1) / -0.16113461655571494025720663927566180D+02 /
  data e11cs( 2) / +0.77940727787426802769272245891741497D+01 /
  data e11cs( 3) / -0.19554058188631419507127283812814491D+01 /
  data e11cs( 4) / +0.37337293866277945611517190865690209D+00 /
  data e11cs( 5) / -0.56925031910929019385263892220051166D-01 /
  data e11cs( 6) / +0.72110777696600918537847724812635813D-02 /
  data e11cs( 7) / -0.78104901449841593997715184089064148D-03 /
  data e11cs( 8) / +0.73880933562621681878974881366177858D-04 /
  data e11cs( 9) / -0.62028618758082045134358133607909712D-05 /
  data e11cs(10) / +0.46816002303176735524405823868362657D-06 /
  data e11cs(11) / -0.32092888533298649524072553027228719D-07 /
  data e11cs(12) / +0.20151997487404533394826262213019548D-08 /
  data e11cs(13) / -0.11673686816697793105356271695015419D-09 /
  data e11cs(14) / +0.62762706672039943397788748379615573D-11 /
  data e11cs(15) / -0.31481541672275441045246781802393600D-12 /
  data e11cs(16) / +0.14799041744493474210894472251733333D-13 /
  data e11cs(17) / -0.65457091583979673774263401588053333D-15 /
  data e11cs(18) / +0.27336872223137291142508012748799999D-16 /
  data e11cs(19) / -0.10813524349754406876721727624533333D-17 /
  data e11cs(20) / +0.40628328040434303295300348586666666D-19 /
  data e11cs(21) / -0.14535539358960455858914372266666666D-20 /
  data e11cs(22) / +0.49632746181648636830198442666666666D-22 /
  data e11cs(23) / -0.16208612696636044604866560000000000D-23 /
  data e11cs(24) / +0.50721448038607422226431999999999999D-25 /
  data e11cs(25) / -0.15235811133372207813973333333333333D-26 /
  data e11cs(26) / +0.44001511256103618696533333333333333D-28 /
  data e11cs(27) / -0.12236141945416231594666666666666666D-29 /
  data e11cs(28) / +0.32809216661066001066666666666666666D-31 /
  data e11cs(29) / -0.84933452268306432000000000000000000D-33 /

  data e12cs( 1) / -0.3739021479220279511668698204827D-01 /
  data e12cs( 2) / +0.4272398606220957726049179176528D-01 /
  data e12cs( 3) / -0.130318207984970054415392055219726D+00 /
  data e12cs( 4) / +0.144191240246988907341095893982137D-01 /
  data e12cs( 5) / -0.134617078051068022116121527983553D-02 /
  data e12cs( 6) / +0.107310292530637799976115850970073D-03 /
  data e12cs( 7) / -0.742999951611943649610283062223163D-05 /
  data e12cs( 8) / +0.453773256907537139386383211511827D-06 /
  data e12cs( 9) / -0.247641721139060131846547423802912D-07 /
  data e12cs(10) / +0.122076581374590953700228167846102D-08 /
  data e12cs(11) / -0.548514148064092393821357398028261D-10 /
  data e12cs(12) / +0.226362142130078799293688162377002D-11 /
  data e12cs(13) / -0.863589727169800979404172916282240D-13 /
  data e12cs(14) / +0.306291553669332997581032894881279D-14 /
  data e12cs(15) / -0.101485718855944147557128906734933D-15 /
  data e12cs(16) / +0.315482174034069877546855328426666D-17 /
  data e12cs(17) / -0.923604240769240954484015923200000D-19 /
  data e12cs(18) / +0.255504267970814002440435029333333D-20 /
  data e12cs(19) / -0.669912805684566847217882453333333D-22 /
  data e12cs(20) / +0.166925405435387319431987199999999D-23 /
  data e12cs(21) / -0.396254925184379641856000000000000D-25 /
  data e12cs(22) / +0.898135896598511332010666666666666D-27 /
  data e12cs(23) / -0.194763366993016433322666666666666D-28 /
  data e12cs(24) / +0.404836019024630033066666666666666D-30 /
  data e12cs(25) / -0.807981567699845120000000000000000D-32 /

  data ae13cs( 1) / -0.60577324664060345999319382737747D+00 /
  data ae13cs( 2) / -0.11253524348366090030649768852718D+00 /
  data ae13cs( 3) / +0.13432266247902779492487859329414D-01 /
  data ae13cs( 4) / -0.19268451873811457249246838991303D-02 /
  data ae13cs( 5) / +0.30911833772060318335586737475368D-03 /
  data ae13cs( 6) / -0.53564132129618418776393559795147D-04 /
  data ae13cs( 7) / +0.98278128802474923952491882717237D-05 /
  data ae13cs( 8) / -0.18853689849165182826902891938910D-05 /
  data ae13cs( 9) / +0.37494319356894735406964042190531D-06 /
  data ae13cs(10) / -0.76823455870552639273733465680556D-07 /
  data ae13cs(11) / +0.16143270567198777552956300060868D-07 /
  data ae13cs(12) / -0.34668022114907354566309060226027D-08 /
  data ae13cs(13) / +0.75875420919036277572889747054114D-09 /
  data ae13cs(14) / -0.16886433329881412573514526636703D-09 /
  data ae13cs(15) / +0.38145706749552265682804250927272D-10 /
  data ae13cs(16) / -0.87330266324446292706851718272334D-11 /
  data ae13cs(17) / +0.20236728645867960961794311064330D-11 /
  data ae13cs(18) / -0.47413283039555834655210340820160D-12 /
  data ae13cs(19) / +0.11221172048389864324731799928920D-12 /
  data ae13cs(20) / -0.26804225434840309912826809093395D-13 /
  data ae13cs(21) / +0.64578514417716530343580369067212D-14 /
  data ae13cs(22) / -0.15682760501666478830305702849194D-14 /
  data ae13cs(23) / +0.38367865399315404861821516441408D-15 /
  data ae13cs(24) / -0.94517173027579130478871048932556D-16 /
  data ae13cs(25) / +0.23434812288949573293896666439133D-16 /
  data ae13cs(26) / -0.58458661580214714576123194419882D-17 /
  data ae13cs(27) / +0.14666229867947778605873617419195D-17 /
  data ae13cs(28) / -0.36993923476444472706592538274474D-18 /
  data ae13cs(29) / +0.93790159936721242136014291817813D-19 /
  data ae13cs(30) / -0.23893673221937873136308224087381D-19 /
  data ae13cs(31) / +0.61150624629497608051934223837866D-20 /
  data ae13cs(32) / -0.15718585327554025507719853288106D-20 /
  data ae13cs(33) / +0.40572387285585397769519294491306D-21 /
  data ae13cs(34) / -0.10514026554738034990566367122773D-21 /
  data ae13cs(35) / +0.27349664930638667785806003131733D-22 /
  data ae13cs(36) / -0.71401604080205796099355574271999D-23 /
  data ae13cs(37) / +0.18705552432235079986756924211199D-23 /
  data ae13cs(38) / -0.49167468166870480520478020949333D-24 /
  data ae13cs(39) / +0.12964988119684031730916087125333D-24 /
  data ae13cs(40) / -0.34292515688362864461623940437333D-25 /
  data ae13cs(41) / +0.90972241643887034329104820906666D-26 /
  data ae13cs(42) / -0.24202112314316856489934847999999D-26 /
  data ae13cs(43) / +0.64563612934639510757670475093333D-27 /
  data ae13cs(44) / -0.17269132735340541122315987626666D-27 /
  data ae13cs(45) / +0.46308611659151500715194231466666D-28 /
  data ae13cs(46) / -0.12448703637214131241755170133333D-28 /
  data ae13cs(47) / +0.33544574090520678532907007999999D-29 /
  data ae13cs(48) / -0.90598868521070774437543935999999D-30 /
  data ae13cs(49) / +0.24524147051474238587273216000000D-30 /
  data ae13cs(50) / -0.66528178733552062817107967999999D-31 /

  data ae14cs( 1) / -0.1892918000753016825495679942820D+00 /
  data ae14cs( 2) / -0.8648117855259871489968817056824D-01 /
  data ae14cs( 3) / +0.7224101543746594747021514839184D-02 /
  data ae14cs( 4) / -0.8097559457557386197159655610181D-03 /
  data ae14cs( 5) / +0.1099913443266138867179251157002D-03 /
  data ae14cs( 6) / -0.1717332998937767371495358814487D-04 /
  data ae14cs( 7) / +0.2985627514479283322825342495003D-05 /
  data ae14cs( 8) / -0.5659649145771930056560167267155D-06 /
  data ae14cs( 9) / +0.1152680839714140019226583501663D-06 /
  data ae14cs(10) / -0.2495030440269338228842128765065D-07 /
  data ae14cs(11) / +0.5692324201833754367039370368140D-08 /
  data ae14cs(12) / -0.1359957664805600338490030939176D-08 /
  data ae14cs(13) / +0.3384662888760884590184512925859D-09 /
  data ae14cs(14) / -0.8737853904474681952350849316580D-10 /
  data ae14cs(15) / +0.2331588663222659718612613400470D-10 /
  data ae14cs(16) / -0.6411481049213785969753165196326D-11 /
  data ae14cs(17) / +0.1812246980204816433384359484682D-11 /
  data ae14cs(18) / -0.5253831761558460688819403840466D-12 /
  data ae14cs(19) / +0.1559218272591925698855028609825D-12 /
  data ae14cs(20) / -0.4729168297080398718476429369466D-13 /
  data ae14cs(21) / +0.1463761864393243502076199493808D-13 /
  data ae14cs(22) / -0.4617388988712924102232173623604D-14 /
  data ae14cs(23) / +0.1482710348289369323789239660371D-14 /
  data ae14cs(24) / -0.4841672496239229146973165734417D-15 /
  data ae14cs(25) / +0.1606215575700290408116571966188D-15 /
  data ae14cs(26) / -0.5408917538957170947895023784252D-16 /
  data ae14cs(27) / +0.1847470159346897881370231402310D-16 /
  data ae14cs(28) / -0.6395830792759094470500610425050D-17 /
  data ae14cs(29) / +0.2242780721699759457250233276170D-17 /
  data ae14cs(30) / -0.7961369173983947552744555308646D-18 /
  data ae14cs(31) / +0.2859308111540197459808619929272D-18 /
  data ae14cs(32) / -0.1038450244701137145900697137446D-18 /
  data ae14cs(33) / +0.3812040607097975780866841008319D-19 /
  data ae14cs(34) / -0.1413795417717200768717562723696D-19 /
  data ae14cs(35) / +0.5295367865182740958305442594815D-20 /
  data ae14cs(36) / -0.2002264245026825902137211131439D-20 /
  data ae14cs(37) / +0.7640262751275196014736848610918D-21 /
  data ae14cs(38) / -0.2941119006868787883311263523362D-21 /
  data ae14cs(39) / +0.1141823539078927193037691483586D-21 /
  data ae14cs(40) / -0.4469308475955298425247020718489D-22 /
  data ae14cs(41) / +0.1763262410571750770630491408520D-22 /
  data ae14cs(42) / -0.7009968187925902356351518262340D-23 /
  data ae14cs(43) / +0.2807573556558378922287757507515D-23 /
  data ae14cs(44) / -0.1132560944981086432141888891562D-23 /
  data ae14cs(45) / +0.4600574684375017946156764233727D-24 /
  data ae14cs(46) / -0.1881448598976133459864609148108D-24 /
  data ae14cs(47) / +0.7744916111507730845444328478037D-25 /
  data ae14cs(48) / -0.3208512760585368926702703826261D-25 /
  data ae14cs(49) / +0.1337445542910839760619930421384D-25 /
  data ae14cs(50) / -0.5608671881802217048894771735210D-26 /
  data ae14cs(51) / +0.2365839716528537483710069473279D-26 /
  data ae14cs(52) / -0.1003656195025305334065834526856D-26 /
  data ae14cs(53) / +0.4281490878094161131286642556927D-27 /
  data ae14cs(54) / -0.1836345261815318199691326958250D-27 /
  data ae14cs(55) / +0.7917798231349540000097468678144D-28 /
  data ae14cs(56) / -0.3431542358742220361025015775231D-28 /
  data ae14cs(57) / +0.1494705493897103237475066008917D-28 /
  data ae14cs(58) / -0.6542620279865705439739042420053D-29 /
  data ae14cs(59) / +0.2877581395199171114340487353685D-29 /
  data ae14cs(60) / -0.1271557211796024711027981200042D-29 /
  data ae14cs(61) / +0.5644615555648722522388044622506D-30 /
  data ae14cs(62) / -0.2516994994284095106080616830293D-30 /
  data ae14cs(63) / +0.1127259818927510206370368804181D-30 /
  data ae14cs(64) / -0.5069814875800460855562584719360D-31 /

  data ntae10 / 0 /
  data ntae11 / 0 /
  data ntae12 / 0 /
  data nte11 / 0 /
  data nte12 / 0 /
  data ntae13 / 0 /
  data ntae14 / 0 /
  data xmax / 0.0D+00 /

  if ( ntae10 == 0 ) then
    eta = 0.1D+00 * r8_mach ( 3 )
    ntae10 = r8_inits ( ae10cs, 50, eta )
    ntae11 = r8_inits ( ae11cs, 60, eta )
    ntae12 = r8_inits ( ae12cs, 41, eta )
    nte11 = r8_inits ( e11cs, 29, eta )
    nte12 = r8_inits ( e12cs, 25, eta )
    ntae13 = r8_inits ( ae13cs, 50, eta )
    ntae14 = r8_inits ( ae14cs, 64, eta )
    xmax = - log ( r8_mach ( 1 ) )
    xmax = xmax - log ( xmax )
  end if

  if ( x <= - 32.0D+00 ) then
    r8_e1 = exp ( - x ) / x * ( 1.0D+00 &
      + r8_csevl ( 64.0D+00 / x + 1.0D+00, ae10cs, ntae10 ) )
  else if ( x <= - 8.0D+00 ) then
    r8_e1 = exp ( - x ) / x * ( 1.0D+00 &
      + r8_csevl ( ( 64.0D+00 / x + 5.0D+00 ) / 3.0D+00, ae11cs, ntae11 ) )
  else if ( x <= - 4.0D+00 ) then
    r8_e1 = exp ( - x ) / x * (1.0D+00 &
      + r8_csevl ( 16.0D+00 / x + 3.0D+00, ae12cs, ntae12 ) )
  else if ( x <= - 1.0D+00 ) then
    r8_e1 = - log ( - x ) &
      + r8_csevl ( ( 2.0D+00 * x + 5.0D+00 ) / 3.0D+00, e11cs, nte11 )
  else if ( x == 0.0D+00 ) then
    !write ( *, '(a)' ) ' '
    !write ( *, '(a)' ) 'R8_E1 - Fatal error!'
    !write ( *, '(a)' ) '  X is zero.'
    !stop
    r8_e1 = 0.0D+00
  else if ( x <= 1.0D+00 ) then
    r8_e1 = ( - log ( abs ( x ) ) - 0.6875D+00 + x ) &
      + r8_csevl ( x, e12cs, nte12 )
  else if ( x <= 4.0D+00 ) then
    r8_e1 = exp ( - x ) / x * ( 1.0D+00 &
      + r8_csevl ( ( 8.0D+00 / x - 5.0D+00 ) / 3.0D+00, ae13cs, ntae13 ) )
  else if ( x <= xmax ) then
    r8_e1 = exp ( - x ) / x * ( 1.0D+00 &
      + r8_csevl ( 8.0D+00 / x - 1.0D+00, ae14cs, ntae14 ) )
  else
    r8_e1 = 0.0D+00
  end if

  return
end function

function r8_mach ( i )
!*****************************************************************************80
!
!! R8_MACH returns real ( kind = 8 ) real machine-dependent constants.
!
!  Discussion:
!
!    R8_MACH can be used to obtain machine-dependent parameters
!    for the local machine environment.  It is a function
!    with one input argument, and can be called as follows:
!
!      D = R8_MACH ( I )
!
!    where I=1,...,5.  The output value of D above is
!    determined by the input value of I:.
!
!    R8_MACH ( 1) = B^(EMIN-1), the smallest positive magnitude.
!    R8_MACH ( 2) = B^EMAX*(1 - B^(-T)), the largest magnitude.
!    R8_MACH ( 3) = B^(-T), the smallest relative spacing.
!    R8_MACH ( 4) = B^(1-T), the largest relative spacing.
!    R8_MACH ( 5) = LOG10(B)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 April 2007
!
!  Author:
!
!    Original FORTRAN77 version by Phyllis Fox, Andrew Hall, Norman Schryer.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Phyllis Fox, Andrew Hall, Norman Schryer,
!    Algorithm 528:
!    Framework for a Portable Library,
!    ACM Transactions on Mathematical Software,
!    Volume 4, Number 2, June 1978, page 176-188.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the index of the desired constant.
!
!    Output, real ( kind = 8 ) R8_MACH, the value of the constant.
!
  implicit none

  real ( kind = 8 ) r8_mach
  integer ( kind = 4 ) i

  if ( i < 1 ) then
    !write ( *, '(a)' ) ' '
    !write ( *, '(a)' ) 'R8_MACH - Fatal error!'
    !write ( *, '(a)' ) '  The input argument I is out of bounds.'
    !write ( *, '(a)' ) '  Legal values satisfy 1 <= I <= 5.'
    !write ( *, '(a,i12)' ) '  I = ', i
    !stop
    r8_mach = 0.0D+00
  else if ( i == 1 ) then
    r8_mach = 4.450147717014403D-308
  else if ( i == 2 ) then
    r8_mach = 8.988465674311579D+307
  else if ( i == 3 ) then
    r8_mach = 1.110223024625157D-016
  else if ( i == 4 ) then
    r8_mach = 2.220446049250313D-016
  else if ( i == 5 ) then
    r8_mach = 0.301029995663981D+000
  else if ( 5 < i ) then
    !write ( *, '(a)' ) ' '
    !write ( *, '(a)' ) 'R8_MACH - Fatal error!'
    !write ( *, '(a)' ) '  The input argument I is out of bounds.'
    !write ( *, '(a)' ) '  Legal values satisfy 1 <= I <= 5.'
    !write ( *, '(a,i12)' ) '  I = ', i
    !stop
    r8_mach = 0.0D+00
  end if

  return
end

function r8_csevl ( x, a, n )
!*****************************************************************************80
!
!! R8_CSEVL evaluates a Chebyshev series.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Roger Broucke,
!    Algorithm 446:
!    Ten Subroutines for the Manipulation of Chebyshev Series,
!    Communications of the ACM,
!    Volume 16, Number 4, April 1973, pages 254-256.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the evaluation point.
!
!    Input, real ( kind = 8 ) A(N), the Chebyshev coefficients.
!
!    Input, integer ( kind = 4 ) N, the number of Chebyshev coefficients.
!
!    Output, real ( kind = 8 ) R8_CSEVL, the Chebyshev series evaluated at X.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) b0
  real ( kind = 8 ) b1
  real ( kind = 8 ) b2
  integer ( kind = 4 ) i
  real ( kind = 8 ) r8_csevl
  real ( kind = 8 ) twox
  real ( kind = 8 ) x

  r8_csevl = 0.0D+00
  if ( n < 1 ) then
    !write ( *, '(a)' ) ' '
    !write ( *, '(a)' ) 'R8_CSEVL - Fatal error!'
    !write ( *, '(a)' ) '  Number of terms <= 0.'
    !stop
  else if ( 1000 < n ) then
    !write ( *, '(a)' ) ' '
    !write ( *, '(a)' ) 'R8_CSEVL - Fatal error!'
    !write ( *, '(a)' ) '  Number of terms > 1000.'
    !stop
  else if ( x < -1.1D+00 .or. 1.1D+00 < x ) then
    !write ( *, '(a)' ) ' '
    !write ( *, '(a)' ) 'R8_CSEVL - Fatal error!'
    !write ( *, '(a)' ) '  X outside (-1,+1)'
    !write ( *, '(a,g14.6)' ) '  X = ', x
    !stop
  else
      twox = 2.0D+00 * x
      b1 = 0.0D+00
      b0 = 0.0D+00
      b2 = 0.0D+00

      do i = n, 1, -1
        b2 = b1
        b1 = b0
        b0 = twox * b1 - b2 + a(i)
      end do

      r8_csevl = 0.5D+00 * ( b0 - b2 )
  end if
  return
end

function r8_inits ( dos, nos, eta )
!*****************************************************************************80
!
!! R8_INITS initializes a Chebyshev series.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Roger Broucke,
!    Algorithm 446:
!    Ten Subroutines for the Manipulation of Chebyshev Series,
!    Communications of the ACM,
!    Volume 16, Number 4, April 1973, pages 254-256.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) DOS(NOS), the Chebyshev coefficients.
!
!    Input, integer ( kind = 4 ) NOS, the number of coefficients.
!
!    Input, real ( kind = 8 ) ETA, the desired accuracy.
!
!    Output, integer ( kind = 4 ) R8_INITS, the number of terms of the
!    series needed to ensure the requested accuracy.
!
  implicit none

  integer ( kind = 4 ) nos

  real ( kind = 8 ) dos(nos)
  real ( kind = 8 ) err
  real ( kind = 8 ) eta
  integer ( kind = 4 ) i
  integer ( kind = 4 ) r8_inits

  if ( nos < 1 ) then
    !write ( *, '(a)' ) ' '
    !write ( *, '(a)' ) 'R8_INITS - Fatal error!'
    !write ( *, '(a)' ) '  Number of coefficients < 1.'
    !stop
  end if

  err = 0.0D+00

  do i = nos, 1, -1
    err = err + abs ( dos(i) )
    if ( eta < err ) then
      r8_inits = i
      return
    end if
  end do

  r8_inits = nos
  !write ( *, '(a)' ) ' '
  !write ( *, '(a)' ) 'R8_INITS - Warning!'
  !write ( *, '(a)' ) '  ETA may be too small.'

  return
end
