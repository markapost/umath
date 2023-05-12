#pragma once
fix SinTable[FIX_RAD_MAX] = {
         0U,        402U,        804U,       1206U,       1608U,       2010U, 
      2412U,       2814U,       3215U,       3617U,       4018U,       4420U, 
      4821U,       5222U,       5622U,       6023U,       6423U,       6823U, 
      7223U,       7623U,       8022U,       8421U,       8819U,       9218U, 
      9616U,      10013U,      10410U,      10807U,      11204U,      11600U, 
     11995U,      12390U,      12785U,      13179U,      13573U,      13966U, 
     14359U,      14751U,      15142U,      15533U,      15924U,      16313U, 
     16702U,      17091U,      17479U,      17866U,      18253U,      18639U, 
     19024U,      19408U,      19792U,      20175U,      20557U,      20938U, 
     21319U,      21699U,      22078U,      22456U,      22834U,      23210U, 
     23586U,      23960U,      24334U,      24707U,      25079U,      25450U, 
     25820U,      26189U,      26557U,      26925U,      27291U,      27656U, 
     28020U,      28383U,      28745U,      29106U,      29465U,      29824U, 
     30181U,      30538U,      30893U,      31247U,      31600U,      31952U, 
     32302U,      32651U,      33000U,      33346U,      33692U,      34036U, 
     34379U,      34721U,      35061U,      35400U,      35738U,      36074U, 
     36409U,      36743U,      37075U,      37406U,      37736U,      38064U, 
     38390U,      38716U,      39039U,      39362U,      39682U,      40002U, 
     40319U,      40636U,      40950U,      41264U,      41575U,      41885U, 
     42194U,      42501U,      42806U,      43110U,      43412U,      43712U, 
     44011U,      44308U,      44604U,      44897U,      45189U,      45480U, 
     45768U,      46055U,      46341U,      46624U,      46906U,      47186U, 
     47464U,      47740U,      48015U,      48288U,      48559U,      48828U, 
     49095U,      49360U,      49624U,      49886U,      50146U,      50404U, 
     50660U,      50914U,      51166U,      51416U,      51665U,      51911U, 
     52156U,      52398U,      52639U,      52877U,      53114U,      53348U, 
     53581U,      53811U,      54040U,      54266U,      54491U,      54713U, 
     54934U,      55152U,      55368U,      55582U,      55794U,      56004U, 
     56212U,      56417U,      56621U,      56822U,      57022U,      57219U, 
     57414U,      57607U,      57797U,      57986U,      58172U,      58356U, 
     58538U,      58718U,      58895U,      59070U,      59243U,      59414U, 
     59583U,      59749U,      59913U,      60075U,      60235U,      60392U, 
     60547U,      60700U,      60850U,      60998U,      61144U,      61288U, 
     61429U,      61568U,      61705U,      61839U,      61971U,      62101U, 
     62228U,      62353U,      62475U,      62596U,      62714U,      62829U, 
     62942U,      63053U,      63162U,      63268U,      63371U,      63473U, 
     63572U,      63668U,      63762U,      63854U,      63943U,      64030U, 
     64115U,      64197U,      64276U,      64354U,      64428U,      64501U, 
     64571U,      64638U,      64703U,      64766U,      64826U,      64884U, 
     64939U,      64992U,      65043U,      65091U,      65136U,      65179U, 
     65220U,      65258U,      65294U,      65327U,      65358U,      65386U, 
     65412U,      65436U,      65457U,      65475U,      65491U,      65505U, 
     65516U,      65524U,      65531U,      65534U,      65535U,      65534U, 
     65531U,      65524U,      65516U,      65505U,      65491U,      65475U, 
     65457U,      65436U,      65412U,      65386U,      65358U,      65327U, 
     65294U,      65258U,      65220U,      65179U,      65136U,      65091U, 
     65043U,      64992U,      64939U,      64884U,      64826U,      64766U, 
     64703U,      64638U,      64571U,      64501U,      64428U,      64353U, 
     64276U,      64196U,      64114U,      64030U,      63943U,      63854U, 
     63762U,      63668U,      63571U,      63472U,      63371U,      63267U, 
     63161U,      63053U,      62942U,      62829U,      62713U,      62595U, 
     62475U,      62353U,      62228U,      62100U,      61971U,      61839U, 
     61704U,      61568U,      61429U,      61288U,      61144U,      60998U, 
     60850U,      60699U,      60547U,      60392U,      60234U,      60075U, 
     59913U,      59749U,      59583U,      59414U,      59243U,      59070U, 
     58895U,      58717U,      58538U,      58356U,      58172U,      57985U, 
     57797U,      57606U,      57413U,      57218U,      57021U,      56822U, 
     56621U,      56417U,      56211U,      56004U,      55794U,      55582U, 
     55368U,      55151U,      54933U,      54713U,      54490U,      54266U, 
     54039U,      53811U,      53580U,      53348U,      53113U,      52877U, 
     52638U,      52398U,      52155U,      51911U,      51664U,      51416U, 
     51166U,      50913U,      50659U,      50403U,      50145U,      49885U, 
     49624U,      49360U,      49094U,      48827U,      48558U,      48287U, 
     48014U,      47740U,      47463U,      47185U,      46905U,      46623U, 
     46340U,      46055U,      45768U,      45479U,      45189U,      44897U, 
     44603U,      44307U,      44010U,      43712U,      43411U,      43109U, 
     42805U,      42500U,      42193U,      41885U,      41575U,      41263U, 
     40950U,      40635U,      40319U,      40001U,      39682U,      39361U, 
     39039U,      38715U,      38390U,      38063U,      37735U,      37406U, 
     37075U,      36742U,      36409U,      36074U,      35737U,      35400U, 
     35061U,      34720U,      34378U,      34035U,      33691U,      33346U, 
     32999U,      32651U,      32301U,      31951U,      31599U,      31246U, 
     30892U,      30537U,      30181U,      29823U,      29465U,      29105U, 
     28744U,      28382U,      28019U,      27655U,      27290U,      26924U, 
     26557U,      26189U,      25819U,      25449U,      25078U,      24706U, 
     24333U,      23960U,      23585U,      23209U,      22833U,      22455U, 
     22077U,      21698U,      21318U,      20938U,      20556U,      20174U, 
     19791U,      19407U,      19023U,      18638U,      18252U,      17865U, 
     17478U,      17090U,      16702U,      16312U,      15923U,      15532U, 
     15141U,      14750U,      14358U,      13965U,      13572U,      13178U, 
     12784U,      12390U,      11994U,      11599U,      11203U,      10806U, 
     10410U,      10012U,       9615U,       9217U,       8819U,       8420U, 
      8021U,       7622U,       7222U,       6822U,       6422U,       6022U, 
      5622U,       5221U,       4820U,       4419U,       4017U,       3616U, 
      3214U,       2813U,       2411U,       2009U,       1607U,       1205U, 
       803U,        401U,          0U, 4294966894U, 4294966491U, 4294966089U, 
4294965687U, 4294965285U, 4294964883U, 4294964482U, 4294964080U, 4294963678U, 
4294963277U, 4294962876U, 4294962475U, 4294962074U, 4294961673U, 4294961272U, 
4294960872U, 4294960472U, 4294960072U, 4294959673U, 4294959273U, 4294958874U, 
4294958476U, 4294958077U, 4294957680U, 4294957282U, 4294956885U, 4294956488U, 
4294956091U, 4294955695U, 4294955300U, 4294954905U, 4294954510U, 4294954116U, 
4294953722U, 4294953329U, 4294952937U, 4294952545U, 4294952153U, 4294951762U, 
4294951372U, 4294950982U, 4294950593U, 4294950204U, 4294949816U, 4294949429U, 
4294949043U, 4294948657U, 4294948272U, 4294947887U, 4294947503U, 4294947120U, 
4294946738U, 4294946357U, 4294945976U, 4294945596U, 4294945217U, 4294944839U, 
4294944462U, 4294944085U, 4294943710U, 4294943335U, 4294942961U, 4294942588U, 
4294942216U, 4294941845U, 4294941475U, 4294941106U, 4294940738U, 4294940371U, 
4294940005U, 4294939639U, 4294939275U, 4294938912U, 4294938550U, 4294938190U, 
4294937830U, 4294937471U, 4294937114U, 4294936757U, 4294936402U, 4294936048U, 
4294935695U, 4294935344U, 4294934993U, 4294934644U, 4294934296U, 4294933949U, 
4294933603U, 4294933259U, 4294932916U, 4294932574U, 4294932234U, 4294931895U, 
4294931557U, 4294931221U, 4294930886U, 4294930552U, 4294930220U, 4294929889U, 
4294929560U, 4294929231U, 4294928905U, 4294928580U, 4294928256U, 4294927934U, 
4294927613U, 4294927294U, 4294926976U, 4294926660U, 4294926345U, 4294926032U, 
4294925720U, 4294925410U, 4294925102U, 4294924795U, 4294924489U, 4294924186U, 
4294923884U, 4294923583U, 4294923284U, 4294922987U, 4294922692U, 4294922398U, 
4294922106U, 4294921816U, 4294921527U, 4294921240U, 4294920955U, 4294920671U, 
4294920390U, 4294920110U, 4294919832U, 4294919555U, 4294919281U, 4294919008U, 
4294918737U, 4294918468U, 4294918200U, 4294917935U, 4294917671U, 4294917410U, 
4294917150U, 4294916892U, 4294916636U, 4294916382U, 4294916129U, 4294915879U, 
4294915631U, 4294915384U, 4294915140U, 4294914897U, 4294914657U, 4294914418U, 
4294914182U, 4294913947U, 4294913715U, 4294913484U, 4294913256U, 4294913029U, 
4294912805U, 4294912582U, 4294912362U, 4294912144U, 4294911928U, 4294911713U, 
4294911502U, 4294911292U, 4294911084U, 4294910878U, 4294910675U, 4294910473U, 
4294910274U, 4294910077U, 4294909882U, 4294909689U, 4294909498U, 4294909310U, 
4294909124U, 4294908940U, 4294908758U, 4294908578U, 4294908400U, 4294908225U, 
4294908052U, 4294907881U, 4294907713U, 4294907546U, 4294907382U, 4294907221U, 
4294907061U, 4294906904U, 4294906749U, 4294906596U, 4294906445U, 4294906297U, 
4294906151U, 4294906008U, 4294905867U, 4294905728U, 4294905591U, 4294905457U, 
4294905325U, 4294905195U, 4294905068U, 4294904943U, 4294904820U, 4294904700U, 
4294904582U, 4294904467U, 4294904353U, 4294904243U, 4294904134U, 4294904028U, 
4294903924U, 4294903823U, 4294903724U, 4294903628U, 4294903534U, 4294903442U, 
4294903353U, 4294903266U, 4294903181U, 4294903099U, 4294903020U, 4294902942U, 
4294902867U, 4294902795U, 4294902725U, 4294902658U, 4294902593U, 4294902530U, 
4294902470U, 4294902412U, 4294902357U, 4294902304U, 4294902253U, 4294902205U, 
4294902160U, 4294902117U, 4294902076U, 4294902038U, 4294902002U, 4294901969U, 
4294901938U, 4294901910U, 4294901884U, 4294901860U, 4294901839U, 4294901821U, 
4294901805U, 4294901791U, 4294901780U, 4294901772U, 4294901765U, 4294901762U, 
4294901761U, 4294901762U, 4294901765U, 4294901772U, 4294901780U, 4294901791U, 
4294901805U, 4294901821U, 4294901840U, 4294901860U, 4294901884U, 4294901910U, 
4294901938U, 4294901969U, 4294902002U, 4294902038U, 4294902076U, 4294902117U, 
4294902160U, 4294902206U, 4294902254U, 4294902304U, 4294902357U, 4294902412U, 
4294902470U, 4294902530U, 4294902593U, 4294902658U, 4294902726U, 4294902796U, 
4294902868U, 4294902943U, 4294903020U, 4294903100U, 4294903182U, 4294903266U, 
4294903353U, 4294903442U, 4294903534U, 4294903628U, 4294903725U, 4294903824U, 
4294903925U, 4294904029U, 4294904135U, 4294904243U, 4294904354U, 4294904467U, 
4294904583U, 4294904701U, 4294904821U, 4294904944U, 4294905069U, 4294905196U, 
4294905326U, 4294905458U, 4294905592U, 4294905729U, 4294905868U, 4294906009U, 
4294906152U, 4294906298U, 4294906446U, 4294906597U, 4294906750U, 4294906905U, 
4294907062U, 4294907222U, 4294907383U, 4294907547U, 4294907714U, 4294907882U, 
4294908053U, 4294908226U, 4294908402U, 4294908579U, 4294908759U, 4294908941U, 
4294909125U, 4294909311U, 4294909500U, 4294909690U, 4294909883U, 4294910078U, 
4294910275U, 4294910474U, 4294910676U, 4294910879U, 4294911085U, 4294911293U, 
4294911503U, 4294911715U, 4294911929U, 4294912145U, 4294912363U, 4294912584U, 
4294912806U, 4294913031U, 4294913257U, 4294913486U, 4294913716U, 4294913949U, 
4294914183U, 4294914420U, 4294914658U, 4294914899U, 4294915141U, 4294915386U, 
4294915632U, 4294915881U, 4294916131U, 4294916383U, 4294916637U, 4294916893U, 
4294917151U, 4294917411U, 4294917673U, 4294917937U, 4294918202U, 4294918469U, 
4294918739U, 4294919009U, 4294919282U, 4294919557U, 4294919833U, 4294920111U, 
4294920391U, 4294920673U, 4294920957U, 4294921242U, 4294921529U, 4294921817U, 
4294922108U, 4294922400U, 4294922694U, 4294922989U, 4294923286U, 4294923585U, 
4294923886U, 4294924188U, 4294924491U, 4294924797U, 4294925103U, 4294925412U, 
4294925722U, 4294926034U, 4294926347U, 4294926662U, 4294926978U, 4294927296U, 
4294927615U, 4294927936U, 4294928258U, 4294928582U, 4294928907U, 4294929234U, 
4294929562U, 4294929891U, 4294930222U, 4294930554U, 4294930888U, 4294931223U, 
4294931559U, 4294931897U, 4294932236U, 4294932577U, 4294932918U, 4294933261U, 
4294933606U, 4294933951U, 4294934298U, 4294934646U, 4294934995U, 4294935346U, 
4294935697U, 4294936050U, 4294936404U, 4294936760U, 4294937116U, 4294937474U, 
4294937832U, 4294938192U, 4294938553U, 4294938915U, 4294939278U, 4294939642U, 
4294940007U, 4294940373U, 4294940740U, 4294941108U, 4294941477U, 4294941847U, 
4294942218U, 4294942590U, 4294942963U, 4294943337U, 4294943712U, 4294944088U, 
4294944464U, 4294944841U, 4294945220U, 4294945599U, 4294945978U, 4294946359U, 
4294946741U, 4294947123U, 4294947506U, 4294947889U, 4294948274U, 4294948659U, 
4294949045U, 4294949432U, 4294949819U, 4294950207U, 4294950595U, 4294950984U, 
4294951374U, 4294951764U, 4294952155U, 4294952547U, 4294952939U, 4294953332U, 
4294953725U, 4294954118U, 4294954513U, 4294954907U, 4294955302U, 4294955698U, 
4294956094U, 4294956490U, 4294956887U, 4294957284U, 4294957682U, 4294958080U, 
4294958478U, 4294958877U, 4294959276U, 4294959675U, 4294960075U, 4294960474U, 
4294960874U, 4294961275U, 4294961675U, 4294962076U, 4294962477U, 4294962878U, 
4294963279U, 4294963681U, 4294964082U, 4294964484U, 4294964886U, 4294965288U, 
4294965690U, 4294966092U, 4294966494U, 4294966896U };

fix CosTable[FIX_RAD_MAX] = {
     65536U,      65534U,      65531U,      65524U,      65516U,      65505U, 
     65491U,      65475U,      65457U,      65436U,      65412U,      65386U, 
     65358U,      65327U,      65294U,      65258U,      65220U,      65179U, 
     65136U,      65091U,      65043U,      64992U,      64939U,      64884U, 
     64826U,      64766U,      64703U,      64638U,      64571U,      64501U, 
     64428U,      64353U,      64276U,      64197U,      64114U,      64030U, 
     63943U,      63854U,      63762U,      63668U,      63571U,      63473U, 
     63371U,      63268U,      63162U,      63053U,      62942U,      62829U, 
     62714U,      62596U,      62475U,      62353U,      62228U,      62100U, 
     61971U,      61839U,      61705U,      61568U,      61429U,      61288U, 
     61144U,      60998U,      60850U,      60700U,      60547U,      60392U, 
     60235U,      60075U,      59913U,      59749U,      59583U,      59414U, 
     59243U,      59070U,      58895U,      58717U,      58538U,      58356U, 
     58172U,      57986U,      57797U,      57606U,      57414U,      57219U, 
     57021U,      56822U,      56621U,      56417U,      56212U,      56004U, 
     55794U,      55582U,      55368U,      55152U,      54933U,      54713U, 
     54491U,      54266U,      54040U,      53811U,      53581U,      53348U, 
     53114U,      52877U,      52638U,      52398U,      52155U,      51911U, 
     51664U,      51416U,      51166U,      50914U,      50659U,      50403U, 
     50145U,      49886U,      49624U,      49360U,      49095U,      48827U, 
     48558U,      48287U,      48015U,      47740U,      47464U,      47185U, 
     46905U,      46624U,      46340U,      46055U,      45768U,      45479U, 
     45189U,      44897U,      44603U,      44308U,      44011U,      43712U, 
     43411U,      43109U,      42806U,      42500U,      42193U,      41885U, 
     41575U,      41263U,      40950U,      40635U,      40319U,      40001U, 
     39682U,      39361U,      39039U,      38715U,      38390U,      38064U, 
     37735U,      37406U,      37075U,      36743U,      36409U,      36074U, 
     35738U,      35400U,      35061U,      34720U,      34379U,      34036U, 
     33692U,      33346U,      32999U,      32651U,      32302U,      31951U, 
     31600U,      31247U,      30893U,      30537U,      30181U,      29824U, 
     29465U,      29105U,      28744U,      28382U,      28019U,      27655U, 
     27290U,      26924U,      26557U,      26189U,      25820U,      25450U, 
     25079U,      24707U,      24334U,      23960U,      23585U,      23210U, 
     22833U,      22456U,      22078U,      21699U,      21319U,      20938U, 
     20557U,      20174U,      19791U,      19408U,      19023U,      18638U, 
     18252U,      17866U,      17478U,      17091U,      16702U,      16313U, 
     15923U,      15533U,      15142U,      14750U,      14358U,      13966U, 
     13572U,      13179U,      12785U,      12390U,      11995U,      11599U, 
     11203U,      10807U,      10410U,      10013U,       9615U,       9217U, 
      8819U,       8420U,       8021U,       7622U,       7223U,       6823U, 
      6423U,       6022U,       5622U,       5221U,       4820U,       4419U, 
      4018U,       3616U,       3215U,       2813U,       2411U,       2009U, 
      1607U,       1205U,        803U,        401U,          0U, 4294966894U, 
4294966492U, 4294966090U, 4294965688U, 4294965286U, 4294964884U, 4294964482U, 
4294964080U, 4294963679U, 4294963277U, 4294962876U, 4294962475U, 4294962074U, 
4294961673U, 4294961273U, 4294960872U, 4294960472U, 4294960073U, 4294959673U, 
4294959274U, 4294958875U, 4294958476U, 4294958078U, 4294957680U, 4294957282U, 
4294956885U, 4294956488U, 4294956092U, 4294955696U, 4294955300U, 4294954905U, 
4294954511U, 4294954116U, 4294953723U, 4294953330U, 4294952937U, 4294952545U, 
4294952153U, 4294951762U, 4294951372U, 4294950982U, 4294950593U, 4294950205U, 
4294949817U, 4294949430U, 4294949043U, 4294948657U, 4294948272U, 4294947887U, 
4294947504U, 4294947121U, 4294946739U, 4294946357U, 4294945976U, 4294945597U, 
4294945218U, 4294944839U, 4294944462U, 4294944086U, 4294943710U, 4294943335U, 
4294942961U, 4294942588U, 4294942216U, 4294941845U, 4294941475U, 4294941106U, 
4294940738U, 4294940371U, 4294940005U, 4294939640U, 4294939276U, 4294938913U, 
4294938551U, 4294938190U, 4294937830U, 4294937472U, 4294937114U, 4294936758U, 
4294936403U, 4294936049U, 4294935696U, 4294935344U, 4294934993U, 4294934644U, 
4294934296U, 4294933949U, 4294933604U, 4294933260U, 4294932917U, 4294932575U, 
4294932234U, 4294931895U, 4294931558U, 4294931221U, 4294930886U, 4294930553U, 
4294930220U, 4294929889U, 4294929560U, 4294929232U, 4294928905U, 4294928580U, 
4294928256U, 4294927934U, 4294927613U, 4294927294U, 4294926976U, 4294926660U, 
4294926345U, 4294926032U, 4294925720U, 4294925410U, 4294925102U, 4294924795U, 
4294924490U, 4294924186U, 4294923884U, 4294923584U, 4294923285U, 4294922988U, 
4294922692U, 4294922398U, 4294922106U, 4294921816U, 4294921527U, 4294921240U, 
4294920955U, 4294920672U, 4294920390U, 4294920110U, 4294919832U, 4294919555U, 
4294919281U, 4294919008U, 4294918737U, 4294918468U, 4294918201U, 4294917935U, 
4294917672U, 4294917410U, 4294917150U, 4294916892U, 4294916636U, 4294916382U, 
4294916130U, 4294915879U, 4294915631U, 4294915385U, 4294915140U, 4294914898U, 
4294914657U, 4294914419U, 4294914182U, 4294913947U, 4294913715U, 4294913484U, 
4294913256U, 4294913029U, 4294912805U, 4294912583U, 4294912362U, 4294912144U, 
4294911928U, 4294911714U, 4294911502U, 4294911292U, 4294911084U, 4294910878U, 
4294910675U, 4294910473U, 4294910274U, 4294910077U, 4294909882U, 4294909689U, 
4294909499U, 4294909310U, 4294909124U, 4294908940U, 4294908758U, 4294908578U, 
4294908401U, 4294908225U, 4294908052U, 4294907882U, 4294907713U, 4294907547U, 
4294907383U, 4294907221U, 4294907061U, 4294906904U, 4294906749U, 4294906596U, 
4294906446U, 4294906297U, 4294906152U, 4294906008U, 4294905867U, 4294905728U, 
4294905591U, 4294905457U, 4294905325U, 4294905195U, 4294905068U, 4294904943U, 
4294904820U, 4294904700U, 4294904582U, 4294904467U, 4294904354U, 4294904243U, 
4294904134U, 4294904028U, 4294903925U, 4294903823U, 4294903724U, 4294903628U, 
4294903534U, 4294903442U, 4294903353U, 4294903266U, 4294903181U, 4294903099U, 
4294903020U, 4294902942U, 4294902868U, 4294902795U, 4294902725U, 4294902658U, 
4294902593U, 4294902530U, 4294902470U, 4294902412U, 4294902357U, 4294902304U, 
4294902253U, 4294902205U, 4294902160U, 4294902117U, 4294902076U, 4294902038U, 
4294902002U, 4294901969U, 4294901938U, 4294901910U, 4294901884U, 4294901860U, 
4294901839U, 4294901821U, 4294901805U, 4294901791U, 4294901780U, 4294901772U, 
4294901765U, 4294901762U, 4294901761U, 4294901762U, 4294901765U, 4294901772U, 
4294901780U, 4294901791U, 4294901805U, 4294901821U, 4294901839U, 4294901860U, 
4294901884U, 4294901910U, 4294901938U, 4294901969U, 4294902002U, 4294902038U, 
4294902076U, 4294902117U, 4294902160U, 4294902205U, 4294902253U, 4294902304U, 
4294902357U, 4294902412U, 4294902470U, 4294902530U, 4294902593U, 4294902658U, 
4294902725U, 4294902795U, 4294902868U, 4294902943U, 4294903020U, 4294903100U, 
4294903182U, 4294903266U, 4294903353U, 4294903442U, 4294903534U, 4294903628U, 
4294903725U, 4294903824U, 4294903925U, 4294904029U, 4294904135U, 4294904243U, 
4294904354U, 4294904467U, 4294904583U, 4294904701U, 4294904821U, 4294904944U, 
4294905069U, 4294905196U, 4294905325U, 4294905457U, 4294905592U, 4294905728U, 
4294905867U, 4294906009U, 4294906152U, 4294906298U, 4294906446U, 4294906597U, 
4294906749U, 4294906905U, 4294907062U, 4294907221U, 4294907383U, 4294907547U, 
4294907714U, 4294907882U, 4294908053U, 4294908226U, 4294908401U, 4294908579U, 
4294908759U, 4294908940U, 4294909125U, 4294909311U, 4294909499U, 4294909690U, 
4294909883U, 4294910078U, 4294910275U, 4294910474U, 4294910676U, 4294910879U, 
4294911085U, 4294911293U, 4294911503U, 4294911715U, 4294911929U, 4294912145U, 
4294912363U, 4294912583U, 4294912806U, 4294913030U, 4294913257U, 4294913485U, 
4294913716U, 4294913948U, 4294914183U, 4294914420U, 4294914658U, 4294914899U, 
4294915141U, 4294915386U, 4294915632U, 4294915880U, 4294916131U, 4294916383U, 
4294916637U, 4294916893U, 4294917151U, 4294917411U, 4294917673U, 4294917936U, 
4294918202U, 4294918469U, 4294918738U, 4294919009U, 4294919282U, 4294919557U, 
4294919833U, 4294920111U, 4294920391U, 4294920673U, 4294920956U, 4294921242U, 
4294921528U, 4294921817U, 4294922108U, 4294922400U, 4294922693U, 4294922989U, 
4294923286U, 4294923585U, 4294923885U, 4294924187U, 4294924491U, 4294924796U, 
4294925103U, 4294925412U, 4294925722U, 4294926033U, 4294926347U, 4294926661U, 
4294926978U, 4294927295U, 4294927615U, 4294927935U, 4294928258U, 4294928581U, 
4294928907U, 4294929233U, 4294929561U, 4294929891U, 4294930222U, 4294930554U, 
4294930888U, 4294931223U, 4294931559U, 4294931897U, 4294932236U, 4294932576U, 
4294932918U, 4294933261U, 4294933605U, 4294933951U, 4294934298U, 4294934646U, 
4294934995U, 4294935345U, 4294935697U, 4294936050U, 4294936404U, 4294936759U, 
4294937116U, 4294937473U, 4294937832U, 4294938192U, 4294938552U, 4294938914U, 
4294939277U, 4294939641U, 4294940006U, 4294940373U, 4294940740U, 4294941108U, 
4294941477U, 4294941847U, 4294942218U, 4294942590U, 4294942963U, 4294943337U, 
4294943711U, 4294944087U, 4294944464U, 4294944841U, 4294945219U, 4294945598U, 
4294945978U, 4294946359U, 4294946740U, 4294947122U, 4294947505U, 4294947889U, 
4294948274U, 4294948659U, 4294949045U, 4294949431U, 4294949818U, 4294950206U, 
4294950595U, 4294950984U, 4294951374U, 4294951764U, 4294952155U, 4294952547U, 
4294952939U, 4294953331U, 4294953724U, 4294954118U, 4294954512U, 4294954907U, 
4294955302U, 4294955698U, 4294956094U, 4294956490U, 4294956887U, 4294957284U, 
4294957682U, 4294958080U, 4294958478U, 4294958876U, 4294959275U, 4294959675U, 
4294960074U, 4294960474U, 4294960874U, 4294961274U, 4294961675U, 4294962076U, 
4294962477U, 4294962878U, 4294963279U, 4294963680U, 4294964082U, 4294964484U, 
4294964886U, 4294965287U, 4294965689U, 4294966091U, 4294966494U, 4294966896U, 
         1U,        403U,        805U,       1207U,       1609U,       2011U, 
      2413U,       2815U,       3216U,       3618U,       4019U,       4421U, 
      4822U,       5223U,       5624U,       6024U,       6424U,       6824U, 
      7224U,       7624U,       8023U,       8422U,       8821U,       9219U, 
      9617U,      10014U,      10412U,      10809U,      11205U,      11601U, 
     11996U,      12392U,      12786U,      13180U,      13574U,      13967U, 
     14360U,      14752U,      15143U,      15534U,      15925U,      16314U, 
     16704U,      17092U,      17480U,      17867U,      18254U,      18640U, 
     19025U,      19409U,      19793U,      20176U,      20558U,      20940U, 
     21320U,      21700U,      22079U,      22457U,      22835U,      23211U, 
     23587U,      23962U,      24335U,      24708U,      25080U,      25451U, 
     25821U,      26191U,      26559U,      26926U,      27292U,      27657U, 
     28021U,      28384U,      28746U,      29107U,      29466U,      29825U, 
     30183U,      30539U,      30894U,      31248U,      31601U,      31953U, 
     32303U,      32653U,      33001U,      33347U,      33693U,      34037U, 
     34380U,      34722U,      35062U,      35401U,      35739U,      36075U, 
     36411U,      36744U,      37076U,      37407U,      37737U,      38065U, 
     38391U,      38717U,      39040U,      39363U,      39683U,      40003U, 
     40320U,      40637U,      40951U,      41265U,      41576U,      41886U, 
     42195U,      42502U,      42807U,      43111U,      43413U,      43713U, 
     44012U,      44309U,      44604U,      44898U,      45190U,      45481U, 
     45769U,      46056U,      46341U,      46625U,      46907U,      47187U, 
     47465U,      47741U,      48016U,      48289U,      48559U,      48829U, 
     49096U,      49361U,      49625U,      49887U,      50146U,      50404U, 
     50660U,      50915U,      51167U,      51417U,      51665U,      51912U, 
     52156U,      52399U,      52639U,      52878U,      53115U,      53349U, 
     53582U,      53812U,      54041U,      54267U,      54492U,      54714U, 
     54934U,      55152U,      55369U,      55583U,      55795U,      56005U, 
     56212U,      56418U,      56622U,      56823U,      57022U,      57219U, 
     57414U,      57607U,      57798U,      57986U,      58173U,      58357U, 
     58539U,      58718U,      58896U,      59071U,      59244U,      59415U, 
     59583U,      59750U,      59914U,      60076U,      60235U,      60392U, 
     60547U,      60700U,      60851U,      60999U,      61145U,      61288U, 
     61429U,      61568U,      61705U,      61839U,      61971U,      62101U, 
     62228U,      62353U,      62476U,      62596U,      62714U,      62830U, 
     62943U,      63054U,      63162U,      63268U,      63372U,      63473U, 
     63572U,      63668U,      63762U,      63854U,      63943U,      64030U, 
     64115U,      64197U,      64277U,      64354U,      64429U,      64501U, 
     64571U,      64638U,      64704U,      64766U,      64826U,      64884U, 
     64940U,      64992U,      65043U,      65091U,      65136U,      65179U, 
     65220U,      65258U,      65294U,      65327U,      65358U,      65386U, 
     65412U,      65436U,      65457U,      65475U,      65491U,      65505U, 
     65516U,      65524U,      65531U,      65534U };
