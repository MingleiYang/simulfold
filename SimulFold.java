import java.io.*;
import java.util.Random;

public class SimulFold{

    // Model constants

    // these ratematrices come from the Pfold program:
    // Pi: stacionary distribution,
    // D: diagonal matrix (eigenvalues), 
    // V and W: eigenvector matrices, 
    // Single: 4x4 rate matrix
    // Double: 16x16 ratematrix.
    static final double[] ePiSingle = {0.364097, 
				 0.273013, 
				 0.211881, 
				 0.151009};
    
    static final double[][] eVSingle = {
	{-0.65540029317563497, 0.25105490546012943, -0.68474016636611168, -0.027451045566647427},
	{-0.65540029317563497, 0.12662289264263851, 0.5756067113246639, 0.47248414840567282},
	{-0.65540029317563497, -0.74587320437216564, -0.11391430399580099, -0.056342506516670117},	
	{-0.65540029317563497, 0.21229414613180295, 0.77015278506218832, -0.70897472897469516}
    };
    
    static final double[] eDSingle = {0, 
				-1.4234675102660919, 
				-0.96349817076723299, 
				-1.9978941673634107};
    
    static final double[][] eWSingle = {
	{-0.5555337765197319, -0.41655916673024379, -0.32328487217081525, -0.23040728173664765},
	{0.60134410527598936, 0.22742217240225829, -1.0396670142421991, 0.21090073656395161},
	{-0.70529165850969244, 0.44456474493367609, -0.068280409630862418, 0.3290073232068787},
	{-0.072532206889486625, 0.93610695167161195, -0.086632944986868063, -0.77694179979525724}
    };
    
    static final double[] ePiDouble =  { 0.001167, 0.177977, 0.001058, 0.001806, 
				   0.177977, 0.002793, 0.049043, 0.000763, 
				   0.001058, 0.049043, 0.000406, 0.266974, 
				   0.001806, 0.000763, 0.266974, 0.000391};
    
    static final double[][] eVDouble = {
	{2.415785562154005e-07, -0.0051677625638404728, 0.31136969613914461, 0.23887614624354481, 3.6136851513674749e-06, 0.50614044046817508, -0.059720838561802943, -0.98002957880161157, -0.14271272282264788, -1.9170628606272861e-06, 1.1412824416526148, 1.1593701753036456, -1.6076000366787746e-05, 0.14516795965333629, 1.4448856662305083e-05, 7.713616688315594e-06},
	{-0.00025368025803975216, 0.00017810989468302569, 0.0047393473292717337, 0.0030968326160297093, -0.0023281300324349214, 0.50614044046817508, -0.0020962975693147165, 0.0025971313227400963, -0.0034473999296903155, -0.38046050876514353, -0.011273114034743724, 1.1932581301624643, 0.00067652646180973964, 0.11210711027189454, 0.35051213659844821, 0.06385201841180721},
	{0.0017070515402485738, 0.0017232488105842316, -0.02324253597689982, -0.031016490986617237, -0.013296558140366911, 0.50614044046817508, 0.0030663288573519842, 0.4054944357785743, 0.052852949582039163, 0.18941091982091052, 1.5792818341373005, 0.8057714613420095, 0.56788066859434239, -0.025443982437077593, 0.10837138686851985, -0.27661323496405332},
	{-0.00070422320707056764, 0.019192279849367656, -0.95244676086143676, -0.56734340930398242, 0.41911341409018682, 0.50614044046817508, 0.045312854740129083, -0.11629699318599682, -0.00018559396772193661, -0.39655639376519392, 0.093382687272942452, 0.47944203153748921, 0.011257076871905729, 0.044576991959171842, 0.1240083873757748, 0.028091792784308584},
	{0.00025376630317049587, 0.00017814608184677668, 0.0047389700527586329, 0.003097009366914167, 0.0023282550141874312, 0.50614044046817508, -0.0020961011887216481, 0.0025971629542114395, -0.0034472923303881268, 0.38046111839066216, -0.01127185373451179, 1.1932511346441479, -0.00067640885851318175, 0.11210880205719324, -0.35051417571159527, -0.063850977538773765},
	{3.2919286874747863e-06, -0.079528183497582339, -0.033327190770807445, 0.12316545085970698, -2.5714774974053924e-06, 0.50614044046817508, -0.46749733846577363, -0.073466286356797639, 0.75181339922680979, 8.1011099245373516e-07, -0.0092334709193971909, 0.7644563194404721, 4.6296350510708923e-07, -0.5591417161299449, 1.2371288999433085e-06, -4.0183424725312456e-06},
	{0.00038011765280826455, -0.009941458926583947, 0.0021925976933923987, 0.00052566963021741621, 0.00058084343143415765, 0.50614044046817508, 0.0048557792043955302, 0.00036092006354416199, -0.01770498618523151, 0.44367173103818636, -0.0045137551818550106, 0.15892413553828438, 0.0064803922694584076, -0.80547190722490225, -0.081821121788630408, 0.44295157807783514},
	{-0.41643091621599387, 0.024849538089776194, -0.078901858118393522, 0.59598718188946354, -0.0015412328598857813, 0.50614044046817508, 1.3861972233918165, -0.058917095126211544, 0.41884157960579504, 0.011907150477574914, 0.014957287945354003, 0.47254741835822028, 0.0036980114301938616, -0.081205667616673666, -0.083967866842570626, 0.021259936137240081},
	{-0.001706870994249751, 0.0017230926833378723, -0.02324467374320539, -0.031017730747214912, 0.013296145535200655, 0.50614044046817508, 0.0030637892234154657, 0.4055075715302191, 0.052854032998023277, -0.18942121429869424, 1.5791729074831335, 0.80579503000357411, -0.56790661888416838, -0.025540368335443594, -0.10832127929182792, 0.27662423812479681},
	{-0.00038010244803936548, -0.0099415146666104806, 0.0021929790279459874, 0.00052550366448363692, -0.00058081712484539473, 0.50614044046817508, 0.0048549722119084534, 0.00036098680439303551, -0.017705381596563239, -0.44367140144997003, -0.004514233462071515, 0.1589278455657428, -0.0064791193476075844, -0.80545964045348328, 0.081820883014264573, -0.44295876579696064},
	{-2.7104241745686485e-06, 10.415658360216749, 0.013538807779062493, 0.022846098635235901, -2.9484005909358016e-07, 0.50614044046817508, -0.039565776441460185, -0.0027907878330858857, 0.030650818236547635, -1.3243787494305306e-07, -0.0035799207344963586, -0.58985017492047875, 1.3780825429192494e-07, -0.1585505526360261, 3.7704486096000654e-07, -7.8411643628183025e-07},
	{0.00017580551632983315, -0.0060133543369858154, 0.0017686689668575631, 0.00186182591323479, -0.0010343986382334515, 0.50614044046817508, -0.0016267572044576511, 0.00014281398281498794, 0.00022440340798533226, -0.48000328773569884, -0.0010410603100483443, -0.83876264542257661, 0.0015681587715277948, 0.076061600257684991, -0.19963218504764393, 0.040884976276425115},
	{0.00070540805982874681, 0.019190719961026271, -0.95238443059778533, -0.56735165625940753, -0.41913475034829289, 0.50614044046817508, 0.045324766858368726, -0.11629673912183509, -0.00018876195070737965, 0.39655428076671567, 0.093381691995513258, 0.47944337890155631, -0.011260026876231026, 0.044581574336225449, -0.12400489442275203, -0.028093555056879979},
	{0.4163490646742628, 0.024892979689051676, -0.078924762081796454, 0.59614914516860584, 0.0015439326043375843, 0.50614044046817508, 1.3864759082374396, -0.058924237297570556, 0.41889828339765006, -0.011867213207417597, 0.014957477252697854, 0.47263673015454216, -0.0036993519839365867, -0.081216401863632021, 0.084000283313671195, -0.02127680557981116},
	{-0.00017568516309773524, -0.0060133765318858695, 0.0017685660364350072, 0.0018616189036447476, 0.0010344326412930586, 0.50614044046817508, -0.0016273732212222598, 0.00014285473769259525, 0.00022419558986025884, 0.4800027595936035, -0.0010418097863957618, -0.83876026113899871, -0.0015682787495091765, 0.076062733720870288, 0.19963319396857354, -0.04088432050568775},
	{9.083888435224922e-06, 0.028377269348334556, 1.2471283399429116, -4.027759166213805, 1.8887958501533901e-05, 0.50614044046817508, 0.62632497917001873, -0.087227714027739128, 0.37631842527198411, 1.5343211517200414e-06, 0.012598827570287288, 0.2953187726563366, -2.8237197807506706e-07, -0.10115062319860327, 1.8989890196374633e-06, -2.2030007632560693e-06}
    };
    
    static final double[] eDDouble = {-7.232211198372914, 
				-6.9378878868913967, 
				-6.2538124817144523, 
				-5.6427104001775943, 
				-6.1176507074310429, 
				6.0000000000000032e-39, 
				-5.1372542947589652, 
				-3.9953894191381081, 
				-4.1171671559439282, 
				-0.53942525998402424,
				-1.9201142732002245, 
				-1.3392093589472744, 
				-2.5828233007199528, 
				-2.7690392475615768, 
				-1.4984178182479098, 
				-2.8332594076881743};
    
    static final double[][] eWDouble = {
	{1.0653029687346831e-06, -0.17060605967789838, 0.0068245842781747123, -0.0048058695546085136, 0.1706639270926672, 3.4755274281453048e-05, 0.070443182303486465, -1.2006361136770121, -0.0068238624772472053, -0.070440364564558738, -4.1582174405671978e-06, 0.17735578803801519, 0.0048139554111097703, 1.2004001227523815, -0.17723437351833879, 1.3421234019045227e-05},
	{-0.00013676909582085311, 0.00071889671165893492, 4.1347401699671562e-05, 0.00078606576445142639, 0.00071904277223053389, -0.0050392050769585859, -0.011057112269893345, 0.00042998908019113539, 4.1343655603384231e-05, -0.011057174265193617, 0.095901925774320385, -0.036408294222150715, 0.00078600187549028404, 0.00043074078081619375, -0.03640842860279854, 0.00025162971635369762},
	{0.090352743949050573, 0.20973773183386701, -0.0061145334754878081, -0.42771314976151364, 0.20972103562769467, -0.023153638944944918, 0.026738074485120431, -0.014969452489291614, -0.0061150958686643176, 0.026742724746191274, 0.0013667862567963718, 0.11741123915594032, -0.42768515924856598, -0.014973797884953429, 0.11740440622751472, 0.12125008539124632},
	{0.034144087543242198, 0.067507727010399532, -0.0040192979516396289, -0.12549766090986106, 0.067511579995670701, 0.042148990162236284, 0.0031576339779650622, 0.055697195575791954, -0.0040194586070546115, 0.0031566370418477092, 0.0011360827148491322, 0.060880738266154213, -0.12549948515528014, 0.055712331640305764, 0.060873969160310126, -0.1928910704649372},
	{6.6160666372766842e-06, -0.65005457435686897, -0.022070064792173362, 1.1874846667294638, 0.65008947144543283, -1.1271658789943236e-05, 0.044690458925015192, -0.0018448943520789518, 0.022069379936539544, -0.044688434879530732, -1.8779811409938125e-07, -0.43324704787429275, -1.1875451192907387, 0.0018481260138354738, 0.43326128970007344, 1.1586185589700386e-05},
	{0.002305684167264991, 0.35163560500198909, 0.0020903289194227595, 0.0035681796110373382, 0.35163560500198909, 0.0055202069951485726, 0.096896031375472966, 0.0015074867348956196, 0.0020903289194227595, 0.096896031375472966, 0.00080214890480684341, 0.52747020126084287, 0.0035681796110373382, 0.0015074867348956196, 0.52747020126084287, 0.00077251286152580246},
	{-0.018763192286605956, -0.10044464514047194, 0.00087340238584703044, 0.022031757661213514, -0.10043523551310225, -0.35165396658892645, 0.064112975838207756, 0.28474725007431367, 0.00087267900540007829, 0.064102320763574458, -0.0043246971730451257, -0.1169234666666577, 0.022037549503369402, 0.28480449643297712, -0.11696774298241364, 0.065930514686319996},
	{-0.74101155813933162, 0.29948338588961637, 0.2779620521341668, -0.13608215394434661, 0.29948703341411986, -0.1329931986626634, 0.011468404331329693, -0.029125997403734661, 0.27797105654991544, 0.011470525053111634, -0.00073412125261013898, 0.02470328119274165, -0.13608185665722489, -0.029129528176323059, 0.024710330776974685, -0.022097655105741722},
	{-0.084652246450890353, -0.31186058546690693, 0.028422340388494231, -0.00017036741813949646, -0.31185085175974769, 1.0676807262384818, -0.44134433876092144, 0.16243487038324719, 0.0284229230091072, -0.44135419545049542, 0.0063251822819019345, 0.030451125987753173, -0.00017327549262362201, 0.16245686121111313, 0.030422925453875524, 0.074788905845750861},
	{-1.1502407143032894e-08, -0.34814085580015464, 0.0010303201823930342, -0.0036821733355316623, 0.34814141363873413, 1.1637305975673383e-08, 0.11187171936734229, 4.6710434414004055e-05, -0.0010303761802638951, -0.11187163626175349, -2.7645201732629787e-10, -0.65886237398491909, 0.0036821537155563521, -4.6553764920167476e-05, 0.6588616490462289, 3.0844274237811596e-09},
	{0.1936650383775759, -0.29173935373335774, 0.24295874857948913, 0.024522872934091024, -0.29170673815135956, -0.0037512725740531808, -0.032188612957993461, 0.0016594517397703781, 0.24294199116291454, -0.032192023682801135, -0.0002113424212904857, -0.040414023751265793, 0.024522611568064377, 0.001659472742668915, -0.040443118468078243, 0.00071629863562531386},
	{0.0015187230553627303, 0.23838772547878406, 0.0009569365774423676, 0.0009719400022034628, 0.23838632792224718, 0.0023975334925359958, 0.0087488806518195728, 0.00040472081288453981, 0.00095696456765474735, 0.0087490848913288136, -0.00026881490907948522, -0.25135854577405448, 0.00097194273362237937, 0.00040479730540445062, -0.25135783125719313, 0.00012961444903676226},
	{-2.724904688259167e-05, 0.17488442078750205, 0.87266027247397349, 0.029528802181252686, -0.17485401993023741, 1.8787786306565514e-06, 0.46161591811715313, 0.0040982197298434012, -0.87270015019519687, -0.46152524444117393, 8.1265033268003813e-08, 0.60808102661161456, -0.029536540432945086, -0.0040997053617569468, -0.60812755017505327, -1.6036175709816705e-07},
	{0.00234870569377881, 0.2766202769770863, -0.00037321382807198469, 0.0011161317868675399, 0.27662445139668179, -0.021658843766687048, -0.54766425736288493, -0.00085900927062601442, -0.00037462762208753593, -0.54765591682090053, -0.00089244357125005131, 0.2815278725581114, 0.0011162465217669614, -0.00085912281956811458, 0.28153206786136692, -0.00054831773358350719},
	{2.5641578315927069e-07, 0.94865294033672543, 0.0017435753236586663, 0.0034057222048077861, -0.94865845914900226, 5.2563196650902014e-08, -0.061021499146121067, -0.00097426718111951415, -0.0017427691483685561, 0.061021321070222782, 2.3278740251895519e-09, -0.81047559844492612, -0.0034056262755894347, 0.00097464330480954672, 0.81047969450686119, 1.1291187673287387e-08},
	{4.137881251322116e-07, 0.52238130446723396, -0.013452646703336372, 0.0023320954175675176, -0.52237278895548467, -5.1608645737037905e-07, 0.99857892285500482, 0.000745651003964035, 0.013453181824637136, -0.99859512666864336, -1.4633752460392855e-08, 0.50174358888137494, -0.0023322417160903568, -0.00074624266692615486, -0.50173554121220121, -3.9595015598063342e-08}
    };
    

    // alignment prior parameters

    static final double gapPenalty = -3.0;
    static final double openingPenalty = -3.0;
    static final double closingPenalty = -3.0;
    static final double bulgePenalty = -1.0;

    // structure prior parameters

    static final int smallLoopTreshold = 8;


    // structure parameters

    static final int minHelixLength = 3;
    static final int minLoopLength = 3;

    // MCMC parameters

    static final double spanForEdgeSampling = 0.1; // works well in practice

    // reasoning for choosing 0.23: if we assume that a basepairing has about three times 
    // more Felsenstein's likelihood than the corresponding single stranded nucleotides and a helix is about
    // 8 basepairs long, and there are 5 such helices in the structure than the probability that none of them
    // will be removed in an MCMC step is 0.5.
    static final double felsRatioPower = 0.023;

    static final double helixRemoveCutoff = 0.3;
    static final double helixRemoveTreshold = 0.9;

    // these parameters are in the alignment proposal and backproposal functions, giving the score for
    // staying in an HMM indel state and jumping out. When exponentiated, they are some kind of 
    // conditional probabilities
    static final double gapExtensionInHMM  = -2.0;
    static final double gapOpeningInHMM = -8.0;

    public static double proposeNewHelix(double score, double temperature){
	
	if(Math.pow(score,temperature) > 11.0){ // works well in practice, 11 is a score of an average helix with length 3-4
	    return 0.9;
	}
	else if(score > 10e-4){
	    return 0.001;
	}
	else{
	    return 1e-20;
	}
	
	//return 0.7;
    }

    static String[] sequences;
    
    static int[] interpret(char cC){
	int[] iTable = new int[4];
	if(cC == 'a' || cC == 'A'){
	    iTable[0] = 1;
	    iTable[1] = 0;
	    iTable[2] = 0;
	    iTable[3] = 0;
	}
	else if(cC == 'u' || cC == 'U' || cC == 't' || cC == 'T'){
	    iTable[0] = 0;
	    iTable[1] = 1;
	    iTable[2] = 0;
	    iTable[3] = 0;
	}
	else if(cC == 'g' || cC == 'G'){
	    iTable[0] = 0;
	    iTable[1] = 0;
	    iTable[2] = 1;
	    iTable[3] = 0;
	}
	else if(cC == 'c' || cC == 'C'){
	    iTable[0] = 0;
	    iTable[1] = 0;
	    iTable[2] = 0;
	    iTable[3] = 1;
	}
	else if(cC == 'b' || cC == 'B'){  /* b = g|c|t */
	    iTable[0] = 0;
	    iTable[1] = 1;
	    iTable[2] = 1;
	    iTable[3] = 1;
	}
	else if(cC == 'd' || cC == 'D'){  /* d = g|a|t */
	    iTable[0] = 1;
	    iTable[1] = 1;
	    iTable[2] = 1;
	    iTable[3] = 0;
	}
	else if(cC == 'h' || cC == 'H'){  /* h = t|c|a */
	    iTable[0] =1;
	    iTable[1] =1;
	    iTable[2] =0;
	    iTable[3] =1;
	}
	else if(cC == 'v' || cC == 'V'){  /* v = g|a|c */
	    iTable[0] = 1;
	    iTable[1] = 0;
	    iTable[2] = 1;
	    iTable[3] = 1;
	}
	else if(cC == 'y' || cC == 'Y'){  /* y = t|c */
	    iTable[0] = 0;
	    iTable[1] = 1;
	    iTable[2] = 0;
	    iTable[3] = 1;
	}
	else if(cC == 'r' || cC == 'R'){  /* r = a|g */
	    iTable[0] = 1;
	    iTable[1] = 0;
	    iTable[2] = 1;
	    iTable[3] = 0;
	}
	else if(cC == 'm' || cC == 'M'){  /* m = a|c */
	    iTable[0] = 1;
	    iTable[1] = 0;
	    iTable[2] = 0;
	    iTable[3] = 1;
	}
	else if(cC == 'k' || cC == 'K'){  /* k = g|t */
	    iTable[0] = 0;
	    iTable[1] = 1;
	    iTable[2] = 1;
	    iTable[3] = 0;
	}
	else if(cC == 'w' || cC == 'W'){  /* w = a|t */
	    iTable[0] = 1;
	    iTable[1] = 1;
	    iTable[2] = 0;
	    iTable[3] = 0;
	}
	else if(cC == 's' || cC == 'S'){  /* s = g|c */
	    iTable[0] = 0;
	    iTable[1] = 0;
	    iTable[2] = 1;
	    iTable[3] = 1;
	}
	else if(cC == 'n' || cC == 'N' || cC == '-'){
	    iTable[0] = 1;
	    iTable[1] = 1;
	    iTable[2] = 1;
	    iTable[3] = 1;
	}
	else{

	    throw new Error("Character cannot be interpreted! "+cC);

	}

	return iTable;

    }

    //                                      a,  u,  g,  c
    static int[][] alignmentScoreTable = {{ 2, -1,  0, -1},  // a
				          {-1,  2, -1,  0},  // u
					  { 0,  1,  2, -1},  // g
					  {-1,  0, -1,  2}}; // c

    public static double alignmentScore(char a, char b){

	// alignmentScore is used for proposing alignments 

	int[] interpret1 = interpret(a);
	int[] interpret2 = interpret(b);
	double returnValue = 0.0;
	for(int i = 0; i < 4; i++){
	    for(int j = 0; j < 4; j++){
		returnValue += interpret1[i] * alignmentScoreTable[i][j] * interpret2[j];
	    }
	}
	return returnValue;
    }

    //                                     a,  u,  g,  c
    static int[][] basepairScoreTable = {{-6,  2, -6, -6},  // a
					 { 2, -6,  1, -6},  // u
					 {-6,  1, -6,  5},  // g
					 {-6, -6,  5, -6}}; // c

    public static int basepairScore(char a, char b){

	// basepairScore is used in helixStatistics and helixTable for proposing new helices

	int[] interpret1 = interpret(a);
	int[] interpret2 = interpret(b);
	int returnValue = 0;
	for(int i = 0; i < 4; i++){
	    for(int j = 0; j < 4; j++){
		returnValue += interpret1[i] * basepairScoreTable[i][j] * interpret2[j];
	    }
	}
	return returnValue;
    }

    public static int[][][] helixTable;

    public static int[][] helixStatistics;

    public static void setHelixStatistics(){

	// for each sequence in the alignment, calculate set of helices

	// allocate memory, use LOWER TRIANGLE (!!!) matrix to save memory

	// helixTable[sequences.length->i][sequences[i].length->j][j+1]
	// helixTable[k][i][j] = score of helix whose outer bp is at j-i (with j<i)
	//
	// helixStatistics[sequences.length->][sequences[i].length]
	// helixStatistics[k][j] = sum of scores of all helices whose outer bp is at j-i* with i* >= j
	//
	// where
	// sequences.length = # of sequences in alignment

	helixTable = new int[sequences.length][][];
	helixStatistics = new int[sequences.length][];
	for(int i = 0; i < sequences.length; i++){
	    helixTable[i] = new int[sequences[i].length()][];
	    helixStatistics[i] = new int[sequences[i].length()];
	    for(int j = 0; j < helixTable[i].length; j++){
		helixTable[i][j] = new int[j+1];
	    }
	}

	// calculate helixTable using dynamic programming for obtaining helices

	for(int k = 0; k < sequences.length; k++){
	    // initialization 7 = minLoop + 2*(minHelix-1)
	    for(int j = 0; j < minLoopLength + 2 * minHelixLength - 1; j++){
		for(int i = j; i < sequences[k].length(); i++){
		    helixTable[k][i][i-j] = 0;
		}
	    }
	    for(int i = minLoopLength + 2 * minHelixLength - 1; i < sequences[k].length(); i++){
		for(int j = i - (minLoopLength + 2 * minHelixLength - 1); j > 0; j--){
		    int score = basepairScore(sequences[k].charAt(i),sequences[k].charAt(j));
		    helixTable[k][i][j] = score < 0 ? 0 : helixTable[k][i-1][j+1] + score;
		    if(helixTable[k][i-1][j+1] < 12){
			helixTable[k][i-1][j+1] = 0;
		    }
		    // we treat the the end of a helix as much better
		    if(helixTable[k][i][j] == 0){
			helixTable[k][i-1][j+1] = helixTable[k][i-1][j+1] * helixTable[k][i-1][j+1];
		    }
		}
	    }
	    // setting the last row
	    for(int j = 1; j < sequences[k].length(); j++){
		if(helixTable[k][sequences[k].length() - 1][j] < 12){
		    helixTable[k][sequences[k].length() - 1][j] = 0;
		}
		helixTable[k][sequences[k].length() - 1][j] = helixTable[k][sequences[k].length() - 1][j] * helixTable[k][sequences[k].length() - 1][j];
	
	    }
	}

	// calculate helixStatistics

	for(int k = 0; k < sequences.length; k++){
	    for(int j = 0; j < sequences[k].length(); j++){
		helixStatistics[k][j] = 0;
		for(int i = j; i < sequences[k].length(); i++){
		    helixStatistics[k][j] += helixTable[k][i][j] * helixTable[k][i][j];
		}
	    }
	}	

	// Printing
	/*
	System.out.println("We have "+sequences.length+" number of sequences\n\n");
	
	for(int k = 0; k < sequences.length; k++){
	    for(int i = 0; i < sequences[k].length(); i++){
		for(int j = 0; j <= i; j++){
		    System.out.print(helixTable[k][i][j]+"\t");
		}
		System.out.println();
	    }
	    System.out.println();
	    for(int j = 0; j < sequences[k].length(); j++){
		System.out.print(helixStatistics[k][j]+"\t");
	    }
	    System.out.println();
	    for(int j = 0; j < sequences[k].length(); j++){
		System.out.print(sequences[k].charAt(j)+"\t");
	    }
	    System.out.println();   
	    System.out.println();   
	}	
	*/
    }

    public static double[][][] dp;

    public static void dpprint(){
	for(int k = 0; k < dp[0][0].length; k++){
	    for(int i = 0; i < dp.length; i++){
		for(int j = 0; j < dp[0].length; j++){
		    System.out.print(dp[i][j][k]+" ");
		}
		System.out.println();
	    }
	    System.out.println();
	    System.out.println();
	}
    }
    
    // tells the number of digits in the edgelengths when printing the tree in Newick format
    public static int rounding = 7;

    //Poisson parameter
    public static final double lambdaForHelixRemove = 3.0;

    // gives a random number, how many helices to be removed
    public static int getNumberofHelicesToBeRemoved(int numberofHelices){
	int i = 0;
	double sum = Math.exp(0.0 - lambdaForHelixRemove);
	double factor = 1.0;
	Random generator = new Random();
	double probability = generator.nextDouble();

	// draw a random number from the truncated Poisson distribution, the last possible 
	// value gets the probability of the entire tail of the distribution.
	while(sum < probability && i < numberofHelices){
	    i++;
	    factor *= lambdaForHelixRemove / (double)i;
	    sum += Math.exp(0.0 - lambdaForHelixRemove) * factor;
	}

	return i;
    }

    // gives the probability of removing a given number of helices from the maximum number
    public static double getProbabilityForRemovedHelices(int numberofRemovedHelices, int maxNumberofHelices){

	// if we removed the maximum number of helices, then the probability is the tail of the Poisson
	// distribution, we calculat the sum fo the other part, the return value is 1.0 - sum.
	if(numberofRemovedHelices == maxNumberofHelices){
	    double sum = 0.0;
	    double factor = 1.0;
	    for(int i = 0; i < numberofRemovedHelices; i++){
		sum += Math.exp(0.0 - lambdaForHelixRemove) * factor;
		factor *= lambdaForHelixRemove / (double) (i + 1);
	    }
	    return 1.0 - sum;

	}
	// otherwise return with the appropriate value from the Poisson distribution
	else{
	    double factor = 1.0;
	    for(int i = 0; i < numberofRemovedHelices; i++){
		factor *= lambdaForHelixRemove / (double) (i + 1);
	    }
	    return Math.exp(0.0 - lambdaForHelixRemove) * factor;
	}
    }

    public static boolean sampleTree = false;
    public static boolean sampleAlignment = false;
    public static boolean sampleStructure = false;

    public static void main(String[] args) throws IOException { // this is the simulfold main
	
	for(int i = 0; i < args.length; i++){
	    if(args[i].equals("-A")){
		sampleAlignment = true;
	    }
	    if(args[i].equals("-S")){
		sampleStructure = true;
	    }
	    if(args[i].equals("-T")){
		sampleTree = true;
	    }
	}

	/*
	for(int i = 0; i < 4; i++){
	    for(int j = 0; j < 4; j++){
		System.out.print((ePiDouble[4* i + j]/(ePiSingle[i]*ePiSingle[j]))+" ");
	    }
	    System.out.println();
	}
	*/
	MCMC mcmc = new MCMC(1,1.0);
	for(int i = 0; i < mcmc.alignment.length; i++){
	    sequences = mcmc.alignment[i].readAlignment(args[0]);
	    mcmc.tree[i] = new Tree(args[1]);
	    mcmc.alignment[i].findLeaves(mcmc.tree[i]);
	}

	// calculating the scores for the sequences
	setHelixStatistics();

	mcmc.doMCMC(10000,2000,100);
    }
    
}
