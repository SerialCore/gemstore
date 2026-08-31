typedef struct{
	double mn;
	double ms;
	double mc;
	double mb;

	double alpha[3];
	double gamma[3];

	double b;
	double c;
	double sigma0;
	double s;
	double f;
	
	double econt;
	double etens;
	double esov;
	double esos;
	double eCoul;
}vargs;

typedef double (*vCenPart)(double,basis_qnum,basis_qnum,vargs,int);

double getsigmaij(double sigma0,double s,double mi,double mj)
{
	return sqrt(sigma0 * sigma0 * (0.5 + 0.5 * pow(4.0 * mi * mj / ((mi + mj) * (mi + mj)), 4.0))+ s * s * pow(2.0 * mi * mj / (mi + mj), 2.0));
}

void getsigmakij(double *sigmakij, double *gammak,double sigma0,double s,double mi,double mj)
{
	int k;
	double sigmaij2=pow(getsigmaij(sigma0,s,mi,mj),2);
	for(k=0;k<3;k++)
	{
		sigmakij[k]=1/sqrt(1/pow(gammak[k],2)+1/sigmaij2);
	}
}

double Ti(double p,basis_qnum qf,basis_qnum qi,vargs varg,int c)
{
	double mc;
	switch(c)
	{
		case 1:
			mc=qf.m1;
			break;
		case 2:
			mc=qf.m2;
			break;
		case 3:
			mc=qf.m3;
			break;
		default:
			printf("error_Ti\n");
			mc=0;
			break;
	}
	return  sqrt(p*p+mc*mc);
}


double InteCenV(double b11,int n,basis_qnum qf,basis_qnum qi,vCenPart v,vargs varg,int c)
{
	#define InteLength 50
	
	const double InteWeight[InteLength]={
		 9.49074370469973154861E-003, 2.20248042432282550986E-002, 3.43744122779471232838E-002, 4.62954423315795426704E-002, 5.74270673524102611643E-002,
		 6.72616743347232737151E-002, 7.51655105427606561651E-002, 8.04487210189303589897E-002, 8.24866653780293862703E-002, 8.08772841479732575406E-002,
		 7.55968131214373958980E-002, 6.70988969096135470989E-002, 5.63037514101928791197E-002, 4.44533647270852744933E-002, 3.28602655529441193691E-002,
		 2.26281047022441036531E-002, 1.44421205638529698147E-002, 8.49988917441698033564E-003, 4.58981582947691516137E-003, 2.26248032155746286264E-003,
		 1.01296345857313010366E-003, 4.09855051025779821045E-004, 1.49103396554347145977E-004, 4.85209725035065972790E-005, 1.40499187954310526217E-005,
		 3.60057445584128798598E-006, 8.12054885660543383779E-007, 1.60236226323028537512E-007, 2.74914498210976715853E-008, 4.07393870688591939297E-009,
		 5.17733909405469218567E-010, 5.59883137228935952164E-011, 5.10833440280922385331E-012, 3.89531180859741328527E-013, 2.45631708318249480239E-014,
		 1.26561330926208323536E-015, 5.25582653443700296954E-017, 1.73149373479240142712E-018, 4.44198243358375823341E-020, 8.68013859892778184465E-022,
		 1.25805074779480909285E-023, 1.30871663315222277155E-025, 9.37690669168311366587E-028, 4.38351942427205493247E-030, 1.25991102592426549802E-032,
		 1.45161317869704793270E-035, 1.44829917427803858346E-037,-2.07564613115271761842E-039, 1.59488572462612746697E-041,-3.26209228139812174140E-044
	};
	const double InteNode[InteLength]={
		 3.69941668941189387078E-003, 1.94655785593001262377E-002, 4.77232575110775569884E-002, 8.82994741423326136810E-002, 1.40944249831512837831E-001,
		 2.05343824732975154536E-001, 2.81130944731155938791E-001, 3.67895701867100816092E-001, 4.65196632382567480466E-001, 5.72571575839739193266E-001,
		 6.89547864293400636314E-001, 8.15651524129035456870E-001, 9.50415293767051842487E-001, 1.09338537085220496150E+000, 1.24412689394047714663E+000,
		 1.40222823252836353434E+000, 1.56730420579590844586E+000, 1.73899837732924138424E+000, 1.91698458425367595630E+000, 2.10096785886658098596E+000,
		 2.29068489295538163123E+000, 2.48590418287233886279E+000, 2.68642597976929227008E+000, 2.89208215616667258410E+000, 3.10273608869143607692E+000,
		 3.31828264842635672287E+000, 3.53864838570112628139E+000, 3.76379199610129929525E+000, 3.99370515987453091261E+000, 4.22841385900133938902E+000,
		 4.46798029678176976527E+000, 4.71250557663335036490E+000, 4.96213334417320913961E+000, 5.21705466625021294844E+000, 5.47751452299329293275E+000,
		 5.74382044124152601111E+000, 6.01635402812180604868E+000, 6.29558651982983683470E+000, 6.58210002638651010927E+000, 6.87661707965657488199E+000,
		 7.18004266516692291795E+000, 7.49352570514240963406E+000, 7.81855215039574067889E+000, 8.15709210450168878089E+000, 8.51184526474371849140E+000,
		 8.88668005924412897894E+000, 9.28749674141648604654E+000, 9.72416586588463146083E+000, 1.02158862585784281522E+001, 1.08129860729453608573E+001
	};
	double r,res;
	int i;
	res=0;
	for(i=0;i<InteLength;i++)
	{
		r=InteNode[i]/sqrt(b11);
		res=res+InteWeight[i]*v(r,qf,qi,varg,c)*pow(r,n);
	}
	res=res/sqrt(b11);
	return res;
}

double Nnl(int l,double nu)
{
	const double Factorial2[300]={
		 1.00000000000000000000E+000, 1.00000000000000000000E+000, 2.00000000000000000000E+000, 3.00000000000000000000E+000, 8.00000000000000000000E+000,
		 1.50000000000000000000E+001, 4.80000000000000000000E+001, 1.05000000000000000000E+002, 3.84000000000000000000E+002, 9.45000000000000000000E+002,
		 3.84000000000000000000E+003, 1.03950000000000000000E+004, 4.60800000000000000000E+004, 1.35135000000000000000E+005, 6.45120000000000000000E+005,
		 2.02702500000000000000E+006, 1.03219200000000000000E+007, 3.44594250000000000000E+007, 1.85794560000000000000E+008, 6.54729075000000000000E+008,
		 3.71589120000000000000E+009, 1.37493105750000000000E+010, 8.17496064000000000000E+010, 3.16234143225000000000E+011, 1.96199055360000000000E+012,
		 7.90585358062500000000E+012, 5.10117543936000000000E+013, 2.13458046676875000000E+014, 1.42832912302080000000E+015, 6.19028335362937500000E+015,
		 4.28498736906240000000E+016, 1.91898783962510625000E+017, 1.37119595809996800000E+018, 6.33265987076285062500E+018, 4.66206625753989120000E+019,
		 2.21643095476699771875E+020, 1.67834385271436083200E+021, 8.20079453263789155937E+021, 6.37770664031457116160E+022, 3.19830986772877770815E+023,
		 2.55108265612582846464E+024, 1.31130704576879886034E+025, 1.07145471557284795514E+026, 5.63862029680583509947E+026, 4.71440074852053100265E+027,
		 2.53737913356262579476E+028, 2.16862434431944426122E+029, 1.19256819277443412353E+030, 1.04093968527333324538E+031, 5.84358414459472720534E+031,
		 5.20469842636666622693E+032, 2.98022791374331087472E+033, 2.70644318171066643800E+034, 1.57952079428395476360E+035, 1.46147931812375987652E+036,
		 8.68736436856175119982E+036, 8.18428418149305530852E+037, 4.95179769008019818390E+038, 4.74688482526597207894E+039, 2.92156063714731692850E+040,
		 2.84813089515958324736E+041, 1.78215198865986332638E+042, 1.76584115499894161336E+043, 1.12275575285571389562E+044, 1.13013833919932263255E+045,
		 7.29791239356214032155E+045, 7.45891303871552937486E+046, 4.88960130368663401543E+047, 5.07206086632655997490E+048, 3.37382489954377747065E+049,
		 3.55044260642859198243E+050, 2.39541567867608200416E+051, 2.55631867662858622735E+052, 1.74865344543353986303E+053, 1.89167582070515380824E+054,
		 1.31149008407515489727E+055, 1.43767362373591689426E+056, 1.00984736473786927090E+057, 1.12138542651401517752E+058, 7.97779418142916724015E+058,
		 8.97108341211212142020E+059, 6.46201328695762546452E+060, 7.35628839793193956456E+061, 5.36347102817482913555E+062, 6.17928225426282923423E+063,
		 4.55895037394860476522E+064, 5.31418273866603314144E+065, 3.96628682533528614574E+066, 4.67648081002610916447E+067, 3.52999527454840466971E+068,
		 4.20883272902349824802E+069, 3.21229569983904824943E+070, 3.87212611070161838818E+071, 2.98743500085031487197E+072, 3.63979854405952128489E+073,
		 2.83806325080779912837E+074, 3.49420660229714043349E+075, 2.75292135328356515452E+076, 3.42432247025119762482E+077, 2.72539213975072950298E+078,
		 3.42432247025119762482E+079, 2.75264606114823679801E+080, 3.49280891965622157732E+081, 2.83522544298268390195E+082, 3.63252127644247044041E+083,
		 2.97698671513181809704E+084, 3.85047255302901866683E+085, 3.18537578519104536384E+086, 4.15851035727134016018E+087, 3.47205960585823944658E+088,
		 4.57436139299847417620E+089, 3.85398616250264578571E+090, 5.12328476015829107734E+091, 4.35500436362798973785E+092, 5.84054462658045182817E+093,
		 5.00825501817218819853E+094, 6.77503176683332412068E+095, 5.85965837126146019228E+096, 7.99453748486332246241E+097, 6.97299346180113762881E+098,
		 9.59344498183598695489E+099, 8.43732208877937653086E+100, 1.17040028778399040849E+102, 1.03779061691986331329E+103, 1.45129635685214810653E+104,
		 1.29723827114982914162E+105, 1.82863340963370661423E+106, 1.64749260436028300985E+107, 2.34065076433114446622E+108, 2.12526545962476508271E+109,
		 3.04284599363048780608E+110, 2.78409775210844225836E+111, 4.01655671159224390403E+112, 3.70285001030422820361E+113, 5.38218599353360683140E+114,
		 4.99884751391070807488E+115, 7.31977295120570529071E+116, 6.84842109405767006259E+117, 1.01012866726638733011E+119, 9.51930532074016138700E+119,
		 1.41418013417294226216E+121, 1.34222205022436275556E+122, 2.00813579052557801227E+123, 1.91937753182083874046E+124, 2.89171553835683233767E+125,
		 2.78309742114021617367E+126, 4.22190468600097521300E+127, 4.09115320907611777529E+128, 6.24841893528144331525E+129, 6.09581828152341548518E+130,
		 9.37262840292216497288E+131, 9.20468560510035738263E+132, 1.42463951724416907587E+134, 1.40831689758035467954E+135, 2.19394485655602037685E+136,
		 2.18289119124954975329E+137, 3.42255397622739178788E+138, 3.42713917026179311266E+139, 5.40763528243927902486E+140, 5.44915128071625104914E+141,
		 8.65221645190284643978E+142, 8.77313356195316418912E+143, 1.40165906520826112324E+145, 1.43002077059836576282E+146, 2.29872086694154824212E+147,
		 2.35953427148730350866E+148, 3.81587663912297008192E+149, 3.94042223338379685946E+150, 6.41067275372658973762E+151, 6.65931357441861669250E+152,
		 1.08981436813352025539E+154, 1.13874262122558345441E+155, 1.87448071318965483928E+156, 1.97002473472025937614E+157, 3.26159644094999942035E+158,
		 3.44754328576045390825E+159, 5.74040973607199897981E+160, 6.10215161579600341760E+161, 1.02179293302081581840E+163, 1.09228513922748461175E+164,
		 1.83922727943746847313E+165, 1.97703610200174714726E+166, 3.34739364857619262110E+167, 3.61797606666319727950E+168, 6.15920431338019442283E+169,
		 6.69325572332691496707E+170, 1.14561200228871616264E+172, 1.25163882026213309884E+173, 2.15375056430278638577E+174, 2.36559737029543155681E+175,
		 4.09212607217529413297E+176, 4.51829097726427427351E+177, 7.85688205857656473530E+178, 8.72030158612004934788E+179, 1.52423511936385355865E+181,
		 1.70045880929340962283E+182, 2.98750083395315297495E+183, 3.34990385430801695699E+184, 5.91525165122724289040E+185, 6.66630867007295374441E+186,
		 1.18305033024544857808E+188, 1.33992804268466370262E+189, 2.38976166709580612772E+190, 2.72005392664986731633E+191, 4.87511380087544450055E+192,
		 5.57611054963222799848E+193, 1.00427344298034156711E+195, 1.15425488377387119568E+196, 2.08888876139911045959E+197, 2.41239270708739079898E+198,
		 4.38666639893813196515E+199, 5.09014861195439458585E+200, 9.29973276574883976613E+201, 1.08420165434628604678E+203, 1.99014281187025170995E+204,
		 2.33103355684451500059E+205, 4.29870847363974369349E+206, 5.05834281835259755128E+207, 9.37118447253464125182E+208, 1.10777707721921886373E+210,
		 2.06166058395762107540E+211, 2.44818734065447368884E+212, 4.57688649638591878739E+213, 5.45945776965947632612E+214, 1.02522257519044580837E+216,
		 1.22837799817338217337E+217, 2.31700301993040752693E+218, 2.78841805585357753356E+219, 5.28276688544132916140E+220, 6.38547734790469255187E+221,
		 1.21503638365150570712E+223, 1.47504526736598397948E+224, 2.81888441007149324052E+225, 3.43685547296274267219E+226, 6.59618951956729418282E+227,
		 8.07661036146244527965E+228, 1.55670072661788142714E+230, 1.91415665566659953127E+231, 3.70494772935055779660E+232, 4.57483440704317287975E+233,
		 8.89187455044133871186E+234, 1.10253509209740466402E+236, 2.15183364120680396827E+237, 2.67916027379669333357E+238, 5.25047408454460168257E+239,
		 6.56394267080189866725E+240, 1.29161662479797201391E+242, 1.62129383968806897081E+243, 3.20320922949897059450E+244, 4.03702166082329173731E+245,
		 8.00802307374742648627E+246, 1.01329243686664622606E+248, 2.01802181458435147454E+249, 2.56362986527261495194E+250, 5.12577540904425274533E+251,
		 6.53725615644516812747E+252, 1.31219850471532870280E+254, 1.68007483220640820876E+255, 3.38547214216554805323E+256, 4.35139381541459726068E+257,
		 8.80222756963042493841E+258, 1.13571378582320988503E+260, 2.30618362324317133386E+261, 2.98692725671504199765E+262, 6.08832476536197232140E+263,
		 7.91535723029486129378E+264, 1.61949438758628463749E+266, 2.11340038048872796544E+267, 4.34024495873124282848E+268, 5.68504702351467822703E+269,
		 1.17186613885743556369E+271, 1.54064774337247779952E+272, 3.18747589769222473323E+273, 4.20596833940686439270E+274, 8.73368395967669576907E+275,
		 1.15664129333688770799E+277, 2.41049677287076803226E+278, 3.20389638254317895114E+279, 6.70118102858073512969E+280, 8.93887090729546927369E+281,
		 1.87633068800260583631E+283, 2.51182272495002686590E+284, 5.29125254016734845840E+285, 7.10845831160857603052E+286, 1.50271572140752696218E+288,
		 2.02591061880844416869E+289, 4.29776696322552711185E+290, 5.81436347598023476416E+291, 1.23775688540895180821E+293, 1.68035104455828784684E+294,
		 3.58949496768596024382E+295, 4.88982153966461763431E+296, 1.04813253056430039119E+298, 1.43271771112173296685E+299, 3.08150963985904315011E+300,
		 4.22651724780911225221E+301, 9.12126853398276772434E+302, 1.25527562259930633890E+304, 2.71813802312686478185E+305, 3.75327411157192595333E+306
	};

	double pi=3.14159265358979324;
        return sqrt(pow(2.,l+2.)*pow(2*nu,l+1.5)/(sqrt(pi)*Factorial2[2*l+1]));
}


void coordinatesTransformation_r_1(double m1,double m2,double m3,int a,int c,double *alpha_ac,double *beta_ac,double *gamma_ac,double *delta_ac)
{       
        double alpha[3][3]={{1,-(m3/(m1+m3)),-(m3/(m2+m3))},{-(m2/(m1+m2)),1,-(m2/(m2+m3))},{-(m1/(m1+m2)),-(m1/(m1+m3)),1}};
        double  beta[3][3]={{0,1,-1},{-1,0,1},{1,-1,0}};
        double gamma[3][3]={{0,-1+(m2*m3)/((m1+m2)*(m1+m3)),1-(m1*m3)/((m1+m2)*(m2+m3))},{1-(m2*m3)/((m1+m2)*(m1+m3)),0,-1+(m1*m2)/((m1+m3)*(m2+m3))},{-1+(m1*m3)/((m1+m2)*(m2+m3)),1-(m1*m2)/((m1+m3)*(m2+m3)),0}};
        double delta[3][3]={{1,-(m2/(m1+m2)),-(m1/(m1+m2))},{-(m3/(m1+m3)),1,-(m1/(m1+m3))},{-(m3/(m2+m3)),-(m2/(m2+m3)),1}};
        *alpha_ac=alpha[a-1][c-1];
         *beta_ac= beta[a-1][c-1];
        *gamma_ac=gamma[a-1][c-1];
        *delta_ac=delta[a-1][c-1];
}                

void coordinatesTransformation_r_2(double m1,double m2,double m3,int a,int c,double *alpha_ac,double *beta_ac,double *gamma_ac,double *delta_ac)
{       
	 double alpha[3][3]={{0,1,-1},{-1,0,1},{1,-1,0}};
	 double  beta[3][3]={{1,-(m3/(m1+m3)),-(m3/(m2+m3))},{-(m2/(m1+m2)),1,-(m2/(m2+m3))},{-(m1/(m1+m2)),-(m1/(m1+m3)),1}};
	 double gamma[3][3]={{1,-(m2/(m1+m2)),-(m1/(m1+m2))},{-(m3/(m1+m3)),1,-(m1/(m1+m3))},{-(m3/(m2+m3)),-(m2/(m2+m3)),1}};
	 double delta[3][3]={{0,-1+(m2*m3)/((m1+m2)*(m1+m3)),1-(m1*m3)/((m1+m2)*(m2+m3))},{1-(m2*m3)/((m1+m2)*(m1+m3)),0,-1+(m1*m2)/((m1+m3)*(m2+m3))},{-1+(m1*m3)/((m1+m2)*(m2+m3)),1-(m1*m2)/((m1+m3)*(m2+m3)),0}};
        *alpha_ac=alpha[a-1][c-1];
         *beta_ac= beta[a-1][c-1];
        *gamma_ac=gamma[a-1][c-1];
        *delta_ac=delta[a-1][c-1];
}                


void coordinatesTransformation_p_1(double m1,double m2,double m3,int a,int c,double *alpha_ac,double *beta_ac,double *gamma_ac,double *delta_ac)
{       
	double alpha[3][3]={{1,-(m2/(m1+m2)),-(m1/(m1+m2))},{-(m3/(m1+m3)),1,-(m1/(m1+m3))},{-(m3/(m2+m3)),-(m2/(m2+m3)),1}};
	double  beta[3][3]={{0,(m1*(m1+m2+m3))/((m1+m2)*(m1+m3)),-((m2*(m1+m2+m3))/((m1+m2)*(m2+m3)))},{-((m1*(m1+m2+m3))/((m1+m2)*(m1+m3))),0,(m3*(m1+m2+m3))/((m1+m3)*(m2+m3))},{(m2*(m1+m2+m3))/((m1+m2)*(m2+m3)),-((m3*(m1+m2+m3))/((m1+m3)*(m2+m3))),0}};
	double gamma[3][3]={{0,-1,1},{1,0,-1},{-1,1,0}};
	double delta[3][3]={{1,-(m3/(m1+m3)),-(m3/(m2+m3))},{-(m2/(m1+m2)),1,-(m2/(m2+m3))},{-(m1/(m1+m2)),-(m1/(m1+m3)),1}};
        *alpha_ac=alpha[a-1][c-1];
         *beta_ac= beta[a-1][c-1];
        *gamma_ac=gamma[a-1][c-1];
        *delta_ac=delta[a-1][c-1];
}                

void coordinatesTransformation_p_2(double m1,double m2,double m3,int a,int c,double *alpha_ac,double *beta_ac,double *gamma_ac,double *delta_ac)
{       
	 double alpha[3][3]={{0,(m1*(m1+m2+m3))/((m1+m2)*(m1+m3)),-((m2*(m1+m2+m3))/((m1+m2)*(m2+m3)))},{-((m1*(m1+m2+m3))/((m1+m2)*(m1+m3))),0,(m3*(m1+m2+m3))/((m1+m3)*(m2+m3))},{(m2*(m1+m2+m3))/((m1+m2)*(m2+m3)),-((m3*(m1+m2+m3))/((m1+m3)*(m2+m3))),0}};
	 double  beta[3][3]={{1,-(m2/(m1+m2)),-(m1/(m1+m2))},{-(m3/(m1+m3)),1,-(m1/(m1+m3))},{-(m3/(m2+m3)),-(m2/(m2+m3)),1}};
	 double gamma[3][3]={{1,-(m3/(m1+m3)),-(m3/(m2+m3))},{-(m2/(m1+m2)),1,-(m2/(m2+m3))},{-(m1/(m1+m2)),-(m1/(m1+m3)),1}};
	 double delta[3][3]={{0,-1,1},{1,0,-1},{-1,1,0}};
        *alpha_ac=alpha[a-1][c-1];
         *beta_ac= beta[a-1][c-1];
        *gamma_ac=gamma[a-1][c-1];
        *delta_ac=delta[a-1][c-1];
}                

void coordinatesTransformation_p_i(double m1,double m2,double m3,int a,int c,double *alpha_ac,double *beta_ac,double *gamma_ac,double *delta_ac)
{       
	double alpha[3][3]={{m2/(m1+m2),-1,m1/(m1+m2)},{-1,m3/(m1+m3),m1/(m1+m3)},{m2/(m2+m3),m3/(m2+m3),-1}};
	double  beta[3][3]={{-(m1/(m1+m2)),-(m2/(m1+m2)),m1/(m1+m2)+m2/(m1+m2)},{-(m1/(m1+m3)),m1/(m1+m3)+m3/(m1+m3),-(m3/(m1+m3))},{m2/(m2+m3)+m3/(m2+m3),-(m2/(m2+m3)),-(m3/(m2+m3))}};
	double gamma[3][3]={{1,0,-1},{0,-1,1},{-1,1,0}};
	double delta[3][3]={{1,-1,0},{-1,0,1},{0,1,-1}};
        *alpha_ac=alpha[a-1][c-1];
         *beta_ac= beta[a-1][c-1];
        *gamma_ac=gamma[a-1][c-1];
        *delta_ac=delta[a-1][c-1];
}                





void getT1r(basis_qnum qf, basis_qnum qi, double* txrp_cent, int c, int vt)
{
	int a = qf.c;
	int b = qi.c;

	int l1 = qf.lrho;
	int l2 = qf.llam;
	int l3 = qi.lrho;
	int l4 = qi.llam;

	double nu1 = qf.nurho;
	double nu2 = qf.nulam;
	double nu3 = qi.nurho;
	double nu4 = qi.nulam;

	double m1 = qf.m1;
	double m2 = qf.m2;
	double m3 = qf.m3;

	double alpha_ac, beta_ac, gamma_ac, delta_ac, alpha_bc, beta_bc, gamma_bc, delta_bc;

	double mi, mj, mk;
	double afrr, airr, afrR, airR, afRR, aiRR;
	double arr, arR, aRR;
	double aDfr[4], aDir[4], aDfR[4], aDiR[4];
	double aDr[4], aDR[4];
	double lrr;
	double lDr[4];
	double coef, coeg;
	double f1, f2, f3, f4;
	double g1, g2, g3, g4;
	double h1, h2, h3, h4;
	double alpha_i, alpha_j, beta_i, beta_j;
	double hi11, hi12, hi13, hi14;
	double hj11, hj12, hj13, hj14;
	double h21, h22, h23, h24;
	double coec, coet, coeso;
	double hrho11, hrho12, hrho13, hrho14;
	double hlambda11, hlambda12, hlambda13, hlambda14;
	double pi = 3.14159265358979324;
	int i;

	coordinatesTransformation_r_1(m1, m2, m3, a, c, &alpha_ac, &beta_ac, &gamma_ac, &delta_ac);
	coordinatesTransformation_r_1(m1, m2, m3, b, c, &alpha_bc, &beta_bc, &gamma_bc, &delta_bc);

	getijk(m1, m2, m3, &mi, &mj, &mk, c);

	afrr = (nu1 * alpha_ac * alpha_ac + nu2 * gamma_ac * gamma_ac);
	airr = (nu3 * alpha_bc * alpha_bc + nu4 * gamma_bc * gamma_bc);
	afrR = 2 * (nu1 * alpha_ac * beta_ac + nu2 * gamma_ac * delta_ac);
	airR = 2 * (nu3 * alpha_bc * beta_bc + nu4 * gamma_bc * delta_bc);
	afRR = (nu1 * beta_ac * beta_ac + nu2 * delta_ac * delta_ac);
	aiRR = (nu3 * beta_bc * beta_bc + nu4 * delta_bc * delta_bc);

	arr = afrr + airr;
	arR = afrR + airR;
	aRR = afRR + aiRR;

	aDfr[0] = 2 * alpha_ac;
	aDfr[1] = 2 * gamma_ac;
	aDfr[2] = 0;
	aDfr[3] = 0;
	aDir[0] = 0;
	aDir[1] = 0;
	aDir[2] = 2 * alpha_bc;
	aDir[3] = 2 * gamma_bc;
	aDfR[0] = 2 * beta_ac;
	aDfR[1] = 2 * delta_ac;
	aDfR[2] = 0;
	aDfR[3] = 0;
	aDiR[0] = 0;
	aDiR[1] = 0;
	aDiR[2] = 2 * beta_bc;
	aDiR[3] = 2 * delta_bc;

	for(i = 0; i < 4; i++)
	{
		aDr[i] = aDfr[i] + aDir[i];
		aDR[i] = aDfR[i] + aDiR[i];
	}

	lrr = arr - (arR * arR) / (4 * aRR);

	for(i = 0; i < 4; i++)
	{
		lDr[i] = aDr[i] - (arR / (2 * aRR)) * aDR[i];
	}

	coef = 1;
	coeg = 1 / (4 * aRR);

	f1 = lDr[0];
	f2 = lDr[1];
	f3 = lDr[2];
	f4 = lDr[3];
	g1 = aDR[0];
	g2 = aDR[1];
	g3 = aDR[2];
	g4 = aDR[3];
	h1 = lDr[0];
	h2 = lDr[1];
	h3 = lDr[2];
	h4 = lDr[3];

	alpha_i = 1;
	alpha_j = -1;
	beta_i = mi / (mi + mj);
	beta_j = mj / (mi + mj);

	hi11 = (alpha_i * aDir[0] + beta_i * aDiR[0]) - ((alpha_i * airR + 2 * beta_i * aiRR) / (2 * aRR)) * aDR[0];
	hi12 = (alpha_i * aDir[1] + beta_i * aDiR[1]) - ((alpha_i * airR + 2 * beta_i * aiRR) / (2 * aRR)) * aDR[1];
	hi13 = (alpha_i * aDir[2] + beta_i * aDiR[2]) - ((alpha_i * airR + 2 * beta_i * aiRR) / (2 * aRR)) * aDR[2];
	hi14 = (alpha_i * aDir[3] + beta_i * aDiR[3]) - ((alpha_i * airR + 2 * beta_i * aiRR) / (2 * aRR)) * aDR[3];

	hj11 = (alpha_j * aDir[0] + beta_j * aDiR[0]) - ((alpha_j * airR + 2 * beta_j * aiRR) / (2 * aRR)) * aDR[0];
	hj12 = (alpha_j * aDir[1] + beta_j * aDiR[1]) - ((alpha_j * airR + 2 * beta_j * aiRR) / (2 * aRR)) * aDR[1];
	hj13 = (alpha_j * aDir[2] + beta_j * aDiR[2]) - ((alpha_j * airR + 2 * beta_j * aiRR) / (2 * aRR)) * aDR[2];
	hj14 = (alpha_j * aDir[3] + beta_j * aDiR[3]) - ((alpha_j * airR + 2 * beta_j * aiRR) / (2 * aRR)) * aDR[3];

	h21 = lDr[0];
	h22 = lDr[1];
	h23 = lDr[2];
	h24 = lDr[3];

	coec = 4 * pi * pow(sqrt(pi / aRR), 3);
	coet = 16 * pi * sqrt(24 * pi) * pow(sqrt(pi / aRR), 3);
	coeso = (4 * pi) / aRR * sqrt(pi / (4 * aRR)) * (-sqrt(3)) * (16 * sqrt(2) * pi * pi) / 3;

	alpha_i = 1;
	alpha_j = 0;
	beta_i = 0;
	beta_j = 1;

	hrho11 = (alpha_i * aDir[0] + beta_i * aDiR[0]) - ((alpha_i * airR + 2 * beta_i * aiRR) / (2 * aRR)) * aDR[0];
	hrho12 = (alpha_i * aDir[1] + beta_i * aDiR[1]) - ((alpha_i * airR + 2 * beta_i * aiRR) / (2 * aRR)) * aDR[1];
	hrho13 = (alpha_i * aDir[2] + beta_i * aDiR[2]) - ((alpha_i * airR + 2 * beta_i * aiRR) / (2 * aRR)) * aDR[2];
	hrho14 = (alpha_i * aDir[3] + beta_i * aDiR[3]) - ((alpha_i * airR + 2 * beta_i * aiRR) / (2 * aRR)) * aDR[3];

	hlambda11 = (alpha_j * aDir[0] + beta_j * aDiR[0]) - ((alpha_j * airR + 2 * beta_j * aiRR) / (2 * aRR)) * aDR[0];
	hlambda12 = (alpha_j * aDir[1] + beta_j * aDiR[1]) - ((alpha_j * airR + 2 * beta_j * aiRR) / (2 * aRR)) * aDR[1];
	hlambda13 = (alpha_j * aDir[2] + beta_j * aDiR[2]) - ((alpha_j * airR + 2 * beta_j * aiRR) / (2 * aRR)) * aDR[2];
	hlambda14 = (alpha_j * aDir[3] + beta_j * aDiR[3]) - ((alpha_j * airR + 2 * beta_j * aiRR) / (2 * aRR)) * aDR[3];

	txrp_cent[ 0] = coef;
	txrp_cent[ 1] = f1;
	txrp_cent[ 2] = f2;
	txrp_cent[ 3] = f3;
	txrp_cent[ 4] = f4;
	txrp_cent[ 5] = coeg;
	txrp_cent[ 6] = g1;
	txrp_cent[ 7] = g2;
	txrp_cent[ 8] = g3;
	txrp_cent[ 9] = g4;

	txrp_cent[18] = lrr;
	txrp_cent[20] = Nnl(l1, nu1) * Nnl(l2, nu2) * Nnl(l3, nu3) * Nnl(l4, nu4);
	


	switch(vt)
	{
		case 1:
			txrp_cent[19] = coec;
			break;

		case 2:
			txrp_cent[10] = h1;
			txrp_cent[11] = h2;
			txrp_cent[12] = h3;
			txrp_cent[13] = h4;
			txrp_cent[19] = coet;
		break;

		case 3:
			txrp_cent[10] = hi11;
			txrp_cent[11] = hi12;
			txrp_cent[12] = hi13;
			txrp_cent[13] = hi14;
			txrp_cent[14] = h21;
			txrp_cent[15] = h22;
			txrp_cent[16] = h23;
			txrp_cent[17] = h24;
			txrp_cent[19] = coeso;
		break;

		case 4:
			txrp_cent[10] = hj11;
			txrp_cent[11] = hj12;
			txrp_cent[12] = hj13;
			txrp_cent[13] = hj14;
			txrp_cent[14] = h21;
			txrp_cent[15] = h22;
			txrp_cent[16] = h23;
			txrp_cent[17] = h24;
			txrp_cent[19] = coeso;
		break;

		case 5:
			txrp_cent[10] = hrho11;
			txrp_cent[11] = hrho12;
			txrp_cent[12] = hrho13;
			txrp_cent[13] = hrho14;
			txrp_cent[14] = h21;
			txrp_cent[15] = h22;
			txrp_cent[16] = h23;
			txrp_cent[17] = h24;
			txrp_cent[19] = coeso;
		break;

		case 6:
			txrp_cent[10] = hlambda11;
			txrp_cent[11] = hlambda12;
			txrp_cent[12] = hlambda13;
			txrp_cent[13] = hlambda14;
			txrp_cent[14] = h21;
			txrp_cent[15] = h22;
			txrp_cent[16] = h23;
			txrp_cent[17] = h24;
			txrp_cent[19] = coeso;
		break;

		default:
		break;
	}
}

void getTir(basis_qnum qf, basis_qnum qi, double* txrp_cent, int c, int vt)
{
	int a = qf.c;
	int b = qi.c;

	int l1 = qf.lrho;
	int l2 = qf.llam;
	int l3 = qi.lrho;
	int l4 = qi.llam;

	double nu1 = qf.nurho;
	double nu2 = qf.nulam;
	double nu3 = qi.nurho;
	double nu4 = qi.nulam;

	double ma1 = qf.m1;
	double ma2 = qf.m2;
	double ma3 = qf.m3;

	double mb1 = qi.m1;
	double mb2 = qi.m2;
	double mb3 = qi.m3;

	double ralpha_ac, rbeta_ac, rgamma_ac, rdelta_ac, ralpha_bc, rbeta_bc, rgamma_bc, rdelta_bc;
	double  alpha_ac,  beta_ac,  gamma_ac,  delta_ac,  alpha_bc,  beta_bc,  gamma_bc,  delta_bc;

	double afrr, airr, afrR, airR, afRR, aiRR;
	double arr, arR, aRR;
	double aDfr[4], aDir[4], aDfR[4], aDiR[4];
	double aDr[4], aDR[4];
	double lrr;
	double lDr[4];
	double coef, coeg;
	double f1, f2, f3, f4;
	double g1, g2, g3, g4;
	double h1, h2, h3, h4;
	double alpha_i, beta_i;
	double h21, h22, h23, h24;
	double coec, coet, coeso;
	double hrho11, hrho12, hrho13, hrho14;
	double pi = 3.14159265358979324;
	int i;

	double mai,maj,mak,mbi,mbj,mbk;
	double alpha_a_i,beta_a_i,alpha_b_i,beta_b_i;

	getijk(ma1, ma2, ma3, &mai, &maj, &mak, c);
	getijk(mb1, mb2, mb3, &mbi, &mbj, &mbk, c);

	alpha_a_i=1; beta_a_i=mai/(mai+maj); alpha_b_i=1; beta_b_i=mbi/(mbi+mbj);

	coordinatesTransformation_r_1(ma1, ma2, ma3, a, c, &ralpha_ac, &rbeta_ac, &rgamma_ac, &rdelta_ac);
	coordinatesTransformation_r_1(mb1, mb2, mb3, b, c, &ralpha_bc, &rbeta_bc, &rgamma_bc, &rdelta_bc);

	alpha_ac=alpha_a_i*ralpha_ac+beta_a_i* rbeta_ac; beta_ac= rbeta_ac;
	gamma_ac=alpha_a_i*rgamma_ac+beta_a_i*rdelta_ac;delta_ac=rdelta_ac;

	alpha_bc=alpha_b_i*ralpha_bc+beta_b_i* rbeta_bc; beta_bc= rbeta_bc;
	gamma_bc=alpha_b_i*rgamma_bc+beta_b_i*rdelta_bc;delta_bc=rdelta_bc;

	afrr = (nu1 * alpha_ac * alpha_ac + nu2 * gamma_ac * gamma_ac);
	airr = (nu3 * alpha_bc * alpha_bc + nu4 * gamma_bc * gamma_bc);
	afrR = 2 * (nu1 * alpha_ac * beta_ac + nu2 * gamma_ac * delta_ac);
	airR = 2 * (nu3 * alpha_bc * beta_bc + nu4 * gamma_bc * delta_bc);
	afRR = (nu1 * beta_ac * beta_ac + nu2 * delta_ac * delta_ac);
	aiRR = (nu3 * beta_bc * beta_bc + nu4 * delta_bc * delta_bc);

	arr = afrr + airr;
	arR = afrR + airR;
	aRR = afRR + aiRR;

	aDfr[0] = 2 * alpha_ac;
	aDfr[1] = 2 * gamma_ac;
	aDfr[2] = 0;
	aDfr[3] = 0;
	aDir[0] = 0;
	aDir[1] = 0;
	aDir[2] = 2 * alpha_bc;
	aDir[3] = 2 * gamma_bc;
	aDfR[0] = 2 * beta_ac;
	aDfR[1] = 2 * delta_ac;
	aDfR[2] = 0;
	aDfR[3] = 0;
	aDiR[0] = 0;
	aDiR[1] = 0;
	aDiR[2] = 2 * beta_bc;
	aDiR[3] = 2 * delta_bc;

	for(i = 0; i < 4; i++)
	{
		aDr[i] = aDfr[i] + aDir[i];
		aDR[i] = aDfR[i] + aDiR[i];
	}

	lrr = arr - (arR * arR) / (4 * aRR);

	for(i = 0; i < 4; i++)
	{
		lDr[i] = aDr[i] - (arR / (2 * aRR)) * aDR[i];
	}

	coef = 1;
	coeg = 1 / (4 * aRR);

	f1 = lDr[0];
	f2 = lDr[1];
	f3 = lDr[2];
	f4 = lDr[3];
	g1 = aDR[0];
	g2 = aDR[1];
	g3 = aDR[2];
	g4 = aDR[3];
	h1 = lDr[0];
	h2 = lDr[1];
	h3 = lDr[2];
	h4 = lDr[3];

	h21 = lDr[0];
	h22 = lDr[1];
	h23 = lDr[2];
	h24 = lDr[3];

	coec = 4 * pi * pow(sqrt(pi / aRR), 3);
	coet = 16 * pi * sqrt(24 * pi) * pow(sqrt(pi / aRR), 3);
	coeso = (4 * pi) / aRR * sqrt(pi / (4 * aRR)) * (-sqrt(3)) * (16 * sqrt(2) * pi * pi) / 3;

	alpha_i = 1;
	beta_i = 0;

	hrho11 = (alpha_i * aDir[0] + beta_i * aDiR[0]) - ((alpha_i * airR + 2 * beta_i * aiRR) / (2 * aRR)) * aDR[0];
	hrho12 = (alpha_i * aDir[1] + beta_i * aDiR[1]) - ((alpha_i * airR + 2 * beta_i * aiRR) / (2 * aRR)) * aDR[1];
	hrho13 = (alpha_i * aDir[2] + beta_i * aDiR[2]) - ((alpha_i * airR + 2 * beta_i * aiRR) / (2 * aRR)) * aDR[2];
	hrho14 = (alpha_i * aDir[3] + beta_i * aDiR[3]) - ((alpha_i * airR + 2 * beta_i * aiRR) / (2 * aRR)) * aDR[3];

	txrp_cent[ 0] = coef;
	txrp_cent[ 1] = f1;
	txrp_cent[ 2] = f2;
	txrp_cent[ 3] = f3;
	txrp_cent[ 4] = f4;
	txrp_cent[ 5] = coeg;
	txrp_cent[ 6] = g1;
	txrp_cent[ 7] = g2;
	txrp_cent[ 8] = g3;
	txrp_cent[ 9] = g4;

	txrp_cent[18] = lrr;
	txrp_cent[20] = Nnl(l1, nu1) * Nnl(l2, nu2) * Nnl(l3, nu3) * Nnl(l4, nu4);
	


	switch(vt)
	{
		case 1:
			txrp_cent[19] = coec;
			break;

		case 2:
			txrp_cent[10] = h1;
			txrp_cent[11] = h2;
			txrp_cent[12] = h3;
			txrp_cent[13] = h4;
			txrp_cent[19] = coet;
		break;
		
		case 5:
			txrp_cent[10] = hrho11;
			txrp_cent[11] = hrho12;
			txrp_cent[12] = hrho13;
			txrp_cent[13] = hrho14;
			txrp_cent[14] = h21;
			txrp_cent[15] = h22;
			txrp_cent[16] = h23;
			txrp_cent[17] = h24;
			txrp_cent[19] = coeso;
		break;
		
		default:
		break;
	}
}

void getTjr(basis_qnum qf, basis_qnum qi, double* txrp_cent, int c, int vt)
{
	int a = qf.c;
	int b = qi.c;

	int l1 = qf.lrho;
	int l2 = qf.llam;
	int l3 = qi.lrho;
	int l4 = qi.llam;

	double nu1 = qf.nurho;
	double nu2 = qf.nulam;
	double nu3 = qi.nurho;
	double nu4 = qi.nulam;

	double ma1 = qf.m1;
	double ma2 = qf.m2;
	double ma3 = qf.m3;

	double mb1 = qi.m1;
	double mb2 = qi.m2;
	double mb3 = qi.m3;

	double ralpha_ac, rbeta_ac, rgamma_ac, rdelta_ac, ralpha_bc, rbeta_bc, rgamma_bc, rdelta_bc;
	double  alpha_ac,  beta_ac,  gamma_ac,  delta_ac,  alpha_bc,  beta_bc,  gamma_bc,  delta_bc;

	double afrr, airr, afrR, airR, afRR, aiRR;
	double arr, arR, aRR;
	double aDfr[4], aDir[4], aDfR[4], aDiR[4];
	double aDr[4], aDR[4];
	double lrr;
	double lDr[4];
	double coef, coeg;
	double f1, f2, f3, f4;
	double g1, g2, g3, g4;
	double h1, h2, h3, h4;
	double alpha_i, beta_i;
	double h21, h22, h23, h24;
	double coec, coet, coeso;
	double hrho11, hrho12, hrho13, hrho14;
	double pi = 3.14159265358979324;
	int i;

	double mai,maj,mak,mbi,mbj,mbk;
	double alpha_a_j,beta_a_j,alpha_b_j,beta_b_j;

	getijk(ma1, ma2, ma3, &mai, &maj, &mak, c);
	getijk(mb1, mb2, mb3, &mbi, &mbj, &mbk, c);

	alpha_a_j=-1; beta_a_j= maj/(mai+maj); alpha_b_j=-1; beta_b_j= mbj/(mbi+mbj);
	//alpha_a_j= 1; beta_a_j=-maj/(mai+maj); alpha_b_j= 1; beta_b_j=-mbj/(mbi+mbj);

	coordinatesTransformation_r_1(ma1, ma2, ma3, a, c, &ralpha_ac, &rbeta_ac, &rgamma_ac, &rdelta_ac);
	coordinatesTransformation_r_1(mb1, mb2, mb3, b, c, &ralpha_bc, &rbeta_bc, &rgamma_bc, &rdelta_bc);

	//alpha_ac=alpha_a_j*ralpha_ac+beta_a_j* rbeta_ac; beta_ac= rbeta_ac;
	//gamma_ac=alpha_a_j*rgamma_ac+beta_a_j*rdelta_ac;delta_ac=rdelta_ac;
	//alpha_bc=alpha_b_j*ralpha_bc+beta_b_j* rbeta_bc; beta_bc= rbeta_bc;
	//gamma_bc=alpha_b_j*rgamma_bc+beta_b_j*rdelta_bc;delta_bc=rdelta_bc;

	alpha_ac=-alpha_a_j*ralpha_ac-beta_a_j* rbeta_ac; beta_ac= rbeta_ac;
	gamma_ac=-alpha_a_j*rgamma_ac-beta_a_j*rdelta_ac;delta_ac=rdelta_ac;
	alpha_bc=-alpha_b_j*ralpha_bc-beta_b_j* rbeta_bc; beta_bc= rbeta_bc;
	gamma_bc=-alpha_b_j*rgamma_bc-beta_b_j*rdelta_bc;delta_bc=rdelta_bc;
	
	afrr = (nu1 * alpha_ac * alpha_ac + nu2 * gamma_ac * gamma_ac);
	airr = (nu3 * alpha_bc * alpha_bc + nu4 * gamma_bc * gamma_bc);
	afrR = 2 * (nu1 * alpha_ac * beta_ac + nu2 * gamma_ac * delta_ac);
	airR = 2 * (nu3 * alpha_bc * beta_bc + nu4 * gamma_bc * delta_bc);
	afRR = (nu1 * beta_ac * beta_ac + nu2 * delta_ac * delta_ac);
	aiRR = (nu3 * beta_bc * beta_bc + nu4 * delta_bc * delta_bc);

	arr = afrr + airr;
	arR = afrR + airR;
	aRR = afRR + aiRR;

	aDfr[0] = 2 * alpha_ac;
	aDfr[1] = 2 * gamma_ac;
	aDfr[2] = 0;
	aDfr[3] = 0;
	aDir[0] = 0;
	aDir[1] = 0;
	aDir[2] = 2 * alpha_bc;
	aDir[3] = 2 * gamma_bc;
	aDfR[0] = 2 * beta_ac;
	aDfR[1] = 2 * delta_ac;
	aDfR[2] = 0;
	aDfR[3] = 0;
	aDiR[0] = 0;
	aDiR[1] = 0;
	aDiR[2] = 2 * beta_bc;
	aDiR[3] = 2 * delta_bc;

	for(i = 0; i < 4; i++)
	{
		aDr[i] = aDfr[i] + aDir[i];
		aDR[i] = aDfR[i] + aDiR[i];
	}

	lrr = arr - (arR * arR) / (4 * aRR);

	for(i = 0; i < 4; i++)
	{
		lDr[i] = aDr[i] - (arR / (2 * aRR)) * aDR[i];
	}

	coef = 1;
	coeg = 1 / (4 * aRR);

	f1 = lDr[0];
	f2 = lDr[1];
	f3 = lDr[2];
	f4 = lDr[3];
	g1 = aDR[0];
	g2 = aDR[1];
	g3 = aDR[2];
	g4 = aDR[3];
	h1 = lDr[0];
	h2 = lDr[1];
	h3 = lDr[2];
	h4 = lDr[3];

	h21 = lDr[0];
	h22 = lDr[1];
	h23 = lDr[2];
	h24 = lDr[3];

	coec = 4 * pi * pow(sqrt(pi / aRR), 3);
	coet = 16 * pi * sqrt(24 * pi) * pow(sqrt(pi / aRR), 3);
	coeso = (4 * pi) / aRR * sqrt(pi / (4 * aRR)) * (-sqrt(3)) * (16 * sqrt(2) * pi * pi) / 3;

	alpha_i = 1;
	beta_i = 0;

	hrho11 = (alpha_i * aDir[0] + beta_i * aDiR[0]) - ((alpha_i * airR + 2 * beta_i * aiRR) / (2 * aRR)) * aDR[0];
	hrho12 = (alpha_i * aDir[1] + beta_i * aDiR[1]) - ((alpha_i * airR + 2 * beta_i * aiRR) / (2 * aRR)) * aDR[1];
	hrho13 = (alpha_i * aDir[2] + beta_i * aDiR[2]) - ((alpha_i * airR + 2 * beta_i * aiRR) / (2 * aRR)) * aDR[2];
	hrho14 = (alpha_i * aDir[3] + beta_i * aDiR[3]) - ((alpha_i * airR + 2 * beta_i * aiRR) / (2 * aRR)) * aDR[3];

	txrp_cent[ 0] = coef;
	txrp_cent[ 1] = f1;
	txrp_cent[ 2] = f2;
	txrp_cent[ 3] = f3;
	txrp_cent[ 4] = f4;
	txrp_cent[ 5] = coeg;
	txrp_cent[ 6] = g1;
	txrp_cent[ 7] = g2;
	txrp_cent[ 8] = g3;
	txrp_cent[ 9] = g4;

	txrp_cent[18] = lrr;
	txrp_cent[20] = Nnl(l1, nu1) * Nnl(l2, nu2) * Nnl(l3, nu3) * Nnl(l4, nu4);
	


	switch(vt)
	{
		case 1:
			txrp_cent[19] = coec;
			break;

		case 2:
			txrp_cent[10] = h1;
			txrp_cent[11] = h2;
			txrp_cent[12] = h3;
			txrp_cent[13] = h4;
			txrp_cent[19] = coet;
		break;
		
		case 5:
			txrp_cent[10] = hrho11;
			txrp_cent[11] = hrho12;
			txrp_cent[12] = hrho13;
			txrp_cent[13] = hrho14;
			txrp_cent[14] = h21;
			txrp_cent[15] = h22;
			txrp_cent[16] = h23;
			txrp_cent[17] = h24;
			txrp_cent[19] =-coeso;
		break;
		
		default:
		break;
	}
}



void getT2r(basis_qnum qf, basis_qnum qi, double* txrp_cent, int c, int vt)
{
	int a = qf.c;
	int b = qi.c;

	int l1 = qf.lrho;
	int l2 = qf.llam;
	int l3 = qi.lrho;
	int l4 = qi.llam;

	double nu1 = qf.nurho;
	double nu2 = qf.nulam;
	double nu3 = qi.nurho;
	double nu4 = qi.nulam;

	double m1 = qf.m1;
	double m2 = qf.m2;
	double m3 = qf.m3;

	double alpha_ac, beta_ac, gamma_ac, delta_ac, alpha_bc, beta_bc, gamma_bc, delta_bc;

	double mi, mj, mk;
	double afrr, airr, afrR, airR, afRR, aiRR;
	double arr, arR, aRR;
	double aDfr[4], aDir[4], aDfR[4], aDiR[4];
	double aDr[4], aDR[4];
	double lrr;
	double lDr[4];
	double coef, coeg;
	double f1, f2, f3, f4;
	double g1, g2, g3, g4;
	double h1, h2, h3, h4;
	double alpha_i, alpha_j, beta_i, beta_j;
	double hi11, hi12, hi13, hi14;
	double hj11, hj12, hj13, hj14;
	double h21, h22, h23, h24;
	double coec, coet, coeso;
	double hrho11, hrho12, hrho13, hrho14;
	double hlambda11, hlambda12, hlambda13, hlambda14;
	double pi = 3.14159265358979324;
	int i;

	coordinatesTransformation_r_2(m1, m2, m3, a, c, &alpha_ac, &beta_ac, &gamma_ac, &delta_ac);
	coordinatesTransformation_r_2(m1, m2, m3, b, c, &alpha_bc, &beta_bc, &gamma_bc, &delta_bc);

	getijk(m1, m2, m3, &mi, &mj, &mk, c);

	afrr = (nu1 * alpha_ac * alpha_ac + nu2 * gamma_ac * gamma_ac);
	airr = (nu3 * alpha_bc * alpha_bc + nu4 * gamma_bc * gamma_bc);
	afrR = 2 * (nu1 * alpha_ac * beta_ac + nu2 * gamma_ac * delta_ac);
	airR = 2 * (nu3 * alpha_bc * beta_bc + nu4 * gamma_bc * delta_bc);
	afRR = (nu1 * beta_ac * beta_ac + nu2 * delta_ac * delta_ac);
	aiRR = (nu3 * beta_bc * beta_bc + nu4 * delta_bc * delta_bc);

	arr = afrr + airr;
	arR = afrR + airR;
	aRR = afRR + aiRR;

	aDfr[0] = 2 * alpha_ac;
	aDfr[1] = 2 * gamma_ac;
	aDfr[2] = 0;
	aDfr[3] = 0;
	aDir[0] = 0;
	aDir[1] = 0;
	aDir[2] = 2 * alpha_bc;
	aDir[3] = 2 * gamma_bc;
	aDfR[0] = 2 * beta_ac;
	aDfR[1] = 2 * delta_ac;
	aDfR[2] = 0;
	aDfR[3] = 0;
	aDiR[0] = 0;
	aDiR[1] = 0;
	aDiR[2] = 2 * beta_bc;
	aDiR[3] = 2 * delta_bc;

	for(i = 0; i < 4; i++)
	{
		aDr[i] = aDfr[i] + aDir[i];
		aDR[i] = aDfR[i] + aDiR[i];
	}

	lrr = arr - (arR * arR) / (4 * aRR);

	for(i = 0; i < 4; i++)
	{
		lDr[i] = aDr[i] - (arR / (2 * aRR)) * aDR[i];
	}

	coef = 1;
	coeg = 1 / (4 * aRR);

	f1 = lDr[0];
	f2 = lDr[1];
	f3 = lDr[2];
	f4 = lDr[3];
	g1 = aDR[0];
	g2 = aDR[1];
	g3 = aDR[2];
	g4 = aDR[3];
	h1 = lDr[0];
	h2 = lDr[1];
	h3 = lDr[2];
	h4 = lDr[3];

	alpha_i = 1;
	alpha_j = -1;
	beta_i = mi / (mi + mj);
	beta_j = mj / (mi + mj);

	hi11 = (alpha_i * aDir[0] + beta_i * aDiR[0]) - ((alpha_i * airR + 2 * beta_i * aiRR) / (2 * aRR)) * aDR[0];
	hi12 = (alpha_i * aDir[1] + beta_i * aDiR[1]) - ((alpha_i * airR + 2 * beta_i * aiRR) / (2 * aRR)) * aDR[1];
	hi13 = (alpha_i * aDir[2] + beta_i * aDiR[2]) - ((alpha_i * airR + 2 * beta_i * aiRR) / (2 * aRR)) * aDR[2];
	hi14 = (alpha_i * aDir[3] + beta_i * aDiR[3]) - ((alpha_i * airR + 2 * beta_i * aiRR) / (2 * aRR)) * aDR[3];

	hj11 = (alpha_j * aDir[0] + beta_j * aDiR[0]) - ((alpha_j * airR + 2 * beta_j * aiRR) / (2 * aRR)) * aDR[0];
	hj12 = (alpha_j * aDir[1] + beta_j * aDiR[1]) - ((alpha_j * airR + 2 * beta_j * aiRR) / (2 * aRR)) * aDR[1];
	hj13 = (alpha_j * aDir[2] + beta_j * aDiR[2]) - ((alpha_j * airR + 2 * beta_j * aiRR) / (2 * aRR)) * aDR[2];
	hj14 = (alpha_j * aDir[3] + beta_j * aDiR[3]) - ((alpha_j * airR + 2 * beta_j * aiRR) / (2 * aRR)) * aDR[3];

	h21 = lDr[0];
	h22 = lDr[1];
	h23 = lDr[2];
	h24 = lDr[3];

	coec = 4 * pi * pow(sqrt(pi / aRR), 3);
	coet = 16 * pi * sqrt(24 * pi) * pow(sqrt(pi / aRR), 3);
	coeso = (4 * pi) / aRR * sqrt(pi / (4 * aRR)) * (-sqrt(3)) * (16 * sqrt(2) * pi * pi) / 3;

	alpha_i = 1;
	alpha_j = 0;
	beta_i = 0;
	beta_j = 1;

	hrho11 = (alpha_i * aDir[0] + beta_i * aDiR[0]) - ((alpha_i * airR + 2 * beta_i * aiRR) / (2 * aRR)) * aDR[0];
	hrho12 = (alpha_i * aDir[1] + beta_i * aDiR[1]) - ((alpha_i * airR + 2 * beta_i * aiRR) / (2 * aRR)) * aDR[1];
	hrho13 = (alpha_i * aDir[2] + beta_i * aDiR[2]) - ((alpha_i * airR + 2 * beta_i * aiRR) / (2 * aRR)) * aDR[2];
	hrho14 = (alpha_i * aDir[3] + beta_i * aDiR[3]) - ((alpha_i * airR + 2 * beta_i * aiRR) / (2 * aRR)) * aDR[3];

	hlambda11 = (alpha_j * aDir[0] + beta_j * aDiR[0]) - ((alpha_j * airR + 2 * beta_j * aiRR) / (2 * aRR)) * aDR[0];
	hlambda12 = (alpha_j * aDir[1] + beta_j * aDiR[1]) - ((alpha_j * airR + 2 * beta_j * aiRR) / (2 * aRR)) * aDR[1];
	hlambda13 = (alpha_j * aDir[2] + beta_j * aDiR[2]) - ((alpha_j * airR + 2 * beta_j * aiRR) / (2 * aRR)) * aDR[2];
	hlambda14 = (alpha_j * aDir[3] + beta_j * aDiR[3]) - ((alpha_j * airR + 2 * beta_j * aiRR) / (2 * aRR)) * aDR[3];

	txrp_cent[ 0] = coef;
	txrp_cent[ 1] = f1;
	txrp_cent[ 2] = f2;
	txrp_cent[ 3] = f3;
	txrp_cent[ 4] = f4;
	txrp_cent[ 5] = coeg;
	txrp_cent[ 6] = g1;
	txrp_cent[ 7] = g2;
	txrp_cent[ 8] = g3;
	txrp_cent[ 9] = g4;

	txrp_cent[18] = lrr;
	txrp_cent[20] = Nnl(l1, nu1) * Nnl(l2, nu2) * Nnl(l3, nu3) * Nnl(l4, nu4);
	


	switch(vt)
	{
		case 1:
			txrp_cent[19] = coec;
			break;

		case 2:
			txrp_cent[10] = h1;
			txrp_cent[11] = h2;
			txrp_cent[12] = h3;
			txrp_cent[13] = h4;
			txrp_cent[19] = coet;
		break;

		case 3:
			txrp_cent[10] = hi11;
			txrp_cent[11] = hi12;
			txrp_cent[12] = hi13;
			txrp_cent[13] = hi14;
			txrp_cent[14] = h21;
			txrp_cent[15] = h22;
			txrp_cent[16] = h23;
			txrp_cent[17] = h24;
			txrp_cent[19] = coeso;
		break;

		case 4:
			txrp_cent[10] = hj11;
			txrp_cent[11] = hj12;
			txrp_cent[12] = hj13;
			txrp_cent[13] = hj14;
			txrp_cent[14] = h21;
			txrp_cent[15] = h22;
			txrp_cent[16] = h23;
			txrp_cent[17] = h24;
			txrp_cent[19] = coeso;
		break;

		case 5:
			txrp_cent[10] = hrho11;
			txrp_cent[11] = hrho12;
			txrp_cent[12] = hrho13;
			txrp_cent[13] = hrho14;
			txrp_cent[14] = h21;
			txrp_cent[15] = h22;
			txrp_cent[16] = h23;
			txrp_cent[17] = h24;
			txrp_cent[19] = coeso;
		break;

		case 6:
			txrp_cent[10] = hlambda11;
			txrp_cent[11] = hlambda12;
			txrp_cent[12] = hlambda13;
			txrp_cent[13] = hlambda14;
			txrp_cent[14] = h21;
			txrp_cent[15] = h22;
			txrp_cent[16] = h23;
			txrp_cent[17] = h24;
			txrp_cent[19] = coeso;
		break;

		default:
		break;
	}
}


void getT1p(basis_qnum qf, basis_qnum qi, double* txrp_cent, int c, int vt)
{
	int a = qf.c;
	int b = qi.c;

	int l1 = qf.lrho;
	int l2 = qf.llam;
	int l3 = qi.lrho;
	int l4 = qi.llam;

	double nu1 = 1/(4*qf.nurho);
	double nu2 = 1/(4*qf.nulam);
	double nu3 = 1/(4*qi.nurho);
	double nu4 = 1/(4*qi.nulam);

	double m1 = qf.m1;
	double m2 = qf.m2;
	double m3 = qf.m3;

	double alpha_ac, beta_ac, gamma_ac, delta_ac, alpha_bc, beta_bc, gamma_bc, delta_bc;

	double mi, mj, mk;
	double afrr, airr, afrR, airR, afRR, aiRR;
	double arr, arR, aRR;
	double aDfr[4], aDir[4], aDfR[4], aDiR[4];
	double aDr[4], aDR[4];
	double lrr;
	double lDr[4];
	double coef, coeg;
	double f1, f2, f3, f4;
	double g1, g2, g3, g4;
	double h1, h2, h3, h4;
	double alpha_i, alpha_j, beta_i, beta_j;
	double hi11, hi12, hi13, hi14;
	double hj11, hj12, hj13, hj14;
	double h21, h22, h23, h24;
	double coec, coet, coeso;
	double hrho11, hrho12, hrho13, hrho14;
	double hlambda11, hlambda12, hlambda13, hlambda14;
	double pi = 3.14159265358979324;
	int i;

	coordinatesTransformation_p_1(m1, m2, m3, a, c, &alpha_ac, &beta_ac, &gamma_ac, &delta_ac);
	coordinatesTransformation_p_1(m1, m2, m3, b, c, &alpha_bc, &beta_bc, &gamma_bc, &delta_bc);

	getijk(m1, m2, m3, &mi, &mj, &mk, c);

	afrr = (nu1 * alpha_ac * alpha_ac + nu2 * gamma_ac * gamma_ac);
	airr = (nu3 * alpha_bc * alpha_bc + nu4 * gamma_bc * gamma_bc);
	afrR = 2 * (nu1 * alpha_ac * beta_ac + nu2 * gamma_ac * delta_ac);
	airR = 2 * (nu3 * alpha_bc * beta_bc + nu4 * gamma_bc * delta_bc);
	afRR = (nu1 * beta_ac * beta_ac + nu2 * delta_ac * delta_ac);
	aiRR = (nu3 * beta_bc * beta_bc + nu4 * delta_bc * delta_bc);

	arr = afrr + airr;
	arR = afrR + airR;
	aRR = afRR + aiRR;

	aDfr[0] = 2 * alpha_ac;
	aDfr[1] = 2 * gamma_ac;
	aDfr[2] = 0;
	aDfr[3] = 0;
	aDir[0] = 0;
	aDir[1] = 0;
	aDir[2] = 2 * alpha_bc;
	aDir[3] = 2 * gamma_bc;
	aDfR[0] = 2 * beta_ac;
	aDfR[1] = 2 * delta_ac;
	aDfR[2] = 0;
	aDfR[3] = 0;
	aDiR[0] = 0;
	aDiR[1] = 0;
	aDiR[2] = 2 * beta_bc;
	aDiR[3] = 2 * delta_bc;

	for(i = 0; i < 4; i++)
	{
		aDr[i] = aDfr[i] + aDir[i];
		aDR[i] = aDfR[i] + aDiR[i];
	}

	lrr = arr - (arR * arR) / (4 * aRR);

	for(i = 0; i < 4; i++)
	{
		lDr[i] = aDr[i] - (arR / (2 * aRR)) * aDR[i];
	}

	coef = 1;
	coeg = 1 / (4 * aRR);

	f1 = lDr[0];
	f2 = lDr[1];
	f3 = lDr[2];
	f4 = lDr[3];
	g1 = aDR[0];
	g2 = aDR[1];
	g3 = aDR[2];
	g4 = aDR[3];
	h1 = lDr[0];
	h2 = lDr[1];
	h3 = lDr[2];
	h4 = lDr[3];

	alpha_i = 1;
	alpha_j = -1;
	beta_i = mi / (mi + mj);
	beta_j = mj / (mi + mj);

	hi11 = (alpha_i * aDir[0] + beta_i * aDiR[0]) - ((alpha_i * airR + 2 * beta_i * aiRR) / (2 * aRR)) * aDR[0];
	hi12 = (alpha_i * aDir[1] + beta_i * aDiR[1]) - ((alpha_i * airR + 2 * beta_i * aiRR) / (2 * aRR)) * aDR[1];
	hi13 = (alpha_i * aDir[2] + beta_i * aDiR[2]) - ((alpha_i * airR + 2 * beta_i * aiRR) / (2 * aRR)) * aDR[2];
	hi14 = (alpha_i * aDir[3] + beta_i * aDiR[3]) - ((alpha_i * airR + 2 * beta_i * aiRR) / (2 * aRR)) * aDR[3];

	hj11 = (alpha_j * aDir[0] + beta_j * aDiR[0]) - ((alpha_j * airR + 2 * beta_j * aiRR) / (2 * aRR)) * aDR[0];
	hj12 = (alpha_j * aDir[1] + beta_j * aDiR[1]) - ((alpha_j * airR + 2 * beta_j * aiRR) / (2 * aRR)) * aDR[1];
	hj13 = (alpha_j * aDir[2] + beta_j * aDiR[2]) - ((alpha_j * airR + 2 * beta_j * aiRR) / (2 * aRR)) * aDR[2];
	hj14 = (alpha_j * aDir[3] + beta_j * aDiR[3]) - ((alpha_j * airR + 2 * beta_j * aiRR) / (2 * aRR)) * aDR[3];

	h21 = lDr[0];
	h22 = lDr[1];
	h23 = lDr[2];
	h24 = lDr[3];

	coec = 4 * pi * pow(sqrt(pi / aRR), 3);
	coet = 16 * pi * sqrt(24 * pi) * pow(sqrt(pi / aRR), 3);
	coeso = (4 * pi) / aRR * sqrt(pi / (4 * aRR)) * (-sqrt(3)) * (16 * sqrt(2) * pi * pi) / 3;

	alpha_i = 1;
	alpha_j = 0;
	beta_i = 0;
	beta_j = 1;

	hrho11 = (alpha_i * aDir[0] + beta_i * aDiR[0]) - ((alpha_i * airR + 2 * beta_i * aiRR) / (2 * aRR)) * aDR[0];
	hrho12 = (alpha_i * aDir[1] + beta_i * aDiR[1]) - ((alpha_i * airR + 2 * beta_i * aiRR) / (2 * aRR)) * aDR[1];
	hrho13 = (alpha_i * aDir[2] + beta_i * aDiR[2]) - ((alpha_i * airR + 2 * beta_i * aiRR) / (2 * aRR)) * aDR[2];
	hrho14 = (alpha_i * aDir[3] + beta_i * aDiR[3]) - ((alpha_i * airR + 2 * beta_i * aiRR) / (2 * aRR)) * aDR[3];

	hlambda11 = (alpha_j * aDir[0] + beta_j * aDiR[0]) - ((alpha_j * airR + 2 * beta_j * aiRR) / (2 * aRR)) * aDR[0];
	hlambda12 = (alpha_j * aDir[1] + beta_j * aDiR[1]) - ((alpha_j * airR + 2 * beta_j * aiRR) / (2 * aRR)) * aDR[1];
	hlambda13 = (alpha_j * aDir[2] + beta_j * aDiR[2]) - ((alpha_j * airR + 2 * beta_j * aiRR) / (2 * aRR)) * aDR[2];
	hlambda14 = (alpha_j * aDir[3] + beta_j * aDiR[3]) - ((alpha_j * airR + 2 * beta_j * aiRR) / (2 * aRR)) * aDR[3];

	txrp_cent[ 0] = coef;
	txrp_cent[ 1] = f1;
	txrp_cent[ 2] = f2;
	txrp_cent[ 3] = f3;
	txrp_cent[ 4] = f4;
	txrp_cent[ 5] = coeg;
	txrp_cent[ 6] = g1;
	txrp_cent[ 7] = g2;
	txrp_cent[ 8] = g3;
	txrp_cent[ 9] = g4;

	txrp_cent[18] = lrr;
	txrp_cent[20] = pow(-1,l3+l4+((l1+l2+l3+l4)%4)/2)*Nnl(l1, nu1) * Nnl(l2, nu2) * Nnl(l3, nu3) * Nnl(l4, nu4);
	


	switch(vt)
	{
		case 1:
			txrp_cent[19] = coec;
			break;

		case 2:
			txrp_cent[10] = h1;
			txrp_cent[11] = h2;
			txrp_cent[12] = h3;
			txrp_cent[13] = h4;
			txrp_cent[19] = coet;
		break;

		case 3:
			txrp_cent[10] = hi11;
			txrp_cent[11] = hi12;
			txrp_cent[12] = hi13;
			txrp_cent[13] = hi14;
			txrp_cent[14] = h21;
			txrp_cent[15] = h22;
			txrp_cent[16] = h23;
			txrp_cent[17] = h24;
			txrp_cent[19] = coeso;
		break;

		case 4:
			txrp_cent[10] = hj11;
			txrp_cent[11] = hj12;
			txrp_cent[12] = hj13;
			txrp_cent[13] = hj14;
			txrp_cent[14] = h21;
			txrp_cent[15] = h22;
			txrp_cent[16] = h23;
			txrp_cent[17] = h24;
			txrp_cent[19] = coeso;
		break;

		case 5:
			txrp_cent[10] = hrho11;
			txrp_cent[11] = hrho12;
			txrp_cent[12] = hrho13;
			txrp_cent[13] = hrho14;
			txrp_cent[14] = h21;
			txrp_cent[15] = h22;
			txrp_cent[16] = h23;
			txrp_cent[17] = h24;
			txrp_cent[19] = coeso;
		break;

		case 6:
			txrp_cent[10] = hlambda11;
			txrp_cent[11] = hlambda12;
			txrp_cent[12] = hlambda13;
			txrp_cent[13] = hlambda14;
			txrp_cent[14] = h21;
			txrp_cent[15] = h22;
			txrp_cent[16] = h23;
			txrp_cent[17] = h24;
			txrp_cent[19] = coeso;
		break;

		default:
		break;
	}
}


void getT2p(basis_qnum qf, basis_qnum qi, double* txrp_cent, int c, int vt)
{
	int a = qf.c;
	int b = qi.c;

	int l1 = qf.lrho;
	int l2 = qf.llam;
	int l3 = qi.lrho;
	int l4 = qi.llam;

	double nu1 = 1/(4*qf.nurho);
	double nu2 = 1/(4*qf.nulam);
	double nu3 = 1/(4*qi.nurho);
	double nu4 = 1/(4*qi.nulam);

	double m1 = qf.m1;
	double m2 = qf.m2;
	double m3 = qf.m3;

	double alpha_ac, beta_ac, gamma_ac, delta_ac, alpha_bc, beta_bc, gamma_bc, delta_bc;

	double mi, mj, mk;
	double afrr, airr, afrR, airR, afRR, aiRR;
	double arr, arR, aRR;
	double aDfr[4], aDir[4], aDfR[4], aDiR[4];
	double aDr[4], aDR[4];
	double lrr;
	double lDr[4];
	double coef, coeg;
	double f1, f2, f3, f4;
	double g1, g2, g3, g4;
	double h1, h2, h3, h4;
	double alpha_i, alpha_j, beta_i, beta_j;
	double hi11, hi12, hi13, hi14;
	double hj11, hj12, hj13, hj14;
	double h21, h22, h23, h24;
	double coec, coet, coeso;
	double hrho11, hrho12, hrho13, hrho14;
	double hlambda11, hlambda12, hlambda13, hlambda14;
	double pi = 3.14159265358979324;
	int i;

	coordinatesTransformation_p_2(m1, m2, m3, a, c, &alpha_ac, &beta_ac, &gamma_ac, &delta_ac);
	coordinatesTransformation_p_2(m1, m2, m3, b, c, &alpha_bc, &beta_bc, &gamma_bc, &delta_bc);

	getijk(m1, m2, m3, &mi, &mj, &mk, c);

	afrr = (nu1 * alpha_ac * alpha_ac + nu2 * gamma_ac * gamma_ac);
	airr = (nu3 * alpha_bc * alpha_bc + nu4 * gamma_bc * gamma_bc);
	afrR = 2 * (nu1 * alpha_ac * beta_ac + nu2 * gamma_ac * delta_ac);
	airR = 2 * (nu3 * alpha_bc * beta_bc + nu4 * gamma_bc * delta_bc);
	afRR = (nu1 * beta_ac * beta_ac + nu2 * delta_ac * delta_ac);
	aiRR = (nu3 * beta_bc * beta_bc + nu4 * delta_bc * delta_bc);

	arr = afrr + airr;
	arR = afrR + airR;
	aRR = afRR + aiRR;

	aDfr[0] = 2 * alpha_ac;
	aDfr[1] = 2 * gamma_ac;
	aDfr[2] = 0;
	aDfr[3] = 0;
	aDir[0] = 0;
	aDir[1] = 0;
	aDir[2] = 2 * alpha_bc;
	aDir[3] = 2 * gamma_bc;
	aDfR[0] = 2 * beta_ac;
	aDfR[1] = 2 * delta_ac;
	aDfR[2] = 0;
	aDfR[3] = 0;
	aDiR[0] = 0;
	aDiR[1] = 0;
	aDiR[2] = 2 * beta_bc;
	aDiR[3] = 2 * delta_bc;

	for(i = 0; i < 4; i++)
	{
		aDr[i] = aDfr[i] + aDir[i];
		aDR[i] = aDfR[i] + aDiR[i];
	}

	lrr = arr - (arR * arR) / (4 * aRR);

	for(i = 0; i < 4; i++)
	{
		lDr[i] = aDr[i] - (arR / (2 * aRR)) * aDR[i];
	}

	coef = 1;
	coeg = 1 / (4 * aRR);

	f1 = lDr[0];
	f2 = lDr[1];
	f3 = lDr[2];
	f4 = lDr[3];
	g1 = aDR[0];
	g2 = aDR[1];
	g3 = aDR[2];
	g4 = aDR[3];
	h1 = lDr[0];
	h2 = lDr[1];
	h3 = lDr[2];
	h4 = lDr[3];

	alpha_i = 1;
	alpha_j = -1;
	beta_i = mi / (mi + mj);
	beta_j = mj / (mi + mj);

	hi11 = (alpha_i * aDir[0] + beta_i * aDiR[0]) - ((alpha_i * airR + 2 * beta_i * aiRR) / (2 * aRR)) * aDR[0];
	hi12 = (alpha_i * aDir[1] + beta_i * aDiR[1]) - ((alpha_i * airR + 2 * beta_i * aiRR) / (2 * aRR)) * aDR[1];
	hi13 = (alpha_i * aDir[2] + beta_i * aDiR[2]) - ((alpha_i * airR + 2 * beta_i * aiRR) / (2 * aRR)) * aDR[2];
	hi14 = (alpha_i * aDir[3] + beta_i * aDiR[3]) - ((alpha_i * airR + 2 * beta_i * aiRR) / (2 * aRR)) * aDR[3];

	hj11 = (alpha_j * aDir[0] + beta_j * aDiR[0]) - ((alpha_j * airR + 2 * beta_j * aiRR) / (2 * aRR)) * aDR[0];
	hj12 = (alpha_j * aDir[1] + beta_j * aDiR[1]) - ((alpha_j * airR + 2 * beta_j * aiRR) / (2 * aRR)) * aDR[1];
	hj13 = (alpha_j * aDir[2] + beta_j * aDiR[2]) - ((alpha_j * airR + 2 * beta_j * aiRR) / (2 * aRR)) * aDR[2];
	hj14 = (alpha_j * aDir[3] + beta_j * aDiR[3]) - ((alpha_j * airR + 2 * beta_j * aiRR) / (2 * aRR)) * aDR[3];

	h21 = lDr[0];
	h22 = lDr[1];
	h23 = lDr[2];
	h24 = lDr[3];

	coec = 4 * pi * pow(sqrt(pi / aRR), 3);
	coet = 16 * pi * sqrt(24 * pi) * pow(sqrt(pi / aRR), 3);
	coeso = (4 * pi) / aRR * sqrt(pi / (4 * aRR)) * (-sqrt(3)) * (16 * sqrt(2) * pi * pi) / 3;

	alpha_i = 1;
	alpha_j = 0;
	beta_i = 0;
	beta_j = 1;

	hrho11 = (alpha_i * aDir[0] + beta_i * aDiR[0]) - ((alpha_i * airR + 2 * beta_i * aiRR) / (2 * aRR)) * aDR[0];
	hrho12 = (alpha_i * aDir[1] + beta_i * aDiR[1]) - ((alpha_i * airR + 2 * beta_i * aiRR) / (2 * aRR)) * aDR[1];
	hrho13 = (alpha_i * aDir[2] + beta_i * aDiR[2]) - ((alpha_i * airR + 2 * beta_i * aiRR) / (2 * aRR)) * aDR[2];
	hrho14 = (alpha_i * aDir[3] + beta_i * aDiR[3]) - ((alpha_i * airR + 2 * beta_i * aiRR) / (2 * aRR)) * aDR[3];

	hlambda11 = (alpha_j * aDir[0] + beta_j * aDiR[0]) - ((alpha_j * airR + 2 * beta_j * aiRR) / (2 * aRR)) * aDR[0];
	hlambda12 = (alpha_j * aDir[1] + beta_j * aDiR[1]) - ((alpha_j * airR + 2 * beta_j * aiRR) / (2 * aRR)) * aDR[1];
	hlambda13 = (alpha_j * aDir[2] + beta_j * aDiR[2]) - ((alpha_j * airR + 2 * beta_j * aiRR) / (2 * aRR)) * aDR[2];
	hlambda14 = (alpha_j * aDir[3] + beta_j * aDiR[3]) - ((alpha_j * airR + 2 * beta_j * aiRR) / (2 * aRR)) * aDR[3];

	txrp_cent[ 0] = coef;
	txrp_cent[ 1] = f1;
	txrp_cent[ 2] = f2;
	txrp_cent[ 3] = f3;
	txrp_cent[ 4] = f4;
	txrp_cent[ 5] = coeg;
	txrp_cent[ 6] = g1;
	txrp_cent[ 7] = g2;
	txrp_cent[ 8] = g3;
	txrp_cent[ 9] = g4;

	txrp_cent[18] = lrr;
	txrp_cent[20] = pow(-1,l3+l4+((l1+l2+l3+l4)%4)/2)*Nnl(l1, nu1) * Nnl(l2, nu2) * Nnl(l3, nu3) * Nnl(l4, nu4);
	


	switch(vt)
	{
		case 1:
			txrp_cent[19] = coec;
			break;

		case 2:
			txrp_cent[10] = h1;
			txrp_cent[11] = h2;
			txrp_cent[12] = h3;
			txrp_cent[13] = h4;
			txrp_cent[19] = coet;
		break;

		case 3:
			txrp_cent[10] = hi11;
			txrp_cent[11] = hi12;
			txrp_cent[12] = hi13;
			txrp_cent[13] = hi14;
			txrp_cent[14] = h21;
			txrp_cent[15] = h22;
			txrp_cent[16] = h23;
			txrp_cent[17] = h24;
			txrp_cent[19] = coeso;
		break;

		case 4:
			txrp_cent[10] = hj11;
			txrp_cent[11] = hj12;
			txrp_cent[12] = hj13;
			txrp_cent[13] = hj14;
			txrp_cent[14] = h21;
			txrp_cent[15] = h22;
			txrp_cent[16] = h23;
			txrp_cent[17] = h24;
			txrp_cent[19] = coeso;
		break;

		case 5:
			txrp_cent[10] = hrho11;
			txrp_cent[11] = hrho12;
			txrp_cent[12] = hrho13;
			txrp_cent[13] = hrho14;
			txrp_cent[14] = h21;
			txrp_cent[15] = h22;
			txrp_cent[16] = h23;
			txrp_cent[17] = h24;
			txrp_cent[19] = coeso;
		break;

		case 6:
			txrp_cent[10] = hlambda11;
			txrp_cent[11] = hlambda12;
			txrp_cent[12] = hlambda13;
			txrp_cent[13] = hlambda14;
			txrp_cent[14] = h21;
			txrp_cent[15] = h22;
			txrp_cent[16] = h23;
			txrp_cent[17] = h24;
			txrp_cent[19] = coeso;
		break;

		default:
		break;
	}
}


void getTpi(basis_qnum qf, basis_qnum qi, double* txrp_cent, int c, int vt)
{
	int a = qf.c;
	int b = qi.c;

	int l1 = qf.lrho;
	int l2 = qf.llam;
	int l3 = qi.lrho;
	int l4 = qi.llam;

	double nu1 = 1/(4*qf.nurho);
	double nu2 = 1/(4*qf.nulam);
	double nu3 = 1/(4*qi.nurho);
	double nu4 = 1/(4*qi.nulam);

	double m1 = qf.m1;
	double m2 = qf.m2;
	double m3 = qf.m3;

	double alpha_ac, beta_ac, gamma_ac, delta_ac, alpha_bc, beta_bc, gamma_bc, delta_bc;

	double mi, mj, mk;
	double afrr, airr, afrR, airR, afRR, aiRR;
	double arr, arR, aRR;
	double aDfr[4], aDir[4], aDfR[4], aDiR[4];
	double aDr[4], aDR[4];
	double lrr;
	double lDr[4];
	double coef, coeg;
	double f1, f2, f3, f4;
	double g1, g2, g3, g4;
	double h1, h2, h3, h4;
	double alpha_i, alpha_j, beta_i, beta_j;
	double hi11, hi12, hi13, hi14;
	double hj11, hj12, hj13, hj14;
	double h21, h22, h23, h24;
	double coec, coet, coeso;
	double hrho11, hrho12, hrho13, hrho14;
	double hlambda11, hlambda12, hlambda13, hlambda14;
	double pi = 3.14159265358979324;
	int i;

	coordinatesTransformation_p_i(m1, m2, m3, a, c, &alpha_ac, &beta_ac, &gamma_ac, &delta_ac);
	coordinatesTransformation_p_i(m1, m2, m3, b, c, &alpha_bc, &beta_bc, &gamma_bc, &delta_bc);

	getijk(m1, m2, m3, &mi, &mj, &mk, c);

	afrr = (nu1 * alpha_ac * alpha_ac + nu2 * gamma_ac * gamma_ac);
	airr = (nu3 * alpha_bc * alpha_bc + nu4 * gamma_bc * gamma_bc);
	afrR = 2 * (nu1 * alpha_ac * beta_ac + nu2 * gamma_ac * delta_ac);
	airR = 2 * (nu3 * alpha_bc * beta_bc + nu4 * gamma_bc * delta_bc);
	afRR = (nu1 * beta_ac * beta_ac + nu2 * delta_ac * delta_ac);
	aiRR = (nu3 * beta_bc * beta_bc + nu4 * delta_bc * delta_bc);

	arr = afrr + airr;
	arR = afrR + airR;
	aRR = afRR + aiRR;

	aDfr[0] = 2 * alpha_ac;
	aDfr[1] = 2 * gamma_ac;
	aDfr[2] = 0;
	aDfr[3] = 0;
	aDir[0] = 0;
	aDir[1] = 0;
	aDir[2] = 2 * alpha_bc;
	aDir[3] = 2 * gamma_bc;
	aDfR[0] = 2 * beta_ac;
	aDfR[1] = 2 * delta_ac;
	aDfR[2] = 0;
	aDfR[3] = 0;
	aDiR[0] = 0;
	aDiR[1] = 0;
	aDiR[2] = 2 * beta_bc;
	aDiR[3] = 2 * delta_bc;

	for(i = 0; i < 4; i++)
	{
		aDr[i] = aDfr[i] + aDir[i];
		aDR[i] = aDfR[i] + aDiR[i];
	}

	lrr = arr - (arR * arR) / (4 * aRR);

	for(i = 0; i < 4; i++)
	{
		lDr[i] = aDr[i] - (arR / (2 * aRR)) * aDR[i];
	}

	coef = 1;
	coeg = 1 / (4 * aRR);

	f1 = lDr[0];
	f2 = lDr[1];
	f3 = lDr[2];
	f4 = lDr[3];
	g1 = aDR[0];
	g2 = aDR[1];
	g3 = aDR[2];
	g4 = aDR[3];
	h1 = lDr[0];
	h2 = lDr[1];
	h3 = lDr[2];
	h4 = lDr[3];

	alpha_i = 1;
	alpha_j = -1;
	beta_i = mi / (mi + mj);
	beta_j = mj / (mi + mj);

	hi11 = (alpha_i * aDir[0] + beta_i * aDiR[0]) - ((alpha_i * airR + 2 * beta_i * aiRR) / (2 * aRR)) * aDR[0];
	hi12 = (alpha_i * aDir[1] + beta_i * aDiR[1]) - ((alpha_i * airR + 2 * beta_i * aiRR) / (2 * aRR)) * aDR[1];
	hi13 = (alpha_i * aDir[2] + beta_i * aDiR[2]) - ((alpha_i * airR + 2 * beta_i * aiRR) / (2 * aRR)) * aDR[2];
	hi14 = (alpha_i * aDir[3] + beta_i * aDiR[3]) - ((alpha_i * airR + 2 * beta_i * aiRR) / (2 * aRR)) * aDR[3];

	hj11 = (alpha_j * aDir[0] + beta_j * aDiR[0]) - ((alpha_j * airR + 2 * beta_j * aiRR) / (2 * aRR)) * aDR[0];
	hj12 = (alpha_j * aDir[1] + beta_j * aDiR[1]) - ((alpha_j * airR + 2 * beta_j * aiRR) / (2 * aRR)) * aDR[1];
	hj13 = (alpha_j * aDir[2] + beta_j * aDiR[2]) - ((alpha_j * airR + 2 * beta_j * aiRR) / (2 * aRR)) * aDR[2];
	hj14 = (alpha_j * aDir[3] + beta_j * aDiR[3]) - ((alpha_j * airR + 2 * beta_j * aiRR) / (2 * aRR)) * aDR[3];

	h21 = lDr[0];
	h22 = lDr[1];
	h23 = lDr[2];
	h24 = lDr[3];

	coec = 4 * pi * pow(sqrt(pi / aRR), 3);
	coet = 16 * pi * sqrt(24 * pi) * pow(sqrt(pi / aRR), 3);
	coeso = (4 * pi) / aRR * sqrt(pi / (4 * aRR)) * (-sqrt(3)) * (16 * sqrt(2) * pi * pi) / 3;

	alpha_i = 1;
	alpha_j = 0;
	beta_i = 0;
	beta_j = 1;

	hrho11 = (alpha_i * aDir[0] + beta_i * aDiR[0]) - ((alpha_i * airR + 2 * beta_i * aiRR) / (2 * aRR)) * aDR[0];
	hrho12 = (alpha_i * aDir[1] + beta_i * aDiR[1]) - ((alpha_i * airR + 2 * beta_i * aiRR) / (2 * aRR)) * aDR[1];
	hrho13 = (alpha_i * aDir[2] + beta_i * aDiR[2]) - ((alpha_i * airR + 2 * beta_i * aiRR) / (2 * aRR)) * aDR[2];
	hrho14 = (alpha_i * aDir[3] + beta_i * aDiR[3]) - ((alpha_i * airR + 2 * beta_i * aiRR) / (2 * aRR)) * aDR[3];

	hlambda11 = (alpha_j * aDir[0] + beta_j * aDiR[0]) - ((alpha_j * airR + 2 * beta_j * aiRR) / (2 * aRR)) * aDR[0];
	hlambda12 = (alpha_j * aDir[1] + beta_j * aDiR[1]) - ((alpha_j * airR + 2 * beta_j * aiRR) / (2 * aRR)) * aDR[1];
	hlambda13 = (alpha_j * aDir[2] + beta_j * aDiR[2]) - ((alpha_j * airR + 2 * beta_j * aiRR) / (2 * aRR)) * aDR[2];
	hlambda14 = (alpha_j * aDir[3] + beta_j * aDiR[3]) - ((alpha_j * airR + 2 * beta_j * aiRR) / (2 * aRR)) * aDR[3];

	txrp_cent[ 0] = coef;
	txrp_cent[ 1] = f1;
	txrp_cent[ 2] = f2;
	txrp_cent[ 3] = f3;
	txrp_cent[ 4] = f4;
	txrp_cent[ 5] = coeg;
	txrp_cent[ 6] = g1;
	txrp_cent[ 7] = g2;
	txrp_cent[ 8] = g3;
	txrp_cent[ 9] = g4;

	txrp_cent[18] = lrr;
	txrp_cent[20] = pow(-1,l3+l4+((l1+l2+l3+l4)%4)/2)*Nnl(l1, nu1) * Nnl(l2, nu2) * Nnl(l3, nu3) * Nnl(l4, nu4);
	


	switch(vt)
	{
		case 1:
			txrp_cent[19] = coec;
			break;

		case 2:
			txrp_cent[10] = h1;
			txrp_cent[11] = h2;
			txrp_cent[12] = h3;
			txrp_cent[13] = h4;
			txrp_cent[19] = coet;
		break;

		case 3:
			txrp_cent[10] = hi11;
			txrp_cent[11] = hi12;
			txrp_cent[12] = hi13;
			txrp_cent[13] = hi14;
			txrp_cent[14] = h21;
			txrp_cent[15] = h22;
			txrp_cent[16] = h23;
			txrp_cent[17] = h24;
			txrp_cent[19] = coeso;
		break;

		case 4:
			txrp_cent[10] = hj11;
			txrp_cent[11] = hj12;
			txrp_cent[12] = hj13;
			txrp_cent[13] = hj14;
			txrp_cent[14] = h21;
			txrp_cent[15] = h22;
			txrp_cent[16] = h23;
			txrp_cent[17] = h24;
			txrp_cent[19] = coeso;
		break;

		case 5:
			txrp_cent[10] = hrho11;
			txrp_cent[11] = hrho12;
			txrp_cent[12] = hrho13;
			txrp_cent[13] = hrho14;
			txrp_cent[14] = h21;
			txrp_cent[15] = h22;
			txrp_cent[16] = h23;
			txrp_cent[17] = h24;
			txrp_cent[19] = coeso;
		break;

		case 6:
			txrp_cent[10] = hlambda11;
			txrp_cent[11] = hlambda12;
			txrp_cent[12] = hlambda13;
			txrp_cent[13] = hlambda14;
			txrp_cent[14] = h21;
			txrp_cent[15] = h22;
			txrp_cent[16] = h23;
			txrp_cent[17] = h24;
			txrp_cent[19] = coeso;
		break;

		default:
		break;
	}
}


typedef void (*txrp_vtype)(basis_qnum,basis_qnum,double*,int*,int);

void t1r_cent(basis_qnum qf,basis_qnum qi,double *t,int *lent,int c)
{
	getT1r(qf,qi,t,c,1);
	*lent=10;
}

void t2r_cent(basis_qnum qf,basis_qnum qi,double *t,int *lent,int c)
{
	getT2r(qf,qi,t,c,1);
	*lent=10;
}

void t1p_cent(basis_qnum qf,basis_qnum qi,double *t,int *lent,int c)
{
	getT1p(qf,qi,t,c,1);
	*lent=10;
}

void t2p_cent(basis_qnum qf,basis_qnum qi,double *t,int *lent,int c)
{
	getT2p(qf,qi,t,c,1);
	*lent=10;
}

void tpi_cent(basis_qnum qf,basis_qnum qi,double *t,int *lent,int c)
{
	getTpi(qf,qi,t,c,1);
	*lent=10;
}

void t1r_tens(basis_qnum qf,basis_qnum qi,double *t,int *lent,int c)
{
	getT1r(qf,qi,t,c,2);
	*lent=14;
}

void t1r_soii(basis_qnum qf,basis_qnum qi,double *t,int *lent,int c)
{
	getT1r(qf,qi,t,c,3);
	*lent=18;
}

void t1r_soij(basis_qnum qf,basis_qnum qi,double *t,int *lent,int c)
{
	getT1r(qf,qi,t,c,3);
	*lent=18;
}

void t1r_soji(basis_qnum qf,basis_qnum qi,double *t,int *lent,int c)
{
	getT1r(qf,qi,t,c,4);
	*lent=18;
}

void t1r_sojj(basis_qnum qf,basis_qnum qi,double *t,int *lent,int c)
{
	getT1r(qf,qi,t,c,4);
	*lent=18;
}

void t1r_sorr(basis_qnum qf,basis_qnum qi,double *t,int *lent,int c)
{
	getT1r(qf,qi,t,c,5);
	*lent=18;
}

void t1r_sorl(basis_qnum qf,basis_qnum qi,double *t,int *lent,int c)
{
	getT1r(qf,qi,t,c,6);
	*lent=18;
}

void tir_cent(basis_qnum qf,basis_qnum qi,double *t,int *lent,int c)
{
	getTir(qf,qi,t,c,1);
	*lent=10;
}

void tir_tens(basis_qnum qf,basis_qnum qi,double *t,int *lent,int c)
{
	getTir(qf,qi,t,c,2);
	*lent=14;
}

void tir_soii(basis_qnum qf,basis_qnum qi,double *t,int *lent,int c)
{
	getTir(qf,qi,t,c,5);
	*lent=18;
}

void tir_soij(basis_qnum qf,basis_qnum qi,double *t,int *lent,int c)
{
	getTir(qf,qi,t,c,5);
	*lent=18;
}

void tjr_cent(basis_qnum qf,basis_qnum qi,double *t,int *lent,int c)
{
	getTjr(qf,qi,t,c,1);
	*lent=10;
}

void tjr_tens(basis_qnum qf,basis_qnum qi,double *t,int *lent,int c)
{
	getTjr(qf,qi,t,c,2);
	*lent=14;
}

void tjr_sojj(basis_qnum qf,basis_qnum qi,double *t,int *lent,int c)
{
	getTjr(qf,qi,t,c,5);
	*lent=18;
}

void tjr_soji(basis_qnum qf,basis_qnum qi,double *t,int *lent,int c)
{
	getTjr(qf,qi,t,c,5);
	*lent=18;
}



double inteVcenPart(txrp_vtype get_trxp_vtype,sumckdk_scdk scdk,basis_qnum qf,basis_qnum qi,vargs varg,vCenPart v,int c)
{
	int i,j;
	double p,sum;
	double t[21];
	int lent;

	get_trxp_vtype(qf,qi,t,&lent,c);

	sum=0;
	for(i=0;i<scdk.len;i++)
	{
		p=scdk.coe[i];
		for(j=0;j<lent;j++)
		{
			p*=pow(t[j],scdk.nfgh[i][j]);
		}
		p*=t[19]*t[20]*InteCenV(t[18],scdk.nfgh[i][18],qf,qi,v,varg,c);
		sum+=p;
	}
	return sum;
}

double InteExprn(double a,int n)
{
	double GammaSemiInt[100]={
 0.00000000000000000000e+000, 1.77245385090551602729e+000, 1.00000000000000000000e+000, 8.86226925452758013649e-001, 1.00000000000000000000e+000, 1.32934038817913702047e+000, 2.00000000000000000000e+000, 3.32335097044784255118e+000, 6.00000000000000000000e+000, 1.16317283965674489291e+001, 2.40000000000000000000e+001, 5.23427777845535201811e+001, 1.20000000000000000000e+002, 2.87885277815044360996e+002, 7.20000000000000000000e+002, 1.87125430579778834647e+003, 5.04000000000000000000e+003, 1.40344072934834125985e+004, 4.03200000000000000000e+004, 1.19292461994609007087e+005,
 3.62880000000000000000e+005, 1.13327838894878556733e+006, 3.62880000000000000000e+006, 1.18994230839622484570e+007, 3.99168000000000000000e+007, 1.36843365465565857255e+008, 4.79001600000000000000e+008, 1.71054206831957321569e+009, 6.22702080000000000000e+009, 2.30923179223142384118e+010, 8.71782912000000000000e+010, 3.34838609873556456972e+011, 1.30767436800000000000e+012, 5.18999845304012508307e+012, 2.09227898880000000000e+013, 8.56349744751620638707e+013, 3.55687428096000000000e+014, 1.49861205331533611773e+015, 6.40237370572800000000e+015, 2.77243229863337181781e+016,
 1.21645100408832000000e+017, 5.40624298233507504473e+017, 2.43290200817664000000e+018, 1.10827981137869038417e+019, 5.10909421717094400000e+019, 2.38280159446418432596e+020, 1.12400072777760768000e+021, 5.36130358754441473342e+021, 2.58520167388849766400e+022, 1.25990634307293746235e+023, 6.20448401733239439360e+023, 3.08677054052869678277e+024, 1.55112100433309859840e+025, 7.87126487834817679606e+025, 4.03291461126605635584e+026, 2.08588519276226685095e+027, 1.08888694504183521607e+028, 5.73618428009623384013e+028, 3.04888344611713860501e+029, 1.63481251982742664443e+030,
 8.84176199373970195454e+030, 4.82269693349090860109e+031, 2.65252859812191058636e+032, 1.47092256471472712333e+033, 8.22283865417792281772e+033, 4.63340607885139043849e+034, 2.63130836933693530167e+035, 1.50585697562670189251e+036, 8.68331761881188649551e+036, 5.04462086834945133991e+037, 2.95232799039604140847e+038, 1.74039419958056071227e+039, 1.03331479663861449296e+040, 6.17839940851099052856e+040, 3.71993326789901217468e+041, 2.25511578410651154292e+042, 1.37637530912263450463e+043, 8.45668419039941828596e+043, 5.23022617466601111760e+044, 3.25582341330377604009e+045,
 2.03978820811974433586e+046, 1.28605024825499153583e+047, 8.15915283247897734345e+047, 5.20850350543271572014e+048, 3.34525266131638071081e+049, 2.16152895475457702386e+050, 1.40500611775287989854e+051, 9.18649805770695235140e+051, 6.04152630633738356373e+052, 3.99612665510252427286e+053, 2.65827157478844876804e+054, 1.77827636152062330142e+055, 1.19622220865480194561e+056, 8.09115744491883602147e+056, 5.50262215981208894985e+057, 3.76238821188725874998e+058, 2.58623241511168180642e+059, 1.78713440064644790624e+060, 1.24139155925360726708e+061, 8.66760184313527234528e+061};
	return 0.5e0*pow(a,-0.5*(1+n))*GammaSemiInt[n+1];
}

double InteExpErfrn(double a, double b, int n)
{
    double result = 0.0;
    
    switch(n)
    {
        case 1:
            result = b / (2.0 * a * sqrt(a + b * b));
            break;
            
        case 3:
            result = (3.0 * a * b + 2.0 * pow(b, 3)) / (4.0 * pow(a, 2) * pow(a + b * b, 1.5));
            break;
            
        case 5:
            result = (15.0 * pow(a, 2) * b + 20.0 * a * pow(b, 3) + 8.0 * pow(b, 5)) / 
                     (8.0 * pow(a, 3) * pow(a + b * b, 2.5));
            break;
            
        case 7:
            result = (3.0 * b * (35.0 * pow(a, 3) + 70.0 * pow(a, 2) * pow(b, 2) + 
                     56.0 * a * pow(b, 4) + 16.0 * pow(b, 6))) / 
                     (16.0 * pow(a, 4) * pow(a + b * b, 3.5));
            break;
            
        case 9:
            result = (3.0 * b * (315.0 * pow(a, 4) + 840.0 * pow(a, 3) * pow(b, 2) + 
                     1008.0 * pow(a, 2) * pow(b, 4) + 576.0 * a * pow(b, 6) + 128.0 * pow(b, 8))) / 
                     (32.0 * pow(a, 5) * pow(a + b * b, 4.5));
            break;
            
        default:
            printf("only n=1,3,5,7,9\n");
            result = 0.0;
            break;
    }
    
    return result;
}

typedef double (*inteCenV)(double b11,int n,basis_qnum qf,basis_qnum qi,vargs varg,int c);

double inteVcenPartA(txrp_vtype get_trxp_vtype,sumckdk_scdk scdk,basis_qnum qf,basis_qnum qi,vargs varg,inteCenV iv,int c)
{
	int i,j;
	double p,sum;
	double t[21];
	int lent;

	get_trxp_vtype(qf,qi,t,&lent,c);

	sum=0;
	for(i=0;i<scdk.len;i++)
	{
		p=scdk.coe[i];
		for(j=0;j<lent;j++)
		{
			p*=pow(t[j],scdk.nfgh[i][j]);
		}
		p*=t[19]*t[20]*iv(t[18],scdk.nfgh[i][18],qf,qi,varg,c);
		sum+=p;
	}
	return sum;
}


double inteNfi(double b11,int n,basis_qnum qf,basis_qnum qi,vargs varg,int c)
{
	return InteExprn(b11,n);
}

double inteVogeG(double b11,int n,basis_qnum qf,basis_qnum qi,vargs varg,int c)
{
	double sigmakij[3];
	int k;
	double mi,mj,mk;
	double res=0;
	getijk(qf.m1,qf.m2,qf.m3,&mi,&mj,&mk,c);
	getsigmakij(sigmakij,varg.gamma,varg.sigma0,varg.s,mi,mj);
	for(k=0;k<3;k++)
	{
		res=res-2*varg.alpha[k]*InteExpErfrn(b11,sigmakij[k],n-1)/3.0;
	}
	return res;
}

double inteVcont(double b11,int n,basis_qnum qf,basis_qnum qi,vargs varg,int c)
{
	double pi=3.1415926535897932385;
	double sigmakij[3];
	int k;
	double mi,mj,mk;
	double res=0;
	getijk(qf.m1,qf.m2,qf.m3,&mi,&mj,&mk,c);
	getsigmakij(sigmakij,varg.gamma,varg.sigma0,varg.s,mi,mj);
	for(k=0;k<3;k++)
	{
		res=res+16*varg.alpha[k]*pow(sigmakij[k],3)*InteExprn(pow(sigmakij[k],2)+b11,n)/(9*mi*mj*sqrt(pi));
	}
	return res;
}

double inteVtens(double b11,int n,basis_qnum qf,basis_qnum qi,vargs varg,int c)
{
	double pi=3.1415926535897932385;
	double sigmakij[3];
	int k;
	double mi,mj,mk;
	double res=0;
	getijk(qf.m1,qf.m2,qf.m3,&mi,&mj,&mk,c);
	getsigmakij(sigmakij,varg.gamma,varg.sigma0,varg.s,mi,mj);
	for(k=0;k<3;k++)
	{
		res=res-4*varg.alpha[k]*sigmakij[k]*InteExprn(pow(sigmakij[k],2)+b11,n-2)/(3*mi*mj*sqrt(pi))-8*varg.alpha[k]*pow(sigmakij[k],3)*InteExprn(pow(sigmakij[k],2)+b11,n)/(9*mi*mj*sqrt(pi))+2*varg.alpha[k]*InteExpErfrn(b11,sigmakij[k],n-3)/(3*mi*mj);

	}
	return res;
}

double inteVsovii(double b11,int n,basis_qnum qf,basis_qnum qi,vargs varg,int c)
{
	double pi=3.1415926535897932385;
	double sigmakij[3];
	int k;
	double mi,mj,mk;
	double res=0;
	getijk(qf.m1,qf.m2,qf.m3,&mi,&mj,&mk,c);
	getsigmakij(sigmakij,varg.gamma,varg.sigma0,varg.s,mi,mj);
	for(k=0;k<3;k++)
	{
		res=res-2*varg.alpha[k]*sigmakij[k]*InteExprn(pow(sigmakij[k],2)+b11,n-2)/(3*sqrt(pi)*mi*mi)+varg.alpha[k]*InteExpErfrn(b11,sigmakij[k],n-3)/(3*mi*mi);
	}
	return res;
}

double inteVsovjj(double b11,int n,basis_qnum qf,basis_qnum qi,vargs varg,int c)
{
	double pi=3.1415926535897932385;
	double sigmakij[3];
	int k;
	double mi,mj,mk;
	double res=0;
	getijk(qf.m1,qf.m2,qf.m3,&mi,&mj,&mk,c);
	getsigmakij(sigmakij,varg.gamma,varg.sigma0,varg.s,mi,mj);
	for(k=0;k<3;k++)
	{
		res=res+2*varg.alpha[k]*sigmakij[k]*InteExprn(pow(sigmakij[k],2)+b11,n-2)/(3*sqrt(pi)*mj*mj)-varg.alpha[k]*InteExpErfrn(b11,sigmakij[k],n-3)/(3*mj*mj);
	}
	return res;
}

double inteVsovij(double b11,int n,basis_qnum qf,basis_qnum qi,vargs varg,int c)
{
	double pi=3.1415926535897932385;
	double sigmakij[3];
	int k;
	double mi,mj,mk;
	double res=0;
	getijk(qf.m1,qf.m2,qf.m3,&mi,&mj,&mk,c);
	getsigmakij(sigmakij,varg.gamma,varg.sigma0,varg.s,mi,mj);
	for(k=0;k<3;k++)
	{
		res=res-4*varg.alpha[k]*sigmakij[k]*InteExprn(pow(sigmakij[k],2)+b11,n-2)/(3*sqrt(pi)*mi*mj)+2*varg.alpha[k]*InteExpErfrn(b11,sigmakij[k],n-3)/(3*mi*mj);
	}
	return res;
}

double inteVsovji(double b11,int n,basis_qnum qf,basis_qnum qi,vargs varg,int c)
{
	double pi=3.1415926535897932385;
	double sigmakij[3];
	int k;
	double mi,mj,mk;
	double res=0;
	getijk(qf.m1,qf.m2,qf.m3,&mi,&mj,&mk,c);
	getsigmakij(sigmakij,varg.gamma,varg.sigma0,varg.s,mi,mj);
	for(k=0;k<3;k++)
	{
		res=res+4*varg.alpha[k]*sigmakij[k]*InteExprn(pow(sigmakij[k],2)+b11,n-2)/(3*sqrt(pi)*mi*mj)-2*varg.alpha[k]*InteExpErfrn(b11,sigmakij[k],n-3)/(3*mi*mj);
	}
	return res;
}

double inteVstring(double b11,int n,basis_qnum qf,basis_qnum qi,vargs varg,int c)
{
	double pi=3.1415926535897932385;
	double sigmaij;
	double mi,mj,mk;
	getijk(qf.m1,qf.m2,qf.m3,&mi,&mj,&mk,c);
	sigmaij=getsigmaij(varg.sigma0,varg.s,mi,mj);
	return varg.b*varg.f*InteExprn(pow(sigmaij,2)+b11,n)/(sqrt(pi)*sigmaij)+varg.b*varg.f*InteExpErfrn(b11,sigmaij,n+1)+varg.b*varg.f*InteExpErfrn(b11,sigmaij,n-1)/(2*pow(sigmaij,2))+varg.c*InteExprn(b11,n)/2;
}

double inteVsosii(double b11,int n,basis_qnum qf,basis_qnum qi,vargs varg,int c)
{
	double pi=3.1415926535897932385;
	double sigmaij;
	double mi,mj,mk;
	getijk(qf.m1,qf.m2,qf.m3,&mi,&mj,&mk,c);
	sigmaij=getsigmaij(varg.sigma0,varg.s,mi,mj);
	return -varg.b*varg.f*InteExprn(pow(sigmaij,2)+b11,n-2)/(2*sqrt(pi)*sigmaij*pow(mi,2))-varg.b*varg.f*InteExpErfrn(b11,sigmaij,n-1)/(2*pow(mi,2))+varg.b*varg.f*InteExpErfrn(b11,sigmaij,n-3)/(4*pow(sigmaij,2)*pow(mi,2));
}

double inteVsosjj(double b11,int n,basis_qnum qf,basis_qnum qi,vargs varg,int c)
{
	double pi=3.1415926535897932385;
	double sigmaij;
	double mi,mj,mk;
	getijk(qf.m1,qf.m2,qf.m3,&mi,&mj,&mk,c);
	sigmaij=getsigmaij(varg.sigma0,varg.s,mi,mj);
	return  varg.b*varg.f*InteExprn(pow(sigmaij,2)+b11,n-2)/(2*sqrt(pi)*sigmaij*pow(mj,2))+varg.b*varg.f*InteExpErfrn(b11,sigmaij,n-1)/(2*pow(mj,2))-varg.b*varg.f*InteExpErfrn(b11,sigmaij,n-3)/(4*pow(sigmaij,2)*pow(mj,2));
}

double betaij(double r,double mi,double mj,double epsilon)
{
	return pow(1+r*r/(sqrt(r*r+mi*mi)*sqrt(r*r+mj*mj)),0.5+epsilon);
}

double deltaij(double r,double mi,double mj,double epsilon)
{
	return pow(mi*mj/(sqrt(r*r+mi*mi)*sqrt(r*r+mj*mj)),0.5+epsilon);
}

double pogeG(double r,basis_qnum qf,basis_qnum qi,vargs varg,int c)
{
	double mi,mj,mk;
	getijk(qf.m1,qf.m2,qf.m3,&mi,&mj,&mk,c);
	return betaij(r,mi,mj,varg.eCoul);
}

double pcont(double r,basis_qnum qf,basis_qnum qi,vargs varg,int c)
{
	double mi,mj,mk;
	getijk(qf.m1,qf.m2,qf.m3,&mi,&mj,&mk,c);
	return deltaij(r,mi,mj,varg.econt);
}

double ptens(double r,basis_qnum qf,basis_qnum qi,vargs varg,int c)
{
	double mi,mj,mk;
	getijk(qf.m1,qf.m2,qf.m3,&mi,&mj,&mk,c);
	return deltaij(r,mi,mj,varg.etens);
}

double psovii(double r,basis_qnum qf,basis_qnum qi,vargs varg,int c)
{
	double mi,mj,mk;
	getijk(qf.m1,qf.m2,qf.m3,&mi,&mj,&mk,c);
	return deltaij(r,mi,mi,varg.esov);
}

double psovjj(double r,basis_qnum qf,basis_qnum qi,vargs varg,int c)
{
	double mi,mj,mk;
	getijk(qf.m1,qf.m2,qf.m3,&mi,&mj,&mk,c);
	return deltaij(r,mj,mj,varg.esov);
}

double psovji(double r,basis_qnum qf,basis_qnum qi,vargs varg,int c)
{
	double mi,mj,mk;
	getijk(qf.m1,qf.m2,qf.m3,&mi,&mj,&mk,c);
	return deltaij(r,mj,mi,varg.esov);
}

double psovij(double r,basis_qnum qf,basis_qnum qi,vargs varg,int c)
{
	double mi,mj,mk;
	getijk(qf.m1,qf.m2,qf.m3,&mi,&mj,&mk,c);
	return deltaij(r,mi,mj,varg.esov);
}

double psosii(double r,basis_qnum qf,basis_qnum qi,vargs varg,int c)
{
	double mi,mj,mk;
	getijk(qf.m1,qf.m2,qf.m3,&mi,&mj,&mk,c);
	return deltaij(r,mi,mi,varg.esos);
}

double psosjj(double r,basis_qnum qf,basis_qnum qi,vargs varg,int c)
{
	double mi,mj,mk;
	getijk(qf.m1,qf.m2,qf.m3,&mi,&mj,&mk,c);
	return deltaij(r,mj,mj,varg.esos);
}

double intepogeG(double b11,int n,basis_qnum qf,basis_qnum qi,vargs varg,int c)
{
	return InteCenV(b11,n,qf,qi,pogeG,varg,c);
}

double intepcont(double b11,int n,basis_qnum qf,basis_qnum qi,vargs varg,int c)
{
	return InteCenV(b11,n,qf,qi,pcont,varg,c);
}

double inteptens(double b11,int n,basis_qnum qf,basis_qnum qi,vargs varg,int c)
{
	return InteCenV(b11,n,qf,qi,ptens,varg,c);
}

double intepsovii(double b11,int n,basis_qnum qf,basis_qnum qi,vargs varg,int c)
{
	return InteCenV(b11,n,qf,qi,psovii,varg,c);
}

double intepsovjj(double b11,int n,basis_qnum qf,basis_qnum qi,vargs varg,int c)
{
	return InteCenV(b11,n,qf,qi,psovjj,varg,c);
}

double intepsovji(double b11,int n,basis_qnum qf,basis_qnum qi,vargs varg,int c)
{
	return InteCenV(b11,n,qf,qi,psovji,varg,c);
}

double intepsovij(double b11,int n,basis_qnum qf,basis_qnum qi,vargs varg,int c)
{
	return InteCenV(b11,n,qf,qi,psovij,varg,c);
}

double intepsosii(double b11,int n,basis_qnum qf,basis_qnum qi,vargs varg,int c)
{
	return InteCenV(b11,n,qf,qi,psosii,varg,c);
}

double intepsosjj(double b11,int n,basis_qnum qf,basis_qnum qi,vargs varg,int c)
{
	return InteCenV(b11,n,qf,qi,psosjj,varg,c);
}

double inteTi(double b11,int n,basis_qnum qf,basis_qnum qi,vargs varg,int c)
{
	return InteCenV(b11,n,qf,qi,Ti,varg,c);
}


double inteRMS(double b11,int n,basis_qnum qf,basis_qnum qi,vargs varg,int c)
{
	return InteExprn(b11,n+2);
}

