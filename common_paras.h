	// parameters used by both the learning and the simulation routines
	const int N=3000;						// number of GCs
	const int M=100;						// number of MFs
	const int PP=1000;						// number of MF patterns for threshold adjustment
	const int MFSYN=4;						// (average) number of MF to GC synapses per GC
	const double PCsp=40.0;					// PC spontaneous activity

	// parameters associated with statistics of the MF input patterns
	double MF_avrg=91.632;
	double MF_std=90.815;
	double MF_var=MF_std*MF_std;

	double MF_avrg_G1=200.0;		// driver in DS model
	double MF_avrg_G2=200.0;
	double MF_avrg_G3=20.0;			// supporter in DS model
	double MF_avrg_G4=20.0;
	double MF_avrg_G5=20.0;

	double MF_std_G1=50.0;			// driver in DS model
	double MF_std_G2=50.0;
	double MF_std_G3=20.0;			// supporter in DS model
	double MF_std_G4=20.0;
	double MF_std_G5=20.0;

	const double pdriver=0.5;				// fraction of D synapses in DS model

	double pact[]={1,1,1,1,1};				// fraction of active MFs (for each group)

	// parameters associated with inhibition
	double cI=10;						// MLI-to-PC synaptic strength
	//~ double cI=1.025;					// MLI-to-PC synaptic strength
	//~ double cI=0.975;					// MLI-to-PC synaptic strength
	double avrgJIE=1;						// GC-to-GoC synaptic strength
	double avrgJEI=1;						// GoC-to-GC synaptic strength

	double gain_global=1.6518;				// equal gain for all GCs for reccurent Golgi setting

	const double Tfac_gc=0.9;				// factor by which GC thresholds are to be reduced when GoC is added

	// parameters associated with GC properties
	const double nfrac=0.2;					// GC coding level
	const double GCtarget=5;				// GC target firing rate
	const double GCnormfac=409.431;			// GC normalisation factor, calibrated on Fig4a parameters
	const double theta_global=750.0;		// global thresholds

	// parameters associated with CF properties [ CFrate=CFsp + CFscale*(hPC-hPC_target) ]
	double CFsp=1;							// spontaneous CF firing rate
	const double CFsp_bayes=5;				// spontaneous CF firing rate for bayesian learning paradigm
	const double CFscale=0.5;				// scaling of the error signal component
	double CFwidth=0.025;					// STD (width) of CF template
	const double errWeight=3.5;				// weight of delay time for delta-learning rule

	// parameters for gradient descend
	const double alpha_learn=0.0025;
	double J2weight=0.0;

	// GC time-constants
	const double tauGC=0.01;

	const int Nsynpara_2pools=12;

	// desired rank correlation between MF- and syn-identities
	double MFUcorr=0.999999;

	const int Ncl=5;

	std::string MF_STP ("on"); 					// MF-to-GC STP: "on", "off", "identical"
	std::string transients ("on");				// transients: "on", "off"
	std::string gaincontrol ("individual");		// GC gain control ? "global", "individual","individual_active", "no"

	std::string threshold ("adjust");			// threshold: identical for all GCs: "fixed"
												// adjusted for each GC: "adjust"
												// adjusted globally: "global"
												// renormalise synapses instead: "renorm"

	std::string normGCs ("off");				// normalise GC activity to [0,1]? "on", "off"

	std::string Groups ("2SYN_2MF");			// 5 MF types and 5 synapse types with corr=1: "5SYN_5MF"						(full model)
												// 2 MF types and 2 synapse types with corr=MFUcorr (DS model): "2SYN_2MF"  	(DS model)
												// MF distrb. with 5 peaks and Ncl MF/SYN types with corr=MFUcorr: "TOY_5MF"	(DS model)

	std::string Driver ("always");				// should every GC receive a driving input? "always", "random"

	std::string Udistrb ("no");					// rel. prob. set to average: "no"
												// mixture of beta distributions with #SYNclusters = #MFclusters: "beta_mixture"
												// single uniform distribution with #SYNclusters = #MFclusters : "uniform"
												// single beta distribution with #SYNclusters = #MFclusters : "beta"
												// use this while setting Groups='2SYN_2MF', '5SYN_5MF', etc. and set correlation bweteen MF- and SYN-identities via 'MFUcorr'

	std::string pattern ("Gauss");		// pattern type: "Gauss", "trunc_Gauss", "lognorm", "exponential", "gamma", "uniform_multi", "uniform_single" (= seamleass uniform)

	std::string Golgi ("on");					// GoC recurrent feedback: "on", "off"

	std::string Tsignal ("delta"); 				// which eye-lid teaching signal ? "Gauss", "delta"
	std::string cuttsyn ("no"); 				// cut away synapses ? "no","bounds","cutS","cutD","cutSLOW","cutFAST","cutD_AND_cutFAST"
	std::string cutGC ("no"); 					// cut away GC transients ? "no","above","below" [DOES NOT WORK WITH RECURRENT GOLGI]

	const double tsyncut_hi=5;				// upper bound for tsyn cut	(tsyncut_hi=+inf: cut everything,  tsyncut_hi=0: cut nothing)
	const double tsyncut_lo=0.125;					// lower bound for tsyn cut	(tsyncut_lo=0: cut everything,  tsyncut_hi=+inf: cut nothing)

	const double GCcut=0.15;						// upper or lower bound for GC transient cut
	//~ const double GCcut=0.175;						// upper or lower bound for GC transient cut
	const double GCcut_onset=0.02;				// lower bound for GC transient onsets
	const double tlate=0.3;						// disregard 'peaks' later than this
	const double ampsmall=1;					// disregard 'peaks' with amplitude smaller than this

	double Tpre=0.1;							// start of signal [s] // time interval before pattern switch in seconds / for eyeblink learning loop

	int ptarget=1;								// substrate pattern for sequence/eyeblink learning
	int pstart=0;								// context pattern for cancellation learning

	const double pi=3.14159265;