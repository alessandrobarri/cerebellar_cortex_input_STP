#include <gsl/gsl_roots.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_linalg.h>
#include <armadillo>
#include "headers/rtnorm.hpp"
#include <sys/stat.h>
using namespace std;

struct presyn {
	int * idx;
	double * strgth;
	int num;
};
struct presyn_2pools {
	int * idx;	int * Gidx;
	double * strgth1; double * strgth2;
	double ** x1p; double ** x2p;
	double * x1; double * u1; double * x2;
	double * x1m; double * x2m;
	double * U1; double * tauD1; double * pf1;
	double * U2; double * tauD2; double * pf2;
	double * alpha1;	double * WU1; double * Upf1;
	double * alpha2;	double * WU2; double * Upf2;
	double ** tsyn1; double ** tsyn2;
	int num;
};
struct cluster_struct {
	std::vector <int> cluster;
	std::vector <int> synID;
	int idxC;
};
struct flags {
	int MF_STP,trans,trans_dum,gain,threshold,pattern,groups;
	int GoC,driver,Udistrb,mode,cutGC,cuttsyn;
};
struct netpara {
	int N,M,MFSYN,NSYN,Ncl,PP;
	double pdriver,MFUcorr,nfrac,GCtarget,theta_global;
	double avrgJEI;
};
struct synpara {
	const double U_dev_slow;

	arma::vec NFAST;
	arma::vec NSLOW;
	arma::vec USLOW;
};

#include "headers/sort_indices.h"
#include "headers/draw_patterns.h"
#include "headers/threshNmom.h"
#include "headers/find_teffMFstats.h"
#include "headers/generate_MF_groups_DS.h"
#include "headers/generate_MF_patterns_DS.h"
#include "headers/generate_weights_DS.h"
#include "headers/generate_relprob_DS.h"
#include "headers/generate_ratePv_correlations.h"
#include "headers/calc_GCpattern_theta_gain_DS.h"
#include "headers/learn_weights_DS.h"
#include "headers/rtnorm.cpp"

int main (int argc, char* argv[]) {

	#include "common_paras.h"
	#include "syn_rate_param_reduced.h"

// 	random variable seed
	int idum= -4;				// default seed
	int Npara=10;				// default number of different parameter sets
	int Npara2=1;				// default number of parameters on second dimension
	int Nreal=2;				// default number of network realisations per parameter sets
	int Ntrials=4000;			// default number of stimulus presentations per network realisation
	string path = "output/";
	string folder,mode;

	if (argc==7){
		// get folder name from command line input
		folder=argv[1];
		// get mode string from command line input
		mode=argv[2];
		// get number of parameter sets
		Npara=atoi(argv[3]);
		// get number of network realisations
		Nreal=atoi(argv[4]);
		// get number of stimulus presentations
		Ntrials=atoi(argv[5]);
		// get random number seed
		idum=atoi(argv[6]);
	}
	else {
		std::cerr << std::endl;
		std::cerr << "1) Please specify folder name." << std::endl;
		std::cerr << "2) Please specify mode string." << std::endl;
		std::cerr << "3) Please choose number of points in MF rate space." << std::endl;
		std::cerr << "4) Please choose number of network realisations." << std::endl;
		std::cerr << "5) Please choose number of learning steps." << std::endl;
		std::cerr << "6) Please choose random number seed (integer)." << std::endl;
		std::cerr << std::endl;
		return(0);
	}
	Ntrials=Ntrials+1;

	// **********************************
	//	set options
	struct flags FLAGS;

    if 		(MF_STP.compare("on") == 0) 		FLAGS.MF_STP=1;
    else if (MF_STP.compare("off") == 0) 		FLAGS.MF_STP=0;
    else if (MF_STP.compare("identical") == 0) 	FLAGS.MF_STP=2;
    else {std::cerr <<"Bad MF_STP string. \n";exit(EXIT_FAILURE);}

    if 		(threshold.compare("adjust") == 0) 	FLAGS.threshold=1;
    else if (threshold.compare("fixed") == 0) 	FLAGS.threshold=0;
    else if (threshold.compare("renorm") == 0) 	FLAGS.threshold=2;
    else if (threshold.compare("global") == 0) 	FLAGS.threshold=4;
    else {std::cerr <<"Bad threshold string. \n";exit(EXIT_FAILURE);}

	if (MFpattern.compare("Gauss") == 0) 				FLAGS.pattern=0;
    else if (MFpattern.compare("lognorm") == 0) 		FLAGS.pattern=1;
    else if (MFpattern.compare("exponential") == 0) 	FLAGS.pattern=2;
    else if (MFpattern.compare("gamma") == 0) 		FLAGS.pattern=3;
    else if (MFpattern.compare("uniform_multi") == 0) FLAGS.pattern=4;
    else if (MFpattern.compare("uniform_single") == 0)FLAGS.pattern=5;
    else if (MFpattern.compare("trunc_Gauss") == 0)	FLAGS.pattern=6;
    else {std::cerr <<"Bad MFpattern string. \n";exit(EXIT_FAILURE);}

    if (gaincontrol.compare("global") == 0) 				FLAGS.gain=0;
    else if (gaincontrol.compare("individual") == 0) 		FLAGS.gain=1;
    else if (gaincontrol.compare("individual_active") == 0) FLAGS.gain=2;
    else if (gaincontrol.compare("no") == 0) 				FLAGS.gain=3;
    else {std::cerr <<"Bad gaincontrol string. \n";exit(EXIT_FAILURE);}

    if (Groups.compare("2SYN_2MF") == 0)	 		FLAGS.groups=5;
    else if (Groups.compare("TOY_5MF") == 0)		FLAGS.groups=55;
    else {std::cerr <<"Bad Groups string: "<< Groups <<" \n";exit(EXIT_FAILURE);}

    if 		(Driver.compare("always") == 0) 	FLAGS.driver=1;
    else if (Driver.compare("random") == 0) 	FLAGS.driver=0;
    else {std::cerr <<"Bad Driver string. \n";exit(EXIT_FAILURE);}

    if 		(transients.compare("on") == 0) 	FLAGS.trans=1;
    else if (transients.compare("off") == 0) 	FLAGS.trans=0;
    else {std::cerr <<"Bad transients string. \n";exit(EXIT_FAILURE);}

	if (Udistrb.compare("uniform") == 0) 			FLAGS.Udistrb=2;
    else if (Udistrb.compare("no") == 0)			FLAGS.Udistrb=0;
    else {std::cerr <<"Bad Udistrb string. \n";exit(EXIT_FAILURE);}

    if (mode.compare("MF_U_corr") == 0)		FLAGS.mode=6;
    else if (mode.compare("cI") == 0)		FLAGS.mode=7;
    else if (mode.compare("Ntrials") == 0)	FLAGS.mode=8;
    else if (mode.compare("bayesian") == 0)	FLAGS.mode=77;
    else if (mode.compare("scan") == 0)		FLAGS.mode=88;
    else if (mode.compare("sample") == 0)	FLAGS.mode=99;
    else {std::cerr <<"Bad mode string. \n";exit(EXIT_FAILURE);}

    if 		(Golgi.compare("on") == 0) 		FLAGS.GoC=1;
    else if (Golgi.compare("off") == 0) 	FLAGS.GoC=0;
    else {std::cerr <<"Bad Golgi string. \n";exit(EXIT_FAILURE);}

    if (cuttsyn.compare("no") == 0)				FLAGS.cuttsyn=0;
    else if (cuttsyn.compare("cutS") == 0)		FLAGS.cuttsyn=2;
    else if (cuttsyn.compare("cutD") == 0)		FLAGS.cuttsyn=3;
    else if (cuttsyn.compare("cutSLOW") == 0)	FLAGS.cuttsyn=4;
    else if (cuttsyn.compare("cutFAST") == 0)	FLAGS.cuttsyn=5;
    else {std::cerr <<"Bad Cuttsyn string. \n";exit(EXIT_FAILURE);}

    if (cutGC.compare("no") == 0)				FLAGS.cutGC=0;
    else if (cutGC.compare("above") == 0) 		FLAGS.cutGC=1;
    else if (cutGC.compare("below") == 0)		FLAGS.cutGC=2;
    else {std::cerr <<"Bad cutGC string. \n";exit(EXIT_FAILURE);}

	// if reduced driver/supporter model is chosen, set driver flag correspondingly
	if(FLAGS.groups==5){
		FLAGS.driver=2;
	}

	if(FLAGS.cuttsyn==2){
		std::cout << std::endl;
		std::cout << "setting supporter synapses to steady state" << std::endl;
	}
	else if(FLAGS.cuttsyn==3){
		std::cout << std::endl;
		std::cout << "setting driver synapses to steady state" << std::endl;
	}
	else if(FLAGS.cuttsyn==4){
		std::cout << std::endl;
		std::cout << "setting slow pools to steady state" << std::endl;
	}
	else if(FLAGS.cuttsyn==5){
		std::cout << std::endl;
		std::cout << "setting fast pools to steady state" << std::endl;
	}
	std::cout << std::endl;

  	const gsl_rng_type * T;
  	gsl_rng * r;

	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc (T);

	gsl_rng_env_setup();
  	gsl_rng_set(r,idum);

	// **********************************
	// 	temporal variables
	double Tlearn=1.4;
	if(Ntrials<2){
		Tlearn=1.1;
	}
	// with GoC, use longer baseline time-interval
	if(FLAGS.GoC==1)	{
		Tpre=Tpre+0.2;
		std::cout<< std::endl;
		std::cout<< "Reccurent Golgi inputs turned on. Using Tpre = "<<Tpre << std::endl;
	}
	else{
		avrgJIE=0;	// otherwise turn off GC-GoC synapses
	}

	const double duration =Tlearn+Tpre;					// duration of total simulation in seconds
	const double signal =Tlearn+0.1;						// duration of relevant simulation in seconds
  	const double dt = 0.0005;
  	const int MAX_data = duration*1.0/dt;
	const double dbint=0.005;							// time step for signal parsing in seconds
	const int BINS=	signal/dbint;						// number of bins for signal learning
	const int ttpre=int(Tpre/dt);

	//	variables
	double err_final,JI;
	double Tavrg,Gavrg,CLavrg;
	double MF_avrg_CL[5],MF_std_CL[5];
	double a,b;
	int i,j,k,m,n,t,idxS,idxD,Uprobidx,inidx;
	int tt,kk,ns,np,tcut;
	bool bool_tfive,bool_tten,bool_tcut;

	const int NSYN=N*MFSYN;							// number of MFGC synapses
	const int NSYN_2=2*NSYN;						// double NSYN for identity generation

	double Mdata[5][Nsynpara_2pools];

	// time-window for GC transient cutting
	const int tcutfrac=10;
	const double tcut_start=Tpre;								// at onset
	const double tcut_end=1.2+Tpre;								// 1.2s after onset
	const int CUT_tsteps=(tcut_end-tcut_start)/(dt*tcutfrac);	// *tcutfrac : downsampling by factor tcutfrac

	const int ttcut_start=int(tcut_start/dt);
	const int ttcut_end=int(tcut_end/dt);

	//	usefull quantities
	const double sqN=sqrt(N);
	const double dtGC=dt/tauGC;

	// learning parameters
	double alpha=alpha_learn/sqN;

	// **********************************
	// structs / arrays
	struct presyn_2pools W[N];
	struct presyn J;
	struct presyn JEI[N];			// GoC-to-GC
	struct cluster_struct CL[5];
	struct netpara NETPARA;
	J.strgth=new double [N];
	double theta[N];
	std::vector <double> dummyD;
	std::vector <size_t> idx_sorted;
	std::vector <int> drivers,supporters;
	std::vector <int> D_synID,S_synID;
	int MFgroup[M];

	double Utoy_slow[NSYN_2];	double Ntoy_fast[NSYN_2];

	arma::mat MFpatterns_adjust(M,PP);
	arma::mat MFpatterns_preCS(M,PP);
	arma::mat MFpatterns_CS(M,PP);
	arma::mat GCpatterns(N,PP);

	arma::vec Time(CUT_tsteps);
	arma::colvec GCmax(N);

	double PC[BINS];
	double GC[N],GCm[N],GCgain[N],hGC[N];
	double MF[M];
	double GoC,GoCm,GCavrg,hGCavrg,hGCvar,MLI;

	double** MF_Tinput=new double* [M];
	for(i= 0; i < M; i++){
		MF_Tinput[i]=new double [MAX_data+1];
	}

	for(i= 0; i < N; i++) {
		JEI[i].strgth=new double [1];
	}

	// **********************************
	// 	prepare output files
	// add 'std::ios::app' for adding new data to existing files instead of clearing content
	path+=folder;

	struct stat sb;
	if (stat(path.c_str(), &sb) == 0 && S_ISDIR(sb.st_mode)) {
		std::cout << "[Folder exists] "<<"\n";
    }
    else{
		string dircommand=string("mkdir ") + path;
		int systemRet = system(dircommand.c_str());
		if(systemRet == -1){
			std::cerr <<"Error: Folder could not be created."<< "\n";
			exit(EXIT_FAILURE);
		}
		else std::cout << "[Creating folder]" <<"\n";
	}

	std::ofstream GC_file(path+string("GC.dat"));						if (!GC_file) exit(EXIT_FAILURE);
	std::ofstream MF_file(path+string("MF.dat"));						if (!MF_file) exit(EXIT_FAILURE);
	std::ofstream STP_file(path+string("MFGC_STP.dat"));				if (!STP_file) exit(EXIT_FAILURE);
	std::ofstream MLI_file(path+string("MLI.dat"));						if (!MLI_file) exit(EXIT_FAILURE);
	std::ofstream GoC_file(path+string("GoC.dat"));						if (!GoC_file) exit(EXIT_FAILURE);
	std::ofstream CF_file(path+string("CF.dat"));						if (!CF_file) exit(EXIT_FAILURE);
	std::ofstream time_file(path+string("time.dat"));					if (!time_file) exit(EXIT_FAILURE);
	std::ofstream time_learn_file(path+string("time_learn.dat"));		if (!time_learn_file) exit(EXIT_FAILURE);
	std::ofstream STPpara_file(path+string("STPpara.dat"));				if (!STPpara_file) exit(EXIT_FAILURE);
	std::ofstream GCtrans_file(path+string("GCtrans.dat"));				if (!GCtrans_file) exit(EXIT_FAILURE);
	std::ofstream timeGCtrans_file(path+string("timeGCtrans.dat"));		if (!timeGCtrans_file) exit(EXIT_FAILURE);

	std::ofstream err_final_file(path+string("err_final.dat"));			if (!err_final_file) exit(EXIT_FAILURE);
	std::ofstream PC_025_file(path+string("PC_final_025.dat"));			if (!PC_025_file) exit(EXIT_FAILURE);
	std::ofstream PC_050_file(path+string("PC_final_050.dat"));			if (!PC_050_file) exit(EXIT_FAILURE);
	std::ofstream PC_100_file(path+string("PC_final_100.dat"));			if (!PC_100_file) exit(EXIT_FAILURE);
	std::ofstream PC_200_file(path+string("PC_final_200.dat"));			if (!PC_200_file) exit(EXIT_FAILURE);
	std::ofstream PC_300_file(path+string("PC_final_300.dat"));			if (!PC_300_file) exit(EXIT_FAILURE);
	std::ofstream PC_500_file(path+string("PC_final_500.dat"));			if (!PC_500_file) exit(EXIT_FAILURE);
	std::ofstream PC_700_file(path+string("PC_final_700.dat"));			if (!PC_700_file) exit(EXIT_FAILURE);

	std::ofstream PC_bayes_1_file(path+string("PC_bayes_025150.dat"));	if (!PC_bayes_1_file) exit(EXIT_FAILURE);
	std::ofstream PC_bayes_2_file(path+string("PC_bayes_050200.dat"));	if (!PC_bayes_2_file) exit(EXIT_FAILURE);
	std::ofstream PC_bayes_3_file(path+string("PC_bayes_100300.dat"));	if (!PC_bayes_3_file) exit(EXIT_FAILURE);
	std::ofstream PC_bayes_4_file(path+string("PC_bayes_200400.dat"));	if (!PC_bayes_4_file) exit(EXIT_FAILURE);
	std::ofstream PC_bayes_5_file(path+string("PC_bayes_300500.dat"));	if (!PC_bayes_5_file) exit(EXIT_FAILURE);

	std::ofstream htarget_BINS_file(path+string("htarget_BINS.dat"));	if (!htarget_BINS_file) exit(EXIT_FAILURE);
	std::ofstream theta_file(path+string("theta.dat"));						if (!theta_file) exit(EXIT_FAILURE);
	std::ofstream J_file(path+string("J.dat"));							if (!J_file) exit(EXIT_FAILURE);
	std::ofstream Jidx_file(path+string("Jidx.dat"));					if (!Jidx_file) exit(EXIT_FAILURE);
	std::ofstream learnparam_file(path+string("learnparam.dat"));		if (!learnparam_file) exit(EXIT_FAILURE);
	std::ofstream stats_file(path+string("stats.dat"));					if (!stats_file) exit(EXIT_FAILURE);
	std::ofstream flags_file(path+string("flags.dat"));					if (!flags_file) exit(EXIT_FAILURE);


	flags_file << 	FLAGS.MF_STP << "\t"<< FLAGS.threshold << "\t"<< FLAGS.pattern << "\t"<< FLAGS.gain << "\t"<< 0 << "\t"<< FLAGS.groups << "\t"<< FLAGS.driver<< "\t"<<
					FLAGS.trans<< "\t"<< FLAGS.Udistrb<< "\t"<< FLAGS.mode << "\t"<< FLAGS.GoC << "\t"<<  1 << "\t"<< 0 << "\t"<<  FLAGS.cutGC <<std::endl;

	arma::vec MF_U_CORR = arma::linspace<arma::vec>(-0.99, 0.99, Npara);
	arma::vec CI_ARR = arma::logspace<arma::vec>(-1, 3, Npara);
	arma::vec NTRIALS_ARR = { 500, 1000, 2000, 4000, 8000, 16000, 32000};
	NETPARA.N=N;
	NETPARA.M=M;
	NETPARA.MFSYN=MFSYN;
	NETPARA.NSYN=NSYN;
	NETPARA.Ncl=Ncl;
	NETPARA.PP=PP;
	NETPARA.pdriver=pdriver;
	NETPARA.MFUcorr=MFUcorr;
	NETPARA.nfrac=nfrac;
	NETPARA.GCtarget=GCtarget;
	NETPARA.theta_global=theta_global;

	if(FLAGS.mode==6){
		std::cout << "Varying nu-p correlation" << "\n";
	}
	else if(FLAGS.mode==7){
		std::cout << "Varying JI" << "\n";
		Npara=CI_ARR.n_elem;
	}
	else if(FLAGS.mode==8){
		std::cout << "Varying Ntrials" << "\n";
		Npara=7;
	}
	else if(FLAGS.mode==77){
		std::cout << "Bayesian interval learning" << "\n";
		//~ Npara=1;	// only one parameter set when in sample mode
		CFsp=CFsp_bayes;
		Ndelays=5;
	}
	else if(FLAGS.mode==99){
		std::cout << "Sample eye-blink learning" << "\n";
		Npara=1;	// only one parameter set when in sample mode
	}
	else if(FLAGS.mode==88){
		std::cout << "Scan over 2D parameter space" << "\n";
		Npara2=Npara;	// scan over a square grid
	}
	// **********************************
	// start Npara loop
	for (int pp = 0; pp < Npara; pp++){

		// loop over different rate-pr identity correlations
		if(FLAGS.mode==6) {
			MFUcorr=MF_U_CORR(pp);
			NETPARA.MFUcorr=MFUcorr;
		}
		// loop over different inhibition strengths
		if(FLAGS.mode==7) {
			cI=CI_ARR(pp);
		}
		if(FLAGS.mode==8) Ntrials=NTRIALS_ARR(pp);
		JI=cI;

	for (int pp2 = 0; pp2 < Npara2; pp2++){

		if(FLAGS.mode==88){
			MF_avrg_D= (nuDmax-nuDmin)/double(Npara) *pp+nuDmin;			// check pp, pp2 !!
			MF_std_D= (sigDmax-sigDmin)/double(Npara) *pp2+sigDmin;
			//~ MF_avrg_S= (nuSmax-nuSmin)/double(Npara) *pp+nuSmin;			// check pp, pp2 !!
			//~ MF_std_S= (sigSmax-sigSmin)/double(Npara) *pp2+sigSmin;
		}

		std::cout << "********************* \n";
		std::cout << "Parameter set number " << pp*Npara2+pp2+1<< "\n";
		if(FLAGS.groups==5){
			std::cout << "nuD = " << MF_avrg_D  << "\t"<< "sigD = " << MF_std_D << "\t"<< "nuS = " << MF_avrg_S<< "\t"<< "sigS = " << MF_std_S << "\t"<< "MFUcorr = " << MFUcorr << "\t "<< "JI = " << JI << "\t"<< "number of steps = " << Ntrials << "\n";
		}
		else{
			if(FLAGS.groups==55){
				if(FLAGS.pattern==5){
					for (i = 0; i < Ncl; i++) {
						std::cout 	<< "nu0 = " << round(100* 0.5*(NUlims[Ncl-i]+NUlims[Ncl-i-1]) )/100.0  << "\t"<< "sig0 = " << round(100* (NUlims[Ncl-i]-NUlims[Ncl-i-1])/sqrt(12) )/100.0 << "\t";
					}
					std::cout	<< "MFUcorr = " << MFUcorr << "\t"<< "JI = " << cI << "\n";
				}
				else{
					std::cout 	<< "nu0 = " << nu0  << "\t"<< "sig0 = " << sig0 << "\t"
							<< "nu1 = " << nu1	<< "\t"<< "sig1 = " << sig1  << "\t"
							<< "nu2 = " << nu2	<< "\t"<< "sig2 = " << sig2 << "\t"
							<< "nu3 = " << nu3	<< "\t"<< "sig3 = " << sig3 << "\t"
							<< "nu4 = " << nu4	<< "\t"<< "sig4 = " << sig4 << "\t"
							<< "MFUcorr = " << MFUcorr << "\t"<< "JI = " << cI << "\n";
				}
			}
			else{
				std::cout << "nuD = " << MF_avrg_D  << "\t"<< "sigD = " << MF_std_D << "\t"<< "nuS = " << MF_avrg_S<< "\t"<< "sigS = " << MF_std_S << "\t"<< "MFUcorr = " << MFUcorr << "\t"<< "JI = " << cI << "\t"<< "number of steps = " << Ntrials << "\n";
			}
		}

	// start Nreal loop
	for (int rr = 0; rr < Nreal; rr++){
		if( ( (rr+1) % 5) ==0) std::cout << "Realisation number " << rr+1 << "\n";
		MF_avrg_CL[0]=nu0;		MF_std_CL[0]=sig0;
		MF_avrg_CL[1]=nu1;		MF_std_CL[1]=sig1;
		MF_avrg_CL[2]=nu2;		MF_std_CL[2]=sig2;
		MF_avrg_CL[3]=nu3;		MF_std_CL[3]=sig3;
		MF_avrg_CL[4]=nu4;		MF_std_CL[4]=sig4;

		for (k = 0; k < Ncl; k++){
			CL[k].cluster.clear();
			CL[k].synID.clear();
		}
		drivers.clear();
		supporters.clear();

	// **********************************
	//generate MF patterns
	// generate group identities
	generate_MF_groups_DS(MFgroup, CL, grouplims, drivers, supporters, NETPARA, FLAGS, r);

	// for threshold and gain adjustment
	generate_MF_patterns_DS(MFpatterns_adjust, MFgroup, MF_avrg_CL, MF_std_CL, CL, NUlims, grouplims, drivers, supporters, MF_avrg_D, MF_std_D, MF_avrg_S, MF_std_S, pact, NETPARA, FLAGS, r);

	// pre CS
	generate_MF_patterns_DS(MFpatterns_preCS, MFgroup, MF_avrg_CL, MF_std_CL, CL, NUlims, grouplims, drivers, supporters, MF_avrg_D, MF_std_D, MF_avrg_S, MF_std_S, pact, NETPARA, FLAGS, r);

	// during CS
	generate_MF_patterns_DS(MFpatterns_CS, MFgroup, MF_avrg_CL, MF_std_CL, CL, NUlims, grouplims, drivers, supporters, MF_avrg_D, MF_std_D, MF_avrg_S, MF_std_S, pact, NETPARA, FLAGS, r);

	// **********************************
	//	generate connectivity
	generate_weights_DS(W, JEI, CL, drivers, supporters, MFgroup, NETPARA, FLAGS, r);

	// -----------------------------------
	// generate uniform or beta release probability distribution
	// -----------------------------------
	if(FLAGS.Udistrb>1){
		for (i = 0; i < NSYN_2; i++) Utoy_slow[i]=gsl_ran_flat(r, 0.1, 0.9);
		for (i = 0; i < NSYN_2; i++) 	dummyD.push_back(Utoy_slow[i]);
		idx_sorted=sort_indices(dummyD);
		for (i = 0; i < NSYN_2; i++) 	Utoy_slow[i]=dummyD[idx_sorted[i]];
		dummyD.clear();
		idx_sorted.clear();
		for (i = 0; i < NSYN_2/5; i++)				Ntoy_fast[i]=NFAST[0];
		for (i = NSYN_2/5; i < NSYN_2*2/5; i++) 	Ntoy_fast[i]=NFAST[1];
		for (i = NSYN_2*2/5; i < NSYN_2*3/5; i++) 	Ntoy_fast[i]=NFAST[2];
		for (i = NSYN_2*3/5; i < NSYN_2*4/5; i++) 	Ntoy_fast[i]=NFAST[3];
		for (i = NSYN_2*4/5; i < NSYN_2; i++) 		Ntoy_fast[i]=NFAST[4];

	}
	// **********************************
	// 	generate correlation between MF- and synapse-types
	generate_ratePv_correlations(CL, D_synID, S_synID, grouplims, NETPARA, FLAGS, r);

	// -----------------------------------
	// prepare hand-made synapses
	// -----------------------------------
	if(FLAGS.groups==55) {
		// 5 groups toymodel
		for (i=0; i<5; i++) {
			// C0
			if(i==0){
				// first pool
				Mdata[i][0]=NSLOW[4];
				Mdata[i][1]=USLOW[4];
				Mdata[i][4]=PSLOW[4];
				// second pool
				Mdata[i][5]=NFAST[4];
				Mdata[i][6]=UFAST[4];
				Mdata[i][9]=PFAST[4];
			}
			// C1
			else if(i==1){
				// first pool
				Mdata[i][0]=NSLOW[3];
				Mdata[i][1]=USLOW[3];
				Mdata[i][4]=PSLOW[3];
				// second pool
				Mdata[i][5]=NFAST[3];
				Mdata[i][6]=UFAST[3];
				Mdata[i][9]=PFAST[3];
			}
			// C2
			else if(i==2){
				// first pool
				Mdata[i][0]=NSLOW[2];
				Mdata[i][1]=USLOW[2];
				Mdata[i][4]=PSLOW[2];
				// second pool
				Mdata[i][5]=NFAST[2];
				Mdata[i][6]=UFAST[2];
				Mdata[i][9]=PFAST[2];
			}
			// C3
			else if(i==3){
				// first pool
				Mdata[i][0]=NSLOW[1];
				Mdata[i][1]=USLOW[1];
				Mdata[i][4]=PSLOW[1];
				// second pool
				Mdata[i][5]=NFAST[1];
				Mdata[i][6]=UFAST[1];
				Mdata[i][9]=PFAST[1];
			}
			// C4
			else if(i==4){
				// first pool
				Mdata[i][0]=NSLOW[0];
				Mdata[i][1]=USLOW[0];
				Mdata[i][4]=PSLOW[0];
				// second pool
				Mdata[i][5]=NFAST[0];
				Mdata[i][6]=UFAST[0];
				Mdata[i][9]=PFAST[0];
			}
			Mdata[i][2]=Tau_slow;
			// facilitation first pool (turned off)
			Mdata[i][3]=0;

			Mdata[i][7]=Tau_fast;
			// facilitation second pool (turned off)
			Mdata[i][8]=0;
		}
	}
	else{
		// 2 groups
		for (i=0; i<=1; i++) {
			// i==0 : driver
			// i==1 : supporter
			//supporter
			if(i==1){
				// first pool (slow)
				Mdata[i][0]=Nslow_S;
				Mdata[i][1]=Uslow_S;
				Mdata[i][4]=pslow_S;
				// second pool (fast)
				Mdata[i][5]=Nfast_S;
				Mdata[i][6]=Ufast_S;
				Mdata[i][9]=pfast_S;
			}
			//driver
			else {
				// first pool (slow)
				Mdata[i][0]=Nslow_D;
				Mdata[i][1]=Uslow_D;
				Mdata[i][4]=pslow_D;
				// second pool (fast)
				Mdata[i][5]=Nfast_D;
				Mdata[i][6]=Ufast_D;
				Mdata[i][9]=pfast_D;
			}
			Mdata[i][2]=Tau_slow;
			// facilitation first pool (turned off)
			Mdata[i][3]=0;

			Mdata[i][7]=Tau_fast;
			// facilitation second pool (turned off)
			Mdata[i][8]=0;
		}
	}

	idxS=idxD=0;
	for (k=0; k<Ncl; k++) {
		CL[k].idxC=0;
	}
	for (i=0; i<N; i++) {
		for (j=0; j<W[i].num; j++) {
			if(FLAGS.groups==5){
				// set up 2 synapse groups according to MF identity
				if(MFgroup[W[i].idx[j]]==0)	{
					m=D_synID[idxD];
					idxD++;
				}
				else{
					m=S_synID[idxS];
					idxS++;
				}
				W[i].Gidx[j]=m;
				if(m==0)	Uprobidx=gsl_rng_uniform_int (r, NSYN)+NSYN;
				else 		Uprobidx=gsl_rng_uniform_int (r, NSYN);
			}
			else if(FLAGS.groups==55){
				// set up Ncl synapse groups according to MF identity
				for (k=0; k<Ncl; k++){
					if(MFgroup[W[i].idx[j]]==k)	{
						m=CL[k].synID[CL[k].idxC];
						CL[k].idxC++;
					}
				}
				W[i].Gidx[j]=m;
				inidx=round((Ncl-1-m)*NSYN_2/Ncl);
				Uprobidx=gsl_rng_uniform_int (r, NSYN_2/Ncl)+inidx;
			}
			else{
				std::cerr <<"Bad groups string, use different routine. \n";
				exit(EXIT_FAILURE);
			}

			W[i].strgth1[j]=Mdata[m][0];
			W[i].U1[j]=Mdata[m][1];
			W[i].tauD1[j]=Mdata[m][2];
			W[i].pf1[j]=Mdata[m][4];
			W[i].strgth2[j]=Mdata[m][5];
			W[i].U2[j]=Mdata[m][6];
			W[i].tauD2[j]=Mdata[m][7];
			W[i].pf2[j]=Mdata[m][9];

			// draw release prob from distribution
			if(FLAGS.Udistrb>1){
				W[i].U1[j]=Utoy_slow[Uprobidx];
				W[i].U2[j]=Utoy_slow[Uprobidx]*2.0/3.0;
				W[i].strgth2[j]=Ntoy_fast[Uprobidx];
			}

			if(FLAGS.MF_STP==0){
				W[i].tauD1[j]=0;	W[i].pf1[j]=0;
				W[i].tauD2[j]=0;	W[i].pf2[j]=0;
			}
			if(FLAGS.MF_STP==2){
				std::cerr <<"Mode not supported, use different routine. \n";
				exit(EXIT_FAILURE);
			}
			// auxillary quantities
			W[i].alpha1[j]=W[i].tauD1[j]*W[i].U1[j]*(1-W[i].pf1[j]);
			W[i].alpha2[j]=W[i].tauD2[j]*W[i].U2[j]*(1-W[i].pf2[j]);
			W[i].Upf1[j]=W[i].U1[j]*(1-W[i].pf1[j]);
			W[i].Upf2[j]=W[i].U2[j]*(1-W[i].pf2[j]);
			W[i].WU1[j]=W[i].strgth1[j]*W[i].U1[j];
			W[i].WU2[j]=W[i].strgth2[j]*W[i].U2[j];
		}
	}

//	**********************************
// 	calculation of GC patterns & adjust thresholds and gains
	calc_GCpattern_theta_gain_DS(W, MFpatterns_adjust, GCpatterns, theta, GCgain, Tavrg, Gavrg, hGCavrg, hGCvar, CLavrg, NETPARA, FLAGS, r);

//------------------------------------------------------------------------
//	if simulation is run with Golgi feedback reset thresholds and calculate weights appropriately
//------------------------------------------------------------------------
	if(FLAGS.GoC==1){
		for (i=0; i<N; i++) 	theta[i]=theta[i]*Tfac_gc;

		gsl_multiroot_function Fgolgi;
		Fgolgi.f = &find_golgi_params;

		struct golgi_para fix_golgipara;
		fix_golgipara.nuE_set=GCtarget;
		fix_golgipara.fc_set=nfrac;
		fix_golgipara.hext=hGCavrg;
		fix_golgipara.hext_var=hGCvar;
		fix_golgipara.theta=Tavrg*Tfac_gc;


		std::vector<double> para_out;
		para_out={0.0,0.0};
		n=idum;
		while (para_out[0]<=0.5 || para_out[0]>4 || para_out[1]<=0) {
			para_out=root_multiD_golgi(fix_golgipara, Fgolgi, 300, n);
			n++;
		}
		gain_global=para_out[0];
		avrgJEI=para_out[1];

		for (i=0; i<N; i++) 	GCgain[i]=gain_global;			// use global gain

		Gavrg=gain_global;
		Tavrg=Tavrg*Tfac_gc;

		// update GoC-to-GC synapses
		for(i= 0; i < N; i++) {
			JEI[i].strgth[0]=avrgJEI;
		}

		std::cout << "JEI = " << avrgJEI << "\t"<< "alpha = " << gain_global << "\n";
	}
//------------------------------------------------------------------------
// 	set up input signals & target patterns
//------------------------------------------------------------------------
	int ts[2];	//	one step
	ts[0]=-ttpre;
	ts[1]=0;
	for (k=0; k<Ndelays; k++) {

		t=0;
		// 	generate MF input signal
		ns=1;
		// pre CS
		for (j=0; j<M; j++) {
			MF_Tinput[j][0]=MFpatterns_preCS.at(j,pstart);
			if(MF_Tinput[j][0]<0) 	MF_Tinput[j][0]=0;
		}

		// during CS
		t++;
		while(t <= MAX_data) {
			for (j=0; j<M; j++) 	MF_Tinput[j][t]=MF_Tinput[j][t-1];
			if((t-ttpre)==ts[ns]){
				if((ns % 2)==0) np=pstart;
				else 			np=ptarget;

				for (j=0; j<M; j++) {
					MF_Tinput[j][t]=MFpatterns_CS.at(j,np);
					if(MF_Tinput[j][t]<0) 	MF_Tinput[j][t]=0;
				}
				ns++;
			}

			t++;
		}
	}

//	*************************
//  set initial conditions
	bool_tten=bool_tfive=0;
	tt=kk=tcut=0;
	bool_tcut=0;
	arma::mat GC_trial(N, BINS);
	arma::mat GCtrans(N, CUT_tsteps);

	for (j=0; j<M; j++) 	MF[j]=MF_Tinput[j][0];
	for (i=0; i<N; i++) 	GC[i]=GCm[i]=GCpatterns.at(i,pstart);

	a=0;
	for (i=0; i<N; i++) 	a+=GC[i];
	GCavrg= a/double(N);

	if(FLAGS.GoC==1){
		GoC=GoCm=GCavrg;
	}
	else{
		GoC=GoCm=0;
	}

	// set MF-GC synapses to steady state
	if(FLAGS.MF_STP>0){
		for (i=0; i<N; i++) {
			for (j=0; j<W[i].num; j++) {
				W[i].x1[j]=1.0/(1.0+W[i].alpha1[j]*MF[W[i].idx[j]]);
				W[i].x2[j]=1.0/(1.0+W[i].alpha2[j]*MF[W[i].idx[j]]);
			}
		}
	}
	FLAGS.trans_dum=FLAGS.trans;

//	*************************
//	compute basis functions
//	*************************
	t=-1;
	while(t < MAX_data) {
		t++;
		bool_tten=( (t % 10)==0);
		bool_tfive=( (t % 5)==0);
		bool_tcut=( (t % tcutfrac)==0);

		// update variables for next timestep
		for (j=0; j<M; j++) {
			MF[j]=MF_Tinput[j][t];
		}
		for (i=0; i<N; i++)	{
			GCm[i]=GC[i];
			for (j=0; j<W[i].num; j++){
				W[i].x1m[j]=W[i].x1[j];
				W[i].x2m[j]=W[i].x2[j];
			}
		}
		GoCm=GoC;

		//	calculate average GC activity
		GCavrg=0;
		for (i=0; i<N; i++) GCavrg+=GCm[i];
		GCavrg=GCavrg/double(N);

		GoC=GCavrg;

		// evolve GCs
		for (i=0; i<N; i++) {
			//	pre-synaptic inputs of GC i
			a=b=0;
			for (j=0; j<W[i].num; j++) {
				// phasic contribution
				a+=MF[W[i].idx[j]]*( W[i].WU1[j]*W[i].x1m[j] + W[i].WU2[j]*W[i].x2m[j] );
			}
			if(FLAGS.GoC==1){
				hGC[i]=a-GoCm*JEI[i].strgth[0];
			}
			else {
				hGC[i]=a;
			}
			GC[i]=GCm[i]+dtGC*(-GCm[i] + GCgain[i]*std::max(hGC[i]-theta[i],0.0));
			// uncomment for instantaneous GCs
			//~ GC[i]=GCgain[i]*std::max(hGC[i]-theta[i],0.0);
		}

		// instantaneous MLIs
		GCavrg=0;
		for (i=0; i<N; i++) GCavrg+=GC[i];
		GCavrg=GCavrg/double(N);
		MLI=GCavrg;

		// MF-to-GC STP
		if(FLAGS.MF_STP>0){
			if(FLAGS.trans_dum==1){
				for(i= 0; i < N; i++) {
					for (j=0; j<W[i].num; j++){
						W[i].x1[j]=W[i].x1m[j] +dt *((1.0-W[i].x1m[j])/W[i].tauD1[j] -W[i].x1m[j]*W[i].Upf1[j]*MF[W[i].idx[j]]);
						W[i].x2[j]=W[i].x2m[j] +dt *((1.0-W[i].x2m[j])/W[i].tauD2[j] -W[i].x2m[j]*W[i].Upf2[j]*MF[W[i].idx[j]]);
					}
				}
			}
			else {	// steady state
				for (i=0; i<N; i++) {
					for (j=0; j<W[i].num; j++) {
						W[i].x1[j]=1.0/(1.0+W[i].alpha1[j]*MF[W[i].idx[j]]);
						W[i].x2[j]=1.0/(1.0+W[i].alpha2[j]*MF[W[i].idx[j]]);
					}
				}
			}
		}

		// set different synapse types to their steady state values
		if(FLAGS.cuttsyn==2){
			// supporter (only works for 2D-2S configuration)
			for (i=0; i<N; i++) {
				for (j=2; j<W[i].num; j++) {
					W[i].x1[j]=1.0/(1.0+W[i].alpha1[j]*MF[W[i].idx[j]]);
					W[i].x2[j]=1.0/(1.0+W[i].alpha2[j]*MF[W[i].idx[j]]);
				}
			}
		}
		else if(FLAGS.cuttsyn==3){
			// driver (only works for 2D-2S configuration)
			for (i=0; i<N; i++) {
				for (j=0; j<2; j++) {
					W[i].x1[j]=1.0/(1.0+W[i].alpha1[j]*MF[W[i].idx[j]]);
					W[i].x2[j]=1.0/(1.0+W[i].alpha2[j]*MF[W[i].idx[j]]);
				}
			}
		}
		else if(FLAGS.cuttsyn==4){
			// slow pools
			for (i=0; i<N; i++) {
				for (j=0; j<W[i].num; j++) {
					W[i].x1[j]=1.0/(1.0+W[i].alpha1[j]*MF[W[i].idx[j]]);
				}
			}
		}
		else if(FLAGS.cuttsyn==5){
			// fast pools
			for (i=0; i<N; i++) {
				for (j=0; j<W[i].num; j++) {
					W[i].x2[j]=1.0/(1.0+W[i].alpha2[j]*MF[W[i].idx[j]]);
				}
			}
		}

		// write basis to file
		if((rr==0) & (FLAGS.mode==99 || FLAGS.mode==77) & ((t-ttpre)*dt>=-0.1) ){
			if (bool_tten){
				time_file << t*dt-Tpre << std::endl;

				for(i= 0; i < N; i++) GC_file << GC[i] << " ";
				GC_file << std::endl;

				MLI_file << MLI <<std::endl;
				if (FLAGS.GoC==1){
					GoC_file << GoC << std::endl;
				}
			}
			if (bool_tfive){

				for(i= 0; i < M; i++) MF_file << MF[i] << " ";
				MF_file << std::endl;

				STP_file <<t*dt-Tpre 	<<"\t ";
				for(i= 0; i < 50; i++) {
					for(j= 0; j < W[i].num; j++) STP_file <<  W[i].WU1[j]*W[i].x1m[j] + W[i].WU2[j]*W[i].x2m[j] <<" ";
				}
				STP_file << std::endl;
			}
		}

		if( (pp==0) & (pp2==0) & (rr==0) & (tt % int(dbint/dt)==0) & (kk<BINS) & ((t-ttpre)*dt>=-0.1) ) {
			time_learn_file << t*dt-Tpre << std::endl;
		}

		if( (tt % int(dbint/dt)==0) & (kk<BINS) & ((t-ttpre)*dt>=-0.1) ) {
			for(i= 0; i < N; i++)  	GC_trial.at(i,kk)=GC[i];
			kk++;
		}

		if (bool_tcut){
			if(t> ttcut_start && t<= ttcut_end){
				for(i= 0; i < N; i++) {
					GCtrans.at(i, tcut)= GC[i];
				}
				Time(tcut)=t*dt-Tpre;
				tcut++;
			}
		}

		tt++;
	}

	//	*************************
	//	remove GCs that are zero for all t
	//	*************************
	//	rows are GCs, cols are time-slices
	arma::uvec Jidx(N);
	for (i = 0; i < N; i++) Jidx(i)=i;
	arma::ucolvec status = arma::all(GCtrans==0,1);
	for(i= (N-1); i >=0; i--) {
		if(status(i)==1) {
			GCtrans.shed_row(i);
			GC_trial.shed_row(i);
			Jidx.shed_row(i);
		}
	}
	int Nred=GCtrans.n_rows;
	//	*************************
	//	cut GC basis (does not work with recurrent GoC)
	//	*************************

	// subtract GC steady state
	for (i = 0; i < Nred; i++) {
		a=GCtrans.at(i, CUT_tsteps-1);
		for (j = 0; j < CUT_tsteps; j++) {
			GCtrans.at(i, j)=GCtrans.at(i, j)-a;
		}
	}

	// find transient maxima based on zero crossings of derivative
	arma::rowvec GCdiff(CUT_tsteps);
	arma::uvec tcr;
	arma::vec npeaks(Nred);
	arma::uvec pkidx(Nred);
	arma::colvec GCmax(Nred);
	//~ int tmpidx;
	for (i = 0; i <Nred; i++) {
		// compute derivative
		GCdiff=arma::diff(GCtrans.row(i));
		// find zero crossings
		tcr=arma::find(GCdiff % arma::shift( GCdiff, -1) <= 0) +1;
		// exclude zero crossings that are too late
		for(j= (tcr.n_rows-1); j >=0; j--) {
			if(Time(tcr(j))>tlate) tcr.shed_row(j);
		}

		// for the DS model this is sufficient
		GCmax(i)=arma::max(arma::abs(GCtrans.row(i)));
		pkidx(i)=arma::index_max(arma::abs(GCtrans.row(i)));
		npeaks(i)=0;
	}

	// normalise GC transients
	for (i = 0; i < Nred; i++) {
		for (j = 0; j < CUT_tsteps; j++) {
			GCtrans.at(i, j)=GCtrans.at(i, j)/std::abs(GCmax(i));
		}
	}

	if(FLAGS.cutGC>0){
		// compute decay times
		arma::vec t_decay(Nred);
		arma::uvec didx;
		for (i = 0; i < Nred; i++) {
			didx=arma::find( (arma::abs(GCtrans(i,arma::span(pkidx(i), CUT_tsteps-1))) -0.1) > 0.0 , 1,"last")+pkidx(i);
			t_decay(i)=Time(didx(0));
		}

		// cut GC transients
		arma::ucolvec status2(Nred);
		if(FLAGS.cutGC==1){
			// cut above
			status2 = t_decay>=GCcut;
		}
		else if(FLAGS.cutGC==2){
			// cut below
			status2 = t_decay<GCcut;
			// also exclude also GC with multi-peaked transients
			didx=npeaks>1;
			status2=arma::sign(didx + status2); // implements logical OR
		}
		if(FLAGS.cutGC>0){
			for(i= (Nred-1); i >=0; i--) {
				if(status2(i)==1) {
					GCtrans.shed_row(i);
					GC_trial.shed_row(i);
					Jidx.shed_row(i);
				}
			}
		}
		Nred=GCtrans.n_rows;
		std::cout << std::endl;
		if(FLAGS.cutGC==1)		std::cout <<"Cutting GC transients longer than "<< GCcut << " s." << std::endl;
		else if(FLAGS.cutGC==2) 	std::cout <<"Cutting GC transients shorter than "<< GCcut << " s." << std::endl;
		std::cout <<"Number of GC transients after cut: "<< GCtrans.n_rows << std::endl;
		std::cout <<"Number of transient time bins: "<< GCtrans.n_cols << std::endl;
	}

	if( (rr==0) & (FLAGS.mode==99 || FLAGS.mode==77) ){
		// write transients to file
		GCtrans.save(GCtrans_file,arma::raw_ascii);
		Time.save(timeGCtrans_file,arma::raw_ascii);
	}

	//	*************************
	//	do learning
	//	*************************
	err_final=0;
	for (k=0; k<Ndelays; k++) {

		err_final+=learn_weights_DS(GC_trial, J, Jidx, PC, delays[k], interv_min[k], interv_max[k], dt, Tpre, PCsp, CFsp, CFscale, errWeight, J2weight, JI, alpha, dbint, N, Nred, Ntrials,
									BINS, MAX_data, ttpre, rr, pp, pp2, FLAGS, CF_file, htarget_BINS_file, r);

		// write final PC firing to files
		if(FLAGS.mode==77){
			if(k==0){
			for (t=0; t<BINS; t++) PC_bayes_1_file << PC[t] << "\t";
				PC_bayes_1_file << endl;
			}
			else if(k==1){
				for (t=0; t<BINS; t++) PC_bayes_2_file << PC[t] << "\t";
				PC_bayes_2_file << endl;
			}
			else if(k==2){
				for (t=0; t<BINS; t++) PC_bayes_3_file << PC[t] << "\t";
				PC_bayes_3_file << endl;
			}
			else if(k==3){
				for (t=0; t<BINS; t++) PC_bayes_4_file << PC[t] << "\t";
				PC_bayes_4_file << endl;
			}
			else if(k==4){
				for (t=0; t<BINS; t++) PC_bayes_5_file << PC[t] << "\t";
				PC_bayes_5_file << endl;
			}
		}
		else{
			if(k==0){
				for (t=0; t<BINS; t++) PC_025_file << PC[t] << "\t";
				PC_025_file << endl;
			}
			else if(k==1){
				for (t=0; t<BINS; t++) PC_050_file << PC[t] << "\t";
				PC_050_file << endl;
			}
			else if(k==2){
				for (t=0; t<BINS; t++) PC_100_file << PC[t] << "\t";
				PC_100_file << endl;
			}
			else if(k==3){
				for (t=0; t<BINS; t++) PC_200_file << PC[t] << "\t";
				PC_200_file << endl;
			}
			else if(k==4){
				for (t=0; t<BINS; t++) PC_300_file << PC[t] << "\t";
				PC_300_file << endl;
			}
			else if(k==5){
				for (t=0; t<BINS; t++) PC_500_file << PC[t] << "\t";
				PC_500_file << endl;
			}
			else if(k==6){
				for (t=0; t<BINS; t++) PC_700_file << PC[t] << "\t";
				PC_700_file << endl;
			}
		}
		// write final PC synaptic weights to file
		if( (FLAGS.mode==77) | (FLAGS.mode==99) ) {
			for (i=0; i<N; i++) J_file << J.strgth[i] << "\t";
			J_file << endl;
			for(i= 0; i < Nred; i++) 	Jidx_file << Jidx(i)+1 << "\t";
			Jidx_file << endl;
		}

	} // end of delays-loop

//	*************************
// 	write stuff to files

	if(( (FLAGS.mode==77) | (FLAGS.mode==99) ) & (rr==0)){
		for(i= 0; i < N; i++) {
			for(j= 0; j < W[i].num; j++) {
				STPpara_file << i+1 << "\t" << W[i].idx[j]+1 <<"\t" << W[i].Gidx[j] <<"\t" << MFgroup[W[i].idx[j]]
								<< "\t" << W[i].strgth1[j]<< "\t" << W[i].U1[j] << "\t" << W[i].tauD1[j]*1000 << "\t" << W[i].pf1[j]
								<< "\t"	<< W[i].strgth2[j]<< "\t" << W[i].U2[j] << "\t" << W[i].tauD2[j]*1000 << "\t" << W[i].pf2[j] << std::endl;
			}
		}
	}

	if(FLAGS.groups==5){
		err_final_file << MF_avrg_D << "\t" << MF_avrg_S << "\t" << MF_std_D << "\t" << MF_std_S << "\t" << err_final/double(Ndelays) << std::endl;
	}
	else{
		for (i = 0; i < Ncl; i++) {
			err_final_file << MF_avrg_CL[i] << "\t" << MF_std_CL[i] << "\t";
		}
		err_final_file  << err_final/double(Ndelays) << std::endl;
	}

	learnparam_file << Ntrials-1<< "\t" << CFsp<< "\t"<< CFscale << "\t"<< NAN << "\t"<< errWeight << "\t"<< alpha_learn << "\t"<< cI
					<< "\t"<< Nreal << "\t"<< Npara << "\t"<< Npara2 << "\t"<< idum << "\t" << avrgJIE<< "\t" << avrgJEI
					<< "\t" << gain_global << "\t" << nfrac << "\t" << GCtarget << "\t" << NAN << "\t" << Tfac_gc
					<< "\t" << MFUcorr << "\t" << Ncl
					<< "\t" << pact[0] << "\t" << pact[1] << "\t" << pact[2] << "\t" << pact[3] << "\t" << pact[4]
					<< "\t" << grouplims[0] << "\t" << grouplims[1] << "\t" << grouplims[2] << "\t" << grouplims[3] << "\t" << grouplims[4] << "\t" << grouplims[5]
					<< "\t" << J2weight << std::endl;

	stats_file << hGCavrg << "\t" << hGCvar << "\t" << Tavrg << "\t" << Gavrg << "\t" << CLavrg <<std::endl;

	for(i= 0; i < N; i++){
		for(j= 0; j < W[i].num; j++) {
			delete [] W[i].x1p[j];		delete [] W[i].x2p[j];
		}
		delete [] W[i].idx; 	delete [] W[i].Gidx;
		delete [] W[i].strgth1;	delete [] W[i].strgth2;
		delete [] W[i].x1p;		delete [] W[i].x2p;
		delete [] W[i].x1;		delete [] W[i].x2;
		delete [] W[i].x1m;		delete [] W[i].x2m;
		delete [] W[i].U1;		delete [] W[i].tauD1;	delete [] W[i].pf1;
		delete [] W[i].U2;		delete [] W[i].tauD2;	delete [] W[i].pf2;

		delete [] W[i].alpha1;	delete [] W[i].WU1; 	delete [] W[i].Upf1;
		delete [] W[i].alpha2;	delete [] W[i].WU2; 	delete [] W[i].Upf2;
	}

	drivers.clear();
	supporters.clear();

	} // end of Nreal-loop
	} // end of Npara-loop (1)
	} // end of Npara-loop (2)

//	*************************
// 	terminate program
	GC_file.close();
	MF_file.close();
	STP_file.close();
	MLI_file.close();
	GoC_file.close();
	CF_file.close();
	time_file.close();
	time_learn_file.close();
	STPpara_file.close();
	GCtrans_file.close();
	timeGCtrans_file.close();

	err_final_file.close();
	PC_025_file.close();
	PC_050_file.close();
	PC_100_file.close();
	PC_200_file.close();
	PC_300_file.close();
	PC_500_file.close();
	PC_700_file.close();
	PC_bayes_1_file.close();
	PC_bayes_2_file.close();
	PC_bayes_3_file.close();
	PC_bayes_4_file.close();
	PC_bayes_5_file.close();
	htarget_BINS_file.close();
	J_file.close();
	Jidx_file.close();
	learnparam_file.close();
	stats_file.close();
	flags_file.close();
	gsl_rng_free (r);

	for(i= 0; i < M; i++){
		delete [] MF_Tinput[i];
	}

	for(j= 0; j < N; j++){
		delete [] JEI[j].strgth;
	}
	delete [] J.strgth;
	delete [] MF_Tinput;

  return(0);

}
