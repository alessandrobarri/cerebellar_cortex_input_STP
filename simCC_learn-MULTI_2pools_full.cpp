#include <sstream>
#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <limits>
#include <float.h>
#include <vector>
#include <algorithm>
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
	double ** x1p; double ** x2p; double ** u1p; double ** u2p; double ** qp;
	double * x1; double * u1; double * x2; double * u2; double * q;
	double * x1m; double * u1m; double * x2m; double * u2m; double * qm;
	double * U1; double * tauD1; double * tauF1; double * pf1;
	double * U2; double * tauD2; double * tauF2; double * pf2;
	double * tauG; double * deltaG;
	double * kp; double * kpp;
	double * spill;
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
	int GoC,driver,Udistrb,mode,cuttsyn,tsignal,cutGC,normGC;
};
struct netpara {
	int N,M,MFSYN,NSYN,Ncl,PP;
	double pdriver,MFUcorr,nfrac,GCtarget,theta_global,taumin;
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
#include "headers/MFGC_STP_2pools_evolve_func.h"
#include "headers/generate_ratePv_correlations.h"
#include "headers/calc_GCpattern_theta_gain_full.h"
#include "headers/learn_weights_DS.h"
#include "headers/calc_Tdim.h"
#include "headers/calc_Ddim.h"
#include "headers/rtnorm.cpp"

int main (int argc, char* argv[]) {

	#include "common_paras.h"
	#include "syn_paras.h"

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

	if (pattern.compare("Gauss") == 0) 				FLAGS.pattern=0;
    else if (pattern.compare("lognorm") == 0) 		FLAGS.pattern=1;
    else if (pattern.compare("exponential") == 0) 	FLAGS.pattern=2;
    else if (pattern.compare("gamma") == 0) 		FLAGS.pattern=3;
    else if (pattern.compare("uniform_multi") == 0) FLAGS.pattern=4;
    else if (pattern.compare("trunc_Gauss") == 0)	FLAGS.pattern=6;
    else {std::cerr <<"Bad pattern string. \n";exit(EXIT_FAILURE);}

    if (gaincontrol.compare("global") == 0) 				FLAGS.gain=0;
    else if (gaincontrol.compare("individual") == 0) 		FLAGS.gain=1;
    else if (gaincontrol.compare("individual_active") == 0) FLAGS.gain=2;
    else if (gaincontrol.compare("no") == 0) 				FLAGS.gain=3;
    else {std::cerr <<"Bad gaincontrol string. \n";exit(EXIT_FAILURE);}

    if 	(Groups.compare("5SYN_5MF") == 0)			FLAGS.groups=2;
    else {std::cerr <<"Bad Groups string: "<< Groups <<" \n";exit(EXIT_FAILURE);}

    if 		(Driver.compare("always") == 0) 	FLAGS.driver=1;
    else if (Driver.compare("random") == 0) 	FLAGS.driver=0;
    else {std::cerr <<"Bad Driver string. \n";exit(EXIT_FAILURE);}

    if 		(transients.compare("on") == 0) 	FLAGS.trans=1;
    else if (transients.compare("off") == 0) 	FLAGS.trans=0;
    else {std::cerr <<"Bad transients string. \n";exit(EXIT_FAILURE);}

	if (Udistrb.compare("no") == 0)			FLAGS.Udistrb=0;
    else {std::cerr <<"Bad Udistrb string. \n";exit(EXIT_FAILURE);}

    if (mode.compare("cI") == 0)			FLAGS.mode=7;
    else if (mode.compare("bayesian") == 0)	FLAGS.mode=77;
    else if (mode.compare("sample") == 0)	FLAGS.mode=99;
    else {std::cerr <<"Bad mode string. \n";exit(EXIT_FAILURE);}

    if 		(Golgi.compare("on") == 0) 		FLAGS.GoC=1;
    else if (Golgi.compare("off") == 0) 	FLAGS.GoC=0;
    else {std::cerr <<"Bad Golgi string. \n";exit(EXIT_FAILURE);}

    if 		(Tsignal.compare("Gauss") == 0) 		FLAGS.tsignal=0;
    else if (Tsignal.compare("delta") == 0) 		FLAGS.tsignal=1;
    else {std::cerr <<"Bad Tsignal string. \n";exit(EXIT_FAILURE);}

    if (cuttsyn.compare("no") == 0)				FLAGS.cuttsyn=0;
    else if (cuttsyn.compare("bounds") == 0) 	FLAGS.cuttsyn=1;
    else if (cuttsyn.compare("cutS") == 0)		FLAGS.cuttsyn=2;
    else if (cuttsyn.compare("cutD") == 0)		FLAGS.cuttsyn=3;
    else if (cuttsyn.compare("cutSLOW") == 0)	FLAGS.cuttsyn=4;
    else if (cuttsyn.compare("cutFAST") == 0)	FLAGS.cuttsyn=5;
    else if (cuttsyn.compare("cutD_AND_cutFAST") == 0)	FLAGS.cuttsyn=35;
    else {std::cerr <<"Bad Cuttsyn string. \n";exit(EXIT_FAILURE);}

    if (cutGC.compare("no") == 0)				FLAGS.cutGC=0;
    else if (cutGC.compare("above") == 0) 		FLAGS.cutGC=1;
    else if (cutGC.compare("below") == 0)		FLAGS.cutGC=2;
    else {std::cerr <<"Bad cutGC string. \n";exit(EXIT_FAILURE);}

    //~ if(FLAGS.mode==6 && FLAGS.Udistrb==0) {std::cerr <<"Set Udistrb= 'yes' in order to run 'MF_U_corr' option. \n";exit(EXIT_FAILURE);}

	if(FLAGS.cuttsyn==1){
		std::cout << std::endl;
		std::cout << "cutting tsyn from " <<  tsyncut_lo << "s to " << tsyncut_hi << "s" << std::endl;
	}
	else if(FLAGS.cuttsyn==2){
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
	else if(FLAGS.cuttsyn==35){
		std::cout << std::endl;
		std::cout << "setting driver synapses AND all fast pools to steady state" << std::endl;
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
	// if GoC-GC synapses exist:
	if(FLAGS.GoC==1)	{
		Tpre=Tpre+0.2;		// use longer baseline time-interval
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
	double kp,kpp,err_final,JI;
	//~ double xtrans1,xtrans2,expfac1,expfac2;
	double Tavrg,Gavrg,CLavrg,hGCpeak,GCpeak,Dim_GCt,Dim_MF,Dim_GCss,Dim_GConset,Dim_GCtrans;
	double MF_avrg_CL[5],MF_std_CL[5];
	double a,b,c,d;
	int i,j,k,m,n,t,idxS,idxD,idxC0,idxC1,idxC2;
	int tt,kk,pp,pp2,ns,np,rr,tsvd,tcut,tchg1,tchg2;
	bool bool_time,bool_tfive,bool_tten,bool_tcut;

	const double ttshift=0.1;					// time-point of transient evaluation

	vector<double> delays = {0.025, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7};
	int Ndelays=7;
	//~ vector<double> delays = {0.025, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.01};
	//~ const int Ndelays=8;

	// for bayesian learning
	vector<double> interv_min = {0.025, 0.05, 0.1, 0.2, 0.3};
	vector<double> interv_max = {0.15, 0.2, 0.3, 0.4, 0.5};
	if(FLAGS.mode==77)	Ndelays=5;

	const int NSYN=N*MFSYN;							// number of MFGC synapses
	const int NSYN_2=2*NSYN;						// double NSYN for identity generation

	double Mdata[5][Nsynpara_2pools];
	double beta_ab[5][4];

	// MF-GC synapse statistics
	const double pG1=0.06;						// MF to GC synapse group probability
	const double pG2=0.16;
	const double pG3=0.38;
	const double pG4=0.24;
	const double pG5=0.16;

	// time-window for SVD
	const int tfrac=5;
	double tsvd_start=0;								// 100ms before onset
	const double tsvd_end=duration;						// 1s after onset
	// use longer baseline time-interval whith recurrent golgi cells
	if(FLAGS.GoC==1){
		tsvd_start=0.2;
	}
	const int SVD_tsteps=(tsvd_end-tsvd_start)/(dt*tfrac);	// *tfrac : downsampling by factor tfrac

	const int ttsvd_start=int(tsvd_start/dt);
	const int ttsvd_end=int(tsvd_end/dt);

	// time-window for GC transient cutting
	const int tcutfrac=10;
	const double tcut_start=Tpre;								// at onset
	const double tcut_end=1.2+Tpre;								// 1.2s after onset
	const int CUT_tsteps=(tcut_end-tcut_start)/(dt*tcutfrac);	// *tcutfrac : downsampling by factor tcutfrac

	const int ttcut_start=int(tcut_start/dt);
	const int ttcut_end=int(tcut_end/dt);

	const double taumin=1.0/1000.0;

	//	usefull quantities
	const double sqN=sqrt(N);
	const double dtGC=dt/tauGC;

	// learning parameters
	const double alpha=alpha_learn/sqN;

	// generate switching times
	const int Nt=2;						//	one step
	int ts[Nt];
	ts[0]=-ttpre;
	ts[1]=0;
	// **********************************
	// arrays
	struct presyn_2pools W[N];
	struct presyn J;
	struct presyn JEI[N];			// GoC-to-GC
	struct cluster_struct CL[5];
	struct netpara NETPARA;
	J.strgth=new double [N];
	double theta[N];
	std::vector <double> dummyD;
	std::vector <int> dummy;
	std::vector <size_t> idx_sorted;
	std::vector <int> drivers,supporters;
	std::vector <int> D_synID,S_synID,C0_synID,C1_synID,C2_synID;
	std::vector<double> para_out;
	int MFgroup[M];

	arma::mat MFpatterns(M,PP);
	arma::mat GCinT(N, SVD_tsteps);
	arma::mat GCpatterns(N, PP);
	arma::mat GCpatterns_onset(N, PP);		arma::mat hGC_onset(N, PP);
	arma::mat GCpatterns_trans(N, PP);		arma::mat hGC_trans(N, PP);

	arma::vec Time(CUT_tsteps);
	arma::colvec GCmax(N);

//	double para1[Npara],para2[Npara];

	double PC[BINS];
	double GC[N],GCm[N],GCgain[N],hGC[N];
	double MF[M],MFm[M];
	double GoC,GoCm;
	double GCavrg,hGCavrg,hGCvar,MLI;
	double tsyn_avrg[Ncl][2],tsyn_std[Ncl][2],foff[Ncl][2],tempDOUBLE[Ncl][2];
	int tempINT[Ncl][2];

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
	//~ std::ofstream STP_NU_file(path+string("STP_NU.dat"));				if (!STP_NU_file) exit(EXIT_FAILURE);
	std::ofstream MLI_file(path+string("MLI.dat"));						if (!MLI_file) exit(EXIT_FAILURE);
	std::ofstream GoC_file(path+string("GoC.dat"));						if (!GoC_file) exit(EXIT_FAILURE);
	std::ofstream CF_file(path+string("CF.dat"));						if (!CF_file) exit(EXIT_FAILURE);
	std::ofstream time_file(path+string("time.dat"));					if (!time_file) exit(EXIT_FAILURE);
	std::ofstream time_learn_file(path+string("time_learn.dat"));		if (!time_learn_file) exit(EXIT_FAILURE);
	std::ofstream STPpara_file(path+string("STPpara.dat"));				if (!STPpara_file) exit(EXIT_FAILURE);
	std::ofstream GCtrans_file(path+string("GCtrans.dat"));				if (!GCtrans_file) exit(EXIT_FAILURE);
	std::ofstream timeGCtrans_file(path+string("timeGCtrans.dat"));		if (!timeGCtrans_file) exit(EXIT_FAILURE);

	std::ofstream err_final_file(path+string("err_final.dat"));			if (!err_final_file) exit(EXIT_FAILURE);
	std::ofstream PC_010_file(path+string("PC_final_010.dat"));			if (!PC_010_file) exit(EXIT_FAILURE);
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

	//~ std::ofstream htarget_file(path+string("htarget.dat"));				if (!htarget_file) exit(EXIT_FAILURE);
	std::ofstream htarget_BINS_file(path+string("htarget_BINS.dat"));	if (!htarget_BINS_file) exit(EXIT_FAILURE);
	std::ofstream theta_file(path+string("theta.dat"));						if (!theta_file) exit(EXIT_FAILURE);
	std::ofstream J_file(path+string("J.dat"));							if (!J_file) exit(EXIT_FAILURE);
	std::ofstream Jidx_file(path+string("Jidx.dat"));					if (!Jidx_file) exit(EXIT_FAILURE);
	std::ofstream trans_file(path+string("trans.dat"));					if (!trans_file) exit(EXIT_FAILURE);
	std::ofstream tsyn_file(path+string("teff.dat"));					if (!tsyn_file) exit(EXIT_FAILURE);
	std::ofstream tsyndistrb_file(path+string("teff_distrb.dat"));		if (!tsyndistrb_file) exit(EXIT_FAILURE);
	std::ofstream Ampdistrb_file(path+string("Amp_distrb.dat"));		if (!Ampdistrb_file) exit(EXIT_FAILURE);
	std::ofstream tsyn_MFsort_file(path+string("teff_MFsort.dat"));		if (!tsyn_MFsort_file) exit(EXIT_FAILURE);
	//~ std::ofstream tsyn_file_D(path+string("teff_D.dat"));				if (!tsyn_file_D) exit(EXIT_FAILURE);
	std::ofstream dim_file(path+string("dimension.dat"));				if (!dim_file) exit(EXIT_FAILURE);
	std::ofstream learnparam_file(path+string("learnparam.dat"));		if (!learnparam_file) exit(EXIT_FAILURE);
	std::ofstream stats_file(path+string("stats.dat"));					if (!stats_file) exit(EXIT_FAILURE);
	std::ofstream flags_file(path+string("flags.dat"));					if (!flags_file) exit(EXIT_FAILURE);


	flags_file << 	FLAGS.MF_STP << "\t"<< FLAGS.threshold << "\t"<< FLAGS.pattern << "\t"<< FLAGS.gain << "\t"<< FLAGS.normGC << "\t"<< FLAGS.groups << "\t"<< FLAGS.driver<< "\t"<<
					FLAGS.trans<< "\t"<< FLAGS.Udistrb<< "\t"<< FLAGS.mode << "\t"<< FLAGS.GoC << "\t"<<  FLAGS.tsignal << "\t"<< FLAGS.cuttsyn << "\t"<<  FLAGS.cutGC <<std::endl;

	arma::vec MF_U_CORR = arma::linspace<arma::vec>(-0.99, 0.99, Npara);
	//~ arma::vec CI_ARR = arma::logspace<arma::vec>(-1, 3, Npara);
	arma::vec CI_ARR = { 0.1, 0.5, 1, 2};
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
	NETPARA.taumin=taumin;
	// **********************************
	// calculate MF parameters that generate desired tsyn distributions

	if(FLAGS.mode==7){
		cout << "Varying JI" << "\n";
		Npara=CI_ARR.n_elem;
	}
	else if(FLAGS.mode==77){
		cout << "Bayesian learning" << "\n";
		//~ Npara=1;	// only one parameter set when in sample mode
		CFsp=CFsp_bayes;
	}
	else if(FLAGS.mode==99){
		cout << "Sample eye-blink learning" << "\n";
		Npara=1;	// only one parameter set when in sample mode
	}
	// **********************************
	// start Npara loop
	for (pp = 0; pp < Npara; pp++){

		// loop over different supporter parameters
		//~ nuS=para1[pp];
		//~ sigS=para2[pp];

		// loop over different rate-pr identity correlations
		//~ MFUcorr=MF_U_CORR(pp);
		MFUcorr=0.999999;
		// loop over different inhibition strengths
		if(FLAGS.mode==7) {
			cI=CI_ARR(pp);
		}
		//~ cI=464.16;
		JI=cI;

	for (pp2 = 0; pp2 < Npara2; pp2++){

		//~ nuD= (nuDmax-nuDmin)/double(Npara) *pp+nuDmin;			// check pp, pp2 !!
		//~ sigD= (sigDmax-sigDmin)/double(Npara) *pp2+sigDmin;
		//~ nuS= (nuSmax-nuSmin)/double(Npara) *pp+nuSmin;			// check pp, pp2 !!
		//~ sigS= (sigSmax-sigSmin)/double(Npara) *pp2+sigSmin;

		cout << "********************* \n";
		cout << "Parameter set number " << pp*Npara2+pp2+1<< "\n";
		cout 	<< "nu0 = " << MF_avrg_G1  << "\t"<< "sig0 = " << MF_std_G1 << "\t"
						<< "nu1 = " << MF_avrg_G2	<< "\t"<< "sig1 = " << MF_std_G2  << "\t"
						<< "nu2 = " << MF_avrg_G3	<< "\t"<< "sig2 = " << MF_std_G3 << "\t"
						<< "nu3 = " << MF_avrg_G4	<< "\t"<< "sig3 = " << MF_std_G4 << "\t"
						<< "nu4 = " << MF_avrg_G5	<< "\t"<< "sig4 = " << MF_std_G5 << "\t"
						<< "MFUcorr = " << MFUcorr << "\t"<< "JI = " << JI << "\t"<< "trials = " << Ntrials << "\n";

	// start Nreal loop
	for (rr = 0; rr < Nreal; rr++){
		if( ( (rr+1) % 5) ==0) cout << "Realisation number " << rr+1 << "\n";

		MF_avrg_CL[0]=MF_avrg_G1;		MF_std_CL[0]=MF_std_G1;
		MF_avrg_CL[1]=MF_avrg_G2;		MF_std_CL[1]=MF_std_G2;
		MF_avrg_CL[2]=MF_avrg_G3;		MF_std_CL[2]=MF_std_G3;
		MF_avrg_CL[3]=MF_avrg_G4;		MF_std_CL[3]=MF_std_G4;
		MF_avrg_CL[4]=MF_avrg_G5;		MF_std_CL[4]=MF_std_G5;

		drivers.clear();
		supporters.clear();

	// **********************************
	//	generate MF patterns

	// calculate mu and sig parameters for thresholded Gaussian
		// calculate mu and sig parameters for thresholded Gaussian
	if(FLAGS.pattern==0){
		struct threshNpara thpara;
		gsl_multiroot_function FthresN;
		FthresN.f = &thrN;

		for (i = 0; i < 5; i++) {
			thpara.mean=MF_avrg_CL[i];
			thpara.std=MF_std_CL[i];
			para_out=find_multiD (thpara, FthresN, 200);
			MF_avrg_CL[i]=para_out[0];
			MF_std_CL[i]=para_out[1];
		}

//		std::cout <<para_out[0] << std::endl;
//		std::cout <<para_out[1] << std::endl;
	}
	// calculate mu and sig parameters for truncated Gaussian
	else if(FLAGS.pattern==6){
		struct threshNpara thpara;
		gsl_multiroot_function FtruncN;
		FtruncN.f = &truncN;

		for (i = 0; i < 5; i++) {
			thpara.mean=MF_avrg_CL[i];
			thpara.std=MF_std_CL[i];
			para_out=find_multiD (thpara, FtruncN, 200);
			MF_avrg_CL[i]=para_out[0];
			MF_std_CL[i]=para_out[1];
		}
	}

	// generate MF patterns according to MF identities
	if(FLAGS.groups==2){
		// 5 groups: different MF types
		for (i = 0; i < M; i++) {
			a=gsl_rng_uniform (r);
			if(a<pG1) {
				MFgroup[i]=0;	drivers.push_back(i);
				for (k = 0; k < PP; k++) MFpatterns.at(i,k)=draw_patterns(MF_avrg_CL[0], MF_std_CL[0], pact[0], FLAGS.pattern,r);
			}
			else if ( (a>=pG1) & (a<pG1+pG2) ) {
				MFgroup[i]=1;	drivers.push_back(i);
				for (k = 0; k < PP; k++) MFpatterns.at(i,k)=draw_patterns(MF_avrg_CL[1], MF_std_CL[1], pact[1], FLAGS.pattern,r);
			}
			else if ( (a>=pG1+pG2) & (a<pG1+pG2+pG3) ){
				MFgroup[i]=2;	supporters.push_back(i);
				for (k = 0; k < PP; k++) MFpatterns.at(i,k)=draw_patterns(MF_avrg_CL[2], MF_std_CL[2], pact[2], FLAGS.pattern,r);
			}
			else if ( (a>=pG1+pG2+pG3) & (a<pG1+pG2+pG3+pG4) ){
				MFgroup[i]=3;	supporters.push_back(i);
				for (k = 0; k < PP; k++) MFpatterns.at(i,k)=draw_patterns(MF_avrg_CL[3], MF_std_CL[3], pact[3], FLAGS.pattern,r);
			}
			else {
				MFgroup[i]=4;	drivers.push_back(i);
				for (k = 0; k < PP; k++) MFpatterns.at(i,k)=draw_patterns(MF_avrg_CL[4], MF_std_CL[4], pact[4], FLAGS.pattern,r);
			}
		}
	}
	else{
		std::cerr <<"Bad groups string, use different routine. \n";
		exit(EXIT_FAILURE);
	}

	// **********************************
	//	generate golgi connectivity

	// set up GoC-to-GC synapses
	for(i= 0; i < N; i++) {
		// log-normal
		//~ JEI[i].strgth[j]=gsl_ran_lognormal (r, zetaJ, sigmaJ);
		JEI[i].strgth[0]=avrgJEI;
	}

	// **********************************
	//	generate MF-to-GC connectivity
	int *shffl_arr=new int[M];
	// set up connectivity
	for(i= 0; i < N; i++) {

		// draw exactly MFSYN distinct random connections
		for (j=0;j<M;j++) shffl_arr[j]=j;
		gsl_ran_shuffle (r, shffl_arr, M, sizeof (int));

		// ---------------------------------------
		// check that GCs receive the correct MF-types
		// ---------------------------------------
		if( FLAGS.groups>0 && FLAGS.driver==1){
			a=0;
			// check if driver is missing; if so, add
			for (j=0;j<MFSYN;j++) {
				if(MFgroup[shffl_arr[j]]==0 || MFgroup[shffl_arr[j]]==1 || MFgroup[shffl_arr[j]]==4) a+=1;
			}
			if(a==0){
				// assign driver MF to first synapse on GC
				m=gsl_rng_uniform_int (r, drivers.size());
				shffl_arr[0]=drivers[m];
			}

		}
		else if( FLAGS.groups>0 && FLAGS.driver==2){
			// draw random distinc drivers & supporters
			int *shffl_drivers=new int[drivers.size()];
			int *shffl_supporters=new int[supporters.size()];

			for (k=0;k<drivers.size();k++) 	shffl_drivers[k]=drivers[k];
			gsl_ran_shuffle (r, shffl_drivers, drivers.size(), sizeof (int));
			for (j=0;j<2;j++)				shffl_arr[j]=shffl_drivers[j];

			for (k=0;k<supporters.size();k++) 	shffl_supporters[k]=supporters[k];
			gsl_ran_shuffle (r, shffl_supporters, supporters.size(), sizeof (int));
			for (j=2;j<MFSYN;j++)				shffl_arr[j]=shffl_supporters[j];

			delete [] shffl_drivers;
			delete [] shffl_supporters;
		}

		for (j=0;j<MFSYN;j++) dummy.push_back (shffl_arr[j]);
		k=MFSYN;

		W[i].num=k;
		W[i].idx=new int [k];			W[i].Gidx=new int [k];
		W[i].strgth1=new double [k];	W[i].strgth2=new double [k];
		W[i].x1p=new double* [k];		W[i].x2p=new double* [k];
		W[i].u1p=new double* [k]; 		W[i].u2p=new double* [k];
		W[i].qp=new double* [k];
		W[i].x1=new double [k];			W[i].x2=new double [k];			W[i].x1m=new double [k];		W[i].x2m=new double [k];
		W[i].u1=new double [k];			W[i].u2=new double [k];			W[i].u1m=new double [k];		W[i].u2m=new double [k];
		W[i].q=new double [k];			W[i].qm=new double [k];
		W[i].U1=new double [k];			W[i].tauD1=new double [k];		W[i].tauF1=new double [k]; 	W[i].pf1=new double [k];
		W[i].U2=new double [k];			W[i].tauD2=new double [k];		W[i].tauF2=new double [k]; 	W[i].pf2=new double [k];

		W[i].tauG=new double [k];		W[i].deltaG=new double [k];
		W[i].kp=new double [k];			W[i].kpp=new double [k];

		//~ W[i].alpha1=new double [k];		W[i].WU1=new double [k];	 	W[i].Upf1=new double [k];
		//~ W[i].alpha2=new double [k];		W[i].WU2=new double [k];	 	W[i].Upf2=new double [k];
		W[i].tsyn1=new double* [k];		W[i].tsyn2=new double* [k];

		for(j= 0; j < k; j++) {
			W[i].idx[j]=dummy[j];
			W[i].x1p[j]=new double [PP];	W[i].x2p[j]=new double [PP];
			W[i].u1p[j]=new double [PP];	W[i].u2p[j]=new double [PP];
			W[i].qp[j]=new double [PP];
			W[i].tsyn1[j]=new double [PP];	W[i].tsyn2[j]=new double [PP];
		}
		dummy.clear();
	}
	delete [] shffl_arr;

	// **********************************
	// 	generate correlation between MF- and synapse-types
	generate_ratePv_correlations(CL, D_synID, S_synID, grouplims, NETPARA, FLAGS, r);

	// -----------------------------------
	// prepare hand-made synapses
	// -----------------------------------
	// 5 groups
	for (i=0; i<5; i++) {

		Mdata[i][2]=Tau_slow;
		// facilitation first pool (turned off)
//			Mdata[i][3]=0;

		Mdata[i][7]=Tau_fast;
		// facilitation second pool (turned off)
//			Mdata[i][8]=0;

		// desensitisation
		//~ Mdata[i][10]=100.0/1000.0;
		//~ Mdata[i][11]=0.1;
		Mdata[i][10]=tauG;
		Mdata[i][11]=deltaG;

		// G1
		if(i==0){
			// first pool
			Mdata[i][0]=Nslow_G1;
			Mdata[i][1]=Uslow_G1;
			Mdata[i][3]=tauF_G1;
			Mdata[i][4]=pslow_G1;
			// second pool
			Mdata[i][5]=Nfast_G1;
			Mdata[i][6]=Ufast_G1;
			Mdata[i][8]=tauF_G1;
			Mdata[i][9]=pfast_G1;
		}
		// G2
		else if(i==1){
			// first pool
			Mdata[i][0]=Nslow_G2;
			Mdata[i][1]=Uslow_G2;
			Mdata[i][3]=tauF_G2;
			Mdata[i][4]=pslow_G2;
			// second pool
			Mdata[i][5]=Nfast_G2;
			Mdata[i][6]=Ufast_G2;
			Mdata[i][8]=tauF_G2;
			Mdata[i][9]=pfast_G2;
		}
		// G3
		else if(i==2){
			// first pool
			Mdata[i][0]=Nslow_G3;
			Mdata[i][1]=Uslow_G3;
			Mdata[i][3]=tauF_G3;
			Mdata[i][4]=pslow_G3;
			// second pool
			Mdata[i][5]=Nfast_G3;
			Mdata[i][6]=Ufast_G3;
			Mdata[i][8]=tauF_G3;
			Mdata[i][9]=pfast_G3;
		}
		// G4
		else if(i==3){
			// first pool
			Mdata[i][0]=Nslow_G4;
			Mdata[i][1]=Uslow_G4;
			Mdata[i][3]=tauF_G4;
			Mdata[i][4]=pslow_G4;
			// second pool
			Mdata[i][5]=Nfast_G4;
			Mdata[i][6]=Ufast_G4;
			Mdata[i][8]=tauF_G4;
			Mdata[i][9]=pfast_G4;
		}
		// G5
		else if(i==4){
			// first pool
			Mdata[i][0]=Nslow_G5;
			Mdata[i][1]=Uslow_G5;
			Mdata[i][3]=tauF_G5;
			Mdata[i][4]=pslow_G5;
			// second pool
			Mdata[i][5]=Nfast_G5;
			Mdata[i][6]=Ufast_G5;
			Mdata[i][8]=tauF_G5;
			Mdata[i][9]=pfast_G5;
		}
	}

	idxS=idxD=idxC0=idxC1=idxC2=0;
	for (i=0; i<N; i++) {
		for (j=0; j<W[i].num; j++) {
			// set up 5 synapse groups according to MF identity
			if(MFgroup[W[i].idx[j]]==0)			{ m=0;	W[i].Gidx[j]=0;	}
			else if(MFgroup[W[i].idx[j]]==1)	{ m=1;	W[i].Gidx[j]=1;	}
			else if(MFgroup[W[i].idx[j]]==2)	{ m=2;	W[i].Gidx[j]=2;	}
			else if(MFgroup[W[i].idx[j]]==3)	{ m=3;	W[i].Gidx[j]=3;	}
			else if(MFgroup[W[i].idx[j]]==4)	{ m=4;	W[i].Gidx[j]=4;	}

			W[i].strgth1[j]=Mdata[m][0];
			W[i].U1[j]=Mdata[m][1];
			W[i].tauD1[j]=Mdata[m][2];
			W[i].tauF1[j]=Mdata[m][3];
			W[i].pf1[j]=Mdata[m][4];
			W[i].strgth2[j]=Mdata[m][5];
			W[i].U2[j]=Mdata[m][6];
			W[i].tauD2[j]=Mdata[m][7];
			W[i].tauF2[j]=Mdata[m][8];
			W[i].pf2[j]=Mdata[m][9];
			W[i].tauG[j]=Mdata[m][10];
			W[i].deltaG[j]=Mdata[m][11];

			// calculate relative amplitude strength
			W[i].kp[j]=W[i].strgth1[j]/(W[i].strgth1[j]+W[i].strgth2[j]);
			W[i].kpp[j]=W[i].strgth2[j]/(W[i].strgth1[j]+W[i].strgth2[j]);

			if(FLAGS.MF_STP==0){
				W[i].tauD1[j]=0;	W[i].tauF1[j]=0; 	W[i].pf1[j]=0;
				W[i].tauD2[j]=0;	W[i].tauF2[j]=0;	W[i].pf2[j]=0;
				W[i].tauG[j]=0;		W[i].deltaG[j]=0;
			}
			if(FLAGS.MF_STP==2){
				std::cerr <<"Mode not supported, use different routine. \n";
				exit(EXIT_FAILURE);
			}
		}
	}

//	**********************************
// 	calculation of GC patterns & adjust thresholds and gains
	calc_GCpattern_theta_gain_full(W, MFpatterns, GCpatterns, theta, GCgain, Tavrg, Gavrg, hGCavrg, hGCvar, CLavrg, NETPARA, FLAGS, r);

//------------------------------------------------------------------------
// 	calculate tau_effs for target patterns
//------------------------------------------------------------------------
	kk=0;
	// calculate tsyn statistics according to SYN identity
	for (n=0; n<Ncl; n++) {
		tempINT[n][0]=0; tempINT[n][1]=0;
		tempDOUBLE[n][0]=0; tempDOUBLE[n][1]=0;
	}
	for (i=0; i<N; i++) {
		for (j=0; j<W[i].num; j++) {
			for (k=0; k<PP; k++) {
				W[i].tsyn1[j][k]=W[i].tauD1[j]*W[i].x1p[j][k];
				W[i].tsyn2[j][k]=W[i].tauD2[j]*W[i].x2p[j][k];

				// slow pools
				for (n=0; n<Ncl; n++) {
					if( W[i].Gidx[j]==n && (W[i].tsyn1[j][k]<W[i].tauD1[j]) ){
						tempINT[n][0]++;
						tempDOUBLE[n][0]+=W[i].tsyn1[j][k];
						tempDOUBLE[n][1]+=W[i].tsyn1[j][k]*W[i].tsyn1[j][k];
					}
					if( W[i].Gidx[j]==n )	tempINT[n][1]++;
				}
			}
		}
	}

	for (n=0; n<Ncl; n++) {
		tsyn_avrg[n][0]=tempDOUBLE[n][0]/double(tempINT[n][0]);
		tsyn_std[n][0]=sqrt(tempDOUBLE[n][1]/double(tempINT[n][0]) - tsyn_avrg[n][0]*tsyn_avrg[n][0]);
		// number of inactive MFs
		foff[n][0]=1.0-double(tempINT[n][0])/double(tempINT[n][1]);
	}

	// calculate tsyn statistics according to MF identity
	for (n=0; n<Ncl; n++) {
		tempINT[n][0]=0; tempINT[n][1]=0;
		tempDOUBLE[n][0]=0; tempDOUBLE[n][1]=0;
	}
	for (i=0; i<N; i++) {
		for (j=0; j<W[i].num; j++) {
			for (k=0; k<PP; k++) {
				// slow pools
				for (n=0; n<Ncl; n++) {
					if( MFgroup[W[i].idx[j]]==n && (W[i].tsyn1[j][k]<W[i].tauD1[j]) ){
						tempINT[n][0]++;
						tempDOUBLE[n][0]+=W[i].tsyn1[j][k];
						tempDOUBLE[n][1]+=W[i].tsyn1[j][k]*W[i].tsyn1[j][k];
					}
					if( MFgroup[W[i].idx[j]]==n )	tempINT[n][1]++;
				}
			}
		}
	}
	for (n=0; n<Ncl; n++) {
		tsyn_avrg[n][1]=tempDOUBLE[n][0]/double(tempINT[n][0]);
		tsyn_std[n][1]=sqrt(tempDOUBLE[n][1]/double(tempINT[n][0]) - tsyn_avrg[n][1]*tsyn_avrg[n][1]);
		// number of inactive MFs
		foff[n][1]=1.0-double(tempINT[n][0])/double(tempINT[n][1]);
	}

//------------------------------------------------------------------------
// 	calculate GC inputs and firing rates at transient onset
//------------------------------------------------------------------------
	for (i=0; i<N; i++) {
		for (k = 0; k < PP; k++) {
			n=k+1;
			if(n> (PP-1)) n=0;
			//	pre-synaptic inputs of GC i
			hGC_onset.at(i,k)=0;
			for (j=0; j<W[i].num; j++) {
				// phasic contribution
				hGC_onset.at(i,k)+=MFpatterns.at(W[i].idx[j],n)*( W[i].strgth1[j]*W[i].x1p[j][k]*W[i].u1p[j][k] + W[i].strgth2[j]*W[i].x2p[j][k]*W[i].u2p[j][k] ) *W[i].qp[j][k];
			}
			GCpatterns_onset.at(i,k)=GCgain[i]*std::max(hGC_onset.at(i,k)-theta[i],0.0);
		}
	}
//------------------------------------------------------------------------
// 	calculate GC inputs and firing rates at different times during transient
// 	only valid for depletion only model
//------------------------------------------------------------------------
	//~ for (i=0; i<N; i++) {
		//~ for (k = 0; k < PP; k++) {
			//~ n=k+1;
			//~ if(n> (PP-1)) n=0;
			//~ //	pre-synaptic inputs of GC i
			//~ hGC_trans(i,k)=0;
			//~ for (j=0; j<W[i].num; j++) {
				//~ expfac1=exp(-ttshift/W[i].tsyn1[j][n]);
				//~ expfac2=exp(-ttshift/W[i].tsyn2[j][n]);
				//~ xtrans1=W[i].x1p[j][n]*(1.0-expfac1) + W[i].x1p[j][k]*expfac1;
				//~ xtrans2=W[i].x2p[j][n]*(1.0-expfac2) + W[i].x2p[j][k]*expfac2;
				//~ // phasic contribution (no desensitisation, no facilitation)
				//~ // hGC_trans(i,k)+=MFpatterns.at(W[i].idx[j],n)*( W[i].strgth1[j]*W[i].U1[j]*xtrans1 + W[i].strgth2[j]*W[i].U2[j]*xtrans2 );
				//~ hGC_trans(i,k)+=MFpatterns.at(W[i].idx[j],n)*( W[i].WU1[j]*xtrans1 + W[i].WU2[j]*xtrans2 );
			//~ }
			//~ GCpatterns_trans(i,k)=GCgain[i]*std::max(hGC_trans(i,k)-theta[i],0.0);
		//~ }
	//~ }
//------------------------------------------------------------------------
// 	set up input signals & target patterns
//------------------------------------------------------------------------
	for (k=0; k<Ndelays; k++) {

		t=0;//=tt=kk=0;
		// 	generate MF input signal
		ns=1;
		for (j=0; j<M; j++) {
			MF_Tinput[j][0]=MFpatterns.at(j,pstart);
			if(MF_Tinput[j][0]<0) 	MF_Tinput[j][0]=0;
		}

		t++;
		while(t <= MAX_data) {
			for (j=0; j<M; j++) 	MF_Tinput[j][t]=MF_Tinput[j][t-1];
			if((t-ttpre)==ts[ns]){
				if((ns % 2)==0) np=pstart;
				else 			np=ptarget;

				for (j=0; j<M; j++) {
					MF_Tinput[j][t]=MFpatterns.at(j,np);
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
	tt=kk=tsvd=tcut=0;
	hGCpeak=GCpeak=0;
	bool_time=bool_tcut=0;
	GCinT.zeros();
	arma::mat GC_trial(N, BINS);
	arma::mat GCtrans(N, CUT_tsteps);

	for (j=0; j<M; j++) 	MF[j]=MFm[j]=MF_Tinput[j][0];
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
	FLAGS.trans_dum=0;
	MFGC_STP_2pools_evolve(W, MF, dt, taumin, N, FLAGS);
	//~ #include "headers/MFGC_STP_2pools_evolve.h"
	FLAGS.trans_dum=FLAGS.trans;

//	*************************
//	compute basis functions
//	*************************
	t=-1;
	//~ bool bool_tten;
	while(t < MAX_data) {
		t++;
		bool_time=( (t % tfrac)==0);
		bool_tten=( (t % 10)==0);
		bool_tfive=( (t % 5)==0);
		bool_tcut=( (t % tcutfrac)==0);

		// update variables for next timestep
		for (j=0; j<M; j++) {
			MF[j]=MF_Tinput[j][t];
		}
//		cout << pp << "\t" << t << endl;
		for (i=0; i<N; i++)	{
			GCm[i]=GC[i];
			for (j=0; j<W[i].num; j++){
				W[i].x1m[j]=W[i].x1[j];
				W[i].u1m[j]=W[i].u1[j];
				W[i].x2m[j]=W[i].x2[j];
				W[i].u2m[j]=W[i].u2[j];
				W[i].qm[j]=W[i].q[j];
			}
		}
		GoCm=GoC;

		//	calculate average GC activity
		GCavrg=0;
		for (i=0; i<N; i++) GCavrg+=GCm[i];
		GCavrg=GCavrg/double(N);

		GoC=GCavrg;

		// evolve GCs
		c=d=0;
		for (i=0; i<N; i++) {
			//	pre-synaptic inputs of GC i
			a=b=0;
			for (j=0; j<W[i].num; j++) {
				// phasic contribution
				a+=MF[W[i].idx[j]]*( W[i].strgth1[j]*W[i].x1m[j]*W[i].u1m[j] + W[i].strgth2[j]*W[i].x2m[j]*W[i].u2m[j] )*W[i].qm[j];
			}
			if(FLAGS.GoC==1){
				hGC[i]=a-GoCm*JEI[i].strgth[0];
			}
			else {
				hGC[i]=a;
			}
			GC[i]=GCm[i]+dtGC*(-GCm[i] + GCgain[i]*std::max(hGC[i]-theta[i],0.0));
			// instantaneous GCs
			//~ GC[i]=GCgain[i]*std::max(hGC[i]-theta[i],0.0);
			c+=hGC[i];
			d+=GC[i];
		}
		// average GC input and rate
		c=c/double(N);
		d=d/double(N);
		if(c>hGCpeak) hGCpeak=c;
		if(d>GCpeak) GCpeak=d;

		// instantaneous MLIs
		GCavrg=0;
		for (i=0; i<N; i++) GCavrg+=GC[i];
		GCavrg=GCavrg/double(N);
		MLI=GCavrg;

		// MF-to-GC STP
		//~ #include "headers/MFGC_STP_2pools_evolve.h"
		MFGC_STP_2pools_evolve(W, MF, dt, taumin, N, FLAGS);

		// set both fast and slow synapse pools with tsyn > tsyncut to their steady state values; tsyn is calculated based on simplified depression-only model
		if(FLAGS.cuttsyn==1){
			for (i=0; i<N; i++) {
				for (j=0; j<W[i].num; j++) {
					tchg1=tchg2=0;
					if( (W[i].tsyn1[j][ptarget]<tsyncut_hi) && (W[i].tsyn1[j][ptarget]>tsyncut_lo) ) {
						tchg1=1;
						if(W[i].tauF1[j]<taumin)	W[i].u1[j]=W[i].U1[j];
						else 						W[i].u1[j]=W[i].U1[j]*(1.0+W[i].tauF1[j]*MF[W[i].idx[j]])/(1.0+W[i].U1[j]*W[i].tauF1[j]*MF[W[i].idx[j]]);
						W[i].x1[j]=1.0/(1.0+W[i].u1[j]*(1-W[i].pf1[j])*W[i].tauD1[j]*MF[W[i].idx[j]]);
					}
					if( (W[i].tsyn2[j][ptarget]<tsyncut_hi) && (W[i].tsyn2[j][ptarget]>tsyncut_lo) ){
						tchg2=1;
						if(W[i].tauF2[j]<taumin)	W[i].u2[j]=W[i].U2[j];
						else 						W[i].u2[j]=W[i].U2[j]*(1.0+W[i].tauF2[j]*MF[W[i].idx[j]])/(1.0+W[i].U2[j]*W[i].tauF2[j]*MF[W[i].idx[j]]);
						W[i].x2[j]=1.0/(1.0+W[i].u2[j]*(1-W[i].pf2[j])*W[i].tauD2[j]*MF[W[i].idx[j]]);
					}
					kp=W[i].strgth1[j]/(W[i].strgth1[j]+W[i].strgth2[j]);
					kpp=W[i].strgth2[j]/(W[i].strgth1[j]+W[i].strgth2[j]);
					if( (tchg1==1) || (tchg2==1) ) {
						if(W[i].tauG[j]<taumin)		W[i].q[j]=1;
						else 						W[i].q[j]=1.0/(1.0+W[i].deltaG[j]*W[i].tauG[j]*( kp*W[i].x1[j]*W[i].u1[j] + kpp*W[i].x2[j]*W[i].u2[j] )*MF[W[i].idx[j]]);
					}
				}
			}
		}
		else if(FLAGS.cuttsyn==2){
			// supporter (only works for 2D-2S configuration)
			for (i=0; i<N; i++) {
				for (j=2; j<W[i].num; j++) {
					W[i].u1[j]=W[i].U1[j];
					W[i].u2[j]=W[i].U2[j];
					W[i].x1[j]=1.0/(1.0+W[i].u1[j]*(1-W[i].pf1[j])*W[i].tauD1[j]*MF[W[i].idx[j]]);
					W[i].x2[j]=1.0/(1.0+W[i].u2[j]*(1-W[i].pf2[j])*W[i].tauD2[j]*MF[W[i].idx[j]]);
					//~ W[i].q[j]=1;
				}
			}
		}
		else if(FLAGS.cuttsyn==3){
			// driver (only works for 2D-2S configuration)
			for (i=0; i<N; i++) {
				for (j=0; j<2; j++) {
					W[i].u1[j]=W[i].U1[j];
					W[i].u2[j]=W[i].U2[j];
					W[i].x1[j]=1.0/(1.0+W[i].u1[j]*(1-W[i].pf1[j])*W[i].tauD1[j]*MF[W[i].idx[j]]);
					W[i].x2[j]=1.0/(1.0+W[i].u2[j]*(1-W[i].pf2[j])*W[i].tauD2[j]*MF[W[i].idx[j]]);
					//~ W[i].q[j]=1;
				}
			}
		}
		else if(FLAGS.cuttsyn==4){
			// slow pools
			for (i=0; i<N; i++) {
				for (j=0; j<W[i].num; j++) {
					if(W[i].tauF1[j]<taumin)	W[i].u1[j]=W[i].U1[j];
					else 						W[i].u1[j]=W[i].U1[j]*(1.0+W[i].tauF1[j]*MF[W[i].idx[j]])/(1.0+W[i].U1[j]*W[i].tauF1[j]*MF[W[i].idx[j]]);
					W[i].x1[j]=1.0/(1.0+W[i].u1[j]*(1-W[i].pf1[j])*W[i].tauD1[j]*MF[W[i].idx[j]]);
					//~ W[i].q[j]=1;
				}
			}
		}
		else if(FLAGS.cuttsyn==5){
			// fast pools
			for (i=0; i<N; i++) {
				for (j=0; j<W[i].num; j++) {
					if(W[i].tauF2[j]<taumin)	W[i].u2[j]=W[i].U2[j];
					else 						W[i].u2[j]=W[i].U2[j]*(1.0+W[i].tauF2[j]*MF[W[i].idx[j]])/(1.0+W[i].U2[j]*W[i].tauF2[j]*MF[W[i].idx[j]]);
					W[i].x2[j]=1.0/(1.0+W[i].u2[j]*(1-W[i].pf2[j])*W[i].tauD2[j]*MF[W[i].idx[j]]);
					//~ W[i].q[j]=1;
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
					for(j= 0; j < W[i].num; j++) STP_file <<  ( W[i].strgth1[j]*W[i].x1m[j]*W[i].u1m[j] + W[i].strgth2[j]*W[i].x2m[j]*W[i].u2m[j] )*W[i].qm[j] <<" ";
					//~ for(j= 0; j < W[i].num; j++) STP_NU_file <<  (W[i].WU1[j]*W[i].x1m[j] + W[i].WU2[j]*W[i].x2m[j])*MF[W[i].idx[j]] <<" ";
				}
				STP_file << std::endl;
				//~ STP_NU_file << std::endl;
			}
		}

		if( (pp==0) & (pp2==0) & (rr==0) & (tt % int(dbint/dt)==0) & (kk<BINS) & ((t-ttpre)*dt>=-0.1) ) {
//			std::cout << tt << "\t" << kk << std::endl;
			time_learn_file << t*dt-Tpre << std::endl;
		}

		if( (tt % int(dbint/dt)==0) & (kk<BINS) & ((t-ttpre)*dt>=-0.1) ) {
			for(i= 0; i < N; i++)  	GC_trial.at(i,kk)=GC[i];
			kk++;
		}

		if (bool_time){
			if(t>= ttsvd_start && t< ttsvd_end){
//				std::cout << t << "\t"<< t*dt << std::endl;
				for(i= 0; i < N; i++) {
					GCinT.at(i, tsvd)= GC[i];
				}
				tsvd++;
			}
		}

		if (bool_tcut){
			if(t> ttcut_start && t<= ttcut_end){
				//~ std::cout << t << "\t"<< t*dt <<  "\t"<< tcut <<std::endl;
				for(i= 0; i < N; i++) {
					GCtrans.at(i, tcut)= GC[i];
				}
				//~ Time(tcut)=t*dt;
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
	//~ std::cout << GCtrans.n_rows << "\t"<< GCtrans.n_cols << std::endl;
	//~ std::cout << GC_trial.n_rows << "\t"<< GC_trial.n_cols << std::endl;

	// subtract GC steady state
	for (i = 0; i < Nred; i++) {
		a=GCtrans.at(i, CUT_tsteps-1);
		for (j = 0; j < CUT_tsteps; j++) {
			GCtrans.at(i, j)=GCtrans.at(i, j)-a;
		}
	}
	// find transient maxima with naive method
	//~ GCmax=arma::max(arma::abs(GCtrans),1);
	//~ arma::ucolvec pkidx=arma::index_max(arma::abs(GCtrans),1);

	// find transient maxima based on zero crossings of derivative
	arma::rowvec GCdiff(CUT_tsteps);
	arma::uvec tcr;
	arma::vec npeaks(Nred);
	arma::uvec pkidx(Nred);
	arma::colvec GCmax(Nred);
	int tmpidx;
	for (i = 0; i <Nred; i++) {
		// compute derivative
		GCdiff=arma::diff(GCtrans.row(i));
		// find zero crossings
		tcr=arma::find(GCdiff % arma::shift( GCdiff, -1) <= 0) +1;
		//~ std::cout<< tcr.n_rows << "\t" << tcr.n_cols<< std::endl;
		// exclude zero crossings that are too late
		for(j= (tcr.n_rows-1); j >=0; j--) {
			if(Time(tcr(j))>tlate) tcr.shed_row(j);
		}
		if(tcr.is_empty()){
			// if there are no peaks simply look for maximum
			GCmax(i)=arma::max(arma::abs(GCtrans.row(i)));
			pkidx(i)=arma::index_max(arma::abs(GCtrans.row(i)));
			npeaks(i)=0;
		}
		else{
			// if there are one or more peaks, choose the largest one
			arma::vec xx(tcr.n_elem);
			for (j = 0; j < tcr.n_elem; j++) {
				// vector of putative maxima
				xx(j)=std::abs( GCtrans(i,tcr(j)) );
			}
			// among putative maxima, find the largest one
			GCmax(i)=arma::max(xx);
			tmpidx=arma::index_max(xx);
			pkidx(i)=tcr(tmpidx);
			GCmax(i)=GCmax(i)*arma::sign(GCtrans(i,pkidx(i)));
			// if amplitude is negative and too small revert to default measure
			if(GCmax(i)<0 && (std::abs(GCmax(i))<ampsmall) ){
				GCmax(i)=arma::max(arma::abs(GCtrans.row(i)));
				pkidx(i)=arma::index_max(arma::abs(GCtrans.row(i)));
			}
			npeaks(i)=tcr.n_elem;
		}
		//~ std::cout<< npeaks(i) << "\t" << GCmax(i) << "\t" << pkidx(i) << std::endl;
	}

	// normalise GC transients
	for (i = 0; i < Nred; i++) {
		for (j = 0; j < CUT_tsteps; j++) {
			GCtrans.at(i, j)=GCtrans.at(i, j)/std::abs(GCmax(i));
		}
	}

	if(FLAGS.normGC==1){
		// compute decay times
		arma::vec t_decay(Nred);
		arma::uvec didx;
		for (i = 0; i < Nred; i++) {
			didx=arma::find( (arma::abs(GCtrans(i,arma::span(pkidx(i), CUT_tsteps-1))) -0.1) > 0.0 , 1,"last")+pkidx(i);
			t_decay(i)=Time(didx(0));
		}
		//~ std::cout << std::endl;
		//~ for(i= 0; i < 22; i++) {
			//~ std::cout << i << "\t" << GCmax(i) << "\t" << t_decay(i) << std::endl;
		//~ }

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
		//~ for(i= 0; i < 22; i++) {
			//~ std::cout << i << "\t" << status(i) << "\t" << t_decay(i) << std::endl;
		//~ }

		//~ if(rr==0 & (FLAGS.mode==99 || FLAGS.mode==77) ){
			//~ for(i= 0; i < Nred; i++) {
				//~ t_decay_file << i << "\t" << GCmax(i) << "\t" << pkidx(i) << "\t" << pkidx(i)*dt << "\t" << t_decay(i) << std::endl;
			//~ }
			std::cout << std::endl;
			if(FLAGS.cutGC==1)		std::cout <<"Cutting GC transients longer than "<< GCcut << " s." << std::endl;
			else if(FLAGS.cutGC==2) 	std::cout <<"Cutting GC transients shorter than "<< GCcut << " s." << std::endl;
			std::cout <<"Number of GC transients after cut: "<< GCtrans.n_rows << std::endl;
			std::cout <<"Number of transient time bins: "<< GCtrans.n_cols << std::endl;
	}

	if(rr==0 & (FLAGS.mode==99 || FLAGS.mode==77) ){
		// write transients to file
		GCtrans.save(GCtrans_file,arma::raw_ascii);
		Time.save(timeGCtrans_file,arma::raw_ascii);
	}

	//	*************************
	//	do learning
	//	*************************
	err_final=0;
	for (k=0; k<Ndelays; k++) {

		err_final+=learn_weights_DS(GC_trial, J, Jidx, PC, delays[k], interv_min[k], interv_max[k], dt, Tpre, PCsp, CFsp, CFwidth, CFscale, errWeight, J2weight, JI, alpha, dbint, N, Nred, Ntrials,
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
			else if(k==7){
				for (t=0; t<BINS; t++) PC_010_file << PC[t] << "\t";
				PC_010_file << endl;
			}
		}
		// write final PC synaptic weights to file
		if(FLAGS.mode==77 | FLAGS.mode==99) {
			for (i=0; i<N; i++) J_file << J.strgth[i] << "\t";
			J_file << endl;
			for(i= 0; i < Nred; i++) 	Jidx_file << Jidx(i)+1 << "\t";
			Jidx_file << endl;
		}

	} // end of delays-loop

	//	*************************
// 	do temporal PCA
	Dim_GCt=calc_Tdim(GCinT, SVD_tsteps);
	//~ Dim_GCt=NAN;

	//	*************************
// 	do pattern PCA
	// fill GC-matrix
	for (i = 0; i < N; i++) {
		for (k = 0; k < PP; k++) {
			GCpatterns_trans.at(i,k)=0;
		}
	}
	calc_Pdim(MFpatterns, GCpatterns, GCpatterns_onset, GCpatterns_trans, PP, N, M, Dim_MF, Dim_GCss, Dim_GConset, Dim_GCtrans);
	//~ Dim_MF=NAN; Dim_GCss=NAN; Dim_GConset=NAN;	Dim_GCtrans=NAN;

//	*************************
// 	write stuff to files

	if((FLAGS.mode==77 | FLAGS.mode==99) & (rr==0)){
		for(i= 0; i < 100; i++) {
			for(j= 0; j < W[i].num; j++) {
				STPpara_file << i+1 << "\t" << W[i].idx[j]+1 <<"\t" << W[i].Gidx[j] <<"\t" << MFgroup[W[i].idx[j]] << "\n";
			}
		}
		// write the distribution of slow time constants corresponding to target pattern to file (driver and supporters are mixed!)
		for (i=0; i<N; i++) {
			for (j=0; j<W[i].num; j++) {
				if(FLAGS.cuttsyn==1){
					if( (W[i].tsyn1[j][ptarget]<tsyncut_hi) && (W[i].tsyn1[j][ptarget]>tsyncut_lo) ) {
						W[i].tsyn1[j][ptarget]=0;
					}
				}
				else if(FLAGS.cuttsyn==2){
					if( j>1) {
						W[i].tsyn1[j][ptarget]=0;
					}
				}
				else if(FLAGS.cuttsyn==3){
					if( j<2) {
						W[i].tsyn1[j][ptarget]=0;
					}
				}
				else if(FLAGS.cuttsyn==4){
					W[i].tsyn1[j][ptarget]=0;
				}
				tsyndistrb_file << W[i].tsyn1[j][ptarget] << "\t";
				//~ temp=(MFpatterns[W[i].idx[j]][ptarget]-MFpatterns[W[i].idx[j]][pstart])/((1+W[i].alpha1[j]*MFpatterns[W[i].idx[j]][ptarget])*(1+W[i].alpha1[j]*MFpatterns[W[i].idx[j]][pstart]));
				//~ B1=W[i].strgth1[j]*W[i].U1[j]*W[i].tauD1[j]*W[i].U1[j]*(1-W[i].pf1[j])*MFpatterns[W[i].idx[j]][ptarget]*temp;
				//~ B1=W[i].WU1[j]*W[i].alpha1[j]*MFpatterns[W[i].idx[j]][ptarget]*(MFpatterns[W[i].idx[j]][ptarget]-MFpatterns[W[i].idx[j]][pstart])/((1+W[i].alpha1[j]*MFpatterns[W[i].idx[j]][ptarget])*(1+W[i].alpha1[j]*MFpatterns[W[i].idx[j]][pstart]));
				//~ Ampdistrb_file << B1 << "\t";
				Ampdistrb_file << NAN << "\t";
			}
			tsyndistrb_file << "\n";
			Ampdistrb_file << "\n";
		}
	}

	// REWRITE !
	err_final_file << MF_avrg_G1 << "\t" << MF_std_G1 << "\t" << MF_avrg_G2 << "\t" << MF_std_G2 << "\t" << MF_avrg_G3 << "\t" << MF_std_G3 << "\t" << MF_avrg_G4 << "\t" << MF_std_G4 << "\t" << MF_avrg_G5 << "\t" << MF_std_G5 << "\t";
	err_final_file 	<< alpha*sqN << "\t" << NAN << "\t" << Tavrg << "\t" << Gavrg << "\t" << CLavrg << "\t" << JI << "\t" << err_final/double(Ndelays) << std::endl;


	for (n = 0; n < Ncl; n++) {
		tsyn_file << tsyn_avrg[n][0] << "\t" << tsyn_std[n][0] << "\t" << foff[n][0] << "\t" ;
	}
	tsyn_file << std::endl;
	for (n = 0; n < Ncl; n++) {
		tsyn_MFsort_file << tsyn_avrg[n][1] << "\t" << tsyn_std[n][1] << "\t" << foff[n][1] << "\t" ;
	}
	tsyn_MFsort_file << std::endl;

	trans_file << hGCpeak << "\t" << GCpeak << std::endl;

	dim_file << Dim_GCt<< "\t"<< Dim_MF << "\t"<< Dim_GCss << "\t"<< NAN << "\t"<< Dim_GConset << std::endl;

	learnparam_file << Ntrials-1<< "\t" << CFsp<< "\t"<< CFscale << "\t"<< CFwidth << "\t"<< errWeight << "\t"<< alpha_learn << "\t"<< JI
					<< "\t"<< Nreal << "\t"<< Npara << "\t"<< Npara2 << "\t"<< idum << "\t" << avrgJIE<< "\t" << avrgJEI << "\t" << NAN
					<< "\t" << gain_global << "\t" << nfrac << "\t" << GCtarget << "\t" << NAN << "\t" << GCnormfac << "\t" << Tfac_gc
					<< "\t" << NAN << "\t" << MFUcorr << "\t" << Ncl
					<< "\t" << pact[0] << "\t" << pact[1] << "\t" << pact[2] << "\t" << pact[3] << "\t" << pact[4]
					<< "\t" << grouplims[0] << "\t" << grouplims[1] << "\t" << grouplims[2] << "\t" << grouplims[3] << "\t" << grouplims[4] << "\t" << grouplims[5]
					<< "\t" << J2weight << std::endl;

	stats_file << hGCavrg << "\t" << hGCvar << "\t" << Tavrg << "\t" << Gavrg << "\t" << CLavrg << "\t" << JI <<std::endl;

	for(i= 0; i < N; i++){
		for(j= 0; j < W[i].num; j++) {
			delete [] W[i].x1p[j];		delete [] W[i].x2p[j];
			delete [] W[i].u1p[j];		delete [] W[i].u2p[j];
			delete [] W[i].qp[j];
			delete [] W[i].tsyn1[j];	delete [] W[i].tsyn2[j];
		}
		delete [] W[i].idx; 	delete [] W[i].Gidx;
		delete [] W[i].strgth1;	delete [] W[i].strgth2;
		delete [] W[i].x1p;		delete [] W[i].x2p;
		delete [] W[i].u1p;		delete [] W[i].u2p;
		delete [] W[i].qp;
		delete [] W[i].x1;		delete [] W[i].x2;		delete [] W[i].x1m;		delete [] W[i].x2m;
		delete [] W[i].u1;		delete [] W[i].u2;		delete [] W[i].u1m;		delete [] W[i].u2m;
		delete [] W[i].q;		delete [] W[i].qm;
		delete [] W[i].U1;		delete [] W[i].tauD1;	delete [] W[i].tauF1;	delete [] W[i].pf1;
		delete [] W[i].U2;		delete [] W[i].tauD2;	delete [] W[i].tauF2;	delete [] W[i].pf2;
		delete [] W[i].tauG;	delete [] W[i].deltaG;
		delete [] W[i].kp;	delete [] W[i].kpp;

		//~ delete [] W[i].alpha1;	delete [] W[i].WU1; 	delete [] W[i].Upf1;
		//~ delete [] W[i].alpha2;	delete [] W[i].WU2; 	delete [] W[i].Upf2;
		delete [] W[i].tsyn1;	delete [] W[i].tsyn2;
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
	//~ STP_NU_file.close();
	MLI_file.close();
	GoC_file.close();
	CF_file.close();
	time_file.close();
	time_learn_file.close();
	STPpara_file.close();
	GCtrans_file.close();
	timeGCtrans_file.close();

	err_final_file.close();
	PC_010_file.close();
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
	trans_file.close();
	tsyn_file.close();
	tsyndistrb_file.close();
	Ampdistrb_file.close();
	tsyn_MFsort_file.close();
	//~ tsyn_file_D.close();
	dim_file.close();
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