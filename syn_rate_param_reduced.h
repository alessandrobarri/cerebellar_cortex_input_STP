	// ----------------------------------
	// recovery from depression
	const double Tau_slow=2.0;					// from Hallermann et al. 2010
	const double Tau_fast=20.0/1000.0;			// from Hallermann et al. 2010

	// ----------------------------------
	// synaptic parameters for 2 group reduced model (D/S)
	const double Nslow_D=3.5;
	const double Uslow_D=0.8;
	const double pslow_D=0.6;

	const double Nfast_D=14;
	const double Ufast_D=0.6;
	const double pfast_D=0;

	const double Nslow_S=4;
	const double Uslow_S=0.4;
	const double pslow_S=0.6;

	const double Nfast_S=6;
	const double Ufast_S=0.2;
	const double pfast_S=0;

	// ----------------------------------
	// synaptic parameters for 5 group reduced model
	arma::vec NSLOW = { 4, 4, 4, 4, 4};
	arma::vec NFAST = { 6, 6, 8, 12, 16};
	arma::vec USLOW = {0.15, 0.325, 0.5, 0.675, 0.85};
	arma::vec UFAST = {0.1, 0.1, 0.25, 0.325, 0.6};
	arma::vec PSLOW = {0.6, 0.6, 0.6, 0.6, 0.6};
	arma::vec PFAST = {0, 0, 0, 0, 0};

	// ----------------------------------
	// MF rate parameters for 2 group reduced model (D/S)
	// Fig. 3 & 4a
	double MF_avrg_D=200.0;
	double MF_std_D=15.0;
	double MF_avrg_S=25.0;
	double MF_std_S=15.0;

	//~ // Fig. 4b
	//~ double MF_avrg_D=200.0;
	//~ double MF_std_D=15.0;
	//~ double MF_avrg_S=70.0;
	//~ double MF_std_S=15.0;

	//~ // Fig. 4c
	//~ double MF_avrg_D=100.0;
	//~ double MF_std_D=50.0;
	//~ double MF_avrg_S=25.0;
	//~ double MF_std_S=15.0;

	//~ // Fig. 4d (set cuttsyn("cutD") )
	//~ double MF_avrg_D=200.0;
	//~ double MF_std_D=15.0;
	//~ double MF_avrg_S=25.0;
	//~ double MF_std_S=15.0;

	//~ // Fig. 6h-k
	//~ double MF_avrg_D=200.0;
	//~ double MF_std_D=10.0;
	//~ double MF_avrg_S=20.0;
	//~ double MF_std_S=15.0;

	// ----------------------------------
	// MF rate parameters for 5 group reduced model
	// for seamless uniform rates
	const double NUmax=270;					// figure 5, should be this to get same average rates as in Gaussian case
	const double NUmin=5;
	arma::vec NUlims = arma::linspace<arma::vec>(NUmin, NUmax, Ncl+1); 			// for equal 'slices' (Fig. 5a-f)
	//~ arma::vec NUlims = {5, 84.5, 150.75, 203.75, 243.5, 270};					// for 'slices' with decreasing SD

	// for uniform mixture
	arma::vec NU_avrg={31.5, 84.5, 137.5, 190.5, 243.5};			// Fig. 5g,h
	arma::vec NU_sd= arma::linspace<arma::vec>(15.3, 5, 5);
	arma::vec grouplims = {0, 0.15, 0.3, 0.45, 0.6, 1.0};

	double nu0=NU_avrg[4]; double sig0=NU_sd[4];

	double nu1=NU_avrg[3]; double sig1=NU_sd[3];

	double nu2=NU_avrg[2]; double sig2=NU_sd[2];

	double nu3=NU_avrg[1]; double sig3=NU_sd[1];

	double nu4=NU_avrg[0]; double sig4=NU_sd[0];

	// ----------------------------------
	// MF rate parameters for 2 group reduced model grid search (Fig. S4)
	const double nuDmax=250,nuDmin=10;
//	const double sigDmax=80,sigDmin=10;
//	const double sigDmax=150,sigDmin=80;
	const double sigDmax=120,sigDmin=10;

	//~ const double nuSmax=70,nuSmin=5;
	//~ const double sigSmax=60,sigSmin=3;
	const double nuSmax=100,nuSmin=5;
	const double sigSmax=70,sigSmin=3;
