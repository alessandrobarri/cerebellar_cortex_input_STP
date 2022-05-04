	// ----------------------------------
	// recovery from depression
	const double Tau_slow=2.0;					// from Hallermann et al. 2010
	const double Tau_fast=20.0/1000.0;			// from Hallermann et al. 2010

	// ----------------------------------
	// desensitisation
	const double tauG=100.0/1000.0;				// from Hallermann et al. 2010
	const double deltaG=0.1;					// corresponds to Hallermann et al. 2010

	// ----------------------------------
	// synaptic parameters D/S
	const double Nslow_D=3.5;
	const double Uslow_D=0.8;
	//~ const double Uslow_D=0.85;
	const double pslow_D=0.6;

	const double Nfast_D=14;
	const double Ufast_D=0.6;
	const double pfast_D=0;

	const double Nslow_S=4;
	const double Uslow_S=0.4;
	//~ const double Uslow_S=0.15;
	//~ const double Uslow_S=Uslow_D;				// for figure 4 (version 8), third row
	const double pslow_S=0.6;

	const double Nfast_S=6;
	const double Ufast_S=0.2;
	//~ const double Ufast_S=0.1;
	//~ const double Ufast_S=Ufast_D;				// for figure 4 (version 8), third row
	const double pfast_S=0;

	const double U_dev_fast=0.025;						// SD of release prob distribution for fast pools
	const double U_dev_slow=0.05;						// SD of release prob distribution for slow pools

	// ----------------------------------
	// synaptic parameters (G1-G5)
	const double Nslow_G1=4;								// Nslow, Nfast, Uslow, Ufast derived from Chabrol et al. 2015
	const double Uslow_G1=0.9;
	const double pslow_G1=0.6;
	const double Nfast_G1=16;
	const double Ufast_G1=0.72;
	const double pfast_G1=0;
	const double tauF_G1=12.0/1000.0;						// tauF from Saviane et al. 2006
//	const double tauF_G1=0/1000.0;

	const double Nslow_G2=3;
	const double Uslow_G2=0.8;
	const double pslow_G2=0.6;
	const double Nfast_G2=12;
	const double Ufast_G2=0.55;
	const double pfast_G2=0;
	const double tauF_G2=tauF_G1;

	const double Nslow_G3=4;
	const double Uslow_G3=0.4;
	const double pslow_G3=0.6;
	const double Nfast_G3=6;
	const double Ufast_G3=0.35;
	const double pfast_G3=0;
	const double tauF_G3=0;

	const double Nslow_G4=0;
	const double Uslow_G4=0.2;
	const double pslow_G4=0.6;
	const double Nfast_G4=10;
	const double Ufast_G4=0.3;
	const double pfast_G4=0;
	const double tauF_G4=tauF_G1;

	const double Nslow_G5=3;
	const double Uslow_G5=0.4;
	const double pslow_G5=0.6;
	const double Nfast_G5=12;
	const double Ufast_G5=0.15;
	const double pfast_G5=0;
	const double tauF_G5=30.0/1000.0;

	// ----------------------------------
	// for seamless uniform rates
	const double NUmax=270;					// figure 5, should be this to get same average rates as in Gaussian case
	const double NUmin=5;
	arma::vec NUlims = arma::linspace<arma::vec>(NUmin, NUmax, Ncl+1); 			// for equal 'slices'
	//~ arma::vec NUlims = {5, 84.5, 150.75, 203.75, 243.5, 270};					// for 'slices' with decreasing SD

	// for uniform mixture
	arma::vec NU_avrg={31.5, 84.5, 137.5, 190.5, 243.5};
	arma::vec NU_sd= arma::linspace<arma::vec>(15.3, 5, 5);
	//~ arma::vec grouplims = {0, 0.2, 0.4, 0.6, 0.8, 1.0};
	arma::vec grouplims = {0, 0.15, 0.3, 0.45, 0.6, 1.0};

	// ----------------------------------
	// synaptic parameters 5 clusters toy model
	arma::vec NSLOW = { 4, 4, 4, 4, 4};
	arma::vec NFAST = { 6, 6, 8, 12, 16};
	//~ arma::vec NFAST = { 6, 6, 6, 6, 6};
	arma::vec USLOW = {0.15, 0.325, 0.5, 0.675, 0.85};
	arma::vec UFAST = {0.1, 0.1, 0.25, 0.325, 0.6};
	arma::vec PSLOW = {0.6, 0.6, 0.6, 0.6, 0.6};
	arma::vec PFAST = {0, 0, 0, 0, 0};

	double nu0=NU_avrg[4]; double sig0=NU_sd[4];

	double nu1=NU_avrg[3]; double sig1=NU_sd[3];

	double nu2=NU_avrg[2]; double sig2=NU_sd[2];

	double nu3=NU_avrg[1]; double sig3=NU_sd[1];

	double nu4=NU_avrg[0]; double  sig4=NU_sd[0];
