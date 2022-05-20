	// ----------------------------------
	// recovery from depression
	const double Tau_slow=2.0;					// from Hallermann et al. 2010
	const double Tau_fast=20.0/1000.0;			// from Hallermann et al. 2010

	// ----------------------------------
	// desensitisation
	const double tauG=100.0/1000.0;				// from Hallermann et al. 2010
	const double deltaG=0.1;					// produces desensitisation similar to the levels observed in Hallermann et al. 2010

	// ----------------------------------
	// synaptic parameters full model (G1-G5)
	const double Nslow_G1=4;								// Nslow, Nfast, Uslow, Ufast derived from Chabrol et al. 2015
	const double Uslow_G1=0.9;
	const double pslow_G1=0.6;
	const double Nfast_G1=16;
	const double Ufast_G1=0.72;
	const double pfast_G1=0;
	const double tauF_G1=12.0/1000.0;						// tauF from Saviane et al. 2006

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
	// MF type probability of occurence
	const double pG1=0.06;
	const double pG2=0.16;
	const double pG3=0.38;
	const double pG4=0.24;
	//~ const double pG5=0.16; (pG5=1-pG1-pG2-pG3-pG4)

	// ----------------------------------
	// MF rate parameters for full (5 group) model
	double MF_avrg_G1=200.0;
	double MF_avrg_G2=200.0;
	double MF_avrg_G3=20.0;
	double MF_avrg_G4=20.0;
	double MF_avrg_G5=20.0;

	double MF_std_G1=50.0;
	double MF_std_G2=50.0;
	double MF_std_G3=20.0;
	double MF_std_G4=20.0;
	double MF_std_G5=20.0;
