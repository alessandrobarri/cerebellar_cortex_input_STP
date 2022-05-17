void generate_relprob_DS(double Ntoy_fast[], double Ntoy_slow[], double Utoy_slow[], int Utoy_slow_ID[], struct netpara NETPARA, struct synpara SYNPARA, struct flags FLAGS, gsl_rng * r) {

	double a,b;
	double USLOW;
	int i;
	std::vector <double> dummyD;

	const int NSYN_2=2*NETPARA.NSYN;
	// -----------------------------------
	// generate uniform or beta release probability distribution
	// -----------------------------------
	std::vector <size_t> idx_sorted;
	if(FLAGS.Udistrb>1){
		if(FLAGS.Udistrb==2){
			for (i = 0; i < NSYN_2; i++) Utoy_slow[i]=gsl_ran_flat(r, 0.1, 0.9);
		}
		else if (FLAGS.Udistrb==3){
			// beta distribution with mean=0.5 and SD=0.25
			a=(1-0.5)*0.5*0.5/(0.25*0.25) - 0.5;
			b=a*( 1.0/0.5 -1);
			for (i = 0; i < NSYN_2; i++) Utoy_slow[i]=gsl_ran_beta(r, a, b);
		}
		for (i = 0; i < NSYN_2; i++) 	dummyD.push_back(Utoy_slow[i]);
		idx_sorted=sort_indices(dummyD);
		for (i = 0; i < NSYN_2; i++) 	Utoy_slow[i]=dummyD[idx_sorted[i]];
		dummyD.clear();
		idx_sorted.clear();
		for (i = 0; i < NSYN_2/5; i++)				Ntoy_fast[i]=SYNPARA.NFAST[0];
		for (i = NSYN_2/5; i < NSYN_2*2/5; i++) 	Ntoy_fast[i]=SYNPARA.NFAST[1];
		for (i = NSYN_2*2/5; i < NSYN_2*3/5; i++) 	Ntoy_fast[i]=SYNPARA.NFAST[2];
		for (i = NSYN_2*3/5; i < NSYN_2*4/5; i++) 	Ntoy_fast[i]=SYNPARA.NFAST[3];
		for (i = NSYN_2*4/5; i < NSYN_2; i++) 		Ntoy_fast[i]=SYNPARA.NFAST[4];

	}
	else if(FLAGS.Udistrb==1 && FLAGS.groups==55){
		// generate beta mixture release probabilies
		double U_dev_slow=SYNPARA.U_dev_slow;
		USLOW=SYNPARA.USLOW[0];
		a=(1-USLOW)*USLOW*USLOW/(U_dev_slow*U_dev_slow) - USLOW;
		b=a*( 1.0/USLOW -1);
		for (i = 0; i < NSYN_2/5; i++){
			Utoy_slow_ID[i]=0;
			Utoy_slow[i]=gsl_ran_beta(r, a, b);
			Ntoy_fast[i]=SYNPARA.NFAST[0];
		}

		USLOW=SYNPARA.USLOW[1];
		a=(1-USLOW)*USLOW*USLOW/(U_dev_slow*U_dev_slow) - USLOW;
		b=a*( 1.0/USLOW -1);
		for (i = NSYN_2/5; i < NSYN_2*2/5; i++) {
			Utoy_slow_ID[i]=1;
			Utoy_slow[i]=gsl_ran_beta(r, a, b);
			Ntoy_fast[i]=SYNPARA.NFAST[1];
		}

		USLOW=SYNPARA.USLOW[2];
		a=(1-USLOW)*USLOW*USLOW/(U_dev_slow*U_dev_slow) - USLOW;
		b=a*( 1.0/USLOW -1);
		for (i = NSYN_2*2/5; i < NSYN_2*3/5; i++) {
			Utoy_slow_ID[i]=2;
			Utoy_slow[i]=gsl_ran_beta(r, a, b);
			Ntoy_fast[i]=SYNPARA.NFAST[2];
		}

		USLOW=SYNPARA.USLOW[3];
		a=(1-USLOW)*USLOW*USLOW/(U_dev_slow*U_dev_slow) - USLOW;
		b=a*( 1.0/USLOW -1);
		for (i = NSYN_2*3/5; i < NSYN_2*4/5; i++) {
			Utoy_slow_ID[i]=3;
			Utoy_slow[i]=gsl_ran_beta(r, a, b);
			Ntoy_fast[i]=SYNPARA.NFAST[3];
		}

		USLOW=SYNPARA.USLOW[4];
		a=(1-USLOW)*USLOW*USLOW/(U_dev_slow*U_dev_slow) - USLOW;
		b=a*( 1.0/USLOW -1);
		for (i = NSYN_2*4/5; i < NSYN_2; i++) {
			Utoy_slow_ID[i]=4;
			Utoy_slow[i]=gsl_ran_beta(r, a, b);
			Ntoy_fast[i]=SYNPARA.NFAST[4];
		}
	}

	return;
}
