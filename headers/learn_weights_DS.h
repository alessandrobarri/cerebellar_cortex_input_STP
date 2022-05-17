double learn_weights_DS(arma::mat& GC_trial, struct presyn& J, arma::uvec& Jidx, double PC[], double delays, double interv_min, double interv_max, double dt, double Tpre, double PCsp, double CFsp, double CFwidth, double CFscale, double errWeight, double J2weight, double JI, double alpha, double dbint, int N, int Nred, int Ntrials, int BINS, int MAX_data, int ttpre, int rr, int pp, int pp2, struct flags FLAGS, std::ofstream& CF_file, std::ofstream& htarget_BINS_file, gsl_rng * r) {

	int i,m,t,tt,kk,jidx;
	double YL[N],YLm[N],LAMBm[N],LAMB[N],Jm[N];
	double Werr2[BINS],PT_BINS[BINS],CFrate[BINS],DELTA[BINS];
	double PT[MAX_data+1];
	double a,shift,delJ,gamma,Err_tot,Err_tot_final,hPC,GCavrg;
	double Werr2_mean;
	std::vector <double> err_vec;

	double sqN=sqrt(N);
	double alphdbint=alpha*dbint;

	//	set initial values
	for (i=0; i<N; i++) {
		Jm[i]=J.strgth[i]=JI;
		LAMB[i]=LAMBm[i]=1.0;
		YL[i]=J.strgth[i];
	}
	for (m=0; m<Ntrials; m++) {

		if(FLAGS.mode==77){
			// draw CF interval from uniform distribution
			shift=std::round(gsl_ran_flat(r, interv_min, interv_max)/dbint)*dbint;		// discretise for delta rule
		}
		else{
			shift =delays;
		}

		t=tt=kk=0;
		// 	generate CF (teaching) signal
		PT[0]=PCsp;
		t++;
		while(t <= MAX_data) {
			// 'delta' target signal
			if( abs(dt*t-(Tpre+shift)) <=0.001)			PT[t]=0.0;
			else 										PT[t]=PCsp;

			if(dt*t<Tpre)	PT[t]=PCsp;

			// downsampling of PT
			if( (tt % int(dbint/dt)==0) & (kk<BINS) & ((t-ttpre)*dt>=-0.1) ) {
				PT_BINS[kk]=PT[t];
				if(PT_BINS[kk]<1.0)	Werr2[kk]=errWeight*errWeight;
				else					Werr2[kk]=1.0;
				if(FLAGS.mode==77){
					// if target function is different on every trial, write 100 samples to file
					if((pp==0) & (pp2==0) & (rr==0) & (m<100)){
						htarget_BINS_file << PT_BINS[kk] << "\t";
					}
				}
				else{
					// if target function is identical on every trial, write it to file once
					if((pp==0) & (pp2==0) & (rr==0) & (m==0)){
						htarget_BINS_file << PT_BINS[kk] << "\t";
					}
				}
				kk++;
			}
			tt++;

			t++;
		}
		if(FLAGS.mode==77){
			// if target function is different on every trial, write 100 samples to file
			if((pp==0) & (pp2==0) & (rr==0) & (m<100)){
				htarget_BINS_file << endl;
			}
		}
		else{
			// if target function is identical on every trial, write it to file once
			if((pp==0) & (pp2==0) & (rr==0) & (m==0)){
				htarget_BINS_file << endl;
			}
		}

		// normalise teaching signal weights
		Werr2_mean=0;
		for (t=0; t<BINS; t++) {
			Werr2_mean+=Werr2[t];
		}
		Werr2_mean=Werr2_mean/double(BINS);
		for (t=0; t<BINS; t++) {
			Werr2[t]=Werr2[t]/Werr2_mean;
		}

		Err_tot=0;
		for (t=0; t<BINS; t++) {

			// evolve PC
			a=0;
			//	pre-synaptic inputs to PC from GCs
			for (i=0; i<Nred; i++) a+=GC_trial.at(i,t)*J.strgth[Jidx(i)];

			//	pre-synaptic inputs from MLI
			GCavrg=0;
			for (i=0; i<Nred; i++) GCavrg+=GC_trial.at(i,t);
			GCavrg=GCavrg/double(N);

			hPC=a/sqN-GCavrg*JI*sqN +PCsp;
			Err_tot+=(hPC-PT_BINS[t])*(hPC-PT_BINS[t])*Werr2[t];
			CFrate[t]=std::max( CFsp + CFscale*(hPC-PT_BINS[t]) ,0.0);
			DELTA[t]=( CFsp -CFrate[t])*Werr2[t];

			// calculate PC activity in final trial
			if(m==Ntrials-1 ) PC[t]=std::max(hPC,0.0);

			if(FLAGS.mode==99){
				if((pp==0) & (pp2==0) & (rr==0) & (m<100)){
					CF_file << CFrate[t] << "\t";
				}
			}

		}
		if(FLAGS.mode==99){
			if((pp==0) & (pp2==0) & (rr==0) & (m<100)){
				CF_file << std::endl;
			}
		}

//		if( ( (m % nout) ==0) || m==Ntrials-1 ) {
		err_vec.push_back(0.5*Err_tot/double(BINS));

		// update weights
		if(m<Ntrials-1 ){
			for (i=0; i<N; i++) {
				YLm[i]=YL[i];
				LAMBm[i]=LAMB[i];
				Jm[i]=J.strgth[i];
			}
			// loop over non-zero GCs only
			for (i=0; i<Nred; i++) {
				jidx=Jidx(i);
				delJ=0;
				for (t=0; t<BINS; t++) {
					// with CF input
					delJ+=DELTA[t]*GC_trial.at(i,t);
				}
				YL[i]=Jm[jidx]+alphdbint*delJ - J2weight*2.0*(Jm[jidx]-JI);
				LAMB[i]=0.5*(1.0+sqrt(1.0+4.0*LAMBm[i]*LAMBm[i]));
				gamma=(1.0-LAMBm[i])/LAMB[i];
				J.strgth[jidx]=std::max( (1.0-gamma)*YL[i]+gamma*YLm[i] ,0.0);
				// additional acceleration by resetting LAMB
				if((J.strgth[jidx]-Jm[jidx])*delJ<0) LAMB[i]=1.0;
			}
		}

	} // end of trial-loop

	Err_tot_final=err_vec.back();

	return Err_tot_final;
}
