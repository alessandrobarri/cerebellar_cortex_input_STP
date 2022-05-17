void calc_GCpattern_theta_gain_full(struct presyn_2pools W[], arma::mat& MFpatterns, arma::mat& GCpatterns, double theta[], double GCgain[], double& Tavrg, double& Gavrg, double& hGCavrg, double& hGCvar, double& CLavrg, struct netpara NETPARA, struct flags FLAGS, gsl_rng * r) {

	double a,b,c,e,d,f,a2,b2,e2,f2;
	double temp,kp,kpp;
	int i,j,j2,j3,k,n;

	const int N=NETPARA.N;
	const int PP=NETPARA.PP;
	const double nfrac=NETPARA.nfrac;
	const double GCtarget=NETPARA.GCtarget;
	const double theta_global=NETPARA.theta_global;
	const double taumin=NETPARA.taumin;

	double avrgGC[N];
	double avrgGC_act[N];

	arma::mat hGCpatterns(N,PP);
	double* temphGC=new double [N*PP];
	double* temphGC2=new double [PP];

	// steady state MF-GC synaptic weights
	for (i=0; i<N; i++) {
		for (j=0; j<W[i].num; j++) {
			kp=W[i].kp[j];
			kpp=W[i].kpp[j];
			for (k = 0; k < PP; k++) {
				if(W[i].tauF1[j]<taumin)	W[i].u1p[j][k]=W[i].U1[j];
				else						W[i].u1p[j][k]=W[i].U1[j]*(1.0+W[i].tauF1[j]*MFpatterns.at(W[i].idx[j],k))/(1.0+W[i].U1[j]*W[i].tauF1[j]*MFpatterns.at(W[i].idx[j],k));
				if(W[i].tauF2[j]<taumin)	W[i].u2p[j][k]=W[i].U2[j];
				else						W[i].u2p[j][k]=W[i].U2[j]*(1.0+W[i].tauF2[j]*MFpatterns.at(W[i].idx[j],k))/(1.0+W[i].U2[j]*W[i].tauF2[j]*MFpatterns.at(W[i].idx[j],k));
				W[i].x1p[j][k]=1.0/(1.0+W[i].u1p[j][k]*(1.0-W[i].pf1[j])*W[i].tauD1[j]*MFpatterns.at(W[i].idx[j],k));
				W[i].x2p[j][k]=1.0/(1.0+W[i].u2p[j][k]*(1.0-W[i].pf2[j])*W[i].tauD2[j]*MFpatterns.at(W[i].idx[j],k));
				if(W[i].tauG[j]<taumin) 	W[i].qp[j][k]=1.0;
				else 						W[i].qp[j][k]=1.0/(1.0+W[i].deltaG[j]*W[i].tauG[j]*( kp*W[i].x1p[j][k]*W[i].u1p[j][k] + kpp*W[i].x2p[j][k]*W[i].u2p[j][k] )*MFpatterns.at(W[i].idx[j],k));
			}
		}
	}

	// excitatory steady state inputs to GCs
	for (i=0; i<N; i++) {
		for (k = 0; k < PP; k++) {
			//	pre-synaptic inputs of GC i
			hGCpatterns.at(i,k)=0;
			for (j=0; j<W[i].num; j++) {
				// phasic contribution
				hGCpatterns.at(i,k)+=MFpatterns.at(W[i].idx[j],k)*( W[i].strgth1[j]*W[i].x1p[j][k]*W[i].u1p[j][k] + W[i].strgth2[j]*W[i].x2p[j][k]*W[i].u2p[j][k] )*W[i].qp[j][k];
			}
			j=k+i*PP;
			temphGC[j]=hGCpatterns.at(i,k);
		}
	}

	// adjust threshold globaly to produce required (steady state) GC coding-level
	if(FLAGS.threshold==4){
		std::sort(temphGC, temphGC + N*PP);
		n=0;
		b=double(N)*double(PP)*(1.0-nfrac);
		for (i=0; i<N; i++) {
			for (k = 0; k < PP; k++) {
				n++;
				if(n >= b) {
					temp=temphGC[k+i*PP];
					i=N;
					break;
				}
			}
		}
		for (i=0; i<N; i++) theta[i]=temp;
	}
	// adjust threshold for each GC to produce required (steady state) GC coding-level
	else if(FLAGS.threshold==1){
		c=double(PP)*(1.0-nfrac);
		for (i=0; i<N; i++) {
			for (k = 0; k < PP; k++) {
				temphGC2[k]=hGCpatterns.at(i,k);
			}
			std::sort(temphGC2, temphGC2 + PP);
			n=0;
			for (k = 0; k < PP; k++) {
				n++;
				if(n >= c) {
					theta[i]=temphGC2[k];
					k=PP;
					break;
				}
			}
		}
	}
	// set all thresholds to common value
	else if(FLAGS.threshold==0){
		for (i=0; i<N; i++) theta[i]=theta_global;
	}
	// adjust synaptic weights instead of thresholds
	else if (FLAGS.threshold==2){
		c=double(PP)*(1.0-nfrac);
		for (i=0; i<N; i++) {
			for (k = 0; k < PP; k++) {
				temphGC2[k]=hGCpatterns.at(i,k);
			}
			std::sort(temphGC2, temphGC2 + PP);
			n=0;
			for (k = 0; k < PP; k++) {
				n++;
				if(n >= c) {
					theta[i]=temphGC2[k];
					k=PP;
					break;
				}
			}
		}
		for (i=0; i<N; i++) {
//			a=gsl_ran_gaussian_ziggurat (r,1.0);
			a=0.0;
			// 	renormalise synaptic weights
			for (j=0; j<W[i].num; j++) {
				W[i].strgth1[j]=W[i].strgth1[j]*abs(theta_global+a)/theta[i];
				W[i].strgth2[j]=W[i].strgth2[j]*abs(theta_global+a)/theta[i];
			}
			//	renormalise pre-synaptic inputs
			for (k = 0; k < PP; k++)	hGCpatterns.at(i,k)=hGCpatterns.at(i,k)*abs(theta_global+a)/theta[i];
			theta[i]=abs(theta_global+a);
		}
	}

	// calculate GC steady-state inputs
	hGCavrg=hGCvar=0;
	for (i=0; i<N; i++) {
		for (k = 0; k < PP; k++) {
			hGCavrg+=hGCpatterns.at(i,k);
			hGCvar+=hGCpatterns.at(i,k)*hGCpatterns.at(i,k);
		}
	}
	hGCavrg=hGCavrg/(double(N)*double(PP));
	hGCvar=hGCvar/(double(N)*double(PP)) -hGCavrg*hGCavrg;

	// calculate GC steady-state firing rates
	a=b=d=e=f=0;
	j=j3=0;
	for (i=0; i<N; i++) {
		a2=b2=e2=f2=0;
		j2=0;
		for (k = 0; k < PP; k++) {
			GCpatterns.at(i,k)=std::max(hGCpatterns.at(i,k)-theta[i],0.0);
			a+=GCpatterns.at(i,k);
			b+=GCpatterns.at(i,k)*GCpatterns.at(i,k);
			a2+=GCpatterns.at(i,k);
			b2+=GCpatterns.at(i,k)*GCpatterns.at(i,k);
			if(GCpatterns.at(i,k)>0){ 			// without silent GCs
				e2+=GCpatterns.at(i,k);
				j2++;
				j3++;
			}
		}
		avrgGC[i]=a2/double(PP);
		avrgGC_act[i]=e2/double(j2);
	}
	double avrgGC_tot=a/double(N*PP);

	// renormalise GC steady state firing rates
	for (i = 0; i < N; i++) {
		if(FLAGS.gain==0)		GCgain[i]=GCtarget/avrgGC_tot;
		else if(FLAGS.gain==1)	GCgain[i]=GCtarget/avrgGC[i];
		else if(FLAGS.gain==2)	GCgain[i]=GCtarget/avrgGC_act[i];
		else if(FLAGS.gain==3)	GCgain[i]=1.0;
		// if gain is infinite turn GC off
		if(isinf(GCgain[i])) GCgain[i]=0;
	}
	for (k = 0; k < PP; k++) {
		for (i = 0; i < N; i++) {
			GCpatterns.at(i,k)=GCpatterns.at(i,k)*GCgain[i];
		}
	}

	Tavrg=Gavrg=0;
	for(i= 0; i < N; i++) {
		Gavrg+=GCgain[i];
		Tavrg+=theta[i];
	}
	Gavrg=Gavrg/double(N);
	Tavrg=Tavrg/double(N);
	CLavrg=double(j3)/double(N*PP);

	// switch off STP
	if(FLAGS.MF_STP==0){
		for (i=0; i<N; i++) {
			for (j=0; j<W[i].num; j++) {
				W[i].x1[j]=W[i].x2[j]=1.0;
			}
		}
	}

	delete [] temphGC;
	delete [] temphGC2;

	return;
}
