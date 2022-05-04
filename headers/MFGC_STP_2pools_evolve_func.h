void MFGC_STP_2pools_evolve(struct presyn_2pools W[], double MF[], double dt, double taumin, int N, struct flags &FLAGS) {

	int i,j;
	double kp,kpp;
	// MF-to-GC STP
	if(FLAGS.MF_STP>0){
		if(FLAGS.trans_dum==1){
			for(i= 0; i < N; i++) {
				for (j=0; j<W[i].num; j++){
					kp=W[i].strgth1[j]/(W[i].strgth1[j]+W[i].strgth2[j]);
					kpp=W[i].strgth2[j]/(W[i].strgth1[j]+W[i].strgth2[j]);

					if(W[i].tauF1[j]<taumin)	W[i].u1[j]=W[i].U1[j];
					else 						W[i].u1[j]=W[i].u1m[j] +dt *((W[i].U1[j]-W[i].u1m[j])/W[i].tauF1[j] +W[i].U1[j]*(1.0-W[i].u1m[j])*MF[W[i].idx[j]]);
					if(W[i].tauF2[j]<taumin)	W[i].u2[j]=W[i].U2[j];
					else 						W[i].u2[j]=W[i].u2m[j] +dt *((W[i].U2[j]-W[i].u2m[j])/W[i].tauF2[j] +W[i].U2[j]*(1.0-W[i].u2m[j])*MF[W[i].idx[j]]);
					W[i].x1[j]=W[i].x1m[j] +dt *((1.0-W[i].x1m[j])/W[i].tauD1[j] -W[i].x1m[j]*W[i].u1m[j]*(1-W[i].pf1[j])*MF[W[i].idx[j]]);
					W[i].x2[j]=W[i].x2m[j] +dt *((1.0-W[i].x2m[j])/W[i].tauD2[j] -W[i].x2m[j]*W[i].u2m[j]*(1-W[i].pf2[j])*MF[W[i].idx[j]]);
					if(W[i].tauG[j]<taumin)		W[i].q[j]=1;
					else 						W[i].q[j]=W[i].qm[j] +dt *((1.0-W[i].qm[j])/W[i].tauG[j] -W[i].deltaG[j]*( kp*W[i].x1m[j]*W[i].u1m[j] + kpp*W[i].x2m[j]*W[i].u2m[j] )*W[i].qm[j]*MF[W[i].idx[j]]);
				}
			}
		}
		else {	// steady state
			for (i=0; i<N; i++) {
				for (j=0; j<W[i].num; j++) {
					kp=W[i].strgth1[j]/(W[i].strgth1[j]+W[i].strgth2[j]);
					kpp=W[i].strgth2[j]/(W[i].strgth1[j]+W[i].strgth2[j]);

					if(W[i].tauF1[j]<taumin)	W[i].u1[j]=W[i].U1[j];
					else 						W[i].u1[j]=W[i].U1[j]*(1.0+W[i].tauF1[j]*MF[W[i].idx[j]])/(1.0+W[i].U1[j]*W[i].tauF1[j]*MF[W[i].idx[j]]);
					if(W[i].tauF2[j]<taumin)	W[i].u2[j]=W[i].U2[j];
					else 						W[i].u2[j]=W[i].U2[j]*(1.0+W[i].tauF2[j]*MF[W[i].idx[j]])/(1.0+W[i].U2[j]*W[i].tauF2[j]*MF[W[i].idx[j]]);
					W[i].x1[j]=1.0/(1.0+W[i].u1[j]*(1-W[i].pf1[j])*W[i].tauD1[j]*MF[W[i].idx[j]]);
					W[i].x2[j]=1.0/(1.0+W[i].u2[j]*(1-W[i].pf2[j])*W[i].tauD2[j]*MF[W[i].idx[j]]);
					if(W[i].tauG[j]<taumin)		W[i].q[j]=1;
					else 						W[i].q[j]=1.0/(1.0+W[i].deltaG[j]*W[i].tauG[j]*( kp*W[i].x1[j]*W[i].u1[j] + kpp*W[i].x2[j]*W[i].u2[j] )*MF[W[i].idx[j]]);
				}
			}
		}
	}

	return;
}
