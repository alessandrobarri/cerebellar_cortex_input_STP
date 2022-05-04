void calc_Pdim(arma::mat& MFin, arma::mat& GCin, arma::mat& GCpatterns_onset, arma::mat& GCpatterns_trans, int PP, int N, int M, double &Dim_MF, double &Dim_GCss, double &Dim_GConset, double &Dim_GCtrans) {
	//	*************************
	double lambda_ss[PP],lambda_onset[PP],lambda_trans[PP];
	double avrgGC[N],avrgGC_onset[N],avrgGC_trans[N];
	double MF_Pavrg[M];

	int i,k;
	double a,b,c,d,e,f;

	// calculate average firing rate for each MF unit
	for (i = 0; i < M; i++) {
		a=0;
		for (k = 0; k < PP; k++) {
			a+=MFin.at(i,k);
		}
		MF_Pavrg[i]=a/double(PP);
	}
	// compute average GC steady-state rate and average GC onset rate
	for (i=0; i<N; i++) {
		a=b=c=0;
		for (k = 0; k < PP; k++) {
			a+=GCin.at(i,k);
			b+=GCpatterns_onset.at(i,k);
			c+=GCpatterns_trans.at(i,k);
		}
		avrgGC[i]=a/double(PP);
		avrgGC_onset[i]=b/double(PP);
		avrgGC_trans[i]=c/double(PP);
	}
	// fill and centre MF-matrix
	for (i = 0; i < M; i++) {
		for (k = 0; k < PP; k++) {
			MFin.at(i,k)= MFin.at(i,k)-MF_Pavrg[i];
		}
	}
	arma::vec Smf = arma::svd(MFin);

	// compute eigenvalues
	for (k = 0; k < M; k++) {
		lambda_ss[k]=Smf(k)*Smf(k)/(double(M)-1);
	}

	a=b=0;
	for (k = 0; k < M; k++) {
		a+=lambda_ss[k];
		b+=lambda_ss[k]*lambda_ss[k];
	}
	Dim_MF=a*a/b;

	// fill and centre GC-matrix
	for (i = 0; i < N; i++) {
		for (k = 0; k < PP; k++) {
			GCin.at(i,k)= GCin.at(i,k)-avrgGC[i];
			GCpatterns_onset.at(i,k)= GCpatterns_onset.at(i,k)-avrgGC_onset[i];
			GCpatterns_trans.at(i,k)= GCpatterns_trans.at(i,k)-avrgGC_trans[i];
		}
	}
	arma::vec Sgc = arma::svd(GCin);
	arma::vec Sgc_onset = arma::svd(GCpatterns_onset);
	arma::vec Sgc_trans = arma::svd(GCpatterns_trans);

	// compute eigenvalues
	for (k = 0; k < PP; k++) {
		lambda_ss[k]=Sgc(k)*Sgc(k)/(double(PP)-1);
		lambda_onset[k]=Sgc_onset(k)*Sgc_onset(k)/(double(PP)-1);
		lambda_trans[k]=Sgc_trans(k)*Sgc_trans(k)/(double(PP)-1);
	}

	a=b=c=d=e=f=0;
	for (k = 0; k < PP; k++) {
		a+=lambda_ss[k];
		b+=lambda_ss[k]*lambda_ss[k];
		c+=lambda_onset[k];
		d+=lambda_onset[k]*lambda_onset[k];
		e+=lambda_trans[k];
		f+=lambda_trans[k]*lambda_trans[k];
	}
	Dim_GCss=a*a/b;
	Dim_GConset=c*c/d;
	Dim_GCtrans=e*e/f;
	return;
}
