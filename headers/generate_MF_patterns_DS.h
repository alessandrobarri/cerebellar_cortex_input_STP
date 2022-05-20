void generate_MF_patterns_DS(arma::mat& MFpatterns, int MFgroup[], double MF_avrg_CL[], double MF_std_CL[], struct cluster_struct CL[], arma::vec& NUlims, arma::vec& grouplims, std::vector <int>& drivers, std::vector <int>& supporters, double MF_avrg_D, double MF_std_D, double MF_avrg_S, double MF_std_S, double pact[], struct netpara NETPARA, struct flags FLAGS, gsl_rng * r) {

	int i,j,k,inidx,finidx;
	double a;
	std::vector <double> para_out;
	std::vector <double> dummyD;
	std::vector <size_t> idx_sorted;
	int M=NETPARA.M;
	int PP=NETPARA.PP;
	int Ncl=NETPARA.Ncl;

	// calculate mu and sig parameters for thresholded Gaussian
	if(FLAGS.pattern==0){
		struct threshNpara thpara;
		gsl_multiroot_function FthresN;
		FthresN.f = &thrN;

		thpara.mean=MF_avrg_D;
		thpara.std=MF_std_D;
		para_out=find_multiD (thpara, FthresN, 200);
		MF_avrg_D=para_out[0];
		MF_std_D=para_out[1];

		thpara.mean=MF_avrg_S;
		thpara.std=MF_std_S;
		para_out=find_multiD (thpara, FthresN, 200);
		MF_avrg_S=para_out[0];
		MF_std_S=para_out[1];

		for (i = 0; i < 5; i++) {
			thpara.mean=MF_avrg_CL[i];
			thpara.std=MF_std_CL[i];
			para_out=find_multiD (thpara, FthresN, 200);
			MF_avrg_CL[i]=para_out[0];
			MF_std_CL[i]=para_out[1];
		}
	}
	// calculate mu and sig parameters for seamless uniform distribution
	else if(FLAGS.pattern==5){
		for (i = 0; i < Ncl; i++) {
			MF_avrg_CL[i]=0.5*(NUlims[Ncl-i]+NUlims[Ncl-i-1]);
			MF_std_CL[i]=(NUlims[Ncl-i]-NUlims[Ncl-i-1])/sqrt(12);
		}
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

	if (FLAGS.groups==5) {
		// 2 groups: driver & supporter MFs
		for (i = 0; i < M; i++) {
			a=gsl_rng_uniform (r);
			if(MFgroup[i]==0) {
				for (k = 0; k < PP; k++) MFpatterns.at(i,k)=draw_patterns(MF_avrg_D, MF_std_D, pact[0], FLAGS.pattern,r);
			}
			else {
				for (k = 0; k < PP; k++) MFpatterns.at(i,k)=draw_patterns(MF_avrg_S, MF_std_S, pact[4], FLAGS.pattern,r);
			}
		}
	}
	else if (FLAGS.groups==55){
		// single uniform distribution cut into Ncl slices
		if(FLAGS.pattern==5){
			for (j = 0; j < Ncl; j++) {
				inidx=round(M*(Ncl-1-j)/Ncl);
				finidx=round(M*(Ncl-j)/Ncl);
				for (i = inidx; i < finidx; i++) {
					MFgroup[i]=j;
					CL[j].cluster.push_back(j);
					for (k = 0; k < PP; k++) MFpatterns.at(i,k)=draw_patterns(MF_avrg_CL[j], MF_std_CL[j], pact[j], FLAGS.pattern,r);
				}
			}
		}
		// 5 peaks with variable identities
		else{
			// uniform mixture
			if(FLAGS.pattern==4){
				for (j = 0; j < 5; j++) {
					inidx=round(M*grouplims(4-j));
					finidx=round(M*grouplims(5-j));
					for (i = inidx; i < finidx; i++) MFgroup[i]=j;
				}
			}
			// other mixtures
			else{
				for (j = 0; j < 5; j++) {
					inidx=round(M*(4-j)/5);
					finidx=round(M*(5-j)/5);
					for (i = inidx; i < finidx; i++) {
						MFgroup[i]=j;
						//~ for (k = 0; k < PP; k++) MFpatterns.at(i,k)=draw_patterns(MF_avrg_CL[j], MF_std_CL[j], pact[j], FLAGS.pattern,r);
					}
				}
			}

			// -----------------------------------
			// distribute MFs into Ncl clusters
			for (i = 0; i < M; i++) 	dummyD.push_back(MFgroup[i]);
			idx_sorted=sort_indices(dummyD);
			for (k = 0; k < Ncl; k++){
				// uniform mixture
				if(FLAGS.pattern==4){
					inidx=round(M*grouplims(k));
					finidx=round(M*grouplims(k+1));
					//~ std::cout << inidx<< "\t"<< finidx<< std::endl;
				}
				// other mixtures
				else{
					inidx=round(k*M/Ncl);
					finidx=round((k+1)*M/Ncl);
				}
				for (i = inidx; i < finidx; i++){
					MFgroup[idx_sorted[i]]=k;
					CL[k].cluster.push_back(idx_sorted[i]);
				}
			}
			dummyD.clear();
			idx_sorted.clear();

			// generate MF rates from true (clustered) Gaussian mixture
			if(Ncl==1){
				for (i = 0; i < M; i++){
					for (k = 0; k < PP; k++){
						a=gsl_rng_uniform (r);
						for (j = 0; j < 5; j++){
							if((a>=j/5.0) & (a<(j+1.0)/5.0))  MFpatterns.at(i,k)=draw_patterns(MF_avrg_CL[j], MF_std_CL[j], pact[j], FLAGS.pattern,r);
						}
					}
				}
			}
			else if(Ncl==2){
				for (i = 0; i < M; i++){
					if(MFgroup[i]==0){
						for (k = 0; k < PP; k++){
							a=gsl_rng_uniform (r);
							if(a<2.0/5.0)  					MFpatterns.at(i,k)=draw_patterns(MF_avrg_CL[0], MF_std_CL[0], pact[0], FLAGS.pattern,r);
							else if((a>=2/5.0) & (a<4/5.0)) MFpatterns.at(i,k)=draw_patterns(MF_avrg_CL[1], MF_std_CL[1], pact[1], FLAGS.pattern,r);
							else  							MFpatterns.at(i,k)=draw_patterns(MF_avrg_CL[2], MF_std_CL[2], pact[2], FLAGS.pattern,r);
						}
					}
					else if(MFgroup[i]==1){
						for (k = 0; k < PP; k++){
							a=gsl_rng_uniform (r);
							if(a<1.0/5.0)  					MFpatterns.at(i,k)=draw_patterns(MF_avrg_CL[2], MF_std_CL[2], pact[2], FLAGS.pattern,r);
							else if((a>=1/5.0) & (a<3/5.0)) MFpatterns.at(i,k)=draw_patterns(MF_avrg_CL[3], MF_std_CL[3], pact[3], FLAGS.pattern,r);
							else  							MFpatterns.at(i,k)=draw_patterns(MF_avrg_CL[4], MF_std_CL[4], pact[4], FLAGS.pattern,r);
						}
					}
				}
			}
			else if(Ncl==3){
				for (i = 0; i < M; i++){
					if(MFgroup[i]==0){
						for (k = 0; k < PP; k++){
							a=gsl_rng_uniform (r);
							if(a<3/5.0)  					MFpatterns.at(i,k)=draw_patterns(MF_avrg_CL[0], MF_std_CL[0], pact[0], FLAGS.pattern,r);
							else							MFpatterns.at(i,k)=draw_patterns(MF_avrg_CL[1], MF_std_CL[1], pact[1], FLAGS.pattern,r);
						}
					}
					else if(MFgroup[i]==1){
						for (k = 0; k < PP; k++){
							a=gsl_rng_uniform (r);
							if(a<1/5.0)  					MFpatterns.at(i,k)=draw_patterns(MF_avrg_CL[1], MF_std_CL[1], pact[1], FLAGS.pattern,r);
							else if((a>=1/5.0) & (a<4/5.0)) MFpatterns.at(i,k)=draw_patterns(MF_avrg_CL[2], MF_std_CL[2], pact[2], FLAGS.pattern,r);
							else  							MFpatterns.at(i,k)=draw_patterns(MF_avrg_CL[3], MF_std_CL[3], pact[3], FLAGS.pattern,r);
						}
					}
					else if(MFgroup[i]==2){
						for (k = 0; k < PP; k++){
							a=gsl_rng_uniform (r);
							if(a<2.0/5.0)  					MFpatterns.at(i,k)=draw_patterns(MF_avrg_CL[3], MF_std_CL[3], pact[3], FLAGS.pattern,r);
							else							MFpatterns.at(i,k)=draw_patterns(MF_avrg_CL[4], MF_std_CL[4], pact[4], FLAGS.pattern,r);
						}
					}
				}
			}
			else if(Ncl==4){
				for (i = 0; i < M; i++){
					if(MFgroup[i]==0){
						for (k = 0; k < PP; k++){
							a=gsl_rng_uniform (r);
							if(a<4.0/5.0)  				MFpatterns.at(i,k)=draw_patterns(MF_avrg_CL[0], MF_std_CL[0], pact[0], FLAGS.pattern,r);
							else						MFpatterns.at(i,k)=draw_patterns(MF_avrg_CL[1], MF_std_CL[1], pact[1], FLAGS.pattern,r);
						}
					}
					else if(MFgroup[i]==1){
						for (k = 0; k < PP; k++){
							a=gsl_rng_uniform (r);
							if(a<3.0/5.0)  				MFpatterns.at(i,k)=draw_patterns(MF_avrg_CL[1], MF_std_CL[1], pact[1], FLAGS.pattern,r);
							else  						MFpatterns.at(i,k)=draw_patterns(MF_avrg_CL[2], MF_std_CL[2], pact[2], FLAGS.pattern,r);
						}
					}
					else if(MFgroup[i]==2){
						for (k = 0; k < PP; k++){
							a=gsl_rng_uniform (r);
							if(a<2.0/5.0)  				MFpatterns.at(i,k)=draw_patterns(MF_avrg_CL[2], MF_std_CL[2], pact[2], FLAGS.pattern,r);
							else						MFpatterns.at(i,k)=draw_patterns(MF_avrg_CL[3], MF_std_CL[3], pact[3], FLAGS.pattern,r);
						}
					}
					else if(MFgroup[i]==3){
						for (k = 0; k < PP; k++){
							a=gsl_rng_uniform (r);
							if(a<1.0/5.0)  				MFpatterns.at(i,k)=draw_patterns(MF_avrg_CL[3], MF_std_CL[3], pact[3], FLAGS.pattern,r);
							else						MFpatterns.at(i,k)=draw_patterns(MF_avrg_CL[4], MF_std_CL[4], pact[4], FLAGS.pattern,r);
						}
					}
				}
			}
			else if(Ncl==5){
				for (i = 0; i < M; i++){
					if(MFgroup[i]==0){
						for (k = 0; k < PP; k++)		MFpatterns.at(i,k)=draw_patterns(MF_avrg_CL[0], MF_std_CL[0], pact[0], FLAGS.pattern,r);
					}
					else if(MFgroup[i]==1){
						for (k = 0; k < PP; k++)		MFpatterns.at(i,k)=draw_patterns(MF_avrg_CL[1], MF_std_CL[1], pact[1], FLAGS.pattern,r);
					}
					else if(MFgroup[i]==2){
						for (k = 0; k < PP; k++)		MFpatterns.at(i,k)=draw_patterns(MF_avrg_CL[2], MF_std_CL[2], pact[2], FLAGS.pattern,r);
					}
					else if(MFgroup[i]==3){
						for (k = 0; k < PP; k++)		MFpatterns.at(i,k)=draw_patterns(MF_avrg_CL[3], MF_std_CL[3], pact[3], FLAGS.pattern,r);
					}
					else if(MFgroup[i]==4){
						for (k = 0; k < PP; k++)		MFpatterns.at(i,k)=draw_patterns(MF_avrg_CL[4], MF_std_CL[4], pact[4], FLAGS.pattern,r);
					}
				}
			}
		}
	}
	else{
		std::cerr <<"Bad groups string, use different routine. \n";
		exit(EXIT_FAILURE);
	}


	return;
}
