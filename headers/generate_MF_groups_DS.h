void generate_MF_groups_DS(int MFgroup[], struct cluster_struct CL[], arma::vec& grouplims, std::vector <int>& drivers, std::vector <int>& supporters, struct netpara NETPARA, struct flags FLAGS, gsl_rng * r) {

	int i,j,k,inidx,finidx;
	double a;
	std::vector <double> dummyD;
	std::vector <size_t> idx_sorted;
	int M=NETPARA.M;
	int Ncl=NETPARA.Ncl;
	double pdriver=NETPARA.pdriver;

	if (FLAGS.groups==5) {
		// 2 groups: driver & supporter MFs
		for (i = 0; i < M; i++) {
			a=gsl_rng_uniform (r);
			if(a<pdriver) {
				MFgroup[i]=0;
				drivers.push_back(i);
			}
			else {
				MFgroup[i]=1;
				supporters.push_back(i);
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
					//~ std::cout << inidx<< "\t"<< finidx<< std::endl;
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
		}
	}
	else{
		std::cerr <<"Bad groups string, use different routine. \n";
		exit(EXIT_FAILURE);
	}


	return;
}
