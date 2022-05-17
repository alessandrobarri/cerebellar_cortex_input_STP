void generate_ratePv_correlations(struct cluster_struct CL[], std::vector <int> &D_synID, std::vector <int> &S_synID, arma::vec grouplims, struct netpara NETPARA, struct flags FLAGS, gsl_rng * r) {

	int i,j,k,inidx,finidx;
	std::vector <double> dummyD;
	std::vector <size_t> idx_sorted;

	const int NSYN_2=2*NETPARA.NSYN;
	const int Ncl=NETPARA.Ncl;
	const double MFUcorr=NETPARA.MFUcorr;

	arma::mat v(NSYN_2,2);
	arma::mat u(NSYN_2,2);
	arma::mat vtemp(NSYN_2,2);
	arma::umat index(NSYN_2,2);
	arma::umat xident(NSYN_2,2);

	// -----------------------------------
	// 	generate correlation between MF- and synapse-types
	// -----------------------------------
	// create Gaussian copula
	double rho[2][2];
	rho[0][0]=1;
	rho[0][1]=sin(0.5*M_PI*MFUcorr);
	rho[1][0]=sin(0.5*M_PI*MFUcorr);
	rho[1][1]=1;

	//	generate covariance matrix
	gsl_matrix * Copcov_simple = gsl_matrix_alloc (2, 2);
	for (i = 0; i < 2; i++) {
		for (j = 0; j < 2; j++) {
			gsl_matrix_set (Copcov_simple, i, j, rho[i][j]);
		}
	}

	// cholesky decomposition of covariance matrix
 	gsl_linalg_cholesky_decomp (Copcov_simple);

 	// generate 2 correlated patterns from unit Gauss distribution of length NSYN_2
	for (i = 0; i < NSYN_2; i++) {
		for (j = 0; j < 2; j++) {
			vtemp(i,j)=gsl_ran_gaussian_ziggurat (r, 1.0);
		}
	}
	for (i = 0; i < NSYN_2; i++) {
		for (j = 0; j < 2; j++) {
			v(i,j)=0;
			for (k = 0; k <= j; k++) {
				v(i,j)+=vtemp(i,k)*gsl_matrix_get (Copcov_simple, j, k);
			}
		}
	}
	// transform into uniformly distributed random numbers on [0,1]
	for (i = 0; i < NSYN_2; i++) {
		for (j = 0; j < 2; j++) {
			u(i,j)= gsl_cdf_ugaussian_P (v(i,j));
		}
	}
	// sort copula random numbers
	for (j = 0; j < 2; j++) {
		for (i = 0; i < NSYN_2; i++) 	dummyD.push_back(u(i,j));
		idx_sorted=sort_indices(dummyD);
		for (i = 0; i < NSYN_2; i++) 	index(i,j)=idx_sorted[i];
		dummyD.clear();
	}
	idx_sorted.clear();

  	// 	generate identity distributions with desired rank-correlations
	for (i=0; i<NSYN_2/2; i++){
		for (j = 0; j < 2; j++) {
			xident(index(i,j),j)=0;
		}
	}
	for (i=NSYN_2/2; i<NSYN_2; i++){
		for (j = 0; j < 2; j++) {
			xident(index(i,j),j)=1;
		}
	}
	// 	get syn-identity vector for drivers and supporters
	S_synID.clear();	D_synID.clear();
	for (i=0; i<NSYN_2; i++){
		if(xident(i,0)==0)		D_synID.push_back(xident(i,1));
		else 					S_synID.push_back(xident(i,1));
	}

	// 5 groups toymodel
	if(FLAGS.groups==55){
		// 	generate identity distributions with desired rank-correlations
		for (k=0; k<Ncl; k++){
			// uniform mixture
			if(FLAGS.pattern==4){
				inidx=round(NSYN_2*grouplims(k));
				finidx=round(NSYN_2*grouplims(k+1));
				//~ inidx=round(NSYN_2*grouplims(4-k));
				//~ finidx=round(NSYN_2*grouplims(5-k));
			}
			// other mixtures
			else{
				inidx=round(k*NSYN_2/Ncl);
				finidx=round((k+1)*NSYN_2/Ncl);
			}
			for (i=inidx; i<finidx; i++){
				for (j = 0; j < 2; j++) {
					xident( index(i,j) ,j)=k;
				}
			}
		}
		for (k=0; k<Ncl; k++){
			CL[k].synID.clear();
			for (i=0; i<NSYN_2; i++){
				if(xident(i,0)==k)		CL[k].synID.push_back(xident(i,1));
			}
		}
	}

  	gsl_matrix_free (Copcov_simple);

	return;
}
