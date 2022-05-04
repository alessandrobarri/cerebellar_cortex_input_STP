double calc_Tdim(arma::mat& RATE, int SVD_tsteps) {
	//	*************************
	// calculate dimensionality
	// 	do temporal PCA
	int i,j,k;
	int N=RATE.n_rows;
	int Nlamb=std::min(SVD_tsteps,N);
	double lambdaT[Nlamb];
	// compute average rate for each unit and centre activity matrix
	double AVRG[N];
	//~ std::cout << SVD_tsteps <<"\t" << RATE.n_rows <<"\t" << RATE.n_cols << std::endl;
	for (i = 0; i < N; i++) {
		AVRG[i]=0;
		for (j = 0; j < SVD_tsteps; j++) {
			AVRG[i]+=RATE.at(i, j);
		}
		AVRG[i]=AVRG[i]/double(SVD_tsteps);
	}
	for (i = 0; i < N; i++) {
		for (j = 0; j < SVD_tsteps; j++) {
			RATE.at(i, j)=RATE.at(i, j)-AVRG[i];
		}
	}
	// compute singular values
	arma::vec ST = arma::svd(RATE);
	//~ std::cout << Nlamb <<"\t" << ST.n_rows <<"\t" << ST.n_cols << std::endl;
	// compute eigenvalues
	for (k = 0; k < Nlamb; k++) {
		lambdaT[k]=ST(k)*ST(k)/(double(SVD_tsteps)-1);
	}
	//~ std::cout << "muh" << std::endl;
	// compute dimensionality
	double a,b;
	a=b=0;
	for (i = 0; i < Nlamb; i++) {
		a+=lambdaT[i] ;
		b+=lambdaT[i] *lambdaT[i] ;
	}
	// return dimensionality
	return a*a/b;
}
