double draw_patterns(double mu, double sigma, double pact, int patternflag, gsl_rng * r){
	double nu;
	std::pair<double, double> s;  // Output argument of rtnorm
	if (gsl_rng_uniform (r)>pact) {
		nu=0;
	}
	else{
		// truncated Gaussian
		if(patternflag==6){
			//~ nu=-1;
			//~ while(nu<0) nu=mu+gsl_ran_gaussian_ziggurat (r, sigma);
			s = rtnorm(r,0,1e7,mu,sigma);
			nu=s.first;
		}
		// uniform
		else if(patternflag==4 || patternflag==5){
			double a=mu-sqrt(3)*sigma;
			double b=mu+sqrt(3)*sigma;
			if(a<0){
				std::cerr <<"Negative firing rates for parameters: \n";
				std::cerr <<"E(nu) = " << mu << ", SD(nu) = " << sigma << "\n";
				exit(EXIT_FAILURE);
			}
			nu=gsl_ran_flat(r, a, b);
		}
		// gamma
		else if(patternflag==3){
			double k=mu*mu/(sigma*sigma);
			double theta=sigma*sigma/mu;
			nu=gsl_ran_gamma (r, k, theta);
		}
		// exponential
		else if(patternflag==2){
			nu=gsl_ran_exponential (r, mu);
		}
		//log-normal
		else if(patternflag==1){
			double muL=log(mu*mu/sqrt(sigma*sigma+mu*mu));
			double sigmaL=sqrt( log( 1.0+sigma*sigma/(mu*mu) ) );
			nu=exp(muL+gsl_ran_gaussian_ziggurat (r, sigmaL));
		}
		// thresholded Gaussian
		else{
			nu=mu+gsl_ran_gaussian_ziggurat (r, sigma);
			if(nu<0) nu=0;
		}

	}
	return nu;
}
