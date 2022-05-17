#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_randist.h>

 struct golgi_para {
	double nuE_set,nuI_set,fc_set,JIE,hext,hext_var,theta;
};

 struct I_para {
	double alpha,mu,sig,theta;
};
inline double I_gauss_TL (double x, void * para) {
	I_para *p = (I_para *) para;;
	double mu=p->mu;
	double sig=p->sig;
	double theta=p->theta;
	//~ return gsl_ran_gaussian_pdf(x, sig)*std::max(mu+x-theta,0.0);
	return gsl_ran_gaussian_pdf(x, sig)*(mu+x-theta);
}


int find_golgi_params (const gsl_vector * x, void *para, gsl_vector * f){
//		condition for finding gain and JEI
  		struct golgi_para *p = (struct golgi_para *) para;
		double nuE_set,nuI_set,fc_set,hext,hext_var,theta;
		double nuE,hE,fc;
		double y0,y1,tp;
		std::vector<double> para_out;
		gsl_function I;
		struct I_para Ipara;

		const int ITER_intgr=100;
		gsl_integration_cquad_workspace * wrksp = gsl_integration_cquad_workspace_alloc(ITER_intgr);

		double gain = abs(gsl_vector_get (x, 0));
		double JEI = abs(gsl_vector_get (x, 1));

		nuE_set=p->nuE_set;
		nuI_set=nuE_set;
		fc_set=p->fc_set;
		hext=p->hext;
		hext_var=p->hext_var;
		theta=p->theta;

		// compute excitatoty input
		hE=hext - JEI*nuI_set;

		// compute average GC rate
		Ipara.mu=hE;
		Ipara.sig=sqrt(hext_var);
		Ipara.alpha=1.0;
		Ipara.theta=theta;
		I.function = &I_gauss_TL;
		I.params = &Ipara;
		gsl_integration_cquad (&I, theta-hE, 1e6, 0, 1e-8, wrksp, &tp, NULL, NULL);

		nuE=gain*tp;

		// compute coding level
		tp=(theta - hext + JEI*nuI_set)/sqrt(hext_var);
		fc=std::erfc(tp/sqrt(2.0))/2.0;

		y0=nuE_set -nuE;
		y1=fc_set -fc;

		gsl_integration_cquad_workspace_free(wrksp);

		gsl_vector_set (f, 0, y0);
		gsl_vector_set (f, 1, y1);

		return GSL_SUCCESS;
}

std::vector<double> root_multiD_golgi(struct golgi_para para, gsl_multiroot_function Fmulti, int Maxitr, int idum) {
		// multidimensional root finding
		int status,iter;
		const int nD=2;
	    std::vector<double> out_para;

	    const gsl_rng_type * RT;
		gsl_rng * r;

		gsl_rng_env_setup();
		RT = gsl_rng_default;
		r = gsl_rng_alloc (RT);

		gsl_rng_env_setup();
		gsl_rng_set(r,idum);

	    double iniguess[nD];
	    iniguess[0]=5*gsl_rng_uniform (r)+0.1;
	    iniguess[1]=5*gsl_rng_uniform (r)+0.1;

		gsl_vector *z = gsl_vector_alloc (nD);

		const gsl_multiroot_fsolver_type *sol_type;
	    gsl_multiroot_fsolver *solver;
		sol_type = gsl_multiroot_fsolver_hybrids;

		Fmulti.n = nD;
		Fmulti.params = &para;

		for (int i=0;i<nD;i++)	gsl_vector_set (z, i, iniguess[i]);

		solver = gsl_multiroot_fsolver_alloc (sol_type, nD);
		gsl_multiroot_fsolver_set (solver, &Fmulti, z);

		iter=0;
		do {
			iter++;
			status = gsl_multiroot_fsolver_iterate (solver);
			// check if solver is stuck
			if (status) {
				//~ std::cout << gsl_vector_get (solver->x, 0) << "\t" << gsl_vector_get (solver->x, 1) << "\n";
				break;
			}
     		status = gsl_multiroot_test_residual (solver->f, 1e-6);
		}
		while (status == GSL_CONTINUE && iter < Maxitr);

		for (int i=0;i<nD;i++)	out_para.push_back(gsl_vector_get (solver->x, i));

		gsl_multiroot_fsolver_free (solver);
		gsl_vector_free (z);

		return out_para;
}
