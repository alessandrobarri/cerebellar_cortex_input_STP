#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

 struct threshNpara {
	double mean,std;
};
double normal_pdf(double x, double m, double s) {
    static const double inv_sqrt_2pi = 0.3989422804014327;
    double a = (x - m)/s;

    return inv_sqrt_2pi/s *std::exp(-0.5*a*a);
}

int thrN (const gsl_vector * x, void *para, gsl_vector * f){
//		conditions for finding mu and sig of thresholded gaussian
//		given mean and std

		struct threshNpara *p = (struct threshNpara *) para;
		double y0,y1;
		double ertmp,phtmp,mean,std,tmp;

		const double mu = gsl_vector_get (x, 0);
		const double sig = gsl_vector_get (x, 1);

		mean=p->mean;
		std=p->std;

		tmp=mu/sig;
		ertmp=erfc(-tmp/sqrt(2));
		phtmp=sig*normal_pdf(tmp,0,1);

		y0=0.5*mu*ertmp+phtmp -mean;
		y1=0.5*(mu*mu+sig*sig)*ertmp + mu*phtmp -mean*mean -std*std;

		gsl_vector_set (f, 0, y0);
		gsl_vector_set (f, 1, y1);

		return GSL_SUCCESS;
}

int truncN (const gsl_vector * x, void *para, gsl_vector * f){
//		conditions for finding mu and sig of truncated gaussian with a=0 and b=+inf
//		given mean and std

		struct threshNpara *p = (struct threshNpara *) para;
		double y0,y1;
		double Z,X,mean,std,alpha;

		const double mu = gsl_vector_get (x, 0);
		const double sig = gsl_vector_get (x, 1);

		mean=p->mean;
		std=p->std;

		alpha=-mu/sig;
		Z=1-0.5*erfc(-alpha/sqrt(2));
		X=sig*normal_pdf(alpha,0,1)/Z;

		y0=mu +X -mean;
		y1=sig*sig -mu*X -X*X -std*std;

		gsl_vector_set (f, 0, y0);
		gsl_vector_set (f, 1, y1);

		return GSL_SUCCESS;
}


 std::vector<double> find_multiD (struct threshNpara para, gsl_multiroot_function Fmulti, int Maxitr) {
		// multidimensional root finding
		int status,iter;
		const int nD=2;
	    std::vector<double> out_para;
	    double iniguess[nD];
	    iniguess[0]=para.mean;
	    iniguess[1]=para.std;

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
//				cout << gsl_vector_get (solver->x, 0) << "\t" << gsl_vector_get (solver->x, 1) << "\n";
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
