void generate_weights_DS(struct presyn_2pools W[], struct presyn JEI[], struct cluster_struct CL[], std::vector<int> &drivers, std::vector<int> &supporters, int MFgroup[], struct netpara &NETPARA, struct flags FLAGS, gsl_rng * r) {

	int i,j,k,m;
	int N=NETPARA.N;
	int M=NETPARA.M;
	int MFSYN=NETPARA.MFSYN;
	int Ncl=NETPARA.Ncl;
	int PP=NETPARA.PP;
	double a;
	std::vector <int> dummy;
	// **********************************
	//	generate golgi connectivity

	// set up GoC-to-GC synapses
	for(i= 0; i < N; i++) {
		JEI[i].strgth[0]=NETPARA.avrgJEI;
	}

	// **********************************
	//	generate MF-to-GC connectivity
	int *shffl_arr=new int[M];
	// set up connectivity
	for(i= 0; i < N; i++) {

		// draw exactly MFSYN distinct random connections
		for (j=0;j<M;j++) shffl_arr[j]=j;
		gsl_ran_shuffle (r, shffl_arr, M, sizeof (int));

		// ---------------------------------------
		// check that GCs receive the correct MF-types
		// ---------------------------------------
		if( FLAGS.groups>0 && FLAGS.driver==1){
			a=0;
			if(FLAGS.groups==55){
				// check if cluster0 is missing; if so, add
				for (j=0;j<MFSYN;j++) {
					if(MFgroup[shffl_arr[j]]==0) a+=1;
				}
				if(a==0){
					// assign cluster0 MF to first synapse on GC
					m=gsl_rng_uniform_int (r, CL[0].cluster.size());
					shffl_arr[0]=CL[0].cluster[m];
				}
			}
			else{
				// check if driver is missing; if so, add
				for (j=0;j<MFSYN;j++) {
					if(MFgroup[shffl_arr[j]]==0 || MFgroup[shffl_arr[j]]==1 || MFgroup[shffl_arr[j]]==4) a+=1;
				}
				if(a==0){
					// assign driver MF to first synapse on GC
					m=gsl_rng_uniform_int (r, drivers.size());
					shffl_arr[0]=drivers[m];
				}
			}
		}
		else if( FLAGS.groups>0 && FLAGS.driver==2){

			if(FLAGS.groups==55){
				// draw one from each cluster
				// generate random permutation of [0,1,..,Ncl]
				int *RINT=new int[Ncl];
				for (k=0;k<Ncl;k++) 	RINT[k]=k;
				gsl_ran_shuffle (r, RINT, Ncl, sizeof (int));

				// force the first Ncl synapses of each GC to come from distinct groups
				for (j=0;j<min(4,Ncl);j++) {
					int *shffl_cluster=new int[CL[RINT[j]].cluster.size()];

					for (k=0;k<CL[RINT[j]].cluster.size();k++) 	shffl_cluster[k]=CL[RINT[j]].cluster[k];
					gsl_ran_shuffle (r, shffl_cluster, CL[RINT[j]].cluster.size(), sizeof (int));
					shffl_arr[j]=shffl_cluster[0];

					delete [] shffl_cluster;
				}
				delete [] RINT;
			}
			else{
				// draw random distinc drivers & supporters
				int *shffl_drivers=new int[drivers.size()];
				int *shffl_supporters=new int[supporters.size()];

				for (k=0;k<drivers.size();k++) 	shffl_drivers[k]=drivers[k];
				gsl_ran_shuffle (r, shffl_drivers, drivers.size(), sizeof (int));
				for (j=0;j<2;j++)				shffl_arr[j]=shffl_drivers[j];

				for (k=0;k<supporters.size();k++) 	shffl_supporters[k]=supporters[k];
				gsl_ran_shuffle (r, shffl_supporters, supporters.size(), sizeof (int));
				for (j=2;j<MFSYN;j++)				shffl_arr[j]=shffl_supporters[j];

				delete [] shffl_drivers;
				delete [] shffl_supporters;
			}
		}

		for (j=0;j<MFSYN;j++) dummy.push_back (shffl_arr[j]);
		k=MFSYN;

		W[i].num=k;
		W[i].idx=new int [k];			W[i].Gidx=new int [k];
		W[i].strgth1=new double [k];	W[i].strgth2=new double [k];
		W[i].x1p=new double* [k];		W[i].x2p=new double* [k];
		W[i].x1=new double [k];			W[i].x2=new double [k];
		W[i].x1m=new double [k];		W[i].x2m=new double [k];
		W[i].U1=new double [k];			W[i].tauD1=new double [k];		W[i].pf1=new double [k];
		W[i].U2=new double [k];			W[i].tauD2=new double [k];	 	W[i].pf2=new double [k];

		W[i].alpha1=new double [k];		W[i].WU1=new double [k];	 	W[i].Upf1=new double [k];
		W[i].alpha2=new double [k];		W[i].WU2=new double [k];	 	W[i].Upf2=new double [k];
		W[i].tsyn1=new double* [k];		W[i].tsyn2=new double* [k];

		for(j= 0; j < k; j++) {
			W[i].idx[j]=dummy[j];
			W[i].x1p[j]=new double [PP];	W[i].x2p[j]=new double [PP];
			W[i].tsyn1[j]=new double [PP];	W[i].tsyn2[j]=new double [PP];
		}
		dummy.clear();
	}
	delete [] shffl_arr;


	return;
}
