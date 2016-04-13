#pragma once
class Factory{
	public:
	static const DensityProfile * const Create(const Input usr){
		if (usr.denstype == 1){
			return new PowerLawDisk(usr.denscoef);
		}else if (usr.denstype == 2){
			return new LogMassDisk(usr.denscoef);
		}else if (usr.denstype == 3){
			return new ExponentialDisk(usr.denscoef);
		}else if (usr.denstype == 4){
			return new PostImpactDisk(usr.denscoef);
		}else if (usr.denstype == 5){
			return new UserDefinedDisk(usr.denscoef);
		}else{
			return NULL;
			//return nullptr;
		}
	}
};

template <class F> class DiskCreator{
	const DensityProfile * const dp;
	const Input usr;
	double mass, psign, rend;
	int number_of_particle;
	int number_of_ring;	
	std::vector<Ring> iring;
	std::vector<Particle> par;

public:
	DiskCreator(const Input usr_) : dp(F::Create(usr_)), usr(usr_){
	}
	~DiskCreator(){
	}

	int make_disk(){
		number_of_particle = 0;
		number_of_ring = 0;
		set_cut_point();
		//from inner edge to outer edge
		rend  = usr.outeredge;
		psign = 1.0; //the vector of placement of rings
		put_iring();
		//from outer edge to inner edge
		rend                      = usr.inneredge;
		psign                     = -1.0;
		iring[0].p               *= -1;
		iring[number_of_ring].r  *= -1.0;
		put_iring();
		output_data();
		return 0;
	}

	void set_cut_point(){
		const double power       = dp->powerindex(usr.cutedge);
		const double power_plus2 = power + 2.0;
		const double preal       = M_PI * power_plus2;
		const double pint        = static_cast<double>(round(preal));
		const int    p_cut       = static_cast<int>(pint);
		int N_cut;
		if(fabs(preal) > 0.5){
			//section 2.4 in article
			const double gamma = 1.0 - (preal / pint);
			const double roff  = dp->roffset(usr);
			N_cut = abs(p_cut * usr.icut);
			mass      = dp->surfacedens(usr.cutedge) * usr.cutedge * usr.cutedge;
			if(roff > usr.tolerance){
				mass -= gamma * dp->surfacedens(roff) * roff * roff;
			}
			mass     *= 4.0 * M_PI / (fabs(power_plus2) * static_cast<double>(N_cut) * usr.icut);
		}else if(fabs(preal) < usr.tolerance){
			//section 2.3 in article
			if(usr.npzero == 0){ 
				fprintf(stderr,"Error: The number of particles in a ring is zero. \n Please change the parameter of the option -n.\n");
				exit(1);
			}
			const double N = static_cast<double>(usr.npzero); 
			N_cut = usr.npzero;
			mass = 2.0 * M_PI * log((N + M_PI) / (N - M_PI)) / N;
		}else{
			//section 2.3 in article
			const double icut_minus05    = static_cast<double>(usr.icut - 0.5);
			const double power_plus2_inv = 1.0 / power_plus2; 
			const double Ain     = pow(1.0 + pow(usr.cutedge / usr.inneredge, power_plus2) / icut_minus05, power_plus2_inv);
			const double Aout    = pow(1.0 + pow(usr.cutedge / usr.outeredge, power_plus2) / icut_minus05, power_plus2_inv);
			const double N       = round(M_PI * sqrt((Ain + 1.0) * (Aout + 1.0) / ((Ain - 1.0) * (Aout - 1.0))));	
			N_cut   = static_cast<int>(N);
			mass = 2.0 * M_PI * dp->surfacedens(usr.cutedge) * usr.cutedge * usr.cutedge / (fabs(power_plus2) * N * icut_minus05);
		}
		//The value of mass is abnormality.
		if(mass < usr.tolerance || isnan(mass)){
			fprintf(stderr, "Error: Please change the cutedge or the inner edge by changing the parameter of the option -C or -I. \n");
			exit(1);
		}
		Ring a = {.r = usr.cutedge, .N = N_cut, .p = static_cast<int>(pint)};
		iring.push_back(a);
	}

	void put_iring(){
		//For searching the first ring when we place from outer to inner
		int near_cut = 0;
		while(iring[number_of_ring].r > usr.tolerance){
			//Pre ring is 0-th ring for the first ring when we place from outer to inner.
			if(near_cut == 0 && psign < 0.0 && number_of_ring != 0){
				iring[number_of_ring] = Ring::iring_parameter(iring[0], psign, rend, mass, usr.tolerance, dp);
				near_cut = 1;
			}else{
				iring.push_back(Ring::iring_parameter(iring[number_of_ring], psign, rend, mass, usr.tolerance, dp)); 
				++number_of_ring;
			}

			//Inner edge = 0
			if(fabs(iring[number_of_ring].r) < usr.tolerance && iring[number_of_ring].N == 0){
				iring[number_of_ring].r = 0.0;
				iring[number_of_ring].N = 1;
				iring[number_of_ring].p = 1;
				number_of_ring += 1;
				break;
			}
			//Parameters of iring are abnormality.
			if(iring[number_of_ring].N == 0 && iring[number_of_ring].r > usr.tolerance){
				fprintf(stderr, "Error: Please increase the number of a ring at the cut edge or decrease the value of the cut edge by changing the parameter of the option -i or -C.\n");
				break;
			}
		}
	}
	
	//Making particle data
	void memorize_pparameter(){
		for(int i = 0; i < number_of_ring; ++i){
			//random rotation
			const double j0 = static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
			const double r  = iring[i].r;
			for(int j = 0; j < iring[i].N; ++j){
				const double angle            = 2.0 * M_PI * (static_cast<double>(j) + j0) / static_cast<double>(iring[i].N);
				const double x = iring[i].r * cos(angle);
				const double y = iring[i].r * sin(angle);	
				double vol;
				if(i == 0){
					vol = fabs(iring[i].r - iring[i+1].r);
				}else if(i == number_of_ring - 1 || (r - iring[0].r) * (iring[i-1].r - iring[0].r) > -usr.tolerance){
					vol = fabs(iring[i].r - iring[i-1].r);
				}else{
					vol = fabs(iring[i].r - iring[0].r);
				}
				Particle a = {x, y, vol, 0.0};
				par.push_back(a);
			}
			number_of_particle += iring[i].N;
		}
	}

	//Deriving the density
	void calculate_density_by_SPH(){
		for(int i = 0; i < number_of_particle; i++){
			const double smooth_length = 4.0 * par[i].volu;
			for(int j = 0; j < number_of_particle; j++){
				const double del_xij  = par[j].x - par[i].x;
				const double del_yij  = par[j].y - par[i].y;
				const double del_rij  = sqrt(del_xij * del_xij + del_yij * del_yij);
				par[i].dens += mass * WF(del_rij, smooth_length);
			}
		}
	}

	//Out put data
	void output_data(){
		memorize_pparameter();
		if(usr.denscalu){
			calculate_density_by_SPH();
		}
		//Out put to file
		FILE *fp;
		char fname[80];
		if(0 == strlen(usr.filename.c_str())){
			fp = stdout;
		}else{
			sprintf(fname,"./%s", usr.filename.c_str());
			fp = fopen(fname, "w");
			if(!fp){
				fprintf(stderr, "Error: Fail to open ./%s for writing Please change the filename by the option -f.\n", usr.filename.c_str());
				exit(1);
			}
		}
		for(int i = 0; i < number_of_particle; ++i){
			const double rad = sqrt(par[i].x * par[i].x + par[i].y * par[i].y); 
			fprintf(fp,"%d %.16e %.16e %.16e %.16e\n",i, par[i].x, par[i].y, rad, par[i].dens);
		}
		fclose(fp);
	}
};


