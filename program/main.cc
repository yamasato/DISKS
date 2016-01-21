//Copyright 2016- Satoko Yamamoto, Natsuki Hosono, Yoko Funato, Junichiro Makino
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <math.h>
#include <iostream>
#include <getopt.h>
#include "class.h"
#include "sampledisk.h"
#include "userdefineddisk.h"
#include "wendland-62.h"
const static int max_number_of_ring = 1000;

//p og p-gon
int pgon(const double r, const double psign, const DensityProfile& dp){
	const double x = M_PI * (dp.powerindex(r) + 2.0);
	return (fabs(x) > 0.5) ? static_cast<int>(psign * round(x)) : 0;
}

//the number of particles at iring
int pnumber_iring(const int N0, const int p0, const int p){
	const int nextn = N0 + 0.5 * (p0 + p);
	//if the number is minums, zero is returned.
	return (nextn < 0) ? 0 : nextn;
}

//calculation parameters of ring i
Ring iring_parameter(Ring *iring, const double mass, const double psign, const double rend, const int ri, const DensityProfile& dp, const double tolerance){

	const int    max_iteration = 100;
	//paramter of a ring i-1.
	const int    pre_p       = iring[ri].p;
	const int    pre_N       = iring[ri].N;
	const double pre_r       = iring[ri].r;
	const double pre_intdens = dp.intdens(pre_r);
	//initial setting of a iteration..
	double r         = pre_r;
	double rmin      = pre_r;
	double rmax      = rend;
	double error     = 1.0;
	double pre_ite_r;
	int    pre_ite_p = pre_p;
	int    N, p;
	int    number_of_ite_r = 0;
	//start iteration
	while(fabs(error) > tolerance){
		if(error > 0.0){
			rmin = r;
			r    = (rmax + r) / 2.0;
		}else{
			rmax = r;
			r    = (rmin + r) / 2.0;
		}
		p = pgon(r, psign, dp);
		N = pnumber_iring(pre_N, pre_p, p);
		error = 1.0 - 2.0 * fabs(dp.intdens(r) - pre_intdens) / (static_cast<double>(N + pre_N) * mass);
		if((psign > 0 && rmax < tolerance) || (psign < 0 && rmin < tolerance) || psign * (1.0 - rmin / rmax) < tolerance){
			//solution of r is the boundary of p
			if(pre_ite_p != p){
				const double term          = M_PI / (psign * (static_cast<double>(p) - 0.5) - 2.0 * M_PI);
				const double pre_ite_error = 1.0 - dp.powerindex(pre_ite_r) * term;
				double error = 1.0 - dp.powerindex(r) * term;
				int    number_of_ite_boundary_p = 0;
				if(pre_ite_p > p) p = pre_ite_p;
				while(fabs(error) > tolerance){
					r     = r - error * (r - pre_ite_r) / (error - pre_ite_error);
					error = 1.0 - dp.powerindex(r) * term;
					++number_of_ite_boundary_p;
					if(number_of_ite_boundary_p > max_iteration){
						fprintf(stderr, "Error: The number of iterations for deriving parameters of the ring is over %d. \n", max_iteration);
						exit(1);
					}
				}
				break;
			}else if(r < tolerance && pre_N + p == 0){
				break;
			}else if((fabs(r - pre_r) < tolerance || fabs(r - rend) < tolerance) && fabs(error) > tolerance){
				pre_ite_r = -1.0;
				break;
			}
		}
		//input pre-parameter of iteration
		pre_ite_p = p;
		pre_ite_r = r;
		++number_of_ite_r;
		if(number_of_ite_r > max_iteration){
			fprintf(stderr, "Error: The number of iterations for deriving parameters of the ring is over %d. \n", max_iteration);
			exit(1);
		}
	}
	//input solutions
	Ring newring = {pre_ite_r, N, pre_ite_p};
	return newring;
}

void put_iring(Ring *iring, const double *mass, const double psign, const double rend, int *number_of_ring, const DensityProfile& dp, const double tolerance){
	
	//For searching the first ring when we place from outer to inner
	int near_cut = 0;
	while(iring[*number_of_ring].r > tolerance){
		//insufficient memory
		if(max_number_of_ring < *number_of_ring){
			fprintf(stderr, "Please increase the parameter of max_number_of_ring (l.12 in main.cc). \n");
			exit(1);
		}
		//Pre ring is 0-th ring for the first ring when we place from outer to inner.
		if(near_cut != 0 || psign > 0.0){
			iring[*number_of_ring + 1] = iring_parameter(iring, *mass, psign, rend, *number_of_ring, dp, tolerance);
			++*number_of_ring;
		}else{
			//If we have place some particles, *number_of_ring-th ring is out of outer edge.
			if(*number_of_ring != 0){
				iring[*number_of_ring] = iring_parameter(iring, *mass, psign, rend, 0, dp, tolerance);
			}else{ 
				iring[*number_of_ring + 1] = iring_parameter(iring, *mass, psign, rend, 0, dp, tolerance);
				++*number_of_ring;
			}
			near_cut = 1;
		}

		//Inner edge = 0
		if(fabs(iring[*number_of_ring].r) < tolerance && iring[*number_of_ring].N == 0){
			iring[*number_of_ring].r = 0.0;
			iring[*number_of_ring].N = iring[*number_of_ring].p  = 1;
			*number_of_ring += 1;
			break;
		}
		//Parameters of iring are abnormality.
		if(iring[*number_of_ring].N == 0 && iring[*number_of_ring].r > tolerance){
			fprintf(stderr, "Error: Please increase the number of a ring at the cut edge or decrease the value of the cut edge by changing the parameter of the option -i or -C.\n");
			break;
		}
	}
}

//Initial setting
void set_cut_point(Ring *iring, const Input *usr, double *mass, const DensityProfile& dp){
	const double power       = dp.powerindex(usr->cutedge);
	const double power_plus2 = power + 2.0;
	const double preal       = M_PI * power_plus2;
	const double pint        = static_cast<double>(round(preal));

	//Input parameter at cutedge
	iring[0].r = usr->cutedge;
	iring[0].p = static_cast<int>(pint);
	
	if(fabs(preal) > 0.5){
	//section 2.4 in article
		const double gamma = 1.0 - (preal / pint);
		const double roff  = dp.roffset(*usr);
		iring[0].N = abs(iring[0].p * usr->icut);
		*mass      = dp.surfacedens(usr->cutedge) * usr->cutedge * usr->cutedge;
		if(roff > usr->tolerance){
			*mass -= gamma * dp.surfacedens(roff) * roff * roff;
		}
		*mass     *= 4.0 * M_PI / (fabs(power_plus2) * static_cast<double>(iring[0].N) * usr->icut);
	}else if(fabs(preal) < usr->tolerance){
	//section 2.3 in article
		if(usr->npzero == 0){ 
			fprintf(stderr,"Error: The number of particles in a ring is zero. \n Please change the parameter of the option -n.\n");
			exit(1);
		}
		const double N = static_cast<double>(usr->npzero); 
		iring[0].N     = usr->npzero;
		*mass = 2.0 * M_PI * log((N + M_PI) / (N - M_PI)) / N;
	}else{
	//section 2.3 in article
		const double icut_minus05    = static_cast<double>(usr->icut - 0.5);
		const double power_plus2_inv = 1.0 / power_plus2; 
		const double Ain  = pow(1.0 + pow(usr->cutedge / usr->inneredge, power_plus2) / icut_minus05, power_plus2_inv);
		const double Aout = pow(1.0 + pow(usr->cutedge / usr->outeredge, power_plus2) / icut_minus05, power_plus2_inv);
		const double N    = round(M_PI * sqrt((Ain + 1.0) * (Aout + 1.0) / ((Ain - 1.0) * (Aout - 1.0))));
		iring[0].N = static_cast<int>(N);
		*mass = 2.0 * M_PI * dp.surfacedens(usr->cutedge) * usr->cutedge * usr->cutedge 
			/ (fabs(power_plus2) * N * icut_minus05);
	}		
	//The value of mass is abnormality.
	if(*mass < usr->tolerance || isnan(*mass)){
		fprintf(stderr, "Error: Please change the cutedge or the inner edge by changing the parameter of the option -C or -I. \n");
		exit(1);
	}
}

//Making ring data
void make_disk(Ring *iring, const Input *usr, double *mass, int *number_of_ring, DensityProfile& dp){
	set_cut_point(iring, usr, mass, dp);
	//from inner edge to outer edge
	double rend  = usr->outeredge;
	double psign = 1.0; //the vector of placement of rings
	put_iring(iring, mass, psign, rend, number_of_ring, dp, usr->tolerance);
	//from outer edge to inner edge
	rend                      = usr->inneredge;
	psign                     = -1.0;
	iring[0].p               *= -1;
	iring[*number_of_ring].r *= -1.0;
	put_iring(iring, mass, psign, rend, number_of_ring, dp, usr->tolerance);
}

//Making particle data
void memorize_pparameter(Particle *par, Ring *iring, const int *number_of_ring){
	int number_of_particle_ati = 0;//summation the number of particles at a ring i.
	for(int i = 0; i < *number_of_ring; ++i){
		//random rotation
		const double j0 = static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
		const double r  = iring[i].r;

		for(int j = 0; j < iring[i].N; ++j){
			const int    particle_number  = j + number_of_particle_ati ;
			const double angle            = 2.0 * M_PI * (static_cast<double>(j) + j0) / static_cast<double>(iring[i].N);
			par[particle_number].carte[0] = iring[i].r * cos(angle);
			par[particle_number].carte[1] = iring[i].r * sin(angle);
			par[particle_number].polar[0] = iring[i].r;
			par[particle_number].polar[1] = angle;
			
			if(i == 0){
				par[particle_number].volu = fabs(iring[i].r - iring[i+1].r);
			}else if(i == *number_of_ring - 1 || (r - iring[0].r) * (iring[i-1].r - iring[0].r) > -1.0e-15){
				par[particle_number].volu = fabs(iring[i].r - iring[i-1].r);
			}else{
				par[particle_number].volu = fabs(iring[i].r - iring[0].r);
			}
		}
		number_of_particle_ati += iring[i].N;
	}
}

//Deriving the density
void calculate_density_by_SPH(struct Particle *par, const double * const mass, const int number_of_particle){
	for(int i = 0; i < number_of_particle; i++){
		const double smooth_length = 4.0 * par[i].volu;
		par[i].dens = 0.0;
		for(int j = 0; j < number_of_particle; j++){
			const double del_xij  = par[j].carte[0] - par[i].carte[0];
			const double del_yij  = par[j].carte[1] - par[i].carte[1];
			const double del_rij  = sqrt(del_xij * del_xij + del_yij * del_yij);
			par[i].dens += *mass * WF(del_rij, smooth_length);
		}
	}
}

//Out put data
void output_data(Ring *iring, const Input *usr, double const * const mass, const int *number_of_ring){
	int number_of_particle = 0; //the summation of particles
	for(int i = 0; i < *number_of_ring; ++i)
		number_of_particle += iring[i].N;

	Particle *par = new Particle[number_of_particle];
	memorize_pparameter(par, iring, number_of_ring);
	if(usr->denscalu){
		calculate_density_by_SPH(par, mass, number_of_particle);
	}else{
		for(int i = 0; i < number_of_particle; ++i)
			par[i].dens = 0.0;
	}

	//Out put to file
	FILE *fp;
	char fname[80];
	if(0 == strlen(usr->filename.c_str())){
		fp = stdout;
	}else{
		sprintf(fname,"./%s", usr->filename.c_str());
		fp = fopen(fname, "w");
		if(!fp){
			fprintf(stderr, "Error: Fail to open ./%s for writing Please change the filename by the option -f.\n", usr->filename.c_str());
			exit(1);
		}
	}
	for(int i = 0; i < number_of_particle; ++i){
		fprintf(fp,"%d %.16e %.16e %.16e %.16e\n",i, par[i].carte[0], par[i].carte[1], par[i].polar[0], par[i].dens);
	}
	fclose(fp);

	delete [] par;
}

void Input::read_options(const int argc, char* argv[]){

	this->denstype    = 1;
	this->denscoef[0] = 1.0;
	this->denscoef[1] = -1.0;
	this->inneredge   = 1.0;
	this->outeredge   = 4.0;
	this->cutedge     = 1.0;
	this->icut        = 20;
	this->npzero      = 0;
	this->tolerance   = 1.0e-8;
	this->denscalu    = 0;
	this->filename    = "";
	struct option longopts[] = {
		{ "denstype"    , required_argument, NULL, 'd' },
        	{ "coefficient" , required_argument, NULL, 'S' },
        	{ "powernumber" , required_argument, NULL, 'P' },
		{ "inneredge"   , required_argument, NULL, 'I' },
        	{ "outeredge"   , required_argument, NULL, 'O' },
        	{ "cutedge"     , required_argument, NULL, 'C' },
        	{ "iofcut"      , required_argument, NULL, 'i' },
        	{ "nofpzero"    , required_argument, NULL, 'n' },
        	{ "output"      , required_argument, NULL, 'o' },
        	{ "tolerance"   , required_argument, NULL, 't' },
        	{ "outputdens"  , no_argument      , NULL, 'p' },
		{ "help"        , no_argument      , NULL, 'h' },
        	{ 0             , 0                , 0   , 0   },
	};
	
	int opt;
	int longindex;
    	while ((opt = getopt_long(argc, argv, "d:S:P:I:O:C:i:n:o:t:ph", longopts, &longindex)) != -1) {
        	switch (opt) {
            	case 'd':
			this->denstype    = atoi(optarg);
                	break; 
		case 'S':
			this->denscoef[0] = atof(optarg);
                	break;
 		case 'P':
			this->denscoef[1] = atof(optarg);
                	break;
		case 'I':
			this->inneredge   = atof(optarg);
                	break;
            	case 'O':
                	this->outeredge   = atof(optarg);
               		break;
            	case 'C':
                	this->cutedge     = atof(optarg);
			break;
           	case 'i':
                	this->icut        = atoi(optarg);
                	break;
            	case 'n':
                	this->npzero      = atoi(optarg);
			break;
            	case 'o':
                	this->filename    = optarg;
			break;
            	case 't':
                	this->tolerance   = atof(optarg);
			break;
            	case 'p':
		 	this->denscalu  = 1;
			break;
		case 'h':
			fprintf(stderr, "usage: ./a.out [option]\n");
			fprintf(stderr, "-d[N]     --denstype    : create the density profile specified by [N]. \n");
			fprintf(stderr, "-P[N]     --powernumber : change the parameter P to [N]. \n");
			fprintf(stderr, "-S[N]     --coefficient : change the parameter S to [N].\n");
			fprintf(stderr, "-I[N]     --inneredge   : change the position of inner edge to [N]. \n");
			fprintf(stderr, "-O[N]     --outeredge   : change the position of outer edge to [N]. \n");
			fprintf(stderr, "-C[N]     --cutedge     : change the position of cut   edge to [N]. \n");
			fprintf(stderr, "-i[N]     --iofcut      : change the number of ring at a cut edge to [N]. \n");
			fprintf(stderr, "-n[N]     --nofpzero    : change the number of particles in a ring to [N]. \n");
			fprintf(stderr, "                        : (for only the power is -2 at cutedge)\n");
			fprintf(stderr, "-o[fname] --output      : write result to FILE instead of standard output [fname].\n");
			fprintf(stderr, "-t[n]     --tolerance   : change the value of tolerance to [N] \n");
			fprintf(stderr, "-p        --outputdens  : calculate the density of particles \n");
			break;
		default:
                	fprintf(stderr, "Error: Invalid option \'%c\'.\n", opt);
			exit(1);
        	}
    	}
	int param_error = 0;
	if(this->denstype > 5 || this->denstype < 1){
		fprintf(stderr, "Error: Invalid type of the density profile. Please change the parameter of the option -d. \n");
		param_error = 1;
	}else if(this->inneredge < -this->tolerance){
		fprintf(stderr, "Error: Invalid value of the inner edge. Please change the parameter of the option -I. \n");
		param_error = 1;
	}else if(this->outeredge <= this->tolerance){
		fprintf(stderr, "Error: Invalid value of the outer edge. Please change the parameter of the option -O. \n");
		param_error = 1;
	}else if(this->cutedge < this->inneredge || this->cutedge > this->outeredge){
		fprintf(stderr, "Error: Invalid value of the cut edge. Please change the parameter of the option -C. \n");
		param_error = 1;
	}else if(this->icut <= 0){
		fprintf(stderr, "Error: Invalid number of a ring at cut edge. Please change the parameter of the option -i. \n");
		param_error = 1;
	}
	if(param_error) exit(1);
}

int main(int argc, char* argv[]){
	int    number_of_ring = 0;  //the summation of the number of ring
	double mass_of_particle;    //the mass of a particle
	Input usr;
	usr.read_options(argc, argv);
	Ring  *iring = new Ring[max_number_of_ring];
	
	DensityProfile *dp;
	if (usr.denstype == 1){
		dp = new PowerLawDisk(usr.denscoef);
	}else if (usr.denstype == 2){
		dp = new LogMassDisk(usr.denscoef);
	}else if (usr.denstype == 3){
		dp = new ExponentialDisk(usr.denscoef);
	}else if (usr.denstype == 4){
		dp = new PostImpactDisk(usr.denscoef);
	}else if (usr.denstype == 5){
		dp = new UserDefinedDisk();
	}
	make_disk(iring, &usr, &mass_of_particle, &number_of_ring, *dp);
	output_data(iring, &usr, &mass_of_particle, &number_of_ring);

	delete [] iring;
	return 0;
}


