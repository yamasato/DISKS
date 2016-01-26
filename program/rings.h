#pragma once
//Ring class.
class Ring{
public:
	double r;
	int    N;
	int    p;

	static Ring iring_parameter(const Ring iring, const double psign, const double rend, const double mass, const double tolerance, const DensityProfile * const dp){

		const int    max_iteration = 100;
		//paramter of a ring i-1.
		const int    pre_p       = iring.p;
		const int    pre_N       = iring.N;
		const double pre_r       = iring.r;
		const double pre_intdens = dp->intdens(pre_r);
		//initial setting of a iteration..
		double r         = pre_r;
		double rmin      = pre_r;
		double rmax      = rend;
		double error     = 1.0;
		double pre_ite_r = pre_r;
		int    pre_ite_p = pre_p;
		int    N(0), p(0);
		int    number_of_ite_r = 0;

		//start iteration
		while(fabs(error) > tolerance){
			if(error > 0.0){
				rmin = r;
				r    = 0.5 * (rmax + r);
			}else{
				rmax = r;
				r    = 0.5 * (rmin + r);
			}
			
			const double x = M_PI * (dp->powerindex(r) + 2.0);
			p = (fabs(x) > 0.5) ? static_cast<int>(psign * round(x)) : 0;
			//if the number is minums, zero is returned.
			const int nextn = static_cast<int>(static_cast<double>(pre_N) + static_cast<double>(pre_p + p) * 0.5);
			N = (nextn < 0) ? 0 : nextn;

			error = 1.0 - 2.0 * fabs(dp->intdens(r) - pre_intdens) / (static_cast<double>(N + pre_N) * mass);
			if((psign > 0 && rmax < tolerance) || (psign < 0 && rmin < tolerance) || psign * (1.0 - rmin / rmax) < tolerance){
				//solution of r is the boundary of p
				if(pre_ite_p != p){
					const double term          = M_PI / (psign * (static_cast<double>(p) - 0.5) - 2.0 * M_PI);
					const double pre_ite_error = 1.0 - dp->powerindex(pre_ite_r) * term;
					double error = 1.0 - dp->powerindex(r) * term;
					int    number_of_ite_boundary_p = 0;
					if(pre_ite_p > p) p = pre_ite_p;
					while(fabs(error) > tolerance){
						r     = r - error * (r - pre_ite_r) / (error - pre_ite_error);
						error = 1.0 - dp->powerindex(r) * term;
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
		return {.r = pre_ite_r, .N = N, .p = pre_ite_p};
	}
};


