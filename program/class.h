#pragma once
//user defined 
class Input{
public:
	double inneredge; 
	double outeredge; 
	double cutedge;   
	int    icut;      
	int    npzero;
	double tolerance;
	int    denscalu;
	int    denstype;
	std::vector<double> denscoef;
	std::string filename;

	Input(const int argc, char* argv[]){
	
		this->denstype    = 1;
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
			{ "denscoef"    , required_argument, NULL, 'P' },
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
		while ((opt = getopt_long(argc, argv, "d:P:I:O:C:i:n:o:t:ph", longopts, &longindex)) != -1) {
			switch (opt) {
			case 'd':
				this->denstype    = atoi(optarg);
				break; 
			case 'P':
				for(char *p, *q=strtok(p=strdup(optarg), ","); q; free((q=strtok(NULL, ",")) ? NULL : p)){
					this->denscoef.push_back(atof(q));
 				}
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
				fprintf(stderr, "-d[N]         --denstype    : create the density profile specified by [N]. \n");
				fprintf(stderr, "-P[N1][N2]... --denscoef    : change the parameters of the density profile to [N1],[N2].... \n");
				fprintf(stderr, "-I[N]         --inneredge   : change the position of inner edge to [N]. \n");
				fprintf(stderr, "-O[N]         --outeredge   : change the position of outer edge to [N]. \n");
				fprintf(stderr, "-C[N]         --cutedge     : change the position of cut   edge to [N]. \n");
				fprintf(stderr, "-i[N]         --iofcut      : change the number of a ring at a cut edge to [N]. \n");
				fprintf(stderr, "-n[N]         --nofpzero    : change the number of particles in a ring to [N]. \n");
				fprintf(stderr, "                            : (for only the power is -2 at cutedge)\n");
				fprintf(stderr, "-o[fname]     --output      : write result to FILE instead of standard output [fname].\n");
				fprintf(stderr, "-t[N]         --tolerance   : change the value of tolerance to [N] \n");
				fprintf(stderr, "-p            --outputdens  : calculate the density of particles \n");
				exit(1);
			default:
                		fprintf(stderr, "Error: Invalid option \'%c\'.\n", opt);
				break;
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
};

//particle class
class Particle{
public:
	double x;
	double y;
	double volu;
	double dens;
};

