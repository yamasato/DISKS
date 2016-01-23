#include <string>

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
	double denscoef[2];
	std::string filename;

	void read_options(const int argc, char* argv[]);
};
//parameters of a ring
class Ring{
public:
	double r;
	int    N;
	int    p;
};

//parameters of a particle
class Particle{
public:
	double volu;
	double dens;
	double carte[2];
	double polar[2];
};


