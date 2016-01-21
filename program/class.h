#include <string>

//ユーザーが入れる情報
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
char ffname[80] = "test.dat" ;
//リングの情報
class Ring{
public:
	double r;
	int    N;
	int    p;
};

//粒子の情報
class Particle{
public:
	double volu;
	double dens;
	double carte[2];
	double polar[2];
};


