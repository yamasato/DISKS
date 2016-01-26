#pragma once
//density profile
class PowerLawDisk;
class LogMassDisk;
class ExponentialDisk;
class PostImpactDisk;

class DensityProfile{
public:
	DensityProfile() {}
	virtual double graddens(double r) const = 0;
	virtual double surfacedens(double r) const = 0;
	virtual double intdens(double r) const = 0;
	double powerindex(const double r) const { return graddens(r) / surfacedens(r) * r ; }
	virtual double roffset(const Input& usr) const { return usr.inneredge; }
	~DensityProfile() {}
};

class PowerLawDisk : public DensityProfile{
	//Ar^(B)
	double A, B; 
public:
	PowerLawDisk(doublevector A_) {
		for(std::size_t i = A_.size(); i < 2; ++i){
			if(i == 0){
        			A_.push_back(1.0);
			}else{
        			A_.push_back(-1.0);
			}
		}
		A = A_[0];
		B = A_[1];
	}
	double graddens(double r) const { return A * B * pow(r, B - 1.0); }
	double surfacedens(double r) const { return A * pow(r, B); }
	double intdens(double r) const { return 2.0 * M_PI * A / (B + 2.0) * pow(r, B + 2.0);}
	double roffset(const Input& usr) const{
		const double rinit_powB2  = pow(usr.inneredge, (B + 2.0));
		const double rout_powB2   = pow(usr.outeredge, (B + 2.0));
		const double sum_rpowB2   = rinit_powB2 + rout_powB2; 
		const double preal        = M_PI * (B + 2.0);
		const double pint         = static_cast<int>(round(preal));
		const double gamma        = 1.0 - (preal / pint);
		const double gamma_plus1  = 1.0 + gamma;

		double roff  = -16.0 * gamma * rinit_powB2 * rout_powB2 / (sum_rpowB2 * sum_rpowB2 * gamma_plus1 * gamma_plus1);
		roff  = 1.0 - sqrt(roff + 1.0);
		roff *= gamma_plus1 * sum_rpowB2 / 4.0 / gamma;
		roff  = pow(roff, 1.0 / (B + 2.0));
		return roff;
	}
	~PowerLawDisk() {}
};


class LogMassDisk : public DensityProfile{
	//Ar^(-2)
	double A;
public:
	LogMassDisk(doublevector A_) {
		for(std::size_t i = A_.size(); i < 1; ++i){
        		A_.push_back(1.0);
		}
		A = A_[0];
	}
	double graddens(double r) const { return -2.0 * A / (r * r * r); }
	double surfacedens(double r) const { return A / (r * r); }
	double intdens(double r) const { return 2.0 * M_PI * A * log(r); }
	~LogMassDisk() {}
};

class ExponentialDisk : public DensityProfile{
	//Aexp(Br)
	double A, B;
public:
	ExponentialDisk(doublevector A_) {
		for(std::size_t i = A_.size(); i < 2; ++i){
			if(i == 0){
        			A_.push_back(1.0);
			}else{
        			A_.push_back(-1.0);
			}
		}
		A = A_[0];
		B = A_[1];
	}
	double graddens(double r) const { return A * B * exp(B * r); }
	double surfacedens(double r) const { return A * exp(B * r); }
	double intdens(double r) const { return 2.0 * M_PI * A / B * (r - 1.0 / B) * exp(B * r); }
	~ExponentialDisk() {}
};

class PostImpactDisk : public DensityProfile{
	//(Ar + B)exp(Cr)
	double A, B, C;
public:
	PostImpactDisk(doublevector A_) {
		for(std::size_t i = A_.size(); i < 3; ++i){
			if(i == 0){
        			A_.push_back(1.0);
			}else{
        			A_.push_back(-1.0);
			}
		}
		A = A_[0]; 
		B = A_[1]; 
		C = A_[2]; 
	}
	double graddens(double r) const { return ((A * r + B) * C + A) * exp(C * r); }
	double surfacedens(double r) const { return (A * r + B) * exp(C * r); }
	double intdens(double r) const { return 2.0 * M_PI * ((2.0 * A * (1.0 / C - r) - B) / C + r * (A * r + B)) * exp(r * C) / C; }
	~PostImpactDisk() {}
};


