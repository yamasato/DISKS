#pragma once
class UserDefinedDisk;

class UserDefinedDisk : public DensityProfile{
public:
	UserDefinedDisk(doublevector A_) {
	}
	//the gradient of the density profile
	double graddens(double r) const { return -1.0 / r / r; }  //<-Please change
	//the density profile
	double surfacedens(double r) const { return 1.0 / r; }    //<-Please change
	//the integral of the density profile
	double intdens(double r) const { return 2.0 * M_PI * r; } //<-Please change
	~UserDefinedDisk() {}
};


/*
class UserDefinedDisk : public DensityProfile{
public:
	UserDefinedDisk(doublevector A_) {
	}
	//the gradient of the density profile
	double graddens(double r) const { return the grad of density; }		   <-Please change
	//the density profile
	double surfacedens(double r) const { return the density; }		   <-Please change
	//the integral of the density profile
	double intdens(double r) const { return the surface integral of density; } <-Please change
	~UserDefinedDisk() {}
};*/

//Examples

//1) density = r^(-1)
/*
 * class UserDefinedDisk : public DensityProfile{
public:
	UserDefinedDisk(doublevector A_) {
	}
	double graddens(double r) const { return -1.0 / r / r; }  //<-Please change
	double surfacedens(double r) const { return 1.0 / r; }	  //<-Please change
	double intdens(double r) const { return 2.0 * M_PI * r; } //<-Please change
	~UserDefinedDisk() {}
};
*/
//2) density = e^(-r)
/*
class UserDefinedDisk : public DensityProfile{
public:
	UserDefinedDisk(doublevector A_) {
	}
	double graddens(double r) const { return -exp(-r); }			     //<-Please change
	double surfacedens(double r) const { return exp(-r); }			     //<-Please change
	double intdens(double r) const { return -2.0 * M_PI * (r + 1.0) * exp(-r); } //<-Please change
	~UserDefinedDisk() {}
};
*/
/*
If you want to use the option -P.

class UserDefinedDisk : public DensityProfile{
	double param[the number of coefficients];
public:
	UserDefinedDisk (doublevector A_) {
		//Set Default values
		const int number_of_param = sizeof(param) / sizeof(param[0]);
		for(std::size_t i = A_.size(); i < number_of_param; ++i){
         		if(i == id of a coefficient)
				A_.push_back(the value of a coeficient);
			else ....
		}

		for(int i = 0; i < number_of_param; ++i){
			param[i] = A_[i];
		}
	}
	
	//the gradient of the density profile
	double graddens(double r) const { return the grad of density; }		   <-Please change
	//the density profile
	double surfacedens(double r) const { return the density; }		   <-Please change
	//the integral of the density profile
	double intdens(double r) const { return the surface integral of density; } <-Please change
	~UserDefinedDisk() {}
};
*/

//Examples
/*
1) density = param[0]r^(-param[1]) Default profile density = 2r^(-1)
class UserDefinedDisk : public DensityProfile{
	double param[2];							<-Please change
public:
	
	UserDefinedDisk(doublevector A_) {
		//Set Default values
		const int number_of_param = sizeof(param) / sizeof(param[0]); 
		for(std::size_t i = A_.size(); i < number_of_param; ++i){
        		if(i == 0)
				A_.push_back(2.0);			       <-Please change
			else if(i == 1)
				A_.push_back(1.0);			       <-Please change

		}
		for(int i = 0; i < number_of_param; ++i){
			param[i] = A_[i];
		}
	}

	double graddens(double r) const { return - param[0] * param[1] * pow(r, -param[1] - 1.0); }			//<-Please change
	double surfacedens(double r) const { return param[0] * pow(r, -param[1]); }				       //<-Please change
	double intdens(double r) const { return 2.0 * M_PI * param[0] / (-param[1] + 2.0) * pow(r, -param[1] + 2.0); } //<-Please change
	~UserDefinedDisk() {}
};
*/
