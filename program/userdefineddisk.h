
class UserDefinedDisk;

class UserDefinedDisk : public DensityProfile{
public:
	UserDefinedDisk() {
	}
	//the gradient of the density profile
	double graddens(double r) const { return -1.0 / r / r; }  //<-Please change
	//the density profile
	double surfacedens(double r) const { return 1.0 / r; }    //<-Please change
	//the integral of the density profile
	double intdens(double r) const { return 2.0 * M_PI * r; } //<-Please change
	~UserDefinedDisk() {}
};


/*class UserDefinedDisk : public DensityProfile{
	 the variable declaration of coefficients  <- Please change
public:
	LogMassDisk() {
			coefficient = value_of_coefficient; <- Please change
	}
	//the gradient of the density profile
	double graddens(double r) const { return the grad of density; }		   <-Please change
	//the density profile
	double surfacedens(double r) const { return the density; }		   <-Please change
	//the integral of the density profile
	double intdens(double r) const { return the surface integral of density; } <-Please change
	~LogMassDisk() {}
};*/

//Examples
/*
1) density = r^(-1)
class UserDefinedDisk : public DensityProfile{
public:
	LogMassDisk() {
	}
	double graddens(double r) const { return -1.0 / r / r; }		   //<-Please change
	double surfacedens(double r) const { return 1.0 / r; }		   //<-Please change
	double intdens(double r) const { return 2.0 * M_PI * r; } //<-Please change
	~LogMassDisk() {}
};
2) density = e^(-r)
class UserDefinedDisk : public DensityProfile{
public:
	LogMassDisk() {
	}
	double graddens(double r) const { return -exp(-r); }		   //<-Please change
	double surfacedens(double r) const { return exp(-r); }		   //<-Please change
	double intdens(double r) const { return -2.0 * M_PI * (r + 1.0) * exp(-r); } //<-Please change
	~LogMassDisk() {}
};


*/
