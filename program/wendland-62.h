#pragma once
template <typename type> inline type pow8(const type x){
	const double x_mult2 = x       * x;
	const double x_mult4 = x_mult2 * x_mult2;
	return x_mult4 * x_mult4;
}

template <typename type> inline type step(const type x){
	return (x > 1.0) ? 0.0 : (1.0 - x);
}

double WF(const double r, const double h){
	const double H = 2.415230 * h;
	const double q = fabs(r) / H;
	return 78.0 / (7.0 * M_PI * H * H) * pow8(step(q)) * (1.0 + q * (8.0 + q * (25.0 + 32.0 * q)));
}

