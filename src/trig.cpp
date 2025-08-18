#include "rt_weekend.h"

inline double atan_approximation(double x) {
	double a1 =  0.99997726;
	double a3 = -0.33262347;
	double a5 =  0.19354346;
	double a7 = -0.11643287;
	double a9 =  0.05265332;
	double a11= -0.01172120;

	double x_sq = x * x;
	return x * (a1 + x_sq * (a3 + x_sq * (a5 + x_sq * (a7 + x_sq * (a9 + x_sq * a11)))));
}

inline double acos_approximation(double x) {
	double a0 =  1.570795207;
	double a1 = -0.214512362;
	double a2 =  0.087876311;
	double a3 = -0.044958884;
	double a4 =  0.019349939;
	double a5 = -0.004337769;

	return sqrt(1 - x) * (a0 + x * (a1 + x * (a2 + x * (a3 + x * (a4 + x * a5)))));
}

double my_acos(double x) {
	bool flip = x < 0;
	x = fabs(x);
	if (x > 1) return NAN;
	double res = acos_approximation(x);
	return flip ? pi - res : res; 
}

double my_atan2(double y, double x) {
	if (x == 0.0 && y == 0.0)
		return 0.0;
	
	// Ensure input is in [-1, +1]
	bool swap = fabs(x) < fabs(y);
	double atan_input = (swap ? x : y) / (swap ? y : x);

	// Approximate atan
	double res = atan_approximation(atan_input);

	// If swapped, adjust atan output
	res = swap ? (atan_input >= 0.0f ? half_pi : -half_pi) - res : res;

	// Adjust the result depending on the input quadrant
	if (x < 0.0) {
		res = (y >= 0.0 ? pi : -pi) + res;
	}

	// Store result
	return res;
}