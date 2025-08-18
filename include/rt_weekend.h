#ifndef RTWEEKEND_H
#define RTWEEKEND_H

#include <cmath>
#include <limits>
#include <omp.h>
#include <random>
#include <tuple>

#include "random.h"

// Constants
const double epsilon = 0.000001;
const double infinity = std::numeric_limits<double>::infinity();
const double pi = 3.1415926535897932385;
const double half_pi = pi / 2.0;
const double inv_pi = 1.0 / pi;
const double two_pi = 2.0 * pi;
const double inv_two_pi = 1.0 / two_pi;
const double four_pi = 4.0 * pi;
const double inv_four_pi = 1.0 / four_pi;

enum AXIS {
    X, Y, Z, W
};

// Utility Functions

inline int rtw_get_num_threads() {
#ifdef NDEBUG
	return omp_get_max_threads();
#else
	return 1;
#endif
}

inline int rtw_get_thread_num() {
#ifdef NDEBUG
	return omp_get_thread_num();
#else
	return 0;
#endif
}

template<typename T = double>
inline T clamp(T x, T min, T max) {
    if (x < min) return min;
    if (x > max) return max;
    return x;
}

inline double degrees_to_radians(double degrees) {
    return degrees * pi / 180.0;
}

template<typename T = double>
inline T random_number() {
    static thread_local random_state state(42 * (1 + (1 << rtw_get_thread_num())));
    return static_cast<T>(my_erand48(state.x));
}

template<typename T = double>
inline T random_number(T min, T max) {
    return min + (max - min) * random_number<T>();
}

template<>
inline int random_number(int min, int max) {
    return static_cast<int>(random_number(1.0 * min, 1.0 * (max + 1)));
}

template<>
inline long random_number(long min, long max) {
    return static_cast<long>(random_number(1.0 * min, 1.0 * (max + 1)));
}

template<>
inline size_t random_number(size_t min, size_t max) {
    return static_cast<size_t>(random_number(1.0 * min, 1.0 * (max + 1)));
}

// Optimized Trig Functions
extern double my_acos(double x);
extern double my_atan2(double x, double y);
extern double my_exp(double x);

#endif
