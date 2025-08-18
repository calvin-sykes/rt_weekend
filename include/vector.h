#ifndef vec_H
#define vec_H

#include "rt_weekend.h"

#include <array>
#include <cassert>
#include <cmath>
#include <iostream>

#ifndef PLATFORM_MACOS
#define SSE_VECTOR
#endif

using std::enable_if;

template <typename...> struct all_arithmetic;

template <> struct all_arithmetic<> : std::true_type {};

template <typename T, typename ...Rest> struct all_arithmetic<T, Rest...>
    : std::integral_constant<bool,
        std::is_arithmetic<T>::value && all_arithmetic<Rest...>::value>
{};

template <typename T, int N>
class vec {
public:
    using this_t = vec<T, N>;
    using elem_t = std::array<T, N>;
    using enabled = typename std::enable_if<std::is_arithmetic<T>::value>::type;

    // Zero-initialise the vector
    vec() { e.fill(static_cast<T>(0)); }

    // Initialise every component of the vector to a single value
    vec(T _e) { e.fill(_e); }

    // Initialise the vector's components to the given values
    template <
        typename... TT,
        typename = typename enable_if<all_arithmetic<TT...>::value>::type,
        typename = typename enable_if<sizeof...(TT) == N>::type
    >
    vec(TT... _ee) : e{ T(_ee)... } { assert(!has_nans()); }

    // Initialise the vector's components from the given array of values
    vec(elem_t ee) : e { ee } { assert (!has_nans()); }

    // Copy constructor
    vec(const vec &v) : e(v.e) { assert(!has_nans()); }

    // Move constructor
    vec(const vec &&v) noexcept : e(v.e) { assert(!has_nans()); }

    // Copy assignment operator
    vec& operator=(const vec &v) { e = v.e; assert(!has_nans());  return *this; }

    // Move assignment operator
    vec& operator=(const vec &&v) noexcept { e = v.e; assert(!has_nans());  return *this; }

    static vec random();
    static vec random(T min, T max);

    T x() const { return e[X]; }
    T y() const { return e[Y]; }
    T z() const { typename enable_if<(N>=3)>::type(); return e[Z]; }
    T w() const { typename enable_if<(N>=4)>::type(); return e[W]; }

    T operator[] (int i) const { return e[i]; }
    T& operator[] (int i) { return e[i]; }
    
    operator T const* () const { return e.data(); }

    vec operator-() const { 
        vec v;
        for (int i = 0; i < N; i++)
            v[i] = -e[i];
        return v;
    }

    vec& operator+=(const vec &v) {
        for (int i = 0; i < N; i++)
            e[i] += v.e[i];
        return *this;
    }

    vec& operator-=(const vec &v) {
        for (int i = 0; i < N; i++)
            e[i] -= v.e[i];
        return *this;
    }

    vec& operator*=(const T t) {
        for (int i = 0; i < N; i++)
            e[i] *= t;
        return *this;
    }

     vec& operator*=(const vec &v) {
        for (int i = 0; i < N; i++)
            e[i] *= v[i];
        return *this;
    }

    vec& operator/=(const T t) {
      return (*this) *= (static_cast<T>(1.0) / t);
    }

    vec& operator/=(const vec& v) {
        for (int i = 0; i < N; i++)
            e[i] /= v[i];
        return *this;
    }

    T dot(const vec& v) const {
        T n = T(0);
        for (int i = 0; i < N; i++)
            n += e[i] * v.e[i];
        return n;
    }

    T mag2() const {
        return dot(*this);
    }

    T mag() const {
        return sqrt(mag2());
    }

    bool near_zero() const {
        const auto s = 1e-8;
        bool z = true;
        for (int i = 0; i < N; i++)
            z &= (fabs(e[i]) < s);
        return z;
    }

    bool has_nans() const {
        constexpr bool is_floating = std::is_floating_point<T>::value;
        if constexpr (is_floating) {
            for (int i = 0; i < N; i++) {
                if (std::isnan(e[i]))
                    return true;
            }
            return false;
        }
        else {
            return false;
        }
    }

    int rank() const { return N; }
    elem_t const& data() const { return e; }

private:
    std::array<T,N> e;
};

#ifdef SSE_VECTOR
#include <immintrin.h>

template <size_t N>
class vec<double, N> {
public:
    using this_t = vec<double, N>;
    using elem_t = __m256d;
    using enabled = typename std::enable_if<(N<=4)>::type;

    // Zero-initialise the vector
    vec() : e(_mm256_setzero_pd()) {}

    // Initialise every component of the vector to a single value
    vec(double _e) : e(_mm256_set1_pd(_e)) {}

    // Initialise the vector's components to the given values
    template <
        typename... TT,
        typename = typename enable_if<all_arithmetic<TT...>::value>::type,
        typename = typename enable_if<sizeof...(TT) == N>::type
    >
    vec(TT... _ee) : e{ double(_ee)... } { assert(!has_nans()); }

    // Initialise the vector's components from the given array of values
    vec(elem_t ee) : e{ ee } { assert(!has_nans()); }

    // Copy constructor
    vec(const vec& v) : e(v.e) { assert(!has_nans()); }

    // Move constructor
    vec(const vec&& v) noexcept : e(v.e) { assert(!has_nans()); }

    // Copy assignment operator
    vec& operator=(const vec& v) { e = v.e; assert(!has_nans());  return *this; }

    // Move assignment operator
    vec& operator=(const vec&& v) noexcept { e = v.e; assert(!has_nans());  return *this; }

    static vec random();
    static vec random(double min, double max);

    double x() const { return e.m256d_f64[X]; }
    double y() const { return e.m256d_f64[Y]; }
    double z() const { typename enable_if<(N >= 3)>::type(); return e.m256d_f64[Z]; }
    double w() const { typename enable_if<(N >= 4)>::type(); return e.m256d_f64[W]; }

    double operator[] (int i) const { return e.m256d_f64[i]; }
    double& operator[] (int i) { return e.m256d_f64[i]; }

    operator double const* () const { return &e.m256d_f64[0]; }
    operator elem_t const& () const { return e; }

    vec operator-() const {
        vec v;
        v.e = _mm256_mul_pd(e, _mm256_set1_pd(-1.0));
        return v;
    }

    vec& operator+=(const vec& v) {
        e = _mm256_add_pd(e, v.e);
        return *this;
    }

    vec& operator-=(const vec& v) {
        e = _mm256_sub_pd(e, v.e);
        return *this;
    }

    vec& operator*=(const double t) {
        __m256d s = _mm256_set1_pd(t);
        e = _mm256_mul_pd(e, s);
        return *this;
    }

    vec& operator*=(const vec& v) {
        e = _mm256_mul_pd(e, v.e);
        return *this;
    }

    vec& operator/=(const double t) {
        return (*this) *= (1.0 / t);
    }

    vec& operator/=(const vec& v) {
        e = _mm256_div_pd(e, v.e);
        return *this;
    }

    double dot(const vec& v) const {
        if constexpr (N == 2) {
            // prod is (i.i. j.j _._ _._)
            __m256d prod = _mm256_mul_pd(e, v.e);
            prod = _mm256_hadd_pd(prod, prod);

#ifndef NDEBUG
            double dp_check = (*this)[0] * v[0] + (*this)[1] * v[1];
            double dp_simd = prod.m256d_f64[0];
            assert(abs((dp_simd - dp_check) / dp_check) < 1e-10 || abs(dp_simd - dp_check) < 1e-10);
#endif
            return prod.m256d_f64[0];
        }
        else if constexpr (N == 3) {
            // prod is (i.i. j.j k.k _._)
            __m256d prod = _mm256_mul_pd(e, v.e);
            __m256d swap = _mm256_permute2f128_pd(prod, prod, 1);
            prod = _mm256_add_pd(prod, swap);
            prod = _mm256_hadd_pd(prod, prod);

#ifndef NDEBUG
            double dp_check = (*this)[0] * v[0] + (*this)[1] * v[1] + (*this)[2] * v[2];
            double dp_simd = prod.m256d_f64[0];
            assert(abs((dp_simd - dp_check) / dp_check) < 1e-10 || abs(dp_simd - dp_check) < 1e-10);
#endif
            return prod.m256d_f64[0];
        }

       /* if constexpr (N == 2) {
            __m256d dot_prod = _mm256_sum_pd(u1, u2);
            return dot_prod.m256d_f64[0];
        }
        else {
            __m256d u3 = _mm256_shuffle_pd(prod, prod, _MM_SHUFFLE4(2, 2, 2, 2));
            __m256d dot_prod = _mm256_add_pd(u1, _mm256_add_pd(u2, u3));
            return dot_prod.m256d_f64[0];
        }*/
    }

    double mag2() const {
        return dot(*this);
    }

    double mag() const {
        return sqrt(mag2());
    }

    bool near_zero() const {
        const auto s = 1e-8;
        bool z = true;
        for (int i = 0; i < N; i++)
            z &= (fabs(e.m256d_f64[i]) < s);
        return z;
    }

    bool has_nans() const {
        for (int i = 0; i < N; i++) {
            if (std::isnan(e.m256d_f64[i]))
                return true;
        }
        return false;
    }

    int rank() const { return N; }
    elem_t const& data() const { return e; }

private:
    __m256d e;
};
#endif

template<size_t N>
using sse_vec = vec<double, N>;

using vec2 = sse_vec<2>;
using vec3 = sse_vec<3>;
using point3 = sse_vec<3>;
using point2 = sse_vec<2>;
using colour = vec<float,3>;

extern float luminance(const colour& c);

template <typename T, int N>
vec<T,N> vec<T,N>::random() {
    typename vec<T, N>::elem_t c{};
    for (int i = 0; i < N; i++)
        c[i] = random_number<T>();
    return vec(c);
}

template <typename T, int N>
vec<T,N> vec<T,N>::random(T min, T max) {
    typename vec<T, N>::elem_t c{};
    for (int i = 0; i < N; i++)
        c[i] = random_number<T>(min, max);
    return vec(c);
}

template <typename T, int N>
std::ostream& operator<<(std::ostream &out, const vec<T,N> &v) {
    for (int i = 0; i < N; i++)
        out << v[i] << ' ';
    return out; 
}

template <typename T, int N>
std::istream& operator>>(std::istream &in, vec<T,N> &v) {
    for (int i = 0; i < N; i++)
        in >> v[i];
    return in; 
}

template <typename T, int N>
vec<T,N> operator+(vec<T,N> u, const vec<T,N> &v) {
    u += v;
    return u;
}

template <typename T, int N>
vec<T,N> operator-(vec<T,N> u, const vec<T,N> &v) {
    u -= v;
    return u;
}

template <typename T, typename Ts, int N>
vec<T,N> operator*(vec<T,N> v, Ts t) {
    v *= static_cast<T>(t);
    return v;
}

template <typename T, typename Ts, int N>
vec<T,N> operator*(Ts t, vec<T,N> v) {
    v *= static_cast<T>(t);
    return v;
}

template <typename T, int N>
vec<T, N> operator*(vec<T, N> u, const vec<T, N>& v) {
    u *= v;
    return u;
}

template <typename T, typename Ts, int N>
vec<T, N> operator/(vec<T,N> v, Ts t) {
    return static_cast<T>(1.0 / t) * v;
}

template <typename T, int N>
vec<T, N> operator/(vec<T, N> u, const vec<T, N>& v) {
    u /= v;
    return u;
}

template <typename T, int N>
bool operator==(const vec<T, N> &v1, const vec<T,N> &v2) {
    return (v1 - v2).near_zero();
}

template <typename T, int N>
bool operator!=(const vec<T, N> &v1, const vec<T,N> &v2) {
    return !(v1 == v2);
}

template <typename T, int N>
T dot(const vec<T, N>& u, const vec<T, N>& v) {
    return u.dot(v);
}

template <typename T, int N>
vec<T,3> cross(const vec<T,N> &u, const vec<T,N> &v) {
    typename enable_if<(N==3)>::type();
    return vec<T,3>(
        u[1] * v[2] - u[2] * v[1],
        u[2] * v[0] - u[0] * v[2],
        u[0] * v[1] - u[1] * v[0]
    );
}

template <typename VT>
VT unit_vector(const VT &v) {
    return v / v.mag();
}

template <typename T, int N>
vec<T,N> min_vector(const vec<T,N> &v1, const vec<T,N> &v2) {
    vec<T,N> min_vec;
    for (int i = 0; i < N; i++)
        min_vec[i] = std::min(v1[i], v2[i]);
    return min_vec;
}

template <typename T, int N>
vec<T,N> max_vector(const vec<T,N> &v1, const vec<T,N> &v2) {
    vec<T,N> max_vec;
    for (int i = 0; i < N; i++)
        max_vec[i] = std::max(v1[i], v2[i]);
    return max_vec;
}

template <typename T, int N>
T min_axis(const vec<T,N> &v) {
    T min_val = static_cast<T>(infinity);
    for (int i = 0; i < N; i++)
        min_val = std::min(min_val, v[i]);
    return min_val;
}

template <typename T, int N>
T max_axis(const vec<T,N> &v) {
    T max_val = static_cast<T>(-infinity);
    for (int i = 0; i < N; i++)
        max_val = std::max(max_val, v[i]);
    return max_val;
}

// Return random vector with magnitude < 1
template<typename VT = vec3>
VT random_in_unit_sphere() {
    while(true) {
        auto p = VT::random(-1.0, 1.0);
        if (p.mag2() >= 1) continue;
        return p;
    }
}

// Return random vector with magnitude < 1, in same hemisphere as input vector
template <typename VT = vec3>
VT random_in_hemisphere(const VT &normal) {
    VT in_unit_sphere = random_in_unit_sphere<VT>();
    if (dot(in_unit_sphere, normal) > 0.0)
        return in_unit_sphere;
    else
        return -in_unit_sphere;
}

// Return random vector with unit magnitude
template <typename VT = vec3>
VT random_unit_vector() {
    return unit_vector(random_in_unit_sphere<VT>());
}

// Return random vector with p(direction) ~ cos(theta) 
extern vec3 random_cosine_direction();

// Map point on unit sphere to (u,v) coordinates
extern std::pair<double, double> sphere_uv(const point3& p);

// Reflect a vector v relative to the surface normal n 
template <typename VT = vec3>
VT reflect(const VT& v, const VT &n) {
    return v - 2 * dot(v, n) * n;
}

// Refract a vector according to Snell's law
template <typename VT = vec3>
VT refract(const VT& uv, const VT &n, double etai_over_etat) {
    auto cos_theta = std::min(dot(-uv, n), 1.0);
    vec3 r_out_perp = etai_over_etat * (uv + cos_theta * n);
    vec3 r_out_parallel = -sqrt(fabs(1.0 - r_out_perp.mag2())) * n;
    return r_out_perp + r_out_parallel;
}

#ifdef SSE_VECTOR

template<size_t N>
vec<double, N> vec<double, N>::random(double min, double max)
{
    typename vec<double, N>::elem_t c{ 0.0 };
    for (int i = 0; i < N; i++)
        c.m256d_f64[i] = random_number(min, max);
    return vec(c);
}

//template<size_t N>
//vec<double, N> min_vector(const vec<double, N>& v1, const vec<double, N>& v2) {
//    __m256d min = _mm256_min_pd(v1, v2);
//    return vec<double, N>{ min };
//}
//
//template<size_t N>
//vec<double, N> max_vector(const vec<double, N>& v1, const vec<double, N>& v2) {
//    __m256d max = _mm256_max_pd(v1, v2);
//    return vec<double, N>{ max };
//}

#endif

#endif // VECTOR_H
