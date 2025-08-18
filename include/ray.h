#ifndef RAY_H
#define RAY_H

#include "vector.h"
#include "rt_weekend.h"

class ray {
public:
    ray() : has_inverse(false), tm(0.0) {}
    ray(const point3 &origin, const vec3 &direction, double time = 0.0)
        : orig(origin), dir(direction), has_inverse(false), tm(time) {}
    
    point3 origin() const { return orig; }
    vec3 direction() const { return dir; }
    double time() const { return tm; }
    
    vec3 direction_inverse() const {
        if (!has_inverse)
            compute_inverse();
        return dir_inv;
    }

    vec<bool,3> direction_signs() const {
        if (!has_inverse)
            compute_inverse();
        return dir_signs;
    }

    point3 at(double t) const {
        return orig + t * dir;
    }

private:
    void compute_inverse() const {
        for (int i = 0; i < 3; i++) {
            dir_inv[i] = 1.0 / dir[i];
            dir_signs[i] = dir[i] < 0;
        }
        has_inverse = true;
    }

    point3 orig;
    vec3 dir;
    mutable vec3 dir_inv;
    mutable vec<bool, 3> dir_signs;
    mutable bool has_inverse;
    double tm;
};

#endif