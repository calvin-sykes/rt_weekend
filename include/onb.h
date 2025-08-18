#ifndef ONB_H
#define ONB_H

#include "rt_weekend.h"

class orthonormal_basis {
public:
    orthonormal_basis() {}
    orthonormal_basis(const vec3& n) {
        axis[2] = unit_vector(n);
        vec3 a = (fabs(axis[2].x()) > 1 - epsilon) ? vec3(0, 1, 0) : vec3(1, 0, 0);
        axis[1] = unit_vector(cross(axis[2], a));
        axis[0] = cross(axis[2], axis[1]);
    }

    inline vec3 operator[](int i) const { return axis[i]; }

    vec3 u() const { return axis[0]; }
    vec3 v() const { return axis[1]; }
    vec3 w() const { return axis[2]; }

    vec3 local_to_world(double a, double b, double c) const {
        return a * u() + b * v() + c * w();
    }

    vec3 local_to_world(const vec3 &a) const {
        return a.x() * u() + a.y() * v() + a.z() * w();
    }

private:
    vec3 axis[3];
};

#endif
