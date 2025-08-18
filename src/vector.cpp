#include "vector.h"

// Return random vector with p(direction) ~ cos(theta) 
vec3 random_cosine_direction() {
    auto r2 = random_number();
    auto phi = random_number(0.0, two_pi);
    auto z = sqrt(1 - r2);

    auto x = cos(phi) * sqrt(r2);
    auto y = sin(phi) * sqrt(r2);

    return vec3(x, y, z);
}

// Map point on unit sphere to (u,v) coordinates
std::pair<double, double> sphere_uv(const point3& p) {
    // p: a given point on the unit sphere centered at the origin.
    // u: returned value [0,1] of angle around the Y axis from X=-1.
    // v: returned value [0,1] of angle from Y=-1 to Y=+1.
    //     <1 0 0> yields <0.50 0.50>       <-1  0  0> yields <0.00 0.50>
    //     <0 1 0> yields <0.50 1.00>       < 0 -1  0> yields <0.50 0.00>
    //     <0 0 1> yields <0.25 0.50>       < 0  0 -1> yields <0.75 0.50>

    auto theta = my_acos(-p.y());
    auto phi = my_atan2(-p.z(), p.x()) + pi;
    return std::make_pair(phi * inv_two_pi, theta * inv_pi);
}

float luminance(const colour& c) {
    const colour Lbasis(0.2126f, 0.7152f, 0.0722f);
    return c.dot(Lbasis);
}