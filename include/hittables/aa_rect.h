#ifndef AARECT_H
#define AARECT_H

#include "rt_weekend.h"

#include "hittable.h"

template<int axis>
class aa_rect : public hittable {
    static_assert(axis >= 0 && axis < 3, "axis must be one of X=0, Y=1, or Z=2.");
public:
    aa_rect() {}

    aa_rect(
        double _a0, double _a1, double _b0, double _b1, double _k, material* mat
    ) : a0(_a0), a1(_a1), b0(_b0), b1(_b1), k(_k), mat(mat) {}

    virtual bool hit(
        const ray &r, double t_min, double t_max, hit_record &rec
    ) const override;
   
    virtual bool bounding_box(
        double time_0, double time_1, aabb &output_box
    ) const override {
        // The bounding box must have non-zero width in each dimension, so
        // pad the Z dimension a small amount.
        if constexpr (axis == X) { // YZ
            output_box = aabb(
                point3(k - epsilon, a0, b0), point3(k + epsilon, a1, b1)
            );
        } else if constexpr (axis == Y) { // XZ
            output_box = aabb(
                point3(a0, k - epsilon, b0), point3(a1, k + epsilon, b1)
            );
        } else { // XY
            output_box = aabb(
                point3(a0, b0, k - epsilon), point3(a1, b1, k + epsilon)
            );
        }
        return true;
    }

    virtual double pdf_value(const point3 &origin, const vec3 &v) const override;
    virtual vec3 random(const point3 &origin) const override;  

private:
    double a0, a1, b0, b1, k;
    material* mat;
};

template <int axis>
bool aa_rect<axis>::hit(const ray &r, double t_min, double t_max, hit_record &rec) const {
    auto t = (k - r.origin()[axis]) * r.direction_inverse()[axis];
    if (t < t_min || t > t_max)
        return false; // not on the desired segment of the ray

    if (r.direction()[axis] == 0.0) // parallel to surface, treat as no hit
        return false;
    
    auto pt = r.at(t);

    // axis = 0 --> a = 1 b = 2
    // axis = 1 --> a = 2 b = 0
    // axis = 2 --> a = 0 b = 1
    auto e0 = (axis + 1) % 3;
    auto e1 = (axis + 2) % 3;

    if constexpr(axis == 1)
        std::swap(e0, e1);

    auto a = pt[e0];
    auto b = pt[e1];

    if (a < a0 || a > a1 || b < b0 || b > b1)
        return false; // intersection outside bounds of the rectangle
    
    rec.u = (a - a0) / (a1 - a0);
    rec.v = (b - b0) / (b1 - b0);
    rec.t = t;

    auto outward_normal = vec3(1.0 * (axis == X), 1.0 * (axis == Y), 1.0 * (axis == Z));
    rec.set_face_normal(r, outward_normal);
    rec.mat = mat;
    rec.p = pt;
    return true;
}

template<int axis>
double aa_rect<axis>::pdf_value(const point3 &origin, const vec3 &v) const {
    hit_record rec;
    if (!this->hit(ray(origin, v), epsilon, infinity, rec))
        return 0;

    auto area = (a1 - a0) * (b1 - b0);
    auto distance_squared = rec.t * rec.t * v.mag2();
    auto cosine = fabs(dot(v, rec.normal) / v.mag());

    return distance_squared / (cosine * area);
}

template<int axis>
vec3 aa_rect<axis>::random(const point3 &origin) const {
    point3 random_point;
    if constexpr (axis == X) { // YZ
        random_point = point3(k, random_number(a0, a1), random_number(b0, b1));
    } else if constexpr (axis == Y) { // XZ
        random_point = point3(random_number(a0, a1), k, random_number(b0, b1));
    } else { // XY
        random_point = point3(random_number(a0, a1), random_number(b0, b1), k);
    }

    return random_point - origin;
}

using yz_rect = aa_rect<X>;
using xz_rect = aa_rect<Y>;
using xy_rect = aa_rect<Z>;

#endif