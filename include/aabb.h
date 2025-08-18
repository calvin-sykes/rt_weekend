#ifndef AABB_H
#define AABB_H

#include "rt_weekend.h"

#include "ray.h"
#include "vector.h"

class aabb {
public:
    aabb() : aabb(infinity, -infinity) {}
    aabb(const point3 &a, const point3 &b)
        : minimum(a), maximum(b) {
        /*for (int i = 0; i < 3; i++) {
            if (minimum[i] > maximum[i])
                std::swap(minimum[i], maximum[i]);
        }*/
    }

    point3 min() const { return minimum; }
    point3 max() const { return maximum; }

    point3& operator[](size_t i) { return i == 0 ? minimum : maximum; }
    const point3& operator[](size_t i) const { return i == 0 ? minimum : maximum; }

    int longest_axis() const {
        auto dim = maximum - minimum;

        auto maximum_axis_length = 0.0;
        auto maximum_axis_idx = -1;

        for (int i = 0; i < 3; i++) {
            auto axis_length = fabs(dim[i]);
            if (axis_length > maximum_axis_length) {
                maximum_axis_idx = i;
                maximum_axis_length = axis_length;
            }
        }

        return maximum_axis_idx;
    }

    double surface_area() const {
        auto dim = maximum - minimum;
        return 2 * (dim[X] * dim[Y] + dim[Y] * dim[Z] + dim[Z] * dim[X]);
    }

    vec3 relative_position(const point3& p) const {
        vec3 o = p - minimum;

        // catch NaN for bounding box with zero thickness
        return o / (maximum - minimum + point3(epsilon));
    }

    aabb& operator|=(aabb const& b) {
        minimum = min_vector(minimum, b.minimum);
        maximum = max_vector(maximum, b.maximum);
        return *this;
    }

    aabb operator|=(point3 const& p) {
        minimum = min_vector(minimum, p);
        maximum = max_vector(maximum, p);
        return *this;
    } 

    bool hit(const ray &r, double t_min, double t_max) const;
    bool hit(const vec3& o, const vec3& inv_dir, const vec<bool, 3>& dir_signs, double tr_min, double tr_max) const;
private:
    point3 minimum, maximum;
};

extern aabb operator|(aabb const& b1, aabb const& b2);
extern aabb operator|(aabb const& b, point3 const& p);

#endif