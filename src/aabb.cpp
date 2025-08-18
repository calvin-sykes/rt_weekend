#include "aabb.h"

bool aabb::hit(const ray &r, double t_min, double t_max) const {
   /* auto invD = r.direction_inverse();
    auto t0 = (minimum - r.origin()) * invD;
    auto t1 = (maximum - r.origin()) * invD;

    for (int a = 0; a < 3; a++) {
        auto t0a = t0[a];
        auto t1a = t1[a];
        if (invD[a] < 0.0)
            std::swap(t0a, t1a);
        t_min = t0a > t_min ? t0a : t_min;
        t_max = t1a < t_max ? t1a : t_max;
        if (t_max <= t_min)
            return false;
    }
    return true;*/
    return hit(r.origin(), r.direction_inverse(), r.direction_signs(), t_min, t_max);
}

bool aabb::hit(const vec3& o, const vec3 &inv_dir, const vec<bool,3>& dir_signs, double tr_min, double tr_max) const {
    const auto& bounds = *this;
    
    double t_min  = (bounds[ dir_signs.x()].x() - o.x()) * inv_dir.x();
    double t_max  = (bounds[!dir_signs.x()].x() - o.x()) * inv_dir.x();
    double ty_min = (bounds[ dir_signs.y()].y() - o.y()) * inv_dir.y();
    double ty_max = (bounds[!dir_signs.y()].y() - o.y()) * inv_dir.y();

    if (t_min > ty_max || ty_min > t_max)
        return false;
    if (ty_min > t_min) t_min = ty_min;
    if (ty_max < t_max) t_max = ty_max;

    double tz_min = (bounds[ dir_signs.z()].z() - o.z()) * inv_dir.z();
    double tz_max = (bounds[!dir_signs.z()].z() - o.z()) * inv_dir.z();

    if (t_min > tz_max || tz_min > t_max)
        return false;
    if (tz_min > t_min) t_min = tz_min;
    if (tz_max < t_max) t_max = tz_max;

    return (t_min < tr_max) && (t_max > tr_min);
}

aabb operator|(aabb const& b1, aabb const& b2) {
    aabb union_box = b1;
    union_box |= b2;
    return union_box;
}

aabb operator|(aabb const& b, point3 const& p) {
    aabb union_box = b;
    union_box |= p;
    return union_box;
}