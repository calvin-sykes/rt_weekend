#include "hittables/moving_sphere.h"

bool moving_sphere::hit(const ray& r, double t_min, double t_max, hit_record& rec) const {
    vec3 oc = r.origin() - centre(r.time());
    auto a = r.direction().mag2();
    auto half_b = dot(oc, r.direction());
    auto c = oc.mag2() - rad * rad;

    auto discriminant = half_b * half_b - a * c;
    if (discriminant < 0)
        return false;
    auto sqrtd = sqrt(discriminant);

    // Find the nearest root that lies in the acceptable range.
    auto root = (-half_b - sqrtd) / a;
    if (root < t_min || t_max < root) {
        root = (-half_b + sqrtd) / a;
        if (root < t_min || t_max < root)
            return false;
    }

    rec.t = root;
    rec.p = r.at(rec.t);
    vec3 outward_normal = (rec.p - centre(r.time())) / rad;
    rec.set_face_normal(r, outward_normal);
    rec.mat = mat;

    return true;
}

bool moving_sphere::bounding_box(double time_0, double time_1, aabb& output_box) const {
    aabb box0(centre(time_0) - vec3(rad), centre(time_0) + vec3(rad));
    aabb box1(centre(time_1) - vec3(rad), centre(time_1) + vec3(rad));
    output_box = box0 | box1;
    return true;
}