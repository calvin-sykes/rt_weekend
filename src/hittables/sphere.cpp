#include "hittables/sphere.h"

#include "onb.h"

bool sphere::hit(const ray& r, double t_min, double t_max, hit_record& rec) const {
    vec3 oc = r.origin() - cen;
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
    vec3 outward_normal = (rec.p - cen) / rad;
    rec.set_face_normal(r, outward_normal);

    std::tie(rec.u, rec.v) = sphere_uv(outward_normal);
    rec.mat = mat;

    return true;
}

bool sphere::bounding_box(double time_0, double time_1, aabb& output_box) const {
    output_box = aabb(cen - vec3(rad), cen + vec3(rad));
    return true;
}

double sphere::pdf_value(const point3 &o, const vec3 &v) const {
    hit_record rec;
    if (!this->hit(ray(o, v), epsilon, infinity, rec))
        return 0;

    // Fake scattering inside the sphere
    // Probably happens if two hittables are overlapping
    if ((cen - o).mag2() < rad * rad)
         return 0;

    auto cos_theta_max = sqrt(1 - rad * rad / (cen - o).mag2());
    auto solid_angle = two_pi * (1 - cos_theta_max);

    return 1.0 / solid_angle;
}

vec3 sphere::random(const point3 &o) const {
    vec3 direction = cen - o;
    auto distance_squared = direction.mag2();
    orthonormal_basis uvw{ direction };
    return uvw.local_to_world(random_to_sphere(rad, distance_squared));
}
