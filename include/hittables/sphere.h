#ifndef SPHERE_H
#define SPHERE_H

#include "rt_weekend.h"

#include "aabb.h"
#include "hittable.h"

class sphere : public hittable {
public:
    sphere() : rad(0), mat(nullptr) {}
    sphere (point3 centre, double radius, material* material)
        : cen(centre), rad(radius), mat(material) {}

    virtual bool hit(
        const ray& r, double t_min, double t_max, hit_record& rec
    ) const override;
    virtual bool bounding_box(
        double time_0, double time_1, aabb& output_box
    ) const override;

    virtual double pdf_value(const point3 &o, const vec3 &v) const override;
    virtual vec3 random (const point3 &o) const override;

private:
    static vec3 random_to_sphere(double radius, double dist2) {
        auto phi = random_number(0.0, two_pi);
        auto r2 = random_number();

        auto root = std::max(0.0, 1 - radius * radius / dist2);
        auto z = 1 + r2 * (sqrt(root) - 1);
        auto x = cos(phi) * sqrt(1 - z * z);
        auto y = sin(phi) * sqrt(1 - z * z);

        return vec3(x, y, z);
    }

    point3 cen;
    double rad;
    material* mat;
};

#endif

