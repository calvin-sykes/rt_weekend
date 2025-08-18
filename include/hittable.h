#ifndef HITTABLE_H
#define HITTABLE_H

#include "rt_weekend.h"

#include "aabb.h"
#include "ray.h"

class material;

struct hit_record {
    point3 p;
    vec3 normal;
    material *mat;
    double t;
    double u, v;
    bool front_face;

    inline void set_face_normal(const ray &r, const vec3 &outward_normal) {
        front_face = dot(r.direction(), outward_normal) < 0;
        normal = front_face ? outward_normal : -outward_normal;
    }
};

class hittable {
public:
    hittable() {}
    virtual ~hittable() {}

    virtual bool hit(const ray& r, double t_min, double t_max, hit_record &rec) const = 0;
    virtual bool bounding_box(double time_0, double time_1, aabb& output_box) const = 0;

    virtual double pdf_value(const point3 &o, const vec3 &v) const {
        std::cerr << "fell through to default implementation for pdf_value()\n";
        return 1.0;
    }

    virtual vec3 random(const vec3 &o) const {
        std::cerr << "fell through to default implementation for random()\n";
        return vec3(1, 0, 0);
    }

    void set_sampling_target(bool is_target) { sampling_target = is_target; }
    bool is_sampling_target() const { return sampling_target; }

private:
    bool sampling_target = false;
};

#endif