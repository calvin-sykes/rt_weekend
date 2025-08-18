#ifndef CONSTANT_MEDIUM_H
#define CONSTANT_MEDIUM_H

#include "rt_weekend.h"

#include "hittable.h"
#include "material.h"
#include "texture.h"

class constant_medium : public hittable {
public:
    constant_medium(mem_arena &arena, hittable* b, double density, texture* a)
        : boundary(b),
            neg_inv_density(-1 / density),
            phase_function(arena.alloc<isotropic>(a))
        {}

    constant_medium(mem_arena &arena, hittable* b, double density, colour c)
        : boundary(b),
            neg_inv_density(-1 / density),
            phase_function(arena.alloc<isotropic>(arena, c))
        {}

    virtual bool hit(
        const ray &r, double t_min, double t_max, hit_record &rec
    ) const override;

    virtual bool bounding_box(
        double time0, double time1, aabb &output_box
    ) const override {
        return boundary->bounding_box(time0, time1, output_box);
    }

private:
    hittable* boundary;
    double neg_inv_density;
    material* phase_function;
};

bool constant_medium::hit(const ray &r, double t_min, double t_max, hit_record &rec) const {
    hit_record rec1, rec2;

    if (!boundary->hit(r, -infinity, infinity, rec1))
        return false; // ray never intersects the medium

    if (!boundary->hit(r, rec1.t + epsilon, infinity, rec2)) {
        return false; // ray never exits the medium (probably grazing collision?)
    }

    // constrain to desired ray segment
    if (rec1.t < t_min) rec1.t = t_min;
    if (rec2.t > t_max) rec2.t = t_max;

    if (rec1.t >= rec2.t)
        return false;

    if (rec1.t < 0)
        rec1.t = 0;

    const auto ray_length = r.direction().mag();
    const auto distance_inside_boundary = (rec2.t - rec1.t) * ray_length;
    const auto hit_distance = neg_inv_density * log(random_number());

    if (hit_distance > distance_inside_boundary)
        return false; // interaction is beyond edge of medium

    rec.t = rec1.t + hit_distance / ray_length;
    rec.p = r.at(rec.t);

    rec.normal = random_in_unit_sphere();
    rec.front_face = true;
    rec.mat = phase_function;
    return true;
}

#endif