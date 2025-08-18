#ifndef BOX_H
#define BOX_H

#include "rt_weekend.h"

#include "aa_rect.h"
#include "arena.h"
#include "hittable.h"
#include "hittable_list.h"

class box : public hittable {
public:
    box(mem_arena &arena, const point3 &p0, const point3 &p1, material* mat)
    : box_min(p0), box_max(p1) {
        for (int i = 0; i < 3; i++) {
            if (box_min[i] > box_max[i])
                std::swap(box_min[i], box_max[i]);
        }

        sides.add(arena.alloc<xy_rect>(p0.x(), p1.x(), p0.y(), p1.y(), p0.z(), mat));
        sides.add(arena.alloc<xy_rect>(p0.x(), p1.x(), p0.y(), p1.y(), p1.z(), mat));

        sides.add(arena.alloc<xz_rect>(p0.x(), p1.x(), p0.z(), p1.z(), p0.y(), mat));
        sides.add(arena.alloc<xz_rect>(p0.x(), p1.x(), p0.z(), p1.z(), p1.y(), mat));

        sides.add(arena.alloc<yz_rect>(p0.y(), p1.y(), p0.z(), p1.z(), p0.x(), mat));
        sides.add(arena.alloc<yz_rect>(p0.y(), p1.y(), p0.z(), p1.z(), p1.x(), mat));
    }

    virtual bool hit(
        const ray &r, double t_min, double t_max, hit_record &rec
    ) const override {
        return sides.hit(r, t_min, t_max, rec);
    }

    virtual bool bounding_box(
        double time_0, double time_1, aabb &output_box
    ) const override {
        output_box = aabb(box_min, box_max);
        return true;
    }

    virtual double pdf_value(const point3 &origin, const vec3 &v) const override {
        auto val = sides.pdf_value(origin, v);
        return val;
    }

    virtual vec3 random(const point3 &origin) const override {
        auto r = sides.random(origin);
        return r;
    }

private:
    point3 box_min, box_max;
    hittable_list sides;
};

#endif