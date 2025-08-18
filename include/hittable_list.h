#ifndef HITTABLE_LIST_H
#define HITTABLE_LIST_H

#include "rt_weekend.h"

#include "hittable.h"

#include <vector>

class hittable_list : public hittable {
public:
    hittable_list() {}
    hittable_list(hittable* object) { add(object); }

    void swap(hittable_list& other) noexcept {
        objects.swap(other.objects);
    }

    void clear() { objects.clear(); }
    void add(hittable* object) { objects.push_back(object); }
    size_t size() const { return objects.size(); }
    const std::vector<hittable*> &get() const { return objects; }
    hittable* operator[](size_t i) { return objects[i]; }
    const hittable* operator[](size_t i) const { return objects[i]; }

    virtual bool hit(
        const ray& r, double t_min, double t_max, hit_record& rec
    ) const override;
    virtual bool bounding_box(
        double time_0, double time_1, aabb& output_box
    ) const override;

    virtual double pdf_value(const point3 &o, const vec3 &v) const override;
    virtual vec3 random(const vec3 &o) const override;

private:
    std::vector<hittable*> objects;
};

#endif