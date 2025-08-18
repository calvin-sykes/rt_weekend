#include "hittable_list.h"

#include "aabb.h"

bool hittable_list::hit(const ray& r, double t_min, double t_max, hit_record& rec) const {
    hit_record temp_rec;
    bool hit_anything = false;
    auto closest_so_far = t_max;

    for (const auto &object : objects) {
        if (object->hit(r, t_min, closest_so_far, temp_rec)) {
            hit_anything = true;
            closest_so_far = temp_rec.t;
            rec = temp_rec;
        }
    }

    return hit_anything;
}

bool hittable_list::bounding_box(double time_0, double time_1, aabb& output_box) const {
    if (objects.empty())
        return false;

    aabb temp_box;
    for (const auto &object : objects) {
        if (!object->bounding_box(time_0, time_1, temp_box))
            return false;
        output_box |= temp_box;
    }

    return true;
}

double hittable_list::pdf_value(const point3 &o, const vec3 &v) const {
    if (objects.size() == 0)
        return 0;

    auto weight = 1.0 / objects.size();
    auto sum = 0.0;

    for (const auto &object : objects)
        sum += weight * object->pdf_value(o, v);
    return sum;

    /*std::vector<double> pdfs(objects.size());
    for (const auto& object : objects) {
        double pdf = object->pdf_value(o, v);
        if (pdf > 0)
            pdfs.push_back(pdf);
    }

    if (pdfs.size() == 0)
        return 0.0;

    double sum = 0.0;
    double weight = 1.0 / pdfs.size();

    for (const auto pdf : pdfs)
        sum += weight * pdf;
    return sum;*/
}

vec3 hittable_list::random(const vec3 &o) const {
    if (objects.size() == 0)
        return vec3(0);

    auto int_size = static_cast<int>(objects.size());
    return objects[random_number(0, int_size - 1)]->random(o);
}
