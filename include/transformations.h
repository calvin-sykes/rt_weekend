#include "hittable.h"

class translate : public hittable {
public:
    translate(hittable *p, const vec3 &displacement)
        : ptr(p), offset(displacement) {}

    virtual bool hit(
        const ray &r, double t_min, double t_max, hit_record& rec
    ) const override {
        ray shifted_ray(r.origin() - offset, r.direction(), r.time());

        if (!ptr->hit(shifted_ray, t_min, t_max, rec))
            return false;

        rec.p += offset;
        return true;
    }

    virtual bool bounding_box(
        double time_0, double time_1, aabb &output_box
    ) const override {
        if (!ptr->bounding_box(time_0, time_1, output_box))
            return false;

        output_box = aabb(output_box.min() + offset, output_box.max() + offset);
        return true;
    }

    virtual double pdf_value(const point3 &o, const vec3 &v) const override {
        return ptr->pdf_value(o - offset, v);
    }

    virtual vec3 random(const vec3 &o) const override {
        return ptr->random(o - offset);
    }

private:
    hittable* ptr;
    vec3 offset;
};

template<int axis>
class rotate : public hittable {
    static_assert(axis >= 0 && axis < 3, "axis must be one of 0, 1, or 2.");
public:
    rotate(hittable* p, double angle);

    virtual bool hit(
        const ray &r, double t_min, double t_max, hit_record& rec
    ) const override;

    virtual bool bounding_box(
        double time_0, double time_1, aabb &output_box
    ) const override {
        output_box = bbox;
        return has_box;
    }

    virtual double pdf_value(const point3 &o, const vec3 &v) const override;
    virtual vec3 random(const vec3 &o) const override;

private:
    hittable* ptr;
    double sin_theta, cos_theta;
    bool has_box;
    aabb bbox;
};

template<int axis>
rotate<axis>::rotate(hittable* p, double angle) : ptr(p) {
    auto radians = degrees_to_radians(angle);
    sin_theta = sin(radians); cos_theta = cos(radians);
    has_box = ptr->bounding_box(0, 1, bbox);

    point3 min(infinity); point3 max(-infinity);

    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            for (int k = 0; k < 2; k++) {
                auto x = i * bbox.max().x() + (1 - i) * bbox.min().x();
                auto y = j * bbox.max().y() + (1 - j) * bbox.min().y();
                auto z = k * bbox.max().z() + (1 - k) * bbox.min().z();

                point3 p(x, y, z), newp;

                auto e0 = (axis + 1) % 3;
                auto e1 = (axis + 2) % 3;

                newp[axis] = p[axis];
                newp[e0] = cos_theta * p[e0] - sin_theta * p[e1];
                newp[e1] = sin_theta * p[e0] + cos_theta * p[e1];

                for (int c = 0; c < 3; c++) {
                    min[c] = std::min(min[c], newp[c]);
                    max[c] = std::max(max[c], newp[c]);
                }
            }
        }
    }
    bbox = aabb(min, max);
}

template<int axis>
bool rotate<axis>::hit(const ray &r, double t_min, double t_max, hit_record &rec) const {
    auto origin = r.origin();
    auto direction = r.direction();

    auto e0 = (axis + 1) % 3;
    auto e1 = (axis + 2) % 3;

    if constexpr(axis == 1)
        std::swap(e0, e1);

    origin[e0] = cos_theta * r.origin()[e0] - sin_theta * r.origin()[e1];
    origin[e1] = sin_theta * r.origin()[e0] + cos_theta * r.origin()[e1];

    direction[e0] = cos_theta * r.direction()[e0] - sin_theta * r.direction()[e1];
    direction[e1] = sin_theta * r.direction()[e0] + cos_theta * r.direction()[e1];

    ray rotated_r(origin, direction, r.time());

    if (!ptr->hit(rotated_r, t_min, t_max, rec))
        return false;

    point3 p = rec.p;
    vec3 normal = rec.normal;

    p[e0] =  cos_theta * rec.p[e0] + sin_theta * rec.p[e1];
    p[e1] = -sin_theta * rec.p[e0] + cos_theta * rec.p[e1];

    normal[e0] =  cos_theta * rec.normal[e0] + sin_theta * rec.normal[e1];
    normal[e1] = -sin_theta * rec.normal[e0] + cos_theta * rec.normal[e1];

    rec.p = p;
    rec.normal = normal;
    return true;
}

template <int axis>
inline double rotate<axis>::pdf_value(const point3 &o, const vec3 &v) const
{
    auto e0 = (axis + 1) % 3;
    auto e1 = (axis + 2) % 3;

    if constexpr(axis == 1)
        std::swap(e0, e1);

    point3 o_rotated = o; 
    o_rotated[e0] =  cos_theta * o[e0] + sin_theta * o[e1];
    o_rotated[e1] = -sin_theta * o[e0] + cos_theta * o[e1];

    vec3 v_rotated = v;
    v_rotated[e0] =  cos_theta * v[e0] + sin_theta * v[e1];
    v_rotated[e1] = -sin_theta * v[e0] + cos_theta * v[e1];

    return ptr->pdf_value(o_rotated, v_rotated);
}

template <int axis>
inline vec3 rotate<axis>::random(const vec3 &o) const
{
    auto e0 = (axis + 1) % 3;
    auto e1 = (axis + 2) % 3;

    if constexpr(axis == 1)
        std::swap(e0, e1);

    point3 o_rotated = o; 
    o_rotated[e0] =  cos_theta * o[e0] + sin_theta * o[e1];
    o_rotated[e1] = -sin_theta * o[e0] + cos_theta * o[e1];

    return ptr->random(o_rotated);
}

using rotate_x = rotate<X>;
using rotate_y = rotate<Y>;
using rotate_z = rotate<Z>;

class flip_face : public hittable {
public:
    flip_face(hittable* p) : ptr(p) {}

    virtual bool hit(
        const ray &r, double t_min, double t_max, hit_record &rec
    ) const override {
        if (!ptr->hit(r, t_min, t_max, rec))
            return false;

        rec.front_face = !rec.front_face;
        return true;
    }

    virtual bool bounding_box(
        double time_0, double time_1, aabb &output_box
    ) const override {
        return ptr->bounding_box(time_0, time_1, output_box);
    }

    virtual double pdf_value(const point3 &o, const vec3 &v) const override {
        return ptr->pdf_value(o, v);
    };

    virtual vec3 random(const vec3 &o) const override {
        return ptr->random(o);
    }

private:
    hittable* ptr;
};

class one_sided : public hittable {
public:
    one_sided(hittable* p) : ptr(p) {}

    virtual bool hit(
        const ray &r, double t_min, double t_max, hit_record &rec
    ) const override {
        if (!ptr->hit(r, t_min, t_max, rec))
            return false;

        return rec.front_face;
    }

    virtual bool bounding_box(
        double time_0, double time_1, aabb &output_box
    ) const override {
        return ptr->bounding_box(time_0, time_1, output_box);
    }

    virtual double pdf_value(const point3 &o, const vec3 &v) const override {
        return ptr->pdf_value(o, v);
    };

    virtual vec3 random(const vec3 &o) const override {
        return ptr->random(o);
    }

private:
    hittable* ptr;
};
