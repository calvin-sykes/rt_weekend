#ifndef MOVING_SPHERE_H
#define MOVING_SPHERE_H

#include "rt_weekend.h"

#include "aabb.h"
#include "hittable.h"
#include "vector.h"

class moving_sphere : public hittable {
public:
    moving_sphere (
        point3 centre_0, point3 centre_1, double time_0, double time_1,
        double radius, material* material
    )
        : cen0(centre_0), cen1(centre_1), t0(time_0), t1(time_1),
          rad(radius), mat(material) {}

    virtual bool hit(
        const ray& r, double t_min, double t_max, hit_record& rec
    ) const override;
    virtual bool bounding_box(
        double time_0, double time_1, aabb& output_box
    ) const override;

    point3 centre(double time) const {
        return cen0 + ((time - t0) / (t1 - t0)) * (cen1 - cen0);
    }

private:
        point3 cen0, cen1;
        double t0, t1;
        double rad;
        material* mat;
};

#endif