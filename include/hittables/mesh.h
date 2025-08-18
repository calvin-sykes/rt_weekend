#ifndef MESH_H
#define MESH_H

#include "rt_weekend.h"

#include "aabb.h"
#include "hittable.h"
#include "hittable_list.h"
#include "vector.h"

//#define OLD_BVH
#ifdef BVH_OLD
#include "bvh_old.h"
#else
#include "bvh.h"
#endif

class mesh : public hittable {
public:
    mesh(mem_arena& arena, const point3& origin, double scale, double rotate_x, const char* filename, material* material);
    mesh(mem_arena& arena, const point3& origin, double scale, const char* filename, material* material) :
        mesh(arena, origin, scale, 0.0, filename, material) {};


    virtual bool hit(
        const ray& r, double t_min, double t_max, hit_record& rec
    ) const override;
    virtual bool bounding_box(
        double time_0, double time_1, aabb& output_box
    ) const override;

private:
    material *mat;
#ifdef BVH_OLD
    bvh_node mesh_bvh;
#else
    bvh_accel mesh_bvh;
#endif
};

#endif