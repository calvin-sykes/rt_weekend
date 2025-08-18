#ifndef BVH_H
#define BVH_H

#include "rt_weekend.h"

#include "arena.h"
#include "hittable.h"
#include "hittable_list.h"

#include <algorithm>

class bvh_node : public hittable
{
public:
    bvh_node() : left(nullptr), right(nullptr) {}

    bvh_node(mem_arena &arena, const hittable_list &list, double time_0, double time_1)
        : bvh_node(arena, list.get(), 0, list.size(), time_0, time_1) {
        std::cerr << "built BVH structure with " << num_nodes << " nodes.\n";
    }

    bvh_node(
        mem_arena &arena, const std::vector<hittable*> &src_objects,
        size_t start, size_t end, double time_0, double time_1);

    virtual bool hit(
        const ray &r, double t_min, double t_max, hit_record &rec) const override;
    virtual bool bounding_box(
        double time_0, double time_1, aabb &output_box) const override;

private:
    static size_t num_nodes;
    hittable *left, *right;
    aabb box;
};

size_t bvh_node::num_nodes = 0;

template <int axis>
inline static bool box_compare(const hittable* a, const hittable* b)
{
    aabb box_a, box_b;

    if (!a->bounding_box(0, 0, box_a) || !b->bounding_box(0, 0, box_b))
        std::cerr << "No bounding box in box_compare.\n";

    return box_a.min()[axis] < box_b.min()[axis];
}

bvh_node::bvh_node(
    mem_arena &arena, const std::vector<hittable*> &src_objects,
    size_t start, size_t end, double time_0, double time_1)
{
    num_nodes++;
    std::vector<hittable*> objects = src_objects;
    size_t object_span = end - start;

#if defined(BVH_SAH) or defined(BVH_LONGEST)
    // Construct bounding box enclosing all objects and find longest axis
    aabb main_box;
    objects[start]->bounding_box(time_0, time_1, main_box);
    for (size_t i = 1; i < object_span; i++)
    {
        aabb new_box;
        objects[start + i]->bounding_box(time_0, time_1, new_box);
        main_box |= new_box;
    }

    int axis = main_box.longest_axis();
    auto comparator = (axis == 0) ? box_compare<0>
                    : (axis == 1) ? box_compare<1>
                                  : box_compare<2>;
#else
    // Choose axis randomly
    auto axis = random_number<int>(0, 2);
    auto comparator = (axis == 0) ? box_compare<0>
                    : (axis == 1) ? box_compare<1>
                                  : box_compare<2>;
#endif

#ifdef BVH_SAH
    switch (object_span)
    {
    case 2:
        left = objects[start];
        right = objects[start + 1];
        break;
    case 3:
        left = arena.alloc<bvh_node>(arena, objects, start, start + 2, time_0, time_1);
        right = objects[start + 2];
        break;
    default:
        std::sort(objects.begin() + start, objects.begin() + end, comparator);

        auto boxes = new aabb[object_span];
        auto left_area = new double[object_span];
        auto right_area = new double[object_span];

        for (size_t i = 0; i < object_span; i++)
            objects[start + i]->bounding_box(time_0, time_1, boxes[i]);

        left_area[0] = boxes[0].surface_area();
        aabb left_box = boxes[0];
        for (size_t i = 1; i < object_span - 1; i++)
        {
            left_box |= boxes[i];
            left_area[i] = left_box.surface_area();
        }

        right_area[object_span - 1] = boxes[object_span - 1].surface_area();
        aabb right_box = boxes[object_span - 1];
        for (size_t i = object_span - 2; i > 0; i--)
        {
            right_box |= boxes[i];
            right_area[i] = right_box.surface_area();
        }

        // Calculate surface area heuristic for each partition and find minimum
        auto min_SAH = infinity;
        size_t min_SAH_idx = 0;
        for (size_t i = 0; i < object_span - 1; i++)
        {
            auto SAH = i * left_area[i] + (object_span - i - 1) * right_area[i + 1];
            if (SAH < min_SAH)
            {
                min_SAH_idx = i;
                min_SAH = SAH;
            }
        }

        if (min_SAH_idx == 0)
            left = objects[start];
        else
            left = arena.alloc<bvh_node>(
                arena, objects, start, start + min_SAH_idx + 1, time_0, time_1);

        if (min_SAH_idx == object_span - 2)
            right = objects[start + object_span - 1];
        else
            right = arena.alloc<bvh_node>(
                arena, objects, start + min_SAH_idx + 1, end, time_0, time_1);

        delete[] boxes;
        delete[] left_area;
        delete[] right_area;

        break;
    }
#else
    switch (object_span)
    {
    case 2:
        left = objects[start];
        right = objects[start + 1];
        break;
    case 3:
        left = arena.alloc<bvh_node>(arena, objects, start, start + 2, time_0, time_1);
        right = objects[start + 2];
        break;
    default:
        std::sort(objects.begin() + start, objects.begin() + end, comparator);

        auto mid = start + object_span / 2;
        left = arena.alloc<bvh_node>(arena, objects, start, mid, time_0, time_1);
        right = arena.alloc<bvh_node>(arena, objects, mid, end, time_0, time_1);
        break;
    }
#endif

    aabb box_left, box_right;

    if (!left->bounding_box(time_0, time_1, box_left) || !right->bounding_box(time_0, time_1, box_right))
        std::cerr << "No bounding box in bvh_node constructor.\n";

    box = box_left | box_right;
}

bool bvh_node::hit(const ray &r, double t_min, double t_max, hit_record &rec) const
{
    /*if (!box.hit(r, t_min, t_max))
        return false;*/
    if (!box.hit(r.origin(), r.direction_inverse(), r.direction_signs(), t_min, t_max))
        return false;

    bool hit_left = left->hit(r, t_min, t_max, rec);
    bool hit_right = right->hit(r, t_min, hit_left ? rec.t : t_max, rec);

    return hit_left || hit_right;
}

bool bvh_node::bounding_box(double time_0, double time_1, aabb &output_box) const
{
    output_box = box;
    return true;
}

#endif