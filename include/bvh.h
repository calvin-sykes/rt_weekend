#ifndef BVH_FAST_H
#define BVH_FAST_H

#include "rt_weekend.h"

#include "arena.h"
#include "hittable.h"
#include "hittable_list.h"

#include <algorithm>

namespace rtw_detail {
    struct bvh_info {
        bvh_info(size_t idx, const aabb& bbox) : index(idx), box(bbox) {
            centre = 0.5 * (box.min() + box.max());
        }
        size_t index;
        aabb box;
        point3 centre;
    };

    struct bucket_info {
        int count = 0;
        aabb box;
    };

    struct bvh_accel_build_node {
        bvh_accel_build_node() = default;

        bvh_accel_build_node(size_t offset_first, size_t n, const aabb& b)
            : box(b), split_axis(-1), first_obj(offset_first), n_obj(n)
        {
            children[0] = children[1] = nullptr;
        }

        bvh_accel_build_node(int axis, bvh_accel_build_node* c0, bvh_accel_build_node* c1)
            : box(c0->box | c1->box), split_axis(axis), first_obj(0), n_obj(0)
        {
            children[0] = c0; children[1] = c1;
        }

        aabb box;
        bvh_accel_build_node* children[2];
        size_t split_axis, first_obj, n_obj;
    };
};

struct bvh_accel_node {
    aabb box;
    union { size_t first_obj, second_child; };
    uint16_t n_obj;
    uint16_t axis;
};
//static_assert(sizeof(bvh_accel_node)==64ULL);

class bvh_accel : public hittable {
public:
    bvh_accel() : leaf_size(1), nodes(nullptr) {};

    bvh_accel(mem_arena &arena, const hittable_list& list, int max_objects_per_leaf, double time_0, double time_1);

    bvh_accel(const bvh_accel&) = delete;
    bvh_accel operator=(const bvh_accel&) = delete;

    bvh_accel(bvh_accel&& other) noexcept
        : leaf_size(other.leaf_size) {
        objects.swap(other.objects);
        nodes = other.nodes;
        other.nodes = nullptr;
        bbox = other.bbox;
    }

    bvh_accel& operator=(bvh_accel&& other) noexcept {
        leaf_size = other.leaf_size;
        objects.swap(other.objects);
        nodes = other.nodes;
        other.nodes = nullptr;
        bbox = other.bbox;
        return *this;
    }

    ~bvh_accel() {
        if (nodes)
            free_aligned(nodes);
    }

    virtual bool hit(
        const ray& r, double t_min, double t_max, hit_record& rec) const override;
    virtual bool bounding_box(
        double time_0, double time_1, aabb& output_box) const override;

    virtual double pdf_value(const point3& o, const vec3& v) const override;
    virtual vec3 random(const vec3& o) const override;

private:
    rtw_detail::bvh_accel_build_node* build_recursive(
        mem_arena &arena, std::vector<rtw_detail::bvh_info> &build_data,
        hittable_list &ordered_objects, size_t start, size_t end, size_t *total_nodes
    );

    size_t build_flattened(rtw_detail::bvh_accel_build_node* node, size_t* offset);

    unsigned int leaf_size;
    hittable_list objects;
    bvh_accel_node* nodes;
    aabb bbox;
};

#endif