#include "bvh.h"

using rtw_detail::bvh_info;
using rtw_detail::bucket_info;
using rtw_detail::bvh_accel_build_node;

bvh_accel::bvh_accel(mem_arena& arena, const hittable_list& list, int max_objects_per_leaf, double time_0, double time_1)
 : leaf_size(std::min(255, max_objects_per_leaf)), objects(list), nodes(nullptr) {
    // Store data about each hittable
    aabb tmp_box;
    aabb union_box;
    std::vector<bvh_info> build_data;
    build_data.reserve(objects.size());
    for (size_t i = 0; i < objects.size(); i++) {
        objects[i]->bounding_box(time_0, time_1, tmp_box);
        build_data.emplace_back(i, tmp_box);
        union_box |= tmp_box;
    }

    mem_arena build_arena(1024 * 1024);
    hittable_list ordered_objects;

    size_t total_nodes = 0;
    bvh_accel_build_node* root = build_recursive(build_arena, build_data, ordered_objects, 0, objects.size(), &total_nodes);
    objects.swap(ordered_objects);

    //nodes = static_cast<bvh_accel_node*>(arena.alloc(total_nodes * sizeof(bvh_accel_node)));
    nodes = alloc_aligned<bvh_accel_node>(total_nodes);
    size_t offset = 0;
    build_flattened(root, &offset);
    bbox = union_box;
    std::cerr << "built BVH structure with " << total_nodes << " nodes.\n";
}

bvh_accel_build_node* bvh_accel::build_recursive(
    mem_arena& arena, std::vector<bvh_info>& build_data,
    hittable_list& ordered_objects, size_t start, size_t end, size_t* total_nodes
) {
    (*total_nodes)++;

    aabb node_bbox = build_data[start].box;
    for (size_t i = start + 1; i < end; i++)
        node_bbox |= build_data[i].box;

    auto object_span = end - start;
    if (object_span <= 2) {
        auto offset_first = ordered_objects.size();
        for (size_t i = start; i < end; i++) {
            auto idx = build_data[i].index;
            ordered_objects.add(objects[idx]);
        }
        return arena.alloc<bvh_accel_build_node>(offset_first, object_span, node_bbox);
    } else {
        aabb centroid_bbox;
        for (size_t i = start; i < end; i++)
            centroid_bbox |= build_data[i].centre;

        int axis = centroid_bbox.longest_axis();
        auto mid = (start + end) / 2;

        // If centroid bounds have zero volume, create a leaf
        if (axis == -1) {
            auto offset_first = ordered_objects.size();
            for (size_t i = start; i < end; i++) {
                auto idx = build_data[i].index;
                ordered_objects.add(objects[idx]);
            }
            return arena.alloc<bvh_accel_build_node>(offset_first, object_span, node_bbox);
        } else {
            if (object_span <= 4) {
                // Partition small numbers of objects into equal sized sets
                std::nth_element(&build_data[start], &build_data[mid], &build_data[end-1] + 1,
                    [axis](const bvh_info& a, const bvh_info& b) {
                        return a.centre[axis] < b.centre[axis];
                    });
            } else {
                // Partition objects into equal-sized buckets and select split minimising SAH
                const int n_buckets = 12;
                bucket_info buckets[n_buckets];

                for (size_t i = start; i < end; i++) {
                    int b = static_cast<int>(n_buckets * centroid_bbox.relative_position(build_data[i].centre)[axis]);
                    if (b == n_buckets) b = n_buckets - 1;
                    buckets[b].count++;
                    buckets[b].box |= build_data[i].box;
                }

                // Estimate SAH cost for each bucket
                double cost[n_buckets - 1];
                for (size_t i = 0; i < n_buckets - 1; i++) {
                    aabb b0, b1;
                    int count0 = 0, count1 = 1;
                    for (size_t j = 0; j <= i; j++) {
                        b0 |= buckets[j].box;
                        count0 += buckets[j].count;
                    }
                    for (size_t j = i + 1; j < n_buckets; j++) {
                        b1 |= buckets[j].box;
                        count1 += buckets[j].count;
                    }
                    cost[i] = 1 + (count0 * b0.surface_area() + count1 * b1.surface_area()) / node_bbox.surface_area();
                }

                double min_cost = cost[0];
                int min_cost_bucket = 0;
                for (int i = 1; i < n_buckets - 1; i++) {
                    if (cost[i] < min_cost) {
                        min_cost = cost[i];
                        min_cost_bucket = i;
                    }
                }

                double leaf_cost = 1.0 * object_span;
                if (object_span > leaf_size || min_cost < leaf_cost) {
                    const bvh_info* pmid = std::partition(&build_data[start], &build_data[end - 1] + 1,
                        [=](const bvh_info& info) {
                            int b = static_cast<int>(n_buckets * centroid_bbox.relative_position(info.centre)[axis]);
                            if (b == n_buckets) b = n_buckets - 1;
                            return b <= min_cost_bucket;
                        });
                    mid = pmid - &build_data[0]; // convert from pointer to index offset
                } else {
                    auto offset_first = ordered_objects.size();
                    for (size_t i = start; i < end; i++) {
                        auto idx = build_data[i].index;
                        ordered_objects.add(objects[idx]);
                    }
                    return arena.alloc<bvh_accel_build_node>(offset_first, object_span, node_bbox);
                }
            }
        }
        bvh_accel_build_node* left = build_recursive(arena, build_data, ordered_objects, start, mid, total_nodes);
        bvh_accel_build_node* right = build_recursive(arena, build_data, ordered_objects, mid, end, total_nodes);
        return arena.alloc<bvh_accel_build_node>(axis, left, right);
    }
}

size_t bvh_accel::build_flattened(bvh_accel_build_node* node, size_t* offset) {
    bvh_accel_node* flat_node = &nodes[*offset];
    flat_node->box = node->box;
    size_t my_offset = (*offset)++;

    if (node->n_obj > 0) {
        flat_node->first_obj = node->first_obj;
        flat_node->n_obj = static_cast<uint16_t>(node->n_obj);
    } else {
        flat_node->axis = static_cast<uint16_t>(node->split_axis);
        flat_node->n_obj = 0;
        build_flattened(node->children[0], offset);
        flat_node->second_child = build_flattened(node->children[1], offset);
    }
    return my_offset;
}

bool bvh_accel::hit(const ray& r, double t_min, double t_max, hit_record& rec) const {
    bool hit = false;
    size_t stack_depth = 0, current_index = 0;
    size_t nodes_to_visit[64];
    
    auto const& origin = r.origin();
    auto const& dir_inv = r.direction_inverse();
    auto const& dir_signs = r.direction_signs();

    while (true) {
        const bvh_accel_node& node = nodes[current_index];
        //if (node.box.hit(r, t_min, t_max)) {
        if (node.box.hit(origin, dir_inv, dir_signs, t_min, t_max)) {
            if (node.n_obj > 0) {
                for (int i = 0; i < node.n_obj; i++) {
                    if (objects[node.first_obj + i]->hit(r, t_min, t_max, rec)) {
                        hit = true;
                        t_max = rec.t;
                    }
                }
                if (stack_depth == 0) break;
                current_index = nodes_to_visit[--stack_depth];
            } else {
                if (dir_signs[node.axis]) {
                    nodes_to_visit[stack_depth++] = current_index + 1;
                    current_index = node.second_child;
                } else {
                    nodes_to_visit[stack_depth++] = node.second_child;
                    current_index = current_index + 1;
                }
            }
        } else {
            if (stack_depth == 0) break;
            current_index = nodes_to_visit[--stack_depth];
        }
    }
    return hit;
}

bool bvh_accel::bounding_box(double time_0, double time_1, aabb& output_box) const {
    output_box = bbox;
    return true;
}

double bvh_accel::pdf_value(const point3& o, const vec3& v) const
{
    return objects.pdf_value(o, v);
}

vec3 bvh_accel::random(const vec3& o) const
{
    return objects.random(o);
}
