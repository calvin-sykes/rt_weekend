#include "hittables/triangle.h"

triangle::triangle(const point3& p0, const point3& p1, const point3& p2
) : p(p0), edge1(p1 - p0), edge2(p2 - p0)
{
    auto n = unit_vector(cross(edge1, edge2));
    normals[0] = normals[1] = normals[2] = n;
    uv_coords[0] = uv_coords[1] = uv_coords[2] = point2(0.5);

    auto minimum = min_vector(min_vector(p0, p1), p2);
    auto maximum = max_vector(max_vector(p0, p1), p2);
    bbox = aabb(minimum, maximum);
}

triangle::triangle(
    const point3& p0, const point3& p1, const point3& p2,
    const vec3& n0, const vec3& n1, const vec3& n2,
    const point2 &t0, const point2 &t1, const point2 &t2
) : normals{n0, n1, n2}, uv_coords{t0, t1, t2} {
    p = p0;
    edge1 = p1 - p0;
    edge2 = p2 - p0;

    auto minimum = min_vector(min_vector(p0, p1), p2);
    auto maximum = max_vector(max_vector(p0, p1), p2);
    bbox = aabb(minimum, maximum);
}

bool triangle::hit(const ray& r, double t_min, double t_max, hit_record& rec) const {
    auto h = cross(r.direction(), edge2);
    auto a = dot(edge1, h);
    if (fabs(a) < epsilon)
        return false; // Parallel to triangle

    auto inv_det = 1.0 / a;
    auto dist = r.origin() - p;
    auto u = inv_det * dot(dist, h);
    if (u < 0.0 || u > 1.0)
        return false; // Misses triangle
    
    auto q = cross(dist, edge1);
    auto v = inv_det * dot(r.direction(), q);
    if (v < 0.0 || u + v > 1.0)
        return false; // Misses triangle

    auto t = inv_det * dot(edge2, q);
    if (t < t_min || t > t_max)
        return false; // Not on ray
   
    rec.t = t;
    rec.p = r.at(rec.t);

    auto norm_interp = (1.0 - u - v) * normals[0] + u * normals[1] + v * normals[2];
    // auto norm_interp =  0.333 * (normals[0] + normals[1] + normals[2]); 
    rec.set_face_normal(r, unit_vector(norm_interp));

    auto uv_interp = (1.0 - u - v) * uv_coords[0] + u * uv_coords[1] + v * uv_coords[2];
    // auto uv_interp = 0.333 * (uv_coords[0] + uv_coords[1] + uv_coords[2]);
    rec.u = uv_interp.x();
    rec.v = uv_interp.y();

    rec.mat = nullptr; // NB: set by owning mesh
    return true;
}

bool triangle::bounding_box(double time_0, double time_1, aabb& output_box) const {
    output_box = bbox;
    return true;
}