#include "hittable.h"
#include "material.h"

class triangle : public hittable {
public:
    triangle(const point3& p0, const point3& p1, const point3& p2);
    triangle(const point3& p0, const point3& p1, const point3& p2,
             const vec3& n0, const vec3& n1, const vec3& n2,
             const point2 &t0, const point2 &t1, const point2 &t2);

    virtual bool hit(
        const ray& r, double t_min, double t_max, hit_record& rec
    ) const override;
    virtual bool bounding_box(
        double time_0, double time_1, aabb& output_box
    ) const override;

private:
    point3 p;
    vec3 normals[3];
    point2 uv_coords[3];
    vec3 edge1, edge2;
    aabb bbox;
};

class textured_triangle : public triangle {
public:
    textured_triangle(const point3& p0, const point3& p1, const point3& p2, material* mat) :
        triangle(p0, p1, p2), mat(mat) {}
    textured_triangle(const point3& p0, const point3& p1, const point3& p2,
        const vec3& n0, const vec3& n1, const vec3& n2,
        const point2& t0, const point2& t1, const point2& t2, material* mat) :
        triangle(p0, p1, p2, n0, n1, n2, t0, t1, t2), mat(mat) {}

    virtual bool hit(
        const ray& r, double t_min, double t_max, hit_record& rec
    ) const override {
        bool hit = triangle::hit(r, t_min, t_max, rec);
        rec.mat = mat;
        return hit;
    }

private:
    material* mat;
};