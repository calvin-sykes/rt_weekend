#ifndef MATERIAL_H
#define MATERIAL_H

#include "rt_weekend.h"

#include "hittable.h"
#include "onb.h"
#include "pdf.h"
#include "texture.h"

struct scatter_record {
    ray specular_ray;
    bool pure_specular;
    colour attenuation;
    probability_density_function* pdf_ptr;
};

class material {
public:
    virtual ~material() {}
    virtual colour emitted(
        const ray &r_in, const hit_record &rec
    ) const {
        return colour(0, 0, 0);
    }

    virtual bool scatter(
        const ray &r_in, const hit_record &rec, scatter_record &srec, mem_arena &arena
    ) const {
        return false;
    }

    virtual colour eval_scattering(
        const ray &r_in, const ray &r_scattered, const hit_record &rec, const scatter_record &srec
    ) const {
        return colour(0.0);
    }
};

// Perfect diffuse material
class lambertian : public material {
public:
    lambertian(mem_arena &arena, const colour &a) : albedo(arena.alloc<solid_colour>(a)) {}
    lambertian(const texture* a) : albedo(a) {}

    virtual bool scatter(
        const ray &r_in, const hit_record &rec, scatter_record &srec, mem_arena &arena
    ) const override {
        srec.pure_specular = false;
        srec.attenuation = albedo->value(rec.u, rec.v, rec.p);
		srec.pdf_ptr = arena.alloc<cosine_pdf>(rec.normal);
        return true;
    }

    colour eval_scattering(
        const ray &r_in, const ray &r_scattered, const hit_record &rec, const scatter_record &srec
    ) const override {
        auto cosine = dot(rec.normal, unit_vector(r_scattered.direction()));
        return srec.attenuation * std::max(0.0, cosine * inv_pi);
    }

private:
    const texture* albedo;
};

// Perfect mirror material
class metal : public material {
public:
    metal(mem_arena &arena, const colour& a, double f)
        : albedo(arena.alloc<solid_colour>(a)), fuzz(f < 1 ? f : 1) {}

    virtual bool scatter(
        const ray& r_in, const hit_record& rec, scatter_record &srec, mem_arena &arena
    ) const override {
        vec3 unit_direction = unit_vector(r_in.direction());
        double cos_theta = std::min(dot(-unit_direction, rec.normal), 1.0);

        srec.attenuation = reflectance(cos_theta, albedo->value(rec.u, rec.v, rec.p));
        srec.pure_specular = true;
        srec.pdf_ptr = nullptr;

        vec3 reflected = reflect(unit_direction, rec.normal);
        srec.specular_ray = ray(rec.p, reflected + fuzz * random_in_unit_sphere(), r_in.time());
        return true;
    }

private:
    const texture* albedo;
    double fuzz;

    static colour reflectance(double cosine, const colour &albedo) {
        const auto white = colour(1);
        return albedo + (white - albedo) * pow(1 - cosine, 5);
    }
};

// Glass material
class dielectric : public material {
public:
    dielectric(double index_of_refraction, double index_scatter = 0.0)
        : eta(index_of_refraction), d_eta(index_scatter) {}

    virtual bool scatter(
        const ray& r_in, const hit_record& rec, scatter_record &srec, mem_arena &arena
    ) const override {
        srec.attenuation = colour(1.0);
        srec.pure_specular = true;
        srec.pdf_ptr = nullptr;
        
        auto refractive_index = random_number(-d_eta, d_eta) + eta;
        auto refraction_ratio = rec.front_face ? (1.0 / refractive_index) : refractive_index;

        vec3 unit_direction = unit_vector(r_in.direction());
        double cos_theta = std::min(dot(-unit_direction, rec.normal), 1.0);
        double sin_theta = sqrt(1.0 - cos_theta * cos_theta);

        bool cannot_refract = refraction_ratio * sin_theta > 1.0;
        vec3 direction;
        if (cannot_refract || reflectance(cos_theta, refraction_ratio) > random_number())
            direction = reflect(unit_direction, rec.normal);
        else
            direction = refract(unit_direction, rec.normal, refraction_ratio);
        
        srec.specular_ray = ray(rec.p, direction, r_in.time());
        return true;
    }

private:
    double eta; // index of refraction
    double d_eta; // scatter in eta (0 = transparent, >0 = translucent)

    static double reflectance(double cosine, double ref_idx) {
        auto r0 = (1 - ref_idx) / (1 + ref_idx);
        r0 *= r0;
        return r0 + (1 - r0) * pow(1 - cosine, 5);
    }
};

class coloured_dielectric : public material {
public:
    coloured_dielectric(mem_arena &arena, const colour &a, double transmittance, double index_of_refraction, double index_scatter = 0.0)
        : albedo(arena.alloc<solid_colour>(a)), abs(1.0 - transmittance), eta(index_of_refraction), d_eta(index_scatter) {}
    
    coloured_dielectric(const texture *a, double transmittance, double index_of_refraction, double index_scatter = 0.0)
        : albedo(a), abs(1.0 - transmittance), eta(index_of_refraction), d_eta(index_scatter) {}

    virtual bool scatter(
        const ray& r_in, const hit_record& rec, scatter_record &srec, mem_arena &arena
    ) const override {
        if (rec.front_face) {
            srec.attenuation = colour(1.0);
        } else {
            auto dist = rec.t * r_in.direction().mag();
            srec.attenuation = albedo->value(rec.u, rec.v, rec.p) * exp(-1 * abs * dist);
        }
        srec.pure_specular = true;
        srec.pdf_ptr = nullptr;
        
        auto refractive_index = random_number(-d_eta, d_eta) + eta;
        auto refraction_ratio = rec.front_face ? (1.0 / refractive_index) : refractive_index;

        vec3 unit_direction = unit_vector(r_in.direction());
        double cos_theta = std::min(dot(-unit_direction, rec.normal), 1.0);
        double sin_theta = sqrt(1.0 - cos_theta * cos_theta);

        bool cannot_refract = refraction_ratio * sin_theta > 1.0;
        vec3 direction;
        if (cannot_refract || reflectance(cos_theta, refraction_ratio) > random_number())
            direction = reflect(unit_direction, rec.normal);
        else
            direction = refract(unit_direction, rec.normal, refraction_ratio);
        
        srec.specular_ray = ray(rec.p, direction, r_in.time());
        return true;
    }

private:
    const texture *albedo;
    double abs; // fractional absorbance per unit length
    double eta, d_eta;

    static double reflectance(double cosine, double ref_idx) {
        auto r0 = (1 - ref_idx) / (1 + ref_idx);
        r0 *= r0;
        return r0 + (1 - r0) * pow(1 - cosine, 5);
    }
};

// Isotropic scattering material
class isotropic : public material {
public:
    isotropic(mem_arena &arena, const colour& c) : albedo(arena.alloc<solid_colour>(c)) {}
    isotropic(const texture* a) : albedo(a) {}

    virtual bool scatter(
        const ray &r_in, const hit_record &rec, scatter_record &srec, mem_arena &arena
    ) const override {
        srec.attenuation = albedo->value(rec.u, rec.v, rec.p);
        srec.pure_specular = false;
		srec.pdf_ptr = arena.alloc<uniform_pdf>();
        return true;
    }

    colour eval_scattering(
        const ray &r_in, const ray &r_scattered, const hit_record &rec, const scatter_record &srec
    ) const override {
         return inv_four_pi * srec.attenuation;
    }

private:
    const texture* albedo;
};

class diffuse_light : public material {
public:
    diffuse_light(mem_arena &arena, const colour &c, bool two_sided = false) : emit(arena.alloc<solid_colour>(c)), two_sided(two_sided) {}
    diffuse_light(const texture* a, bool two_sided = true) : emit(a), two_sided(two_sided) {}

    virtual bool scatter(
        const ray &r_in, const hit_record &rec, scatter_record &srec, mem_arena &arena
    ) const override {
        return false;
    }

    virtual colour emitted(const ray &r_in, const hit_record &rec) const override {
        if (two_sided || rec.front_face)
            return emit->value(rec.u, rec.v, rec.p);
        else
            return colour(0);
    }

private:
    const texture* emit;
    bool two_sided;
};

// Phong specular highlight material
class phong : public material {
public:
    phong(mem_arena &arena, const colour &c, double specular_fraction, double shininess)
        : albedo(arena.alloc<solid_colour>(c)), fspec(specular_fraction), shine(shininess) {};
    phong(const texture *a, double specular_fraction, double shininess)
        : albedo(a), fspec(specular_fraction), shine(shininess) {};

    virtual bool scatter(
        const ray& r_in, const hit_record& rec, scatter_record &srec, mem_arena &arena
    ) const override {
        if (random_number() < fspec) {
            srec.attenuation = colour(1.0);
        } else {
            srec.attenuation = albedo->value(rec.u, rec.v, rec.p);
        }
        srec.pure_specular = false;
        srec.pdf_ptr = arena.alloc<blinn_phong_pdf>(rec.normal, r_in.direction(), shine, fspec);
        return true;
    }

    colour eval_scattering(
        const ray &r_in, const ray &r_scattered, const hit_record &rec, const scatter_record &srec
    ) const override {
        auto cosine = std::max(0.0, dot(rec.normal, unit_vector(r_scattered.direction())));
        auto half_vector = unit_vector(unit_vector(r_scattered.direction()) - unit_vector(r_in.direction()));
        auto cos_alpha = std::max(0.0, dot(half_vector, rec.normal));
        auto specular = (shine + 1) * inv_two_pi * pow(cos_alpha, shine);
        return cosine * ((1 - fspec) * srec.attenuation * inv_pi + fspec * colour(1.0) * specular);
        // if (srec.type == DIFFUSE) {
        //     return cosine * srec.attenuation * inv_pi;
        // } else {
        //     auto half_vector = unit_vector(unit_vector(r_scattered.direction()) - unit_vector(r_in.direction()));
        //     auto cos_alpha = std::max(0.0, dot(half_vector, rec.normal));
        //     auto specular = (shine + 1) * inv_two_pi * pow(cos_alpha, shine);
        //     return cosine * colour(1.0) * specular;
        // }
    };

private:
    const texture *albedo;
    double fspec, shine;
};

class ashikhmin_shirley : public material {
    public:
    ashikhmin_shirley(mem_arena &arena, const colour &c_diffuse, const colour &c_specular,
        double specular_fraction, double nu, double nv
    ) : r_d(arena.alloc<solid_colour>(c_diffuse)), r_s(arena.alloc<solid_colour>(c_specular)),
        fspec(specular_fraction), nu(nu), nv(nv) {}

    ashikhmin_shirley(const texture *a_diffuse, const texture *a_specular,
        double specular_fraction, double nu, double nv
    ) : r_d(a_diffuse), r_s(a_specular), fspec(specular_fraction), nu(nu), nv(nv) {}

    virtual bool scatter(
        const ray& r_in, const hit_record& rec, scatter_record &srec, mem_arena &arena
    ) const override {
        srec.attenuation = fspec * r_s->value(rec.u, rec.v, rec.p) + (1 - fspec) * r_d->value(rec.u, rec.v, rec.p);
        srec.pure_specular = false;
        srec.pdf_ptr = arena.alloc<ashikhmin_shirley_pdf>(rec.normal, r_in.direction(), nu, nv, fspec);
        return true;
    }

    colour eval_scattering(
        const ray &r_in, const ray &r_scattered, const hit_record &rec, const scatter_record &srec
    ) const override {
        auto v = -unit_vector(r_in.direction());
        auto l = unit_vector(r_scattered.direction());

        if (dot(l, rec.normal) < 0)
            return colour(0.0);
        
        auto h = unit_vector(v + l);
        auto rs_corr = r_s->value(rec.u, rec.v, rec.p) * fspec;
        auto rd_corr = r_d->value(rec.u, rec.v, rec.p) * (1 - fspec);
        orthonormal_basis onb_normal{ rec.normal };

        auto exponent = nu * pow(dot(h, onb_normal.u()), 2) + nv * pow(dot(h, onb_normal.v()), 2);
        exponent /= 1.0 - pow(dot(h, onb_normal.w()), 2);

        auto denominator = dot(h, v) * std::max(dot(rec.normal, v), dot(rec.normal, l));
        auto fresnel = rs_corr + (colour(1.0) - rs_corr) * pow(1 - dot(v, h), 5);

        auto specular_brdf = (1.0 / (8 * pi)) * fresnel * sqrt((nu + 1) * (nv + 1)) \
            * pow(dot(rec.normal, h), exponent) / denominator;
        
        auto diff_const = rd_corr * (colour(1.0) - rs_corr) * 28 / (23.0 * pi);
        auto diffuse_brdf = diff_const \
            * (1.0 - pow(1.0 - dot(rec.normal, v) / 2.0, 5)) \
            * (1.0 - pow(1.0 - dot(rec.normal, l) / 2.0, 5));
        
        return dot(rec.normal, l) * (diffuse_brdf + specular_brdf);
    };

private:
    const texture *r_d, *r_s;
    double fspec, nu, nv;
};

#endif