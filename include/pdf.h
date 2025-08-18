#ifndef PDF_H
#define PDF_H

#include "rt_weekend.h"

#include "hittable_list.h"
#include "onb.h"

/**
 * A probability density function.
 **/
class probability_density_function {
public:
    virtual ~probability_density_function() {}

    /**
     * @brief Return the value of the PDF for a given direction vector.
     * @param direction The direction vector to evaluate the PDF for.
     **/
    virtual double value(const vec3 &direction) const = 0;

    /**
     * @brief Return a random vector distributed according to the PDF.
     **/
    virtual vec3 generate() const = 0;
};

class cosine_pdf : public probability_density_function {
public:
    cosine_pdf(const vec3& normal) : uvw(normal) {}

    virtual double value(const vec3 &direction) const override {
        auto cosine = dot(unit_vector(direction), uvw.w());
        return std::max(0.0, cosine * inv_pi);
    }

    virtual vec3 generate() const override {
        return uvw.local_to_world(random_cosine_direction());
    }

private:
    orthonormal_basis uvw;
};

class uniform_pdf : public probability_density_function {
public:
    uniform_pdf() {}

    virtual double value(const vec3 &direction) const override {
       return inv_four_pi;
    }

    virtual vec3 generate() const override {
        return random_unit_vector();
    }
};

class blinn_phong_pdf : public probability_density_function {
public:
    blinn_phong_pdf(const vec3& normal, const vec3& in_dir, double shininess, double fspec)
        : shine(shininess), fspec(fspec), rin_direction(in_dir),
        onb_normal(normal), onb_reflected(reflect(in_dir, normal)) {};

    virtual double value(const vec3 &direction) const override {
        auto half_vector = unit_vector(unit_vector(direction) - unit_vector(rin_direction));
        auto cos_alpha = std::max(0.0, dot(half_vector, onb_normal.w()));

        auto cosine = std::max(0.0, dot(unit_vector(direction), onb_normal.w()));
        return cosine * ((1 - fspec) * inv_pi) + fspec * (shine + 1) * inv_two_pi * pow(cos_alpha, shine);
    }

    virtual vec3 generate() const override {
        if (random_number() < fspec) {
            vec3 direction;
            int iters = 0;
            while (true) {
                double u = random_number();
                double cos_theta = pow(1 - u, 1.0 / (2.0 + shine));
                double sin_theta = sqrt(1 - cos_theta * cos_theta);
                double phi = random_number(0.0, two_pi);
                auto v = vec3(cos(phi) * sin_theta, sin(phi) * sin_theta, cos_theta);

                direction = onb_reflected.local_to_world(v);

                if ((dot(direction, onb_normal.w()) > 0) || ++iters > 100)
                    break;
            }
            return direction;
        } else {
            return onb_normal.local_to_world(random_cosine_direction());   
        }
    }

private:
    double shine, fspec;
    vec3 rin_direction;
    orthonormal_basis onb_normal, onb_reflected;
};

class ashikhmin_shirley_pdf : public probability_density_function {
public:
    ashikhmin_shirley_pdf(const vec3& normal, const vec3& in_direction, double nu, double nv, double fspec)
        : rin_direction(in_direction), onb_normal(normal), nu(nu), nv(nv), fspec(fspec) {}

    virtual double value(const vec3 &direction) const override {
        auto v = -unit_vector(rin_direction);
        auto l = unit_vector(direction);
        auto h = unit_vector(v + l);

        if (dot(l, onb_normal.w()) < 0)
            return 0.0;

        auto exponent = nu * pow(dot(h, onb_normal.u()), 2) + nv * pow(dot(h, onb_normal.v()), 2);
        exponent /= 1.0 - pow(dot(h, onb_normal.w()), 2);

        auto pdf_h = sqrt((nu + 1) * (nv + 1)) * inv_two_pi * pow(dot(h, onb_normal.w()), exponent);
        auto cosine = std::max(0.0, dot(l, onb_normal.w()));

        return ((1.0 - fspec) * cosine * inv_pi) + (fspec * pdf_h / (4 * dot(v, h)));
    }

    virtual vec3 generate() const override {
        if (random_number() < fspec) {
            auto h = onb_normal.local_to_world(ashikhmin_shirley_random_vector());
            auto v = unit_vector(rin_direction);
            return unit_vector(reflect(v, h));
        } else {
            return onb_normal.local_to_world(random_cosine_direction());
        }
    }

private:
    vec3 ashikhmin_shirley_random_vector() const {
        double r1 = random_number();
        double r1_corr;
        int correction;

        if (r1 < 0.25) {
            r1_corr = 1.0 - 4.0 * (0.25 - r1);
            correction = 0;
        } else if (r1 < 0.5) {
            r1_corr = 1.0 - 4.0 * (0.5 - r1);
            correction = 1;
        } else if (r1 < 0.75) {
            r1_corr = 1.0 - 4.0 * (0.75 - r1);
            correction = 2;
        } else {
            r1_corr = 1.0 - 4.0 * (1.0 - r1);
            correction = 3;
        }

        double phi = atan(sqrt((nu + 1) / (nv + 1)) * tan(pi * r1_corr * 0.5));
        double phi_corr;

        switch(correction) {
        case 0:
            phi_corr = phi; break;
        case 1:
            phi_corr = pi - phi; break;
        case 2:
            phi_corr = pi + phi; break;
        case 3:
            phi_corr = two_pi - phi; break;
        }

        double r2 = random_number();
        
        double exponent = 1.0 / (nu * pow(cos(phi_corr), 2) + nv * pow(sin(phi_corr), 2) + 1.0);
        double cos_theta = pow(1.0 - r2, exponent);
        double sin_theta = sqrt(1.0 - cos_theta * cos_theta);

        return vec3(cos(phi_corr) * sin_theta, sin(phi_corr) * sin_theta, cos_theta);
    }

    vec3 rin_direction;
    orthonormal_basis onb_normal;
    double nu, nv, fspec;
};

class hittable_pdf : public probability_density_function {
public:
    hittable_pdf(const hittable_list &_objects, const point3 &_origin)
        : objects(_objects), origin(_origin) {}

    virtual double value(const vec3 &direction) const override {
        return objects.pdf_value(origin, direction);
    }

    virtual vec3 generate() const override {
        return objects.random(origin);
    }

private:
    const hittable_list &objects;
    point3 origin;
};

class mixture_pdf : public probability_density_function {
public:
    mixture_pdf(const probability_density_function* p0, const probability_density_function* p1, double chance = 0.5)
        : chance(chance), p{p0, p1} {};

    virtual double value(const vec3 &direction) const override {
        return chance * p[0]->value(direction) + (1 - chance) * p[1]->value(direction);
    }

    virtual vec3 generate() const override {
        if (random_number() < chance)
            return p[0]->generate();
        else
            return p[1]->generate();
    }

private:
    double chance;
    const probability_density_function* p[2];
};

#endif