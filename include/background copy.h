#ifndef ENV_MAP_H
#define ENV_MAP_H

#include "rt_weekend.h"

#include "hittable.h"
#include "texture.h"

template<typename Tval, typename Tweight>
class walker_sampler {
public:
    walker_sampler() : N(0) {}
    walker_sampler(const std::vector<Tval> &v, const std::vector<Tweight> &w)
        : N(v.size()), values(v), aliases(N), P(N)
    {
        for (int i = 0; i < N; i++)
            weights.emplace_back(w[i], i);

        for (int k = N-1; k >= 0; k--) {
            weights.sort(
                [](const auto& a, const auto& b) {
                    return a.first < b.first;
                });
            auto& [q1, a1] = weights.front();
            auto& [qk, ak] = weights.back();
            P[a1] = N * q1;
            aliases[a1] = values[ak];
            qk -= (1.0f / N - q1);
            weights.pop_front();
        }
    }

    Tval sample() const {
        double v = random_number<Tweight>();
        int j = random_number(0, N - 1);

        Tval val;
        if (v < P[j])
            return values[j];
        else
            return aliases[j];
    }

private:
    int N;
    std::vector<Tval> values, aliases;
    std::list<std::pair<Tweight, size_t>> weights;
    std::vector<size_t> P;
};

class background : public hittable {
public:
    background(mem_arena &arena, const colour& c) : background_texture(arena.alloc<solid_colour>(c)) {}
    background(const texture* t);

    virtual bool hit(const ray& r, double t_min, double t_max, hit_record& rec) const override {
        // hacky, we don't actually need to hit against this, but just want to be able to use PDFs
        return false;
    }

    virtual bool bounding_box(double time_0, double time_1, aabb& output_box) const override {
        return false;
    }

    virtual double pdf_value(const point3& o, const vec3& vec) const override {
        double u, v;
        get_sphere_uv(vec, u, v);

        int i = round(u * importance_resolution - 0.5);
        int j = round(v * importance_resolution - 0.5);

        auto theta = v * pi;
        // auto phi = u * two_pi;
        auto sin_theta = sin(theta);

        // const auto d_phi = 0.5 * two_pi / importance_resolution;
        // const auto d_theta = 0.5 * pi / importance_resolution * pi;
        
        return pdf[i * importance_resolution + j] * importance_resolution * importance_resolution / sin_theta;
    }

    virtual vec3 random(const vec3& o) const override {
        return pdf_sampler.sample();
    }
    
    colour value(const ray& r) const {
        auto unit_direction = unit_vector(r.direction());
        
        double u, v;
        get_sphere_uv(unit_direction, u, v);
        return background_texture->value(u, v, unit_direction);

        // int i = round(u * importance_resolution - 0.5);
        // int j = round(v * importance_resolution - 0.5);
        // auto dir = directions[i * importance_resolution + j];
        // get_sphere_uv(dir, u, v);
        // return background_texture->value(u, v, dir);
    }

private:
    static void get_sphere_uv(const point3& p, double& u, double& v) {
        auto theta = my_acos(-p.y());
        auto phi = my_atan2(-p.z(), p.x()) + pi;

        u = phi * inv_two_pi;
        v = theta * inv_pi;
    }

    static float luminance(const colour& c) {
        return 0.2126f * c.x() + 0.7152f * c.y() + 0.0722f * c.z();
    }

    float sample_luminance(double phi, double theta) {
        const auto d_phi = 0.5 * two_pi / importance_resolution;
        const auto d_theta = 0.5 * pi / importance_resolution * pi;
        const int importance_samples = 100;
        
        float L = 0;
        for (int i = 0; i < importance_samples; i++) {
            auto sample_phi = phi + random_number(-d_phi, d_phi);
            auto sample_theta = theta + random_number(-d_theta, d_theta);
            auto sin_phi = sin(sample_phi); auto cos_phi = cos(sample_phi);
            auto sin_theta = sin(sample_theta); auto cos_theta = cos(sample_theta);
            auto dir = vec3(-cos_phi * sin_theta, -cos_theta, sin_phi * sin_theta);

            double u, v;
            get_sphere_uv(dir, u, v);
            L += luminance(background_texture->value(u, v, vec3()));
        }
        return L / importance_samples;
    }

    const texture* background_texture;
    const int importance_resolution = 150;
    const int N = importance_resolution * importance_resolution;
    std::vector<float> pdf;
    walker_sampler<vec3,float> pdf_sampler;
    std::vector<vec3> directions;
};

background::background(const texture* t) : background_texture(t) {
    directions.reserve(N);
    pdf.reserve(N);

    double total_lum = 0.0;
    for (int i = 0; i < importance_resolution; i++) {
        for (int j = 0; j < importance_resolution; j++) {
            auto phi = (i + 0.5) * two_pi / importance_resolution;
            auto theta = (j + 0.5) * pi / importance_resolution;
            auto sin_phi = sin(phi); auto cos_phi = cos(phi);
            auto sin_theta = sin(theta); auto cos_theta = cos(theta);

            auto lum = static_cast<double>(sample_luminance(phi, theta)) * sin_theta;
            total_lum += lum;
            pdf.push_back(lum);

            auto dir = vec3(-cos_phi * sin_theta, -cos_theta, sin_phi * sin_theta);
            directions.push_back(dir);
        }
    }

    // Normalise the probabilities
    const auto luminance_norm = static_cast<float>(1.0 / total_lum);
    for (auto &prob : pdf)
        prob *= luminance_norm;

    // Construct sampler
    pdf_sampler = walker_sampler(directions, pdf);
}
#endif