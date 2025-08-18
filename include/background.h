#ifndef ENV_MAP_H
#define ENV_MAP_H

#include "rt_weekend.h"

#include "hittable.h"
#include "pdf.h"
#include "texture.h"

#include <fstream>

namespace rtw_detail {
    class piecewise_distribution_1d {
    public:
        friend class piecewise_distribution_2d;
        piecewise_distribution_1d() : N(0), integral(0.0f) {};
        piecewise_distribution_1d(const std::vector<float> &pdf)
            : N(pdf.size()), func(pdf), cdf(N + 1) {
                cdf[0] = 0;
                for (size_t i = 1; i < N + 1; i++)
                    cdf[i] = cdf[i - 1] + func[i - 1] / N;
                integral = cdf[N];

                if (integral == 0) {
                    for (size_t i = 1; i < N + 1; i++)
                        cdf[i] = static_cast<float>(i) / N;
                } else {
                    for (size_t i = 1; i < N + 1; i++)
                        cdf[i] /= integral;
                }
            }

        float random(float u, float *pdf, size_t *o = nullptr) const {
            size_t offset = std::lower_bound(cdf.begin(), cdf.end(), u) - 1 - cdf.begin();
            offset = std::min(offset, N - 1);
            
            float du = u - cdf[offset];
            if ((cdf[offset + 1] - cdf[offset]) > 0)
                du /= (cdf[offset + 1] - cdf[offset]);
            if (pdf)
                *pdf = func[offset] / integral;
            if (o)
                *o = offset;

            return (offset + du) / N;
        }
        
    private:
        size_t N;
        std::vector<float> func, cdf;
        float integral;
    };

    class piecewise_distribution_2d : public probability_density_function {
    public:
        piecewise_distribution_2d() = default;
        piecewise_distribution_2d(const std::vector<float> &pdf, size_t nu, size_t nv) {
            std::vector<float> marginals;
            marginals.reserve(nv);
            conditional_pdfs.reserve(nv);
            for (size_t v = 0; v < nv; v++) {
                auto pdf_line = std::vector<float>(&pdf[v * nu], &pdf[v * nu] + nu);
                conditional_pdfs.emplace_back(pdf_line);
                marginals.push_back(conditional_pdfs[v].integral);
            }
            marginal_pdf = piecewise_distribution_1d(marginals);

            // v_counts.resize(nv);
            // u_counts.resize(nv);
            // for (auto& row_count : u_counts)
            //     row_count.resize(nu);
        }

        ~piecewise_distribution_2d() {
            // {
            //     std::ofstream vstream("marginal.txt");
            //     for (auto c : v_counts)
            //         vstream << c << ' ';
            // }

            // {
            //     std::ofstream ustream("conditional.txt");
            //     for (auto const& row : u_counts) {
            //         for (auto c : row) 
            //             ustream << c << ' ';
            //         ustream << '\n';
            //     }
            // }

            // {
            //     std::ofstream dstream("directions.txt");
            //     for (auto const& dir : directions) {
            //         dstream << dir << '\n';
            //     }
            // }
        }

        virtual vec3 generate() const override {
            float u = random_number<float>();
            float v = random_number<float>();

            size_t iv, iu;
            float d1 = marginal_pdf.random(v, nullptr, &iv);
            float d0 = conditional_pdfs[iv].random(u, nullptr, &iu);

            // v_counts[iv]++;
            // u_counts[iv][iu]++;

            auto phi = d0 * two_pi;
            auto theta = d1 * pi;
            auto sin_phi = sin(phi); auto cos_phi = cos(phi);
            auto sin_theta = sin(theta); auto cos_theta = cos(theta);
            return vec3(-cos_phi * sin_theta, -cos_theta, sin_phi * sin_theta);
        }

        virtual double value(const vec3 &direction) const override {
            // static int every = 0;
            auto unit_direction = unit_vector(direction);
            auto [u, v] = sphere_uv(unit_direction);

            // if (rtw_get_thread_num() == 0) {
            //     if (every++ % 1000 == 0)
            //         directions.push_back(unit_direction);
            // }

            size_t iu = clamp(size_t(u * conditional_pdfs[0].N), 0ULL, conditional_pdfs[0].N - 1);
            size_t iv = clamp(size_t(v * marginal_pdf.N), 0ULL, marginal_pdf.N - 1);

            auto theta = v * pi;
            auto sin_theta = sin(theta);
            //double pdf_norm = (conditional_pdfs[0].N * marginal_pdf.N) / (2 * pi * pi * sin_theta);
            double pdf_norm = 2 * pi * pi * sin_theta;

            return conditional_pdfs[iv].func[iu] / (marginal_pdf.integral * pdf_norm);
        }

    private:
        std::vector<piecewise_distribution_1d> conditional_pdfs;
        piecewise_distribution_1d marginal_pdf;

        // mutable std::vector<std::vector<int>> u_counts;
        // mutable std::vector<int> v_counts;
        // mutable std::list<vec3> directions;
    };
};

class background : public hittable {
public:
    background(mem_arena &arena, const colour& c) : background(arena.alloc<solid_colour>(c)) {}
    background(const texture* t) : background_texture(t) {
        std::vector<float> luminances;
        luminances.reserve(N);

        // std::ofstream lstream("luminance.txt");

        double total_lum = 0.0;
        for (int j = 0; j < importance_resolution; j++) {
            for (int i = 0; i < importance_resolution; i++) {
                auto phi = (i + 0.5) * two_pi / importance_resolution;
                auto theta = (j + 0.5) * pi / importance_resolution;
                auto lum = sample_luminance(phi, theta);
                // lstream << lum << ' ';
                total_lum += lum;
                luminances.push_back(lum);
            }
            // lstream << '\n';
        }

        // Normalise lumiances to get pdf
        for (auto &l: luminances)
            l /= static_cast<float>(total_lum);

        dist = rtw_detail::piecewise_distribution_2d(luminances, importance_resolution, importance_resolution);
        std::cerr << "sampled environment map along " << N << " directions.\n";
    }

    virtual bool hit(const ray& r, double t_min, double t_max, hit_record& rec) const override {
        // we don't actually need to hit against the background, but just want to be able to use PDFs
        return false;
    }

    virtual bool bounding_box(double time_0, double time_1, aabb& output_box) const override {
        return false;
    }

    virtual double pdf_value(const point3& o, const vec3& vec) const override {
        return dist.value(vec);
    }

    virtual vec3 random(const vec3& o) const override {
        return dist.generate();
    }
    
    colour value(const ray& r) const {
        auto unit_direction = unit_vector(r.direction());
        auto [u, v] = sphere_uv(unit_direction);
        return background_texture->value(u, v, unit_direction);
    }

private:
    float sample_luminance(double phi, double theta) {
        const auto d_phi = 0.5 * two_pi / importance_resolution;
        const auto d_theta = 0.5 * pi / importance_resolution;
        const int importance_samples = 20;
        
        float L = 0;
        for (int i = 0; i < importance_samples; i++) {
            auto sample_phi = phi + random_number(-d_phi, d_phi);
            auto sample_theta = theta + random_number(-d_theta, d_theta);
            auto sin_phi = sin(sample_phi), cos_phi = cos(sample_phi);
            auto sin_theta = sin(sample_theta), cos_theta = cos(sample_theta);
            auto dir = vec3(-cos_phi * sin_theta, -cos_theta, sin_phi * sin_theta);

            auto [u, v] = sphere_uv(dir);
            L += luminance(background_texture->value(u, v, vec3())) * static_cast<float>(sin_theta);
        }
        return L / importance_samples;
    }

    const texture* background_texture;
    const int importance_resolution = 1000;
    const int N = importance_resolution * importance_resolution;
    rtw_detail::piecewise_distribution_2d dist;
};

#endif