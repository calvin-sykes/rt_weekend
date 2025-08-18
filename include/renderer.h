#ifndef RENDERER_H
#define RENDERER_H

#include "background.h"
#include "hittable_list.h"
#include "scene.h"
#include "writer.h"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

enum FILE_TYPE {
        FILE_TYPE_PPM,
        FILE_TYPE_PNG,
        FILE_TYPE_HDR,
        FILE_TYPE_TEXT,
        N_FILE_TYPES,
        FILE_TYPE_STDOUT = static_cast<unsigned>(-1)
};
const char *const _filetypes[N_FILE_TYPES] = {"ppm", "png", "hdr", "txt"};
#define FILE_TYPE_STRING(x) (_filetypes[x])

template<bool use_roulette = true>
extern colour ray_colour(
    mem_arena &arena, ray &r, const background &bkg,
    const hittable_list &world, const hittable_list &lights, int depth
);

template <typename T>
class online_stats {
public:
    online_stats() : count(0), mu(T(0.0)), M2(T(0.0)) {}

    void update(const T& new_value) {
        ++count;
        T delta = new_value - mu;
        mu += delta / count;
        T delta2 = new_value - mu;
        M2 += delta * delta2;

        // samples.push_back(new_value);
    }

    T mean() { return mu; }
    // T median() {
    //     size_t n = samples.size() / 2;
    //     std::nth_element(
    //         samples.begin(), samples.begin() + n, samples.end(),
    //         [](const T& a, const T& b) { return luminance(a) < luminance(b);}
    //     );
    //     return samples[n];
    // }
    T variance(size_t k = 1) {
        if (count < 2)
            return 0.0;
        else
            return M2 / (count - k);
    }

private:
    size_t count;
    T mu, M2;
    // std::vector<T> samples;
};

class renderer {
public:
    renderer(
        int image_width,
        int image_height,
        int samples_per_pixel,
        int resample_factor,
        int max_depth
    ) : w(image_width * resample_factor), h(image_height * resample_factor),
        ow(image_width), oh(image_height),
        num_samples(samples_per_pixel),
        resample_factor(resample_factor),
        max_depth(max_depth) {
        buffer.reserve(w * h);
        for (int i = 0; i < w * h; i++)
            buffer.emplace_back(0.0, 0.0, 0.0);
        variance_buffer.reserve(w * h);
        for (int i = 0; i < w * h; i++)
            variance_buffer.emplace_back(0.0, 0.0, 0.0);
    }

    void render(const scene &scene);
    std::vector<colour> resample() const;
    void dump(FILE_TYPE type, const char *filename) const;
    void dump_variance(FILE_TYPE type, const char *filename) const;
    void clear() {
        buffer.assign(buffer.size(), 0.0);
        variance_buffer.assign(variance_buffer.size(), 0.0);
    }

private:
    const int w, h, ow, oh;
    const int num_samples, resample_factor, max_depth;
    std::vector<colour> buffer;
    std::vector<colour> variance_buffer;

    static colour sanitise_pixel(const colour &c_in) {
        colour c_out;
        for (int i = 0; i < 3; i++) {
            if (c_in[i] != c_in[i])
                c_out[i] = (i == 1 ? 0.0f : 1.0f);
            else if (c_in[i] == infinity)
                c_out[i] = (i == 0 ? 0.0f : 1.0f);
            else
                c_out[i] = c_in[i];
        }
        return c_out;
    }

    static void tone_mapping(std::vector<colour>& buffer) {
        auto max = std::max_element(
            buffer.begin(), buffer.end(),
            [](colour const& a, colour const& b) {
                return luminance(a) < luminance(b);
            }
        );
        auto white_point = luminance(*max);

        for (size_t i = 0; i < buffer.size(); i++) {
            auto l_old = luminance(buffer[i]);
            auto numerator = l_old * (1.0f + (l_old / (white_point * white_point)));
            auto l_new = numerator / (1.0f + l_old);
            buffer[i] *= l_new / l_old;

            auto r = buffer[i].x();
            auto g = buffer[i].y();
            auto b = buffer[i].z();

            // Gamma correction
            constexpr float gamma = 1.0f / 2.2f;
            r = powf(r, gamma);
            g = powf(g, gamma);
            b = powf(b, gamma);

            buffer[i] = colour(r, g, b);
        }
    }

    colour resample_kernel(int i, int j) const {
        auto c = colour(0.0);
        for (int jj = 0; jj < resample_factor; ++jj) {
            for (int ii = 0; ii < resample_factor; ++ii) {
                auto idx = (i * resample_factor + ii) + w * (j * resample_factor + jj);
                c += buffer[idx];
            }
        }
        return c / float(resample_factor * resample_factor);        
    }
};

#endif