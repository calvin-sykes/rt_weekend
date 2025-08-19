#include "rt_weekend.h"

#include "renderer.h"

#include "camera.h"
#include "material.h"
#include "pdf.h"

#include <sstream>

template<bool use_roulette>
colour ray_colour(
    mem_arena &arena, ray &r, const background &bkg,
    const hittable_list &world, const hittable_list &lights, int depth
) {
    auto c = colour(0.0);
    auto throughput = colour(1.0);

    for (int i = 0; i < depth; i++) {
        hit_record rec;

        // If the ray hits nothing, add the background colour and terminate.
        if (!world.hit(r, epsilon, infinity, rec)) {
            c += throughput * bkg.value(r);
            break;
        }

        // Add the emitted light
        c += throughput * rec.mat->emitted(r, rec);

        // If the ray did not scatter, terminate
        scatter_record srec;
        if (!rec.mat->scatter(r, rec, srec, arena))
            break;

        // If the scattering was pure specular, follow the specular ray
        if (srec.pure_specular) {
            r = srec.specular_ray;
            throughput *= srec.attenuation;
            continue;
        } else {
            ray scattered;
            double pdf_val;

            if (lights.size() > 0) {
                // If we have targets for importance sampling, construct a
                // mixture pdf for the implicit and directly sampled rays
                auto light_pdf = hittable_pdf(lights, rec.p);
                mixture_pdf mix_pdf(&light_pdf, srec.pdf_ptr);
                scattered = ray(rec.p, mix_pdf.generate(), r.time());
                pdf_val = mix_pdf.value(scattered.direction());
                assert(!std::isnan(pdf_val));
            } else {
                // Otherwise, just sample from the implicit ray
                scattered = ray(rec.p, srec.pdf_ptr->generate(), r.time());
                pdf_val = srec.pdf_ptr->value(scattered.direction());
                assert(!std::isnan(pdf_val));
            }

            // If the scattering is impossible, terminate.
            if (pdf_val == 0)
                break;

            // auto scattering_pdf = rec.mat->scattering_pdf(r, rec, srec, scattered);
            // throughput *= srec.attenuation * (scattering_pdf / pdf_val);
            throughput *= rec.mat->eval_scattering(r, scattered, rec, srec) / pdf_val;

            if constexpr(use_roulette) {
                if (luminance(throughput) < 0.1) {
                    auto p = max_axis(throughput);
                    if (random_number() < p)
                        break;
    
                    // Add energy from terminated paths
                    throughput *= 1.0f / (1 - p);
                }
            }
            r = scattered;
        }
    }

    return c;
}

void renderer::render(const scene &scn) {
    auto const& world = scn.world;
    auto const& lights = scn.lights;
    auto const& cam = scn.cam;
    auto const& background = scn.bkg;
    
    size_t total_size = static_cast<size_t>(w) * static_cast<size_t>(h);
    size_t step_size = std::max(total_size / 1000ULL, 100ULL);
    int n_threads = rtw_get_num_threads();
    size_t local_count_max = step_size / n_threads;

    volatile size_t count = 0;
    volatile int threads_finished = 0;

    clock_t t0 = clock();

    std::cerr << std::setprecision(3);
    #pragma omp parallel num_threads(n_threads)
    {
        long long reported_count = 0;
        unsigned int local_count = 0;
        size_t msg_length = 0;
        int tid = rtw_get_thread_num();
        mem_arena arena;

        #pragma omp for schedule(dynamic, 100) nowait
        for (int n = 0; n < h * w; n++) {
            int i = n % w;
            int j = n / w;

            online_stats<colour> stddev;
            online_stats<float> stddev_lum;
            for (int s = 0; s < num_samples; s++) {
                ray r = cam.get_ray(i, j);
                auto pix_colour = ray_colour<false>(arena, r, background, world, lights, max_depth);
                stddev.update(pix_colour);
                stddev_lum.update(luminance(pix_colour));
            }
            auto idx = i + w * j;
            buffer[idx] = sanitise_pixel(stddev.mean());
            variance_buffer[idx] = stddev_lum.variance();
			arena.reset();

            if (++local_count >= local_count_max) {
                #pragma omp atomic
                count += local_count;
                local_count = 0;
            }
            if (tid == 0 && (count - reported_count >= step_size)) {
                auto elapsed = static_cast<long long>(clock()) - t0;
                auto est_remaining = (elapsed * static_cast<long long>(total_size)) / count - elapsed;

                std::cerr << "\r" << std::string(msg_length, ' ');
                std::ostringstream ss;
#if defined(PLATFORM_MACOS)
                ss << "\rTraced " << count * num_samples << '/' << static_cast<long long>(total_size) * num_samples << " rays ("
                    << (1.0 * elapsed) / (n_threads * CLOCKS_PER_SEC) << "s elapsed / "
                    << (1.0 * est_remaining) / (n_threads * CLOCKS_PER_SEC) << "s remaining, "
                    << (count * num_samples * n_threads * CLOCKS_PER_SEC) / elapsed << " rays/sec, "
                    << (count * n_threads * CLOCKS_PER_SEC) / elapsed << " px/sec)";
                    #elif defined(PLATFORM_WINDOWS)
                ss << "\rTraced " << count * num_samples << '/' << static_cast<long long>(total_size) * num_samples << " rays ("
                    << (1.0 * elapsed) / CLOCKS_PER_SEC << "s elapsed / "
                    << (1.0 * est_remaining) / CLOCKS_PER_SEC << "s remaining, "
                    << (count * num_samples * CLOCKS_PER_SEC) / elapsed << " rays/sec, "
                    << (count * CLOCKS_PER_SEC) / elapsed << " px/sec)";
#endif
                auto msg = ss.str();
                msg_length = msg.size();
                std::cerr << msg;
                reported_count = count;
            }
        }
        
        #pragma omp atomic
        threads_finished++;

        if (tid == 0) {
            std::cerr << '\n' << std::flush;
            while (true) {
                if (count - reported_count >= step_size) {
                    std::cerr << "\r" << threads_finished << '/' << n_threads << " threads completed tracing";
                    reported_count = count;
                }
                if(threads_finished == n_threads) break;
            }
        }
    }
    std::cerr << "\nDone.\n";
}

std::vector<colour> renderer::resample() const {
    if (resample_factor == 1) {
        return buffer;
    } else {
        std::vector<colour> resampled_image(ow * oh);
        for (int j = oh - 1; j >= 0; --j) {
            for (int i = 0; i < ow; i++) {
                auto idx = i + ow * j;
                resampled_image[idx] = resample_kernel(i, j);
            }
        }
        return resampled_image;
    }
}

void renderer::dump(FILE_TYPE ftype, const char *filename) const {
    auto resampled_buffer = resample();
    
    if (ftype != FILE_TYPE_HDR) {
        tone_mapping(resampled_buffer);
    }

    std::ofstream os;
    switch (ftype) {
    default:
        std::cerr << "invalid file type.\n";
    case FILE_TYPE_STDOUT:
        ppm_writer(nullptr).write(resampled_buffer, ow, oh);
        break;
    case FILE_TYPE_PPM:
        ppm_writer(filename).write(resampled_buffer, ow, oh);
        break;
    case FILE_TYPE_PNG:
        png_writer(filename).write(resampled_buffer, ow, oh);
        break;
    case FILE_TYPE_HDR:
        hdr_writer(filename).write(resampled_buffer, ow, oh);
        break;
    case FILE_TYPE_TEXT:
        text_writer(filename).write(resampled_buffer, ow, oh);
        break;
    }
}

void renderer::dump_variance(FILE_TYPE ftype, const char *filename) const {
   std::ofstream os;
    switch (ftype) {
    default:
        std::cerr << "invalid file type.\n";
    case FILE_TYPE_STDOUT:
        ppm_writer(nullptr).write(variance_buffer, w, h);
        break;
    case FILE_TYPE_PPM:
        ppm_writer(filename).write(variance_buffer, w, h);
        break;
    case FILE_TYPE_PNG:
        png_writer(filename).write(variance_buffer, w, h);
        break;
    case FILE_TYPE_HDR:
        hdr_writer(filename).write(variance_buffer, w, h);
        break;
    case FILE_TYPE_TEXT:
        text_writer(filename).write(variance_buffer, w, h);
        break;
    }
}
