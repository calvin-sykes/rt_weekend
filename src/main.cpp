#include "rt_weekend.h"

#include "scene.h"
#include "renderer.h"

#include <iostream>
#include <string>

int main(int argc, char **argv) {
    std::string filename;
    FILE_TYPE filetype;

    // int size = 400;
    // uint8_t *buf = new uint8_t[size * size * 3];

    // colour c1(1.0);
    // colour c2(0.8);

    // for (int i = 0; i < size; i++) {
    //     for (int j = 0; j < size; j++) {
    //         int dx = i - size / 2;
    //         int dy = j - size / 2;
    //         float dist = std::max(abs(dx), abs(dy)) * 1.0f / (size / 2);
    //         float lerp = pow(dist, 5);
    //         colour clerp = lerp * c2 + (1 - lerp) * c1;
    //         buf[3 * (i + size * j)] = static_cast<uint8_t>(clerp[0] * 255);
    //         buf[3 * (i + size * j) + 1] = static_cast<uint8_t>(clerp[1] * 255);
    //         buf[3 * (i + size * j) + 2] = static_cast<uint8_t>(clerp[2] * 255);
    //     }
    // }

    // stbi_write_png("media/box_gradient.png", size, size, 3, buf, 0);
    // delete[] buf;

    if (argc < 2) {
        std::cout << "Filename? (blank for stdout)\n";
        std::getline(std::cin, filename);
    }
    else {
        filename.assign(argv[1]);
    }

    if (filename.empty()) {
        filetype = FILE_TYPE_STDOUT;
    } else {
        auto idx = filename.find(".");
        if (idx != std::string::npos) {
            auto ext = filename.substr(idx + 1);
            bool found = false;
            for (auto i = 0U; i < N_FILE_TYPES; i++) {
                if (ext == FILE_TYPE_STRING(i)) {
                    filetype = static_cast<FILE_TYPE>(i);
                    found = true;
                    break;
                }
            }
            if (!found) {
                std::cerr << "unknown file extension '" << ext << "'\n";
                exit(1);
            }
        } else {
            std::cerr << "filename must be of the form <name>.<ext>\n";
            exit(1);
        }
    }

    // Image
    //const auto aspect_ratio = 4.0 / 3.0;
    //const auto aspect_ratio = 16.0 / 9.0;
    const auto aspect_ratio = 1.0;
    
    const int image_width = 400;
    const int image_height = static_cast<int>(image_width / aspect_ratio);
    
    const int samples_per_pixel = 1000;
    const int supersample = 1;
    const int max_depth = 10;
    
    mem_arena arena;
    scene scn = get_scene(arena, SCENE_CORNELL_BOX_GLASS_SPHERES,
        static_cast<size_t>(image_width) * supersample,
        static_cast<size_t>(image_height) * supersample);

    std::cerr << "Scene contains " << scn.world.size() << \
        " objects using " << arena.total_allocated() << " bytes.\n";

    // Render
    renderer rendrr(image_width, image_height, samples_per_pixel, supersample, max_depth);
    rendrr.render(scn);
    rendrr.dump(filetype, filename.c_str());

    auto idx = filename.find(".");
    std::string filename_hdr = filename.replace(idx+1, 3, "hdr");
    rendrr.dump(FILE_TYPE_HDR, filename_hdr.c_str());

    // std::string variance_filename("variance_" + filename.replace(idx+1, 3, "png"));
    // rendrr.dump_variance(FILE_TYPE_PNG, variance_filename.c_str());

    arena.reset();
}
