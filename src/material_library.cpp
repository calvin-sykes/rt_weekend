#include "material_library.h"
#include "texture.h"

#include <fstream>
#include <sstream>
#include <string.h>

std::optional<materials_map> load_material_library(mem_arena& arena, const char* filename) {
    std::string path(filename);
    size_t last_sep = path.find_last_of("/");
    path = path.substr(0, last_sep + 1);

    std::ifstream mtl_file(filename);
    if (!mtl_file.is_open()) {
        return std::nullopt;
    }

    materials_map materials;
    
    double Ns, d;
    colour /*Ka,*/ Kd, Ks;
    std::optional<std::string> map_Kd;

    char line[128];
    std::string header;
    std::string name;
    bool havemat = false;
    while (mtl_file.getline(line, 128)) {
        std::stringstream ss(line);
        ss >> std::skipws >> header;
        if (header == "newmtl") {
            ss >> name;
            havemat = true;
        }
        else if (header == "Ns") {
            ss >> Ns;
        }
        else if (header == "d") {
            ss >> d;
        }
        /*} else if (!strcmp(header, "Ka")) {
                ss >> Ka;*/
        else if (header == "Kd") {
            ss >> Kd;
        }
        else if (header == "Ks") {
            ss >> Ks;
        }
        else if (header == "map_Kd") {
            std::string tmp;
            ss >> tmp;
            map_Kd = path + tmp;
        }
        else if (havemat && header.empty()) {
            texture* tex_diffuse;
            texture* tex_specular = arena.alloc<solid_colour>(Ks);

            if (map_Kd)
                tex_diffuse = arena.alloc<image_texture>(map_Kd->c_str());
            else
                tex_diffuse = arena.alloc<solid_colour>(Kd);
            float fspec = Ks.mag() / Kd.mag();
            fspec = std::min(std::max(0.0f, fspec), 1.0f);
            materials[name] = arena.alloc<ashikhmin_shirley>(tex_diffuse, tex_specular, fspec, Ns, Ns);
            //materials[name] = arena.alloc<lambertian>(tex_diffuse);
            map_Kd.reset();
            havemat = false;
        }
        header.clear();
    }

    // Parse the final material
    if (havemat) {
        texture* tex_diffuse;
        texture* tex_specular = arena.alloc<solid_colour>(Ks);

        if (map_Kd)
            tex_diffuse = arena.alloc<image_texture>(map_Kd->c_str());
        else
            tex_diffuse = arena.alloc<solid_colour>(Kd);
        float fspec = Ks.mag() / Kd.mag();
        fspec = std::min(std::max(0.0f, fspec), 1.0f);
        materials[name] = arena.alloc<ashikhmin_shirley>(tex_diffuse, tex_specular, fspec, Ns, Ns);
        //materials[name] = arena.alloc<lambertian>(tex_diffuse);
        map_Kd.reset();
        havemat = false;
    }

    return materials;
}
