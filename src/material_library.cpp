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
    colour Kd, Ks;
    std::optional<std::string> map_Kd, map_Ks;

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
        } else if (header == "Ns") {
            ss >> Ns;
        } else if (header == "d") {
            ss >> d;
        } else if (header == "Kd") {
            ss >> Kd;
        } else if (header == "Ks") {
            ss >> Ks;
        } else if (header == "map_Kd") {
            std::string tmp;
            ss >> tmp;
            map_Kd = path + tmp;
        } else if (header == "map_Ks") {
            std::string tmp;
            ss >> tmp;
            map_Ks = path + tmp;
        } else if (havemat && header.empty()) {
            texture *tex_diffuse, *tex_specular;
            
            if (map_Kd)
                tex_diffuse = arena.alloc<image_texture>(map_Kd->c_str(), Kd);
            else
                tex_diffuse = arena.alloc<solid_colour>(Kd);
            
            if (map_Ks)
                tex_specular = arena.alloc<image_texture>(map_Ks->c_str(), Ks);
            else
                tex_specular = arena.alloc<solid_colour>(Ks);
            
            float fspec = clamp(Ks.mag() / Kd.mag(), 0.0f, 1.0f);
            //materials[name] = arena.alloc<lambertian>(tex_diffuse);
            materials[name] = arena.alloc<phong>(tex_diffuse, tex_specular, fspec, Ns);
            // materials[name] = arena.alloc<ashikhmin_shirley>(tex_diffuse, tex_specular, fspec, Ns, Ns);

            Ns = d = 0.0;
            Kd = Ks = colour();
            map_Kd.reset();
            map_Ks.reset();
            havemat = false;
        }
        header.clear();
    }

    // Parse the final material
    if (havemat) {
        texture *tex_diffuse, *tex_specular;
            
        if (map_Kd)
            tex_diffuse = arena.alloc<image_texture>(map_Kd->c_str(), Kd);
        else
            tex_diffuse = arena.alloc<solid_colour>(Kd);
        
        if (map_Ks)
            tex_specular = arena.alloc<image_texture>(map_Ks->c_str(), Ks);
        else
            tex_specular = arena.alloc<solid_colour>(Ks);
        
        float fspec = clamp(Ks.mag() / Kd.mag(), 0.0f, 1.0f);
        //materials[name] = arena.alloc<lambertian>(tex_diffuse);
        materials[name] = arena.alloc<phong>(tex_diffuse, tex_specular, fspec, Ns);
        // materials[name] = arena.alloc<ashikhmin_shirley>(tex_diffuse, tex_specular, fspec, Ns, Ns);

        Ns = d = 0.0;
        Kd = Ks = colour();
        map_Kd.reset();
        map_Ks.reset();
        havemat = false;
    }

    return materials;
}
