#include "hittables/mesh.h"
#include "hittables/triangle.h"
#include "material_library.h"

#include <fstream>
#include <regex>
#include <sstream>
#include <string.h>
#include <vector>
//#include <unistd.h>

// Disable strict warnings for this header from the Microsoft Visual C++ compiler.
#ifdef _MSC_VER
#pragma warning (push, 0)
#define _CRT_SECURE_NO_WARNINGS
#pragma warning(disable : 4996)
#endif

mesh::mesh(mem_arena &arena, const point3& origin, double scale, double rotate_x, const char *filename, material* material) {
    mat = material;

    std::vector<point3> temp_vertices;
    std::vector<vec3> temp_normals;
    std::vector<point2> temp_uvs;

    std::vector<size_t> index_vertices;
    std::vector<size_t> index_normals;
    std::vector<size_t> index_uvs;

    std::vector<size_t> temp_index_vertices;
    std::vector<size_t> temp_index_normals;
    std::vector<size_t> temp_index_uvs;

    std::ifstream obj_file(filename);
    if (!obj_file.is_open()) {
        std::cerr << "ERROR: couldn't load obj file '" << filename << "'.\n";
        std::exit(1);
    }

    char* mtl_filename = _strdup(filename);
    char* ext_pos = strstr(mtl_filename, "obj");
    strncpy(ext_pos, "mtl", 3);
    std::optional<materials_map> mtl_lib = load_material_library(arena, mtl_filename);
    free(mtl_filename);

    if (!mtl_lib && !material) {
        std::cerr << "ERROR: failed to load mtl library and no default material provided.\n";
        std::exit(1);
    }

    std::regex face_line_regex(R"(((?:\d+\/*)+)\s*)", std::regex_constants::optimize);
    std::regex face_regex(R"((\d+)(?:\/(\d+)(?:\/?(\d+))?)?)", std::regex_constants::optimize);
    std::smatch line_match, face_match;

    std::vector<std::pair<size_t, ::material*>> material_defs;

    char line[128];
    char header[8];
    while(obj_file.getline(line, 128)) {
        std::stringstream ss(line);
        ss.getline(header, 8, ' ');
        if (!strcmp(header, "v")) { // Vertex
            point3 vtx;
            ss >> vtx;
            temp_vertices.push_back(vtx);
        } else if (!strcmp(header, "vn")) { // Normal
            vec3 norm;
            ss >> norm;
            temp_normals.push_back(norm);
        } else if (!strcmp(header, "vt")) { // UVs
            vec2 uv;
            ss >> uv;
            temp_uvs.push_back(uv);
        } else if (mtl_lib && !strcmp(header, "usemtl")) { // Material definition
            std::string material_name;
            ss >> material_name;
            if (mtl_lib->find(material_name) != mtl_lib->end()) {
                material_defs.emplace_back(index_vertices.size(), mtl_lib->at(material_name));
            }
            else {
                std::cerr << "ERROR: material '" << material_name <<"' not found in mtl library.\n";
                std::exit(1);
            }
        } else if (!strcmp(header, "f")) { // Face
            std::string pos(line);
            while (std::regex_search(pos, line_match, face_line_regex)) {
                auto face_def = line_match[1].str();
                if (std::regex_search(face_def, face_match, face_regex)) {
                    if(face_match[0].matched)
                        temp_index_vertices.push_back(stol(face_match[0].str()));
                    if(face_match[2].matched)
                        temp_index_uvs.push_back(stol(face_match[2].str()));
                    if(face_match[3].matched)
                        temp_index_normals.push_back(stol(face_match[3].str()));
                } else {
                    std::cerr << "found a face vertex, but couldn't interpret it:\n";
                    std::cerr << face_def << '\n';
                }
                pos = line_match.suffix();
            }
            size_t n_vertices = temp_index_vertices.size();

            switch(n_vertices) {
            case 3: // can add immediately
                for (auto iv : temp_index_vertices)
                    index_vertices.push_back(iv);
                for (auto iuv : temp_index_uvs)
                    index_uvs.push_back(iuv);
                for (auto in : temp_index_normals)
                    index_normals.push_back(in);
                break;
            case 4: // need to triangulate quad
                // triangle 1: 0 1 2
                // triangle 2: 2 3 0
                for (size_t idx : {0, 1, 2, 2, 3, 0}) {
                    index_vertices.push_back(temp_index_vertices[idx]);
                    if (!temp_index_uvs.empty()) index_uvs.push_back(temp_index_uvs[idx]);
                    if (!temp_index_normals.empty()) index_normals.push_back(temp_index_normals[idx]);
                } 
                break;
            default:
                std::cerr << "found a face with " << n_vertices << " vertices, skipping.\n"; 
            }

            temp_index_vertices.clear();
            temp_index_normals.clear();
            temp_index_uvs.clear();
        // } else {
        //     std::cerr << "found a line, but couldn't interpret it:\n";
        //     std::cerr << line << '\n';
        }
    }
    
    std::vector<point3> vertices;
    std::vector<vec3> normals;
    std::vector<point2> uvs;

    // Indexing
    for (size_t i = 0; i < index_vertices.size(); i++) {
        size_t vertex_index = index_vertices[i];
        vertices.push_back(temp_vertices[vertex_index - 1]);
    }

    if (index_normals.size() > 0) {
        for (size_t i = 0; i < index_normals.size(); i++) {
            size_t normal_index = index_normals[i];
            normals.push_back(temp_normals[normal_index - 1]);
        }
    }

    if (index_uvs.size() > 0) {
        for (size_t i = 0; i < index_uvs.size(); i++) {
            size_t uv_index = index_uvs[i];
            uvs.push_back(temp_uvs[uv_index - 1]);
        }
    }

    // Shift to average position in xz plane, and to minimum y = 0
    // Scale so that model fits in unit cube
    double avg_x = 0.0, avg_z = 0.0, min_y = infinity;
    point3 min_pt(infinity), max_pt(-infinity);

    double angle = -rotate_x * pi / 180.0;
    double sin_angle = sin(angle), cos_angle = cos(angle);

    for (size_t i = 0; i < vertices.size(); i++) {
        auto const& v = vertices[i];

        // Rotation about x is done first, so that model is correctly oriented
        // on y axis before shifting its base to the origin
        if (rotate_x != 0.0) {
            point3 tmp_vtx;
            tmp_vtx[0] = v[0];
            tmp_vtx[1] = cos_angle * v[1] - sin_angle * v[2];
            tmp_vtx[2] = sin_angle * v[1] + cos_angle * v[2];
            vertices[i] = tmp_vtx;
        }

        avg_x += v.x();
        avg_z += v.z();
        min_y = std::min(min_y, v.y());
        min_pt = min_vector(min_pt, v);
        max_pt = max_vector(max_pt, v);
    }

    avg_x /= vertices.size();
    avg_z /= vertices.size();
    
    auto shift = point3(avg_x, min_y, avg_z);
    auto norm = 1.0 / max_axis(max_pt - min_pt);

    for (size_t i = 0; i < vertices.size(); i++) {
        vertices[i] = (vertices[i] - shift) * norm;
        vertices[i] = vertices[i] * scale + origin;
    }

    // Construct mesh and BVH
    hittable_list triangles;

    if (material_defs.empty() || material) {
        if (normals.size() > 0) {
            for (size_t i = 0; i < vertices.size(); i += 3) {
                triangles.add(arena.alloc<triangle>(
                    vertices[i], vertices[i + 1], vertices[i + 2],
                    normals[i], normals[i + 1], normals[i + 2],
                    uvs[i], uvs[i + 1], uvs[i + 2]
                ));
            }
        }
        else {
            for (size_t i = 0; i < vertices.size(); i += 3) {
                triangles.add(arena.alloc<triangle>(
                    vertices[i], vertices[i + 1], vertices[i + 2]
                ));
            }
        }
    }
    else {
        auto current_mat = material_defs[0].second;
        size_t idx_material = 1;
        if (normals.size() > 0) {
            for (size_t i = 0; i < vertices.size(); i += 3) {
                if (i >= material_defs[idx_material].first && idx_material < material_defs.size() - 1) {
                    current_mat = material_defs[idx_material].second;
                    ++idx_material;
                }
                triangles.add(arena.alloc<textured_triangle>(
                    vertices[i], vertices[i + 1], vertices[i + 2],
                    normals[i], normals[i + 1], normals[i + 2],
                    uvs[i], uvs[i + 1], uvs[i + 2],
                    current_mat
                ));
            }
        }
        else {
            for (size_t i = 0; i < vertices.size(); i += 3) {
                if (i >= material_defs[idx_material].first && idx_material < material_defs.size() - 1) {
                    current_mat = material_defs[idx_material].second;
                    ++idx_material;
                }
                triangles.add(arena.alloc<textured_triangle>(
                    vertices[i], vertices[i + 1], vertices[i + 2],
                    current_mat
                ));
            }
        }
    }

#ifdef OLD_BVH
    mesh_bvh = bvh_node(arena, triangles, 0.0, 1.0);
#else
    mesh_bvh = bvh_accel(arena, triangles, 4, 0.0, 1.0);
#endif
}

bool mesh::hit(const ray& r, double t_min, double t_max, hit_record& rec) const {
    bool hit;
    if ((hit = mesh_bvh.hit(r, t_min, t_max, rec)) && rec.mat == nullptr)
        rec.mat = mat;
    
    return hit;
}

bool mesh::bounding_box(double time_0, double time_1, aabb& output_box) const {
    return mesh_bvh.bounding_box(time_0, time_1, output_box);
}

// Restore MSVC compiler warnings
#ifdef _MSC_VER
#pragma warning (pop)
#endif