#ifndef MATERIAL_LIBRARY_H
#define MATERIAL_LIBRARY_H

#include "material.h"

#include <optional>
#include <string>
#include <unordered_map>

using materials_map = std::unordered_map<std::string, material*>;

std::optional<materials_map> load_material_library(mem_arena& arena, const char* filename);

#endif // MATERIAL_LIBRARY_H