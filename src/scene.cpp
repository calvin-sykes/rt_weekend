#include "scene.h"

#include "background.h"
#include "bvh.h"
#include "camera.h"
#include "material.h"
#include "transformations.h"

#include "hittables/aa_rect.h"
#include "hittables/box.h"
#include "hittables/constant_medium.h"
#include "hittables/mesh.h"
#include "hittables/moving_sphere.h"
#include "hittables/sphere.h"

hittable_list *get_axes_marker(mem_arena &arena, const point3 &origin, double scale)
{
    auto red = arena.alloc<lambertian>(arena, colour(1.0, 0.0, 0.0));
    auto green = arena.alloc<lambertian>(arena, colour(0.0, 1.0, 0.0));
    auto blue = arena.alloc<lambertian>(arena, colour(0.0, 0.0, 1.0));

    double l = scale;
    double w = 0.05 * scale;

    hittable_list *axes = arena.alloc<hittable_list>();
    axes->add(arena.alloc<translate>(arena.alloc<box>(arena, point3( 0, -w, -w), point3(l, w, w), red), origin));   // x
    axes->add(arena.alloc<translate>(arena.alloc<box>(arena, point3(-w,  0, -w), point3(w, l, w), green), origin)); // y
    axes->add(arena.alloc<translate>(arena.alloc<box>(arena, point3(-w, -w,  0), point3(w, w, l), blue), origin));  // z

    return axes;
}

hittable_list world_book1(mem_arena &arena) {
    hittable_list world;

    auto checker = arena.alloc<checker_texture>(arena, colour(0.2, 0.3, 0.1), colour(0.9f));
    world.add(arena.alloc<sphere>(point3(0,-1000,0), 1000, arena.alloc<lambertian>(checker)));

    hittable_list small_spheres;
    for (int a = -11; a < 11; a++) {
        for (int b = -11; b < 11; b++) {
            auto choose_mat = random_number();
            point3 centre(a + 0.9 * random_number(), 0.2, b + 0.9 * random_number());

            if ((centre - point3(4, 0.2, 0)).mag() > 0.9) {
                material* sphere_material;

                if (choose_mat < 0.8) {
                    // diffuse
                    auto albedo = colour::random() * colour::random();
                    //sphere_material = arena.alloc<lambertian>(arena, albedo);
                    sphere_material = arena.alloc<phong>(arena, albedo, 0.05, pow(10, random_number(-2, 4)));
                    //auto centre2 = centre + vec3(0, random_number(0.0, 0.5), 0);
                    //small_spheres.add(arena.alloc<moving_sphere>(
                    //    centre, centre2, 0.0, 1.0, 0.2, sphere_material));
                    small_spheres.add(arena.alloc<sphere>(centre, 0.2, sphere_material));
                } else if (choose_mat < 0.95) {
                    // metal
                    auto albedo = colour::random(0.5, 1);
                    auto fuzz = random_number(0.0, 0.3);
                    sphere_material = arena.alloc<metal>(arena, albedo, fuzz);
                    small_spheres.add(arena.alloc<sphere>(centre, 0.2, sphere_material));
                } else {
                    // glass
                    auto translucency = random_number(0.0, 0.5);
                    sphere_material = arena.alloc<dielectric>(1.5, translucency);
                    small_spheres.add(arena.alloc<sphere>(centre, 0.2, sphere_material));
                }
            }
        }
    }
    world.add(arena.alloc<bvh_accel>(arena, small_spheres, 1, 0, 1));

    auto material1 = arena.alloc<dielectric>(1.5);
    world.add(arena.alloc<sphere>(point3(0, 1, 0), 1.0, material1));

    auto material2 = arena.alloc<lambertian>(arena, colour(0.4, 0.2, 0.1));
    world.add(arena.alloc<sphere>(point3(-4, 1, 0), 1.0, material2));

    auto material3 = arena.alloc<metal>(arena, colour(0.7, 0.6, 0.5), 0.0);
    world.add(arena.alloc<sphere>(point3(4, 1, 0), 1.0, material3));

    return world;
}

hittable_list world_two_spheres(mem_arena &arena) {
    /*hittable_list objects;

    auto checker = arena.alloc<checker_texture>(arena, colour(0.2, 0.3, 0.1), colour(0.9, 0.9, 0.9));
    objects.add(arena.alloc<sphere>(point3(0, -10, 0), 10, arena.alloc<lambertian>(checker)));
    objects.add(arena.alloc<sphere>(point3(0,  10, 0), 10, arena.alloc<lambertian>(checker)));

    return objects;*/

    hittable_list objects;

    auto grey = arena.alloc<lambertian>(arena, colour(0.3f));
    objects.add(arena.alloc<xz_rect>(-5000, 5000, -5000, 5000, 0, grey));

    auto light = arena.alloc<diffuse_light>(arena, colour(100.0));

    auto c = colour(0.2, 0.4, 0.9);

    // Subsurface scattering sphere
    auto boundary = arena.alloc<sphere>(point3(-150, 100, 0), 100, arena.alloc<dielectric>(1.5));
    objects.add(boundary);
    objects.add(arena.alloc<constant_medium>(arena, boundary, 0.2, c));

    // Coloured sphere
    auto coloured_glass = arena.alloc<coloured_dielectric>(arena, c, 1.00, 1.5);
    objects.add(arena.alloc<sphere>(point3(150, 100, 0), 100, coloured_glass));

    return objects;
}

hittable_list world_perlin_spheres(mem_arena &arena) {
    hittable_list objects;

    auto per_tex = arena.alloc<noise_texture>(4);
    objects.add(arena.alloc<sphere>(point3(0, -1000, 0), 1000, arena.alloc<lambertian>(per_tex)));
    objects.add(arena.alloc<sphere>(point3(0, 2, 0), 2, arena.alloc<lambertian>(per_tex)));

    return objects;
}

hittable_list world_earth(mem_arena &arena) {
    auto earth_texture = arena.alloc<image_texture>("media/earthmap.jpg");
    auto earth_surface = arena.alloc<lambertian>(earth_texture);
    auto globe = arena.alloc<sphere>(point3(0), 2, earth_surface);

    return hittable_list(globe);
}

hittable_list world_simple_light(mem_arena &arena) {
    hittable_list objects;

    auto shiny = arena.alloc<metal>(arena, colour(0.8, 0.85, 0.8), 0.0);
    auto per_tex = arena.alloc<noise_texture>(4);
    auto diff_light = arena.alloc<diffuse_light>(arena, colour(1));
    
    auto shiny_sphere = arena.alloc<sphere>(point3(0, -1000, 0), 1000, shiny);
    auto perlin_sphere = arena.alloc<sphere>(point3(0, 2, 0), 2, arena.alloc<lambertian>(per_tex));
    auto light_sphere = arena.alloc<sphere>(point3(0, 1000, 0), 990, diff_light);
    light_sphere->set_sampling_target(true);
    
    auto rotated_box = arena.alloc<translate>(
        arena.alloc<rotate_z>(
            arena.alloc<box>(arena, point3(0), point3(2, 2, 1), shiny),
            22.5
        ), vec3(3, 1, -2.5)
    );

    objects.add(shiny_sphere);
    objects.add(perlin_sphere);
    objects.add(light_sphere);
    objects.add(rotated_box);
    //objects.add(arena.alloc<xy_rect>(3, 5, 1, 3, -2, diff_light));

    return objects;
}

hittable_list world_cornell_box(mem_arena &arena, bool big_light) {
    hittable_list objects;

    auto red = arena.alloc<lambertian>(arena, colour(0.65, 0.05, 0.05));
    // auto green = arena.alloc<lambertian>(arena, colour(0.12, 0.45, 0.15));
    auto green = arena.alloc<lambertian>(arena.alloc<checker_texture>(
        arena, colour(0.12, 0.45, 0.15), 0.6 * colour(0.12, 0.45, 0.15), 0.03));
    auto blue = arena.alloc<lambertian>(arena, colour(0.05, 0.05, 0.65));
    auto white = arena.alloc<lambertian>(arena, colour(0.73, 0.73, 0.73));
    auto light = arena.alloc<diffuse_light>(arena, colour(15, 15, 15) * (big_light ? 0.5 : 1.0));

    objects.add(arena.alloc<yz_rect>(0, 555, 0, 555, 555, blue));
    objects.add(arena.alloc<yz_rect>(0, 555, 0, 555, 0, red));
    objects.add(arena.alloc<xz_rect>(0, 555, 0, 555, 0, green));
    objects.add(arena.alloc<xz_rect>(0, 555, 0, 555, 555, white));
    objects.add(arena.alloc<xy_rect>(0, 555, 0, 555, 555, white));
    objects.add(arena.alloc<one_sided>(arena.alloc<xy_rect>(0, 555, 0, 555, 0, white)));

    if (big_light) {
        auto light_obj = arena.alloc<flip_face>(arena.alloc<xz_rect>(113, 443, 127, 432, 554, light));
        light_obj->set_sampling_target(true);
        objects.add(light_obj);
    } else {
        auto light_obj = arena.alloc<flip_face>(arena.alloc<xz_rect>(213, 343, 227, 332, 554, light));
        light_obj->set_sampling_target(true);
        objects.add(light_obj);
    }

    auto box1 = arena.alloc<translate>(
        arena.alloc<rotate_y>(
            arena.alloc<box>(arena, point3(0, 0, 0), point3(165, 330, 165), white),
            15),
        vec3(265, 0, 295)
    );
    objects.add(box1);

    auto box2 = arena.alloc<translate>(
        arena.alloc<rotate_y>(
            arena.alloc<box>(arena, point3(0, 0, 0), point3(165, 165, 165), white),
            -18),
        vec3(130, 0, 65)
    );
    objects.add(box2);

    return objects;
}

hittable_list world_cornell_box_smoke(mem_arena &arena, bool big_light) {
    hittable_list objects;

    auto red = arena.alloc<lambertian>(arena, colour(0.65, 0.05, 0.05));
    auto white = arena.alloc<lambertian>(arena, colour(0.73, 0.73, 0.73));
    auto green = arena.alloc<lambertian>(arena, colour(0.12, 0.45, 0.15));
    auto light = arena.alloc<diffuse_light>(arena, colour(15, 15, 15) * (big_light ? 0.5 : 1.0));

    objects.add(arena.alloc<yz_rect>(0, 555, 0, 555, 555, green));
    objects.add(arena.alloc<yz_rect>(0, 555, 0, 555, 0, red));
    objects.add(arena.alloc<xz_rect>(0, 555, 0, 555, 0, white));
    objects.add(arena.alloc<xz_rect>(0, 555, 0, 555, 555, white));
    objects.add(arena.alloc<xy_rect>(0, 555, 0, 555, 555, white));

    if (big_light) {
        auto light_obj = arena.alloc<flip_face>(arena.alloc<xz_rect>(113, 443, 127, 432, 554, light));
        light_obj->set_sampling_target(true);
        objects.add(light_obj);
    } else {
        auto light_obj = arena.alloc<flip_face>(arena.alloc<xz_rect>(213, 343, 227, 332, 554, light));
        light_obj->set_sampling_target(true);
        objects.add(light_obj);
    }

    auto box1 = arena.alloc<translate>(
        arena.alloc<rotate_y>(
            arena.alloc<box>(arena, point3(0, 0, 0), point3(165, 330, 165), white),
            15),
        vec3(265, 0, 295)
    );

    auto box2 = arena.alloc<translate>(
        arena.alloc<rotate_y>(
            arena.alloc<box>(arena, point3(0, 0, 0), point3(165, 165, 165), white),
            -18),
        vec3(130, 0, 65)
    );
    
    objects.add(arena.alloc<constant_medium>(arena, box1, 0.01, colour(0.0)));
    objects.add(arena.alloc<constant_medium>(arena, box2, 0.01, colour(1.0)));

    return objects;
}

hittable_list world_cornell_box_specular(mem_arena &arena, bool big_light) {
    hittable_list objects;

    auto red   = arena.alloc<lambertian>(arena, colour(0.65, 0.05, 0.05));
    auto white = arena.alloc<lambertian>(arena, colour(0.73, 0.73, 0.73));
    auto green = arena.alloc<lambertian>(arena, colour(0.12, 0.45, 0.15));
    auto light = arena.alloc<diffuse_light>(arena, colour(15, 15, 15) * (big_light ? 0.5 : 1.0));

    objects.add(arena.alloc<yz_rect>(0, 555, 0, 555, 555, green));
    objects.add(arena.alloc<yz_rect>(0, 555, 0, 555, 0, red));
    objects.add(arena.alloc<xz_rect>(0, 555, 0, 555, 0, white));
    objects.add(arena.alloc<xz_rect>(0, 555, 0, 555, 555, white));
    objects.add(arena.alloc<xy_rect>(0, 555, 0, 555, 555, white));

    if (big_light) {
        auto light_obj = arena.alloc<flip_face>(arena.alloc<xz_rect>(113, 443, 127, 432, 554, light));
        light_obj->set_sampling_target(true);
        objects.add(light_obj);
    } else {
        auto light_obj = arena.alloc<flip_face>(arena.alloc<xz_rect>(213, 343, 227, 332, 554, light));
        light_obj->set_sampling_target(true);
        objects.add(light_obj);
    }

    auto shiny = arena.alloc<metal>(arena, colour(0.8, 0.85, 0.88), 0.3);
    auto glass = arena.alloc<dielectric>(1.5);
    auto phongy = arena.alloc<phong>(arena, colour(0.2, 0.4, 0.9), 0.05, 100);

    // auto box1 = arena.alloc<translate>(
    //     arena.alloc<rotate_y>(
    //         arena.alloc<box>(point3(0, 0, 0), point3(165, 330, 165), shiny),
    //         15),
    //     vec3(265, 0, 295)
    // );
    // box1->set_sampling_target(true);
    // objects.add(box1);

    auto shiny_sphere = arena.alloc<sphere>(point3(310, 135, 340), 135, shiny);
    //shiny_sphere->set_sampling_target(true);
    objects.add(shiny_sphere);

    auto glass_sphere = arena.alloc<sphere>(point3(420, 90, 75), 90, glass);
    //glass_sphere->set_sampling_target(true);
    objects.add(glass_sphere);

    auto glass_sphere2 = arena.alloc<sphere>(point3(420, 90, 75), -75, glass);
    //glass_sphere2->set_sampling_target(true);
    objects.add(glass_sphere2);

    auto phong_sphere = arena.alloc<sphere>(point3(135, 90, 190), 90, phongy);
    objects.add(phong_sphere);

    // auto ss_scat_bdy = arena.alloc<sphere>(point3(420, 90, 75), 90, glass);
    // objects.add(ss_scat_bdy);
    // objects.add(arena.alloc<constant_medium>(arena, ss_scat_bdy, 0.1, colour(0.2, 0.4, 0.9)));
    
    return objects;
}

hittable_list world_cornell_box_bunny(mem_arena &arena, bool big_light) {
    hittable_list objects;

    auto red = arena.alloc<lambertian>(arena, colour(0.65, 0.05, 0.05));
    // auto green = arena.alloc<lambertian>(arena, colour(0.12, 0.45, 0.15));
    auto green = arena.alloc<lambertian>(arena.alloc<checker_texture>(
        arena, colour(0.12, 0.45, 0.15), 0.6 * colour(0.12, 0.45, 0.15), 0.03));
    auto blue = arena.alloc<lambertian>(arena, colour(0.05, 0.05, 0.65));
    auto white = arena.alloc<lambertian>(arena, colour(0.73, 0.73, 0.73));
    auto light = arena.alloc<diffuse_light>(arena, colour(15, 15, 15) * (big_light ? 0.5 : 1.0));

    objects.add(arena.alloc<yz_rect>(0, 555, 0, 555, 555, blue));
    objects.add(arena.alloc<yz_rect>(0, 555, 0, 555, 0, red));
    objects.add(arena.alloc<xz_rect>(0, 555, 0, 555, 0, green));
    objects.add(arena.alloc<xz_rect>(0, 555, 0, 555, 555, white));
    objects.add(arena.alloc<xy_rect>(0, 555, 0, 555, 555, white));
    objects.add(arena.alloc<one_sided>(arena.alloc<xy_rect>(0, 555, 0, 555, 0, white)));

    if (big_light) {
        auto light_obj = arena.alloc<flip_face>(arena.alloc<xz_rect>(113, 443, 127, 432, 554, light));
        light_obj->set_sampling_target(true);
        objects.add(light_obj);
    } else {
        auto light_obj = arena.alloc<flip_face>(arena.alloc<xz_rect>(213, 343, 227, 332, 554, light));
        light_obj->set_sampling_target(true);
        objects.add(light_obj);
    }

    auto phong_mat = arena.alloc<phong>(arena, colour(1.0, 0.918, 0.0), 0.05, 100.0);
    // auto gold_mat = arena.alloc<metal>(arena, colour(1.0, 0.71, 0.29), 0.5);
    // auto glass_mat = arena.alloc<dielectric>(1.5);
    auto bunny = arena.alloc<translate>(
        arena.alloc<rotate_y>(
            arena.alloc<mesh>(arena, point3(0, 0, 0), 400.0, "media/bunny.obj", phong_mat),
            // arena.alloc<mesh>(arena, point3(0, 0, 0), 400.0, "media/bunny.obj", gold_mat),
            // arena.alloc<mesh>(arena, point3(0, 0, 0), 400.0, "media/bunny.obj", glass_mat),
            180),
        vec3(277, -10, 277)
    );
    // auto medium = arena.alloc<constant_medium>(arena, bunny, 0.02, colour(1.0, 0.918, 0.0));

    objects.add(bunny);
    // objects.add(medium);

    return objects;
}

hittable_list world_cornell_box_dragon(mem_arena &arena, bool big_light) {
    hittable_list objects;

    auto red = arena.alloc<lambertian>(arena, colour(0.65, 0.05, 0.05));
    auto white = arena.alloc<lambertian>(arena, colour(0.73, 0.73, 0.73));
    auto green = arena.alloc<lambertian>(arena, colour(0.12, 0.45, 0.15));
    auto light = arena.alloc<diffuse_light>(arena, colour(15, 15, 15) * (big_light ? 0.5 : 1.0));

    objects.add(arena.alloc<yz_rect>(0, 555, 0, 555, 555, green));
    objects.add(arena.alloc<yz_rect>(0, 555, 0, 555, 0, red));
    objects.add(arena.alloc<xz_rect>(0, 555, 0, 555, 0, white));
    objects.add(arena.alloc<xz_rect>(0, 555, 0, 555, 555, white));
    objects.add(arena.alloc<xy_rect>(0, 555, 0, 555, 555, white));

    if (big_light) {
        auto light_obj = arena.alloc<flip_face>(arena.alloc<xz_rect>(113, 443, 127, 432, 554, light));
        light_obj->set_sampling_target(true);
        objects.add(light_obj);
    } else {
        auto light_obj = arena.alloc<flip_face>(arena.alloc<xz_rect>(213, 343, 227, 332, 554, light));
        light_obj->set_sampling_target(true);
        objects.add(light_obj);
    }

    auto phong_mat = arena.alloc<phong>(arena, colour(0.2, 0.4, 0.9), 0.05, 100);
    // auto glass_mat = arena.alloc<dielectric>(1.5);
    auto dragon = arena.alloc<translate>(
        arena.alloc<rotate_y>(
            arena.alloc<mesh>(arena, point3(0, 0, 0), 500.0, "media/dragon.obj", phong_mat),
            // arena.alloc<mesh>(arena, point3(0, 0, 0), 500.0, "media/dragon.obj", glass_mat),
            205),
        vec3(277, 0, 277)
    );
    // auto medium = arena.alloc<constant_medium>(arena, dragon, 0.2, colour(0.2, 0.4, 0.9));

    objects.add(dragon);
    //objects.add(medium);

    return objects;
}

hittable_list world_cornell_box_glass(mem_arena& arena, bool big_light) {
    hittable_list objects;

    auto red = arena.alloc<lambertian>(arena, colour(0.65, 0.05, 0.05));
    auto white = arena.alloc<lambertian>(arena, colour(0.73, 0.73, 0.73));
    auto green = arena.alloc<lambertian>(arena, colour(0.12, 0.45, 0.15));
    auto light = arena.alloc<diffuse_light>(arena, colour(15, 15, 15) * (big_light ? 0.5 : 1.0));
    //auto glass = arena.alloc<dielectric>(1.5);
    auto blue_glass = arena.alloc<coloured_dielectric>(arena, colour(0.2, 0.4, 0.9), 0.999, 1.5, 0.0);

    objects.add(arena.alloc<yz_rect>(0, 555, 0, 555, 555, green));
    objects.add(arena.alloc<yz_rect>(0, 555, 0, 555, 0, red));
    objects.add(arena.alloc<xz_rect>(0, 555, 0, 555, 0, white));
    objects.add(arena.alloc<xz_rect>(0, 555, 0, 555, 555, white));
    objects.add(arena.alloc<xy_rect>(0, 555, 0, 555, 555, white));
    objects.add(arena.alloc<one_sided>(arena.alloc<xy_rect>(0, 555, 0, 555, 0, white)));

    if (big_light) {
        auto light_obj = arena.alloc<xz_rect>(113, 443, 127, 432, 1, light);
        light_obj->set_sampling_target(true);
        objects.add(light_obj);
    }
    else {
        auto light_obj = arena.alloc<xz_rect>(213, 343, 227, 332, 1, light);
        light_obj->set_sampling_target(true);
        objects.add(light_obj);
    }

    constexpr int cube_depth = 5;
    constexpr double box_size = 50;
    constexpr double box_spacing = 50;
    constexpr double box_offset = 45;

    hittable_list cubes;
    for (int z = 0; z < cube_depth; z++) {
        for (int y = 0; y < cube_depth; y++) {
            for (int x = 0; x < cube_depth; x++) {
                const vec3 loc(
                    box_offset + ((box_size + box_spacing) * x),
                    box_offset + ((box_size + box_spacing) * y),
                    -42.5 + ((box_size + box_spacing) * z)
                );
                const double rot = (23 * x) + (13 * y) + (3 * z);

                cubes.add(arena.alloc<translate>(
                    arena.alloc<rotate_y>(
                        arena.alloc<box>(arena, vec3(0), vec3(box_size), blue_glass),
                        rot),
                    loc
                ));
            }
        }
    }

    objects.add(arena.alloc<bvh_accel>(arena, cubes, 1, 0.0, 1.0));
    return objects;
}

hittable_list world_cornell_box_glass_spheres(mem_arena& arena, bool big_light) {
    hittable_list objects;

    auto red = arena.alloc<lambertian>(arena, colour(0.65, 0.05, 0.05));
    auto white = arena.alloc<lambertian>(arena, colour(0.73, 0.73, 0.73));
    auto green = arena.alloc<lambertian>(arena, colour(0.12, 0.45, 0.15));
    auto light = arena.alloc<diffuse_light>(arena, colour(15, 15, 15) * (big_light ? 0.5 : 8.0));
    
    objects.add(arena.alloc<yz_rect>(0, 555, 0, 555, 555, green));
    objects.add(arena.alloc<yz_rect>(0, 555, 0, 555, 0, red));
    objects.add(arena.alloc<xz_rect>(0, 555, 0, 555, 0, white));
    objects.add(arena.alloc<xz_rect>(0, 555, 0, 555, 555, white));
    objects.add(arena.alloc<xy_rect>(0, 555, 0, 555, 555, white));
    objects.add(arena.alloc<one_sided>(arena.alloc<xy_rect>(0, 555, 0, 555, 0, white)));

    if (big_light) {
        auto light_obj = arena.alloc<flip_face>(arena.alloc<xz_rect>(113, 443, 127, 432, 554, light));
        light_obj->set_sampling_target(true);
        objects.add(light_obj);
    }
    else {
        // auto light_obj = arena.alloc<flip_face>(arena.alloc<xz_rect>(213, 343, 227, 332, 554, light));
        auto light_obj = arena.alloc<flip_face>(arena.alloc<xz_rect>(250, 310, 250, 310, 554, light));
        light_obj->set_sampling_target(true);
        objects.add(light_obj);
    }

    hittable_list spheres;

    constexpr int grid = 3;
    constexpr int size = 555 / grid;
    for (int i = 0; i < grid; i++) {
        for (int j = 0; j < grid; j++) {
            for (int k = 0; k < grid; k++) {
                int r = random_number(0.167 * size, 0.333 * size);
                int x = random_number(i * size + 1.1 * r, (i + 1) * size - 1.1 * r);
                int y = random_number(j * size + 1.1 * r, (j + 1) * size - 1.1 * r);
                int z = random_number(k * size + 1.1 * r, (k + 1) * size - 1.1 * r);

                // auto translucency = random_number(0.0, 0.25);
                auto glass = arena.alloc<dielectric>(1.5);

                auto s = arena.alloc<sphere>(point3(x, y, z), r, glass);
                s->set_sampling_target(true);
                spheres.add(s);
            }
        }
    }

    auto sphere_bvh = arena.alloc<bvh_accel>(arena, spheres, 1, 0.0, 1.0);
    sphere_bvh->set_sampling_target(true);
    objects.add(sphere_bvh);

    // Mist/fog medium
    auto boundary = arena.alloc<box>(arena, point3(0), point3(555), arena.alloc<dielectric>(1.5));
    objects.add(arena.alloc<constant_medium>(arena, boundary, 0.001, colour(1)));
    
    return objects;
}

hittable_list world_caustics(mem_arena& arena) {
    hittable_list objects;

    auto light = arena.alloc<diffuse_light>(arena, colour(100));
    auto light_obj = arena.alloc<yz_rect>(0, 100, -10, 10, -300, light);
    light_obj->set_sampling_target(true);
    objects.add(light_obj);
    
    // for (int i = 0; i < 3; i++) {
    //         auto light = arena.alloc<diffuse_light>(arena, 100 * colour(i == 0, i == 1, i == 2));
    //         auto light_obj = arena.alloc<yz_rect>(0, 100, -15 + i*10, -15+(i+1)*10, -300, light);
    //         light_obj->set_sampling_target(true);
    //         objects.add(light_obj);
    // }
    
    auto glass = arena.alloc<dielectric>(1.5);
    auto sphere_obj = arena.alloc<sphere>(point3(-200, 50, 0), 50, glass);
    sphere_obj->set_sampling_target(true);
    objects.add(sphere_obj);
    
    // Mist/fog medium
    auto boundary = arena.alloc<box>(arena, point3(-500), point3(500), arena.alloc<dielectric>(1.5));
    objects.add(arena.alloc<constant_medium>(arena, boundary, 0.001, colour(1)));
    
    return objects;
}

hittable_list world_phong(mem_arena &arena) {
    hittable_list objects;

    auto grey = arena.alloc<lambertian>(arena, colour(0.3f));
    objects.add(arena.alloc<xz_rect>(-5000, 5000, -5000, 5000, 0, grey));

    auto light = arena.alloc<diffuse_light>(arena, colour(100.0));

    hittable_list spheres;
    int num_spheres = 10;
    double min_shine = -4, max_shine = 4;
    double spacing = 1000.0 / (num_spheres + 1);
    double xpos = -500 + 1.5 * spacing;
    double radius = 0.4 * spacing;

    auto ball_colour = colour(1.0, 0.0, 0.0);
    for (int i = 0; i < num_spheres; i++) {
        auto shininess = pow(10, min_shine + i * ((max_shine - min_shine) / (num_spheres - 1)));
        //auto fspec = (0.05 * i) / num_spheres + (0.001 * (num_spheres - i)) / num_spheres;
        //auto phong_mat = arena.alloc<phong>(arena, ball_colour, fspec, shininess);
        auto phong_mat = arena.alloc<ashikhmin_shirley>(arena, ball_colour, colour(1.0), 0.05, shininess, shininess);
        
        auto ball = arena.alloc<sphere>(point3(xpos, radius, 0), radius, phong_mat);
        spheres.add(ball);
        xpos += spacing;
    }

    auto red = arena.alloc<lambertian>(arena, ball_colour);
    auto ball = arena.alloc<sphere>(point3(-500 + 0.5 * spacing, radius, 0), radius, red);
    spheres.add(ball);

    auto bvh = arena.alloc<bvh_accel>(arena, spheres, 1, 0, 1);
    objects.add(bvh);

    auto light_obj = arena.alloc<flip_face>(arena.alloc<xz_rect>(-100, 100, -100, 100, 1000, light));
    light_obj->set_sampling_target(true);
    objects.add(light_obj);
    
    return objects;
}

hittable_list world_book2(mem_arena &arena, bool rgb_light) {
    // Grid of random height green boxes on the ground
    hittable_list boxes1;
    auto ground = arena.alloc<lambertian>(arena, colour(0.48, 0.83, 0.53));

    const int boxes_per_side = 20;
    auto w = 2000 / boxes_per_side;
    for (int i = 0; i < boxes_per_side; i++) {
        for (int j = 0; j < boxes_per_side; j++) {
            auto x0 = -1000.0 + i * w;
            auto z0 = -1000.0 + j * w;
            auto y0 = 0.0;

            auto x1 = x0 + w;
            auto y1 = random_number(0.0, 101.0);
            auto z1 = z0 + w;

            boxes1.add(arena.alloc<box>(arena, point3(x0, y0, z0), point3(x1, y1, z1), ground));
        }
    }

    hittable_list objects;

    objects.add(arena.alloc<bvh_accel>(arena, boxes1, 4, 0.0, 1.0));

    if (rgb_light) {
        // R, G, B lights
    auto red_light = arena.alloc<flip_face>(
        arena.alloc<xz_rect>(123, 223, 147, 412, 554, arena.alloc<diffuse_light>(arena, colour(7, 0, 0)))
    );
    auto green_light = arena.alloc<flip_face>(
        arena.alloc<xz_rect>(223, 323, 147, 412, 554, arena.alloc<diffuse_light>(arena, colour(0, 7, 0)))
    );
    auto blue_light = arena.alloc<flip_face>(
        arena.alloc<xz_rect>(323, 423, 147, 412, 554, arena.alloc<diffuse_light>(arena, colour(0, 0, 7)))
    );
    red_light->set_sampling_target(true);
    green_light->set_sampling_target(true);
    blue_light->set_sampling_target(true);
    objects.add(red_light);
    objects.add(green_light);
    objects.add(blue_light);
    } else {
        // Light
        auto light = arena.alloc<flip_face>(
            arena.alloc<xz_rect>(123, 423, 147, 412, 554, arena.alloc<diffuse_light>(arena, colour(7, 5.41, 3.93)))
        );
        light->set_sampling_target(true);
        objects.add(light);
    }

    // Blurred orange sphere
    auto centre1 = point3(400, 400, 200);
    auto centre2 = centre1 + vec3(30, 0, 0);
    auto moving_sphere_material = arena.alloc<lambertian>(arena, colour(0.7, 0.3, 0.1));
    objects.add(arena.alloc<moving_sphere>(centre1, centre2, 0, 1, 50, moving_sphere_material));

    // Glass sphere
    objects.add(arena.alloc<sphere>(point3(260, 150, 45), 50, arena.alloc<dielectric>(1.5)));
    
    // Silver metal sphere
    objects.add(arena.alloc<sphere>(
        point3(0, 150, 145), 50, arena.alloc<metal>(arena, colour(0.8, 0.8, 0.9), 0.1)
    ));

    //  Sub-surface scattering sphere
    auto boundary = arena.alloc<sphere>(point3(360, 150, 145), 70, arena.alloc<dielectric>(1.5));
    objects.add(boundary);
    objects.add(arena.alloc<constant_medium>(arena, boundary, 0.2, colour(0.2, 0.4, 0.9)));
    
    // Mist/fog medium
    boundary = arena.alloc<sphere>(point3(0, 0, 0), 5000, arena.alloc<dielectric>(1.5));
    objects.add(arena.alloc<constant_medium>(arena, boundary, 0.0001, colour(1)));

    // Earth texture sphere
    auto emat = arena.alloc<lambertian>(arena.alloc<image_texture>("media/earthmap.jpg"));
    //auto emat = arena.alloc<phong>(arena.alloc<image_texture>("media/earthmap.jpg"), 0.05, 50);
    objects.add(arena.alloc<sphere>(point3(400,200,400), 100, emat));
    
    // Perlin noise textured sphere
    auto pertex = arena.alloc<noise_texture>(0.1);
    //auto permat = arena.alloc<lambertian>(pertex);
    auto permat = arena.alloc<phong>(pertex, 0.05, 50);
    objects.add(arena.alloc<sphere>(point3(220, 280, 300), 80, permat));

    // Cube full of random white spheres
    hittable_list boxes2;
    auto white = arena.alloc<lambertian>(arena, colour(0.73, 0.73, 0.73));
    int ns = 1000;
    for (int j = 0; j < ns; j++) {
        boxes2.add(arena.alloc<sphere>(point3::random(0.0, 165.0), 10, white));
    }

    objects.add(arena.alloc<translate>(
        arena.alloc<rotate_y>(
            arena.alloc<bvh_accel>(arena, boxes2, 4, 0.0, 1.0), 15),
        vec3(-100, 270, 395)
    ));

    return objects;
}

hittable_list world_tokyo(mem_arena& arena) {
    auto glass_mat = arena.alloc<dielectric>(1.5);
    auto metal_mat = arena.alloc<metal>(arena, colour(1.0, 0.86, 0.57), 0.0);
    auto marble_tex = arena.alloc<noise_texture>(2.0);
    auto marble_mat = arena.alloc<phong>(marble_tex, 0.001, 10);
    auto dragon_mat = arena.alloc<ashikhmin_shirley>(arena, 2 * colour(0.104, 0.0629, 0.0225), 2 * colour(0.0558, 0.037, 0.015), 0.999, 4.53e3, 4.53e3);
    auto light = arena.alloc<diffuse_light>(
        arena.alloc<blend_texture>(
            arena.alloc<solid_colour>(6.0f * colour(0.38, 0.08, 0.54)),
            arena.alloc<image_texture>("media/box_gradient.png"),
            blend_mode::MULTIPLY
        ), true
    );

    auto light_box = arena.alloc<translate>(
        arena.alloc<rotate_y>(
            arena.alloc<box>(arena, point3(-0.33), point3(0.33), light),
            25),
        point3(-1.35, -0.65, -1.1)
    );
    light_box->set_sampling_target(true);
    
    auto metal_sphere = arena.alloc<sphere>(point3(0.5, 0, 1), 0.98, metal_mat);

    auto glass_sphere = arena.alloc<sphere>(point3(0.5, 0, -1), 0.98, glass_mat);
    glass_sphere->set_sampling_target(true);
    
    auto ground = arena.alloc<xz_rect>(-5.0, 5.0, -5.0, 5.0, -0.98, marble_mat);
    //auto ground = arena.alloc<xz_rect>(-5.0, 5.0, -5.0, 5.0, -0.98, arena.alloc<lambertian>(arena, colour(0.0001f)));
    //auto ground = arena.alloc<xz_rect>(-5.0, 5.0, -5.0, 5.0, -0.98, metal_mat);

    auto dragon = arena.alloc<translate>(
        arena.alloc<rotate_y>(
            arena.alloc<mesh>(arena, vec3(0), 2.0, "media/dragon.obj", dragon_mat),
            -130.0),
        vec3(-1, -0.98, 0.35)
    );

    hittable_list objects;
    objects.add(light_box);
    objects.add(metal_sphere);
    objects.add(glass_sphere);
    objects.add(dragon);
    objects.add(ground);

    return objects;
}

hittable_list world_studio(mem_arena& arena) {
    auto marble_tex = arena.alloc<noise_texture>(2.0);
    auto marble_mat = arena.alloc<phong>(marble_tex, 0.01, 10);

    // auto teapot_mat = arena.alloc<coloured_dielectric>(arena, colour(0.0, 1.0, 1.0), 10.3, 1.5);
    auto teapot_mat = arena.alloc<ashikhmin_shirley>(arena, colour(0.0), colour(0.983, 0.991, 0.995), 1.0, 1000, 10);
    // auto bunny_mat = arena.alloc<phong>(arena, colour(0.85, 0.0, 0.1), 0.1, 100);
    auto bunny_mat = arena.alloc<ashikhmin_shirley>(arena, colour(0.95, 0.0, 0.1), colour(1.0), 0.1, 100, 100);
    // auto sphere_mat = arena.alloc<metal>(arena, colour(0.542, 0.497, 0.449), 0.0);
    auto sphere_mat = arena.alloc<coloured_dielectric>(arena, colour(0.0, 0.45, 0.7), 1.1 , 1.5);

    hittable_list objects;

    auto teapot = arena.alloc<translate>(
        arena.alloc<rotate_y>(
            arena.alloc<mesh>(arena, vec3(0), 2.0, "media/teapot_highpoly.obj", teapot_mat),
            60.0),
        vec3(-2.1, -0.98, -0.2)
    );

    auto bunny = arena.alloc<translate>(
        arena.alloc<rotate_y>(
            arena.alloc<mesh>(arena, vec3(0), 2.0, "media/bunny.obj", bunny_mat),
            -70.0),
        vec3(0.0, -0.98, -1.3)
    );

    auto dragon = arena.alloc<translate>(
        arena.alloc<rotate_y>(
            arena.alloc<mesh>(arena, vec3(0), 2.5, "media/dragon.obj", sphere_mat),
            -140.0),
        vec3(-0.2, -0.98, 0.7)
    );
 
    //auto sfere = arena.alloc<sphere>(vec3(0.0, 0.0, 0.98), 0.98, sphere_mat);
    auto ground = arena.alloc<xz_rect>(-5.0, 5.0, -5.0, 5.0, -0.98, marble_mat);

    objects.add(teapot);
    objects.add(bunny);
    //objects.add(sfere);
    objects.add(dragon);
    objects.add(ground);

    return objects;
}

hittable_list world_park(mem_arena& arena) {
    auto marble_tex = arena.alloc<noise_texture>(2.0);
    auto marble_mat = arena.alloc<phong>(marble_tex, 0.01, 10);

    auto teapot_mat = arena.alloc<phong>(arena, colour(0.12, 0.45, 0.15), 0.08, 50);
    // auto dragon_mat = arena.alloc<metal>(arena, colour(0.983, 0.991, 0.995), 0.45);
    // auto dragon_mat = arena.alloc<phong>(arena, colour(0.6, 0.2, 0.1), 0.1, 100);
    auto dragon_mat = arena.alloc<ashikhmin_shirley>(arena, colour(0.0), colour(0.983, 0.991, 0.995), 1.0, 1000, 10);
    auto bunny_mat = arena.alloc<phong>(arena, colour(0.6, 0.2, 0.1), 0.1, 100);
    //auto sphere_mat = arena.alloc<metal>(arena, colour(0.542, 0.497, 0.449), 0.0);

    hittable_list objects;

    auto teapot = arena.alloc<translate>(
        arena.alloc<rotate_y>(
            arena.alloc<mesh>(arena, vec3(0), 2.0, "media/teapot_highpoly.obj", teapot_mat),
            60.0),
        vec3(-2.1, -0.98, 0.0)
    );

    auto bunny = arena.alloc<translate>(
        arena.alloc<rotate_y>(
            arena.alloc<mesh>(arena, vec3(0), 2.0, "media/bunny.obj", bunny_mat),
            -70.0),
        vec3(0.0, -0.98, -1.3)
    );

    auto dragon = arena.alloc<translate>(
        arena.alloc<rotate_y>(
            arena.alloc<mesh>(arena, vec3(0), 2.5, "media/dragon.obj", dragon_mat),
            -140.0),
        vec3(-0.2, -0.98, 0.7)
    );
 
    //auto sfere = arena.alloc<sphere>(vec3(0.0, 0.0, 0.98), 0.98, sphere_mat);
    auto ground = arena.alloc<xz_rect>(-5.0, 5.0, -5.0, 5.0, -0.98, marble_mat);

    objects.add(teapot);
     objects.add(dragon);
    objects.add(bunny);
    //objects.add(sfere);
    objects.add(ground);

    return objects;
}

hittable_list world_xmas(mem_arena& arena) {
    hittable_list objects;

    auto white = arena.alloc<lambertian>(arena, colour(0.73, 0.73, 0.73));
    auto floor = arena.alloc<xz_rect>(-5000, 5000, -5000, 5000, 0, white);
    objects.add(floor);

    /*auto red = arena.alloc<lambertian>(arena, colour(1.0, 0.0, 0.0));
    auto green = arena.alloc<lambertian>(arena, colour(0.0, 1.0, 0.0));
    auto blue = arena.alloc<lambertian>(arena, colour(0.0, 0.0, 1.0));

    objects.add(arena.alloc<box>(arena, point3(200 - 50, -5, -5), point3(200 + 50, 5, 5), red));
    objects.add(arena.alloc<box>(arena, point3(200 - 5, -50, -5), point3(200 + 5, 50, 5), green));
    objects.add(arena.alloc<box>(arena, point3(200 - 5, -5, -50), point3(200 + 5, 5, 50), blue));*/

    // auto tree_mat = arena.alloc<lambertian>(arena, colour(0.2, 0.8, 0.2));

    aabb tree_box;
    auto tree = arena.alloc<mesh>(arena, point3(0, 0, -400), 300.0, 90.0, "media/12150_Christmas_Tree_V2_L2.obj", nullptr);
    tree->bounding_box(0, 1, tree_box);
    objects.add(tree);

    //objects.add(arena.alloc<box>(arena, tree_box.min(), tree_box.max(), red));

    auto y_range = tree_box.max().y() - tree_box.min().y() - 10;
    /*auto tree_centre = 0.5 * point2(
        tree_box.max().x() + tree_box.min().x(),
        tree_box.max().z() + tree_box.min().z()
    );*/
    // auto tree_radius = 0.25 * (
    //     (tree_box.max().x() - tree_box.min().x()) +
    //     (tree_box.max().z() - tree_box.min().z())
    // );

    point2 tree_centre{ 0, -400 };
    //double tree_radius = 200;

    auto light_intensity = 20.0;
    material* light_colours[] = {
        arena.alloc<diffuse_light>(arena, light_intensity * colour(1.0, 0.0, 0.0)),
        arena.alloc<diffuse_light>(arena, light_intensity * colour(0.0, 1.0, 0.0)),
        arena.alloc<diffuse_light>(arena, light_intensity * colour(1.0, 1.0, 0.0)),
        arena.alloc<diffuse_light>(arena, light_intensity * colour(0.0, 0.0, 1.0)),
        arena.alloc<diffuse_light>(arena, light_intensity * colour(1.0))
    };
    int n_colours = sizeof(light_colours) / sizeof(material*);

    int num_lights = 400;
    hittable_list fairy_lights;
    for (int i = 0; i < num_lights; i++) {
        while (true) {
            auto y = 50 + y_range - sqrt(random_number(50.0 * 50.0, y_range * y_range));
            auto xz = tree_centre + random_unit_vector<point2>() * 100;
            point3 origin{ xz.x(), y, xz.y() };
            vec3 dir{ -xz.x(), 0.0, -xz.y() };
            ray r{ origin, unit_vector(dir) };

            hit_record rec;
            if (!tree->hit(r, epsilon, infinity, rec))
                continue;
            auto p = r.at(rec.t - epsilon);
            auto c = random_number(0, n_colours - 1);
            fairy_lights.add(arena.alloc<sphere>(p, 0.75, light_colours[c]));
            break;
        }
    }
    auto lights_bvh = arena.alloc<bvh_accel>(arena, fairy_lights, 1, 0.0, 1.0);
    //lights_bvh->set_sampling_target(true);
    objects.add(lights_bvh);

    // Mist/fog medium
    auto boundary = arena.alloc<sphere>(point3(0, 0, 0), 5000, arena.alloc<dielectric>(1.5));
    objects.add(arena.alloc<constant_medium>(arena, boundary, 0.00005, colour(1)));

    auto moonlight = arena.alloc<diffuse_light>(arena, 200 * colour(0.75, 0.75, 0.9));
    auto moon = arena.alloc<sphere>(point3(1000, 10000, 1000), 100, moonlight);
    moon->set_sampling_target(true);
    objects.add(moon);

    auto fill_light = arena.alloc<diffuse_light>(arena, 1.2 * colour(0.9, 0.9, 0.75));
    auto fill = arena.alloc<xy_rect>(-1000, 1000, 0, 1000, -1000, fill_light);
    fill->set_sampling_target(true);
    objects.add(fill);

    return objects;
}

camera cam_book1(size_t image_width, size_t image_height) {
    point3 look_from(13, 2, 3);
    point3 look_at(0, 0, 0);
    vec3 vup(0, 1, 0);
    auto vfov = 20.0;
    auto aperture = 0.1;
    auto dist_to_focus = 10.0;

    return camera(
        look_from, look_at, vup, vfov, aperture, dist_to_focus, image_width, image_height
    );
}

camera cam_two_spheres(size_t image_width, size_t image_height) {
    point3 look_from(13, 2, 3);
    point3 look_at(0, 0, 0);
    vec3 vup(0, 1, 0);
    auto vfov = 20.0;
    auto aperture = 0.0;
    auto dist_to_focus = 10.0;

    return camera(
        look_from, look_at, vup, vfov, aperture, dist_to_focus, image_width, image_height
    );
}

camera cam_earth(size_t image_width, size_t image_height) {
    point3 look_from(13, 2, 3);
    point3 look_at(0, 0, 0);
    vec3 vup(0, 1, 0);
    auto vfov = 20.0;
    auto aperture = 0.0;
    auto dist_to_focus = 10.0;

    return camera(
        look_from, look_at, vup, vfov, aperture, dist_to_focus, image_width, image_height
    );
}

camera cam_simple_light(size_t image_width, size_t image_height) {
    point3 look_from(26, 3, 6);
    point3 look_at(0, 2, 0);
    vec3 vup(0, 1, 0);
    auto vfov = 20.0;
    auto aperture = 0.0;
    auto dist_to_focus = 10.0;

    return camera(
        look_from, look_at, vup, vfov, aperture, dist_to_focus, image_width, image_height
    );
}

camera cam_cornell_box(size_t image_width, size_t image_height) {
    point3 look_from(278, 278, -760);
    point3 look_at(278, 278, 0);
    vec3 vup(0, 1, 0);
    auto vfov = 40.0;
    auto aperture = 0.0;
    auto dist_to_focus = 10.0;

    return camera(
        look_from, look_at, vup, vfov, aperture, dist_to_focus, image_width, image_height
    );
}

camera cam_caustics(size_t image_width, size_t image_height) {
    point3 look_from(0, 500, 0);
    point3 look_at(0, 0, 0);
    vec3 vup(0, 0, -1);
    auto vfov = 40.0;
    auto aperture = 0.0;
    auto dist_to_focus = 10.0;

    return camera(
        look_from, look_at, vup, vfov, aperture, dist_to_focus, image_width, image_height
    );
}

camera cam_phong(size_t image_width, size_t image_height) {
    point3 look_from(0, 200, -1000);
    point3 look_at(0, 100, 0);
    vec3 vup(0, 1, 0);

    auto vfov = 20.0;
    auto aperture = 0.0;
    auto dist_to_focus = 10.0;

    double r = 1000;
    double theta = degrees_to_radians(0);
    double phi = degrees_to_radians(20);

    return camera(
        r, theta, phi, look_at, vup, vfov, aperture, dist_to_focus, image_width, image_height
    );
}

camera cam_book2(size_t image_width, size_t image_height) {
    point3 look_from(478, 278, -600);
    point3 look_at(278, 278, 0);
    vec3 vup(0, 1, 0);
    auto vfov = 40.0;
    auto aperture = 0.0;
    auto dist_to_focus = 10.0;

    return camera(
        look_from, look_at, vup, vfov, aperture, dist_to_focus, image_width, image_height
    );
}

camera cam_tokyo(size_t image_width, size_t image_height) {
    point3 look_from(-6.0, 1.0, 0.0);
    point3 look_at(0.0, 0.0, 0.0);
    vec3 vup(0.0, 1.0, 0.0);
    auto vfov = 40.0;
    auto aperture = 0.05;
    auto dist_to_focus = (look_from - look_at).mag() - 1;

    return camera(
        look_from, look_at, vup, vfov, aperture, dist_to_focus, image_width, image_height
    );
}

camera cam_studio(size_t image_width, size_t image_height) {
    point3 look_from(-6.0, 0.8, -0.1);
    point3 look_at(0.0, 0.0, -0.1);
    vec3 vup(0.0, 1.0, 0.0);
    auto vfov = 40.0;
    auto aperture = 0.0;
    auto dist_to_focus = (look_from - look_at).mag();

    return camera(
        look_from, look_at, vup, vfov, aperture, dist_to_focus, image_width, image_height
    );
}

camera cam_park(size_t image_width, size_t image_height) {
    point3 look_from(-6.0, 1.0, -0.1);
    point3 look_at(0.0, 0.0, -0.1);
    vec3 vup(0.0, 1.0, 0.0);
    auto vfov = 40.0;
    auto aperture = 0.05;
    auto dist_to_focus = (look_from - look_at).mag();

    return camera(
        look_from, look_at, vup, vfov, aperture, dist_to_focus, image_width, image_height
    );
}

camera cam_xmas(size_t image_width, size_t image_height) {
    point3 look_from(0, 250, -1000);
    point3 look_at(0, 0, 0);
    vec3 vup(0, 1, 0);
    auto vfov = 40.0;
    auto aperture = 0.1;
    auto dist_to_focus = 800;

    return camera(
        look_from, look_at, vup, vfov, aperture, dist_to_focus, image_width, image_height
    );
}

scene get_scene(mem_arena &arena, int s, size_t image_width, size_t image_height) {
    // auto white = arena.alloc<solid_colour>(1.0f);
    auto black = arena.alloc<solid_colour>(0.0f);
    auto grey = arena.alloc<solid_colour>(0.05f);
    auto sky_blue = arena.alloc<solid_colour>(colour(0.7f, 0.8f, 1.0f));

    switch(s) {
    default:
        std::cerr << "Unknown scene\n.";
    case SCENE_BOOK_1:
        return scene(world_book1(arena), cam_book1(image_width, image_height), sky_blue);
        break;
    case SCENE_TWO_SPHERES:
        return scene(world_two_spheres(arena), cam_phong(image_width, image_height), sky_blue);
        break;
    case SCENE_PERLIN_SPHERES:
        return scene(world_perlin_spheres(arena), cam_two_spheres(image_width, image_height), sky_blue);
        break;
    case SCENE_EARTH:
        return scene(world_earth(arena), cam_earth(image_width, image_height), sky_blue);
        break;
    case SCENE_SIMPLE_LIGHT:
        return scene(world_simple_light(arena), cam_simple_light(image_width, image_height), black);
        break;
    case SCENE_CORNELL_BOX:
        return scene(world_cornell_box(arena, /* big_light */ true), cam_cornell_box(image_width, image_height), black);
        break;
    case SCENE_CORNELL_BOX_SMOKE:
        return scene(world_cornell_box_smoke(arena, /* big_light */ true), cam_cornell_box(image_width, image_height), black);
        break;
    case SCENE_CORNELL_BOX_SPECULAR:
        return scene(world_cornell_box_specular(arena, /* big_light */ true), cam_cornell_box(image_width, image_height), black);
        break;
    case SCENE_CORNELL_BOX_BUNNY:
        return scene(world_cornell_box_bunny(arena, /* big_light */ true), cam_cornell_box(image_width, image_height), black);
        break;
    case SCENE_CORNELL_BOX_DRAGON:
        return scene(world_cornell_box_dragon(arena, /* big_light */ true), cam_cornell_box(image_width, image_height), black);
        break;
    case SCENE_CORNELL_BOX_GLASS:
        return scene(world_cornell_box_glass(arena, /* big_light */ true), cam_cornell_box(image_width, image_height), sky_blue);
        break;
    case SCENE_CORNELL_BOX_GLASS_SPHERES:
        return scene(world_cornell_box_glass_spheres(arena, /* big_light */ false), cam_cornell_box(image_width, image_height), sky_blue);
        break;
    case SCENE_CAUSTICS:
        return scene(world_caustics(arena), cam_caustics(image_width, image_height), grey);
        break;
    case SCENE_PHONG:
        return scene(world_phong(arena), cam_phong(image_width, image_height), black);
        break;
    case SCENE_BOOK_2:
        return scene(world_book2(arena, /* rgb_light */ false), cam_book2(image_width, image_height), black);
        break;
    case SCENE_TOKYO:
        return scene(world_tokyo(arena), cam_tokyo(image_width, image_height), arena.alloc<image_texture>("media/tokyo.hdr"), true);
        //return scene(world_tokyo(arena), cam_tokyo(image_width, image_height), arena.alloc<solid_colour>(colour(1.0)), true);
        break;
    case SCENE_STUDIO:
        return scene(world_studio(arena), cam_studio(image_width, image_height), arena.alloc<image_texture>("media/room.hdr"), true);
        // return scene(world_studio(arena), cam_studio(image_width, image_height), arena.alloc<shifted_texture>(arena.alloc<image_texture>("media/studio.hdr"), 0.25, 0.0), true);
        break;
    case SCENE_PARK:
        return scene(world_park(arena), cam_park(image_width, image_height), arena.alloc<image_texture>("media/path.hdr"), true);
        break;
    case SCENE_XMAS:
        return scene(world_xmas(arena), cam_xmas(image_width, image_height), black);
        break;
    }
}
