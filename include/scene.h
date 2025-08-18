#ifndef SCENE_H
#define SCENE_H

#include "rt_weekend.h"

#include "background.h"
#include "camera.h"
#include "hittable_list.h"
#include "texture.h"

struct scene {
    scene(
        const hittable_list& _world, const camera& _cam, const texture *bkg_texture, bool sample_bkg = false
    ) : world(_world), cam(_cam), bkg(bkg_texture) {
        auto const& objects = world.get();

        for (size_t i = 0 ; i < objects.size(); i++) {
            if (objects[i]->is_sampling_target())
                lights.add(objects[i]);
        }
        if (sample_bkg)
            lights.add(&bkg);
    }

    hittable_list world, lights;
    camera cam;
    background bkg;
};

enum SCENES {
    SCENE_BOOK_1,
    SCENE_TWO_SPHERES,
    SCENE_PERLIN_SPHERES,
    SCENE_EARTH,
    SCENE_SIMPLE_LIGHT,
    SCENE_CORNELL_BOX,
    SCENE_CORNELL_BOX_SMOKE,
    SCENE_CORNELL_BOX_SPECULAR,
    SCENE_CORNELL_BOX_BUNNY,
    SCENE_CORNELL_BOX_DRAGON,
    SCENE_CORNELL_BOX_GLASS,
    SCENE_PHONG,
    SCENE_BOOK_2,
    SCENE_TOKYO,
    SCENE_STUDIO,
    SCENE_PARK,
    SCENE_XMAS,
    NUM_SCENES
};

scene get_scene(mem_arena &arena, int s, size_t image_width, size_t image_height);

#endif