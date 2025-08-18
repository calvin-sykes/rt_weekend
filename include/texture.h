#ifndef TEXTURE_H
#define TEXTURE_H

#include "rt_weekend.h"

#include "arena.h"
#include "perlin.h"
#include "rtw_stb_image.h"

#include <iostream>

class texture {
public:
    virtual ~texture() {}
    virtual colour value(double u, double v, const point3 &p) const = 0;
};

class solid_colour : public texture {
public:
    solid_colour() {}
    solid_colour(const colour& c) : colour_value(c) {}

    solid_colour(double red, double green, double blue)
        : solid_colour(colour(red, green, blue)) {}

    virtual colour value(double u, double v, const vec3& p) const override {
        return colour_value;
    }

private:
    colour colour_value;
};

class checker_texture : public texture {
public:
    checker_texture(const texture* even, const texture* odd, double scale=10)
        : tex_even(even), tex_odd(odd), scl(scale) {}
    checker_texture(mem_arena &arena, const colour& c1, const colour& c2, double scale=10)
        : tex_even(arena.alloc<solid_colour>(c1)),
          tex_odd (arena.alloc<solid_colour>(c2)), scl(scale) {}

    virtual colour value(double u, double v, const point3 &p) const override {
        auto sines = sin(scl * p.x()) * /* sin(scl * p.y()) * */ sin(scl * p.z());
        return sines < 0 ? tex_odd->value(u, v, p) : tex_even->value(u, v, p);
    }

private:
    const texture *tex_even, *tex_odd;
    double scl;
};

class noise_texture : public texture {
public:
    noise_texture() : scl(0) {}
    noise_texture(double scale) : scl(scale) {}

    virtual colour value(double u, double v, const point3 &p) const override {
        //return colour(1, 1, 1) * 0.5 * (1.0 + noise.noise(scl * p));
        //return colour(1, 1, 1) * noise.turb(scl * p);
        return colour(1, 1, 1) * 0.5 * (1 + sin(scl * p.z() + 10 * noise.turb(scl * p)));
    }
    
private:
    perlin noise;
    double scl;
};

class image_texture : public texture {
public:
    image_texture() {}

    image_texture(const char* filename) : image(filename) {}

    virtual colour value(double u, double v, const vec3 &p) const override {
        // Clamp input texture coordinates to [0, 1] x [1, 0]
        u = clamp(u, 0.0, 1.0);
        v = 1.0 - clamp(v, 0.0, 1.0);

        auto i = static_cast<int>(u * image.width());
        auto j = static_cast<int>(v * image.height());

        auto pixel = image.pixel_data(i, j);

        //constexpr float colour_scale = 1.0f / 255.0f;
        return colour(pixel[0], pixel[1], pixel[2]);
    }

private:
    rtw_image image;
};

class shifted_texture : public texture {
public:
    shifted_texture(const texture *source_texture, double du, double dv) 
    : t(source_texture), du(du), dv(dv) {}

    virtual colour value(double u, double v, const vec3 &p) const override {
        // Clamp input texture coordinates to [0, 1] x [1, 0]
        u = fmod(u + du, 1.0);
        v = fmod(v + dv, 1.0);
        return t->value(u, v, p);
    }

private:
    const texture *t;
    double du, dv;
};

enum blend_mode {
    MULTIPLY,
    ADD,
    AVERAGE
};

class blend_texture : public texture {
public:
    blend_texture(const texture *t1, const texture *t2, blend_mode mode) 
    : t1(t1), t2(t2), mode(mode) {}

    virtual colour value(double u, double v, const vec3 &p) const override {
        switch(mode) {
        case MULTIPLY:
            return t1->value(u, v, p) * t2->value(u, v, p);
        case ADD:
            return t1->value(u, v, p) + t2->value(u, v, p);
        case AVERAGE:
            return 0.5 * (t1->value(u, v, p) * t2->value(u, v, p));
        default:
            return colour(0);
        }
    }

private:
    const texture *t1, *t2;
    blend_mode mode;
};

#endif
