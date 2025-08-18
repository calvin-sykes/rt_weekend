#ifndef COLOUR_H
#define COLOUR_H

#include "vector.h"

#include "stb_image_write.h"

#include <iostream>
#include <fstream>
#include <vector>

class image_writer {
public:
    virtual void write(std::vector<colour> const& buffer, size_t width, size_t height) const = 0;
};
class text_writer : public image_writer {
public:
    text_writer() : file(nullptr) {}
    text_writer(const char *filename) : file(filename) {}

    virtual void write(std::vector<colour> const& buffer, size_t width, size_t height) const override {
        if (file) {
            std::ofstream stream(file);
            write_stream(stream, buffer, width, height);
        } else {
            write_stream(std::cout, buffer, width, height);
        }
    }

private:
    void write_stream(std::ostream &os, std::vector<colour> buffer, size_t width, size_t height) const {
        for (size_t j = height - 1; j != static_cast<size_t>(-1); j--) {
            for (size_t i = 0; i < width; i++) {
                auto idx = i + width * j;
                os << buffer[idx] << '\n';
            }
        }
    }

    const char *file;
};

class ppm_writer : public image_writer {
public:
    ppm_writer() : file(nullptr) {}
    ppm_writer(const char *filename) : file(filename) {}

    virtual void write(std::vector<colour> const& buffer, size_t width, size_t height) const override {
        if (file) {
            std::ofstream stream(file);
            write_stream(stream, buffer, width, height);
        } else {
            write_stream(std::cout, buffer, width, height);
        }
    }

private:
    void write_stream(std::ostream &os, std::vector<colour> buffer, size_t width, size_t height) const {
        os << "P3\n" << width << ' ' << height << "\n255\n";
        for (size_t j = height - 1; j != static_cast<size_t>(-1); j--) {
            for (size_t i = 0; i < width; i++) {
                auto idx = i + width * j;
                auto const& c = buffer[idx];
                os << static_cast<uint8_t>(255 * clamp(c[0], 0.0f, 1.0f)) << ' '
                   << static_cast<uint8_t>(255 * clamp(c[1], 0.0f, 1.0f)) << ' '
                   << static_cast<uint8_t>(255 * clamp(c[2], 0.0f, 1.0f)) << '\n';
            }
        }
    }

    const char *file;
};

class png_writer : public image_writer {
public:
    png_writer(const char *filename) : file(filename) {}

    virtual void write(std::vector<colour> const& buffer, size_t width, size_t height) const override {
        // write image to 8-bit buffer
        const int bytes_per_pixel = 3;
        std::vector<uint8_t> uc_buffer(width * height * bytes_per_pixel);

        for (size_t j = 0; j < height; j++) {
            for (size_t i = 0; i < width; i++) {
                auto idx = i + width * j;
                auto const& c = buffer[idx];
                uc_buffer[3 * idx    ] = static_cast<uint8_t>(255 * clamp(c[0], 0.0f, 1.0f));
                uc_buffer[3 * idx + 1] = static_cast<uint8_t>(255 * clamp(c[1], 0.0f, 1.0f));
                uc_buffer[3 * idx + 2] = static_cast<uint8_t>(255 * clamp(c[2], 0.0f, 1.0f));
            }
        }
        stbi_flip_vertically_on_write(1);
        stbi_write_png(file, static_cast<int>(width), static_cast<int>(height), bytes_per_pixel, uc_buffer.data(), 0);
    }

private:
    const char *file;
};

class hdr_writer : public image_writer {
public:
    hdr_writer(const char *filename) : file(filename) {}

    virtual void write(std::vector<colour> const& buffer, size_t width, size_t height) const override {
        const float *buffer_view = *buffer.data();
        stbi_flip_vertically_on_write(1);
        stbi_write_hdr(file, static_cast<int>(width), static_cast<int>(height), 3, buffer_view);
    }

private:
    const char *file;
};

extern void tone_mapping(std::vector<colour>& buffer, float exposure);
extern void write_colour(std::ostream &out, colour const& pixel_colour);
extern void write_colour(std::vector<unsigned char> &buffer, size_t pos, colour const& pixel_colour);

#endif