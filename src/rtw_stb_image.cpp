#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION

#include "rtw_stb_image.h"

rtw_image::~rtw_image() {
    STBI_FREE(data);
}

bool rtw_image::load(const std::string filename) {
    // Loads image data from the given file name. Returns true if the load succeeded.
    auto n = bytes_per_pixel; // Dummy out parameter: original components per pixel
    data = stbi_loadf(filename.c_str(), &image_width, &image_height, &n, bytes_per_pixel);
    bytes_per_scanline = image_width * bytes_per_pixel;
 
    //float max_lum = -1;
    //float min_lum = INFINITY;
    //for (int i = 0; i < image_width * image_height * bytes_per_pixel; i+= bytes_per_pixel) {
    //    float lum = data[i] * 0.2126f + data[i+1] * 0.7152f + data[i+2] * 0.0722f;
    //    max_lum = std::max(max_lum, lum);
    //    min_lum = std::min(min_lum, lum);
    //}
    //max_lum = 0; //log10f(max_lum);
    //min_lum = -1; // log10f(min_lum);

    //for (int i = 0; i < image_width * image_height * bytes_per_pixel; i+= bytes_per_pixel) {
    //    float lum = log10f(data[i] * 0.2126f + data[i + 1] * 0.7152f + data[i + 2] * 0.0722f);
    //    float scl = (lum - min_lum) / (max_lum - min_lum);
    //    scl = scl < 0 ? 0 : (scl > 1 ? 1 : scl);

    //    data[i] *= scl;
    //    data[i+1] *= scl;
    //    data[i+2] *= scl;
    //}

    return data != nullptr;
}