#ifndef CAMERA_H
#define CAMERA_H

#include "rt_weekend.h"

class camera {
public:
    camera(
        point3 pos, point3 look_at, vec3 vup,
        double vfov, double aperture, double focus_dist,
        size_t image_width, size_t image_height
    ) {
        pixel_delta_x = 1.0 / image_width;
        pixel_delta_y = 1.0 / image_height;
        auto aspect_ratio = static_cast<double>(image_width) / image_height;

        auto theta = degrees_to_radians(vfov);
        auto h = tan(theta / 2);
        auto viewport_height = 2.0 * h;
        auto viewport_width = aspect_ratio * viewport_height;

        w = unit_vector(pos - look_at);
        u = unit_vector(cross(vup, w));
        v = cross(w, u);
        
        origin = pos;
        horizontal = focus_dist * viewport_width * u;
        vertical = focus_dist * viewport_height * v;
        lower_left_corner = origin - horizontal/2 - vertical/2 - focus_dist * w;

        lens_radius = aperture / 2;
        t0 = 0.0;
        t1 = 1.0;
    }

    camera(
        double r, double theta, double phi,
        point3 look_at, vec3 vup,
        double vfov, double aperture, double focus_dist,
        size_t image_width, size_t image_height
    ) : camera(orbit_pos(r, theta, phi) + look_at, look_at, vup,
        vfov, aperture, focus_dist, image_width, image_height) {};

    ray get_ray(size_t i, size_t j) const {
        // Pixel position and sample jitter
        double ds = random_number();
        double dt = random_number();
        
        // double ux = 2 * random_number();
        // double uy = 2 * random_number();

        // double ds, dt;
        // if (ux < 1)
        //     ds = sqrt(ux) / 2;
        // else
        //     ds = (2 - sqrt(2 - ux)) / 2;
        // if (uy < 1)
        //     dt = sqrt(uy) / 2;
        // else
        //     dt = 2 - sqrt(2 - uy) / 2;

        double s = (i + ds) * pixel_delta_x;
        double t = (j + dt) * pixel_delta_y;

        // Depth-of-field jitter
        vec3 offset;
        if (lens_radius > 0.0) {
            auto lens_pos = random_in_unit_sphere<vec2>();
            vec3 rd = lens_radius * vec3(lens_pos.x(), lens_pos.y(), 0.0);
            offset = u * rd.x() + v * rd.y();
        }

        return ray(
            origin + offset,
            lower_left_corner + s * horizontal + t * vertical - origin - offset,
            random_number(t0, t1)
        );
    }

private:
    static point3 orbit_pos(double r, double theta, double phi) {
        double sin_theta = sin(theta), cos_theta = cos(theta);
        double cos_phi = cos(phi), sin_phi = sin(phi);
        return point3(
            r * cos_phi * sin_theta,
            r * sin_phi,
            r * cos_phi * cos_theta
        );
    };

    point3 origin;
    point3 lower_left_corner;
    vec3 horizontal;
    vec3 vertical;
    vec3 u, v, w;
    double lens_radius;
    double t0, t1; // shutter open/close times
    double pixel_delta_x, pixel_delta_y;
};

#endif