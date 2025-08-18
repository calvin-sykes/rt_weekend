#ifndef PERLIN_H
#define PERLIN_H

#include "rt_weekend.h"

#include <algorithm>
#include <random>

class perlin {
public:
    perlin() {
        random_vecs = new vec3[point_count];

        for (int i = 0; i < point_count; i++)
            random_vecs[i] = random_unit_vector();
        
        perm_x = perlin_generate_permutation();
        perm_y = perlin_generate_permutation();
        perm_z = perlin_generate_permutation();
    }

    ~perlin() {
        delete[] random_vecs;
        delete[] perm_x;
        delete[] perm_y;
        delete[] perm_z;
    }

    double turb(const point3 &p, int depth=7) const {
        auto accum = 0.0;
        auto ptmp = p;
        auto weight = 1.0;

        for (int i = 0; i < depth; i++) {
            accum += weight * noise(ptmp);
            weight *= 0.5;
            ptmp *= 2;
        }

        return fabs(accum);
    }

    double noise(const point3 &p) const {
        auto u = p.x() - floor(p.x());
        auto v = p.y() - floor(p.y());
        auto w = p.z() - floor(p.z());

        auto i = static_cast<int>(floor(p.x())) & 255;
        auto j = static_cast<int>(floor(p.y())) & 255;
        auto k = static_cast<int>(floor(p.z())) & 255;
        vec3 c[2][2][2];

        for (int di = 0; di < 2; di++)
            for (int dj = 0; dj < 2; dj++)
                for (int dk = 0; dk < 2; dk++)
                    c[di][dj][dk] = random_vecs[
                        perm_x[(i + di) & 255] ^
                        perm_y[(j + dj) & 255] ^
                        perm_z[(k + dk) & 255]
                    ];

        return perlin_interp(c, u, v, w);
    }

private:
    static const int point_count = 256;
    vec3 *random_vecs;
    int *perm_x, *perm_y, *perm_z;

    static int *perlin_generate_permutation() {
        auto p = new int[point_count];

        for (int i = 0; i < point_count; i++)
            p[i] = i;

        std::shuffle(&p[0], &p[point_count - 1], std::mt19937());

        return p;
    }

    static double perlin_interp(vec3 c[2][2][2], double u, double v, double w) {
        auto uu = u * u * (3 - 2 * u);
        auto vv = v * v * (3 - 2 * v);
        auto ww = w * w* (3 - 2 * w);
        
        auto val = 0.0;
        for (int i = 0; i < 2; i++)
            for (int j = 0; j < 2; j++)
                for (int k = 0; k < 2; k++) {
                    vec3 weight_v(u - i, v - j, w - k);
                    val += (i * uu + (1 - i) * (1 - uu)) * 
                           (j * vv + (1 - j) * (1 - vv)) * 
                           (k * ww + (1 - k) * (1 - ww)) * 
                           dot(c[i][j][k], weight_v);
                }

        return val;
    }
};

#endif