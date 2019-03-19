#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "qbmp/qdbmp.h"

int limit = 64;
float fPow = 8;
float shiftValue = 0.95;
float epsilon = 0.0001;
int itLimit = 64;

/**********************************************************************************************************************/

typedef struct Vector {
    float x, y, z;
} Vector;

typedef struct Ray {
    Vector pos, dir;
} Ray;

typedef struct Pair {
    float x, y;
} Pair;

typedef struct Color {
    unsigned char r, g, b;
} Color;

typedef struct Hit {
    float distance;
    int depth;
} Hit;

typedef struct Camera {

    Vector pos, dir, up, u, w, v;

    float view_plane_distance, ratio, shift_multiplier;
    int width, height;

} Camera;

/**********************************************************************************************************************/

float len(Vector * vector) {
    return sqrtf(vector->x * vector->x + vector->y * vector->y + vector->z * vector->z);
}

void normalize(Vector * vector) {
    float l = len(vector);

    if (l == 0.0f) {
        vector->x = 0.0;
        vector->y = 0.0;
        vector->z = 0.0;
    }

    else {
        vector->x = vector->x / l;
        vector->y = vector->y / l;
        vector->z = vector->z / l;
    }
}

void add(Vector * a, Vector * b, Vector * result) {
    result->x = a->x + b->x;
    result->y = a->y + b->y;
    result->z = a->z + b->z;
}

void sub(Vector * a, Vector * b, Vector * result) {
    result->x = a->x - b->x;
    result->y = a->y - b->y;
    result->z = a->z - b->z;
}

void cross(Vector * v1, Vector * v2, Vector * result) {
    float a = v2->x, b = v2->y, c = v2->z;
    float x = v1->x, y = v1->y, z = v1->z;

    result->x = y * c - z * b;
    result->y = x * c - a * z;
    result->z = x * b - y * a;
}

float dot(Vector * v1, Vector * v2) {
    return v1->x * v2->x + v1->y * v2->y + v1->z * v2->z;
}

void scalar_mul(Vector *v, float s, Vector *result) {
    result->x = v->x * s;
    result->y = v->y * s;
    result->z = v->z * s;
}

/**********************************************************************************************************************/

void pow_vec(Vector *v, Vector *result) {
    float ph = atanf(v->y / v->x);
    float th = acosf(v->z / len(v));

    result->x = sinf(fPow * th) * cosf(fPow * ph);
    result->y = sinf(fPow * th) * sinf(fPow * ph);
    result->z = cosf(fPow * th);

    scalar_mul(result, powf(len(v), fPow), result);
}

Pair iterate_z(float dr, Vector z, Vector *c) {
    Pair pair;
    pair.x = 0.0f;
    pair.y = dr;

    for (int i = 0;; i++) {

        float r = len(&z);
        Vector zn;
        pow_vec(&z, &zn);
        add(&zn, c, &zn);

        if (i > limit || r > 2.0f) {

            pair.x = r;
            pair.y = dr;

            break;

        } else {
            dr = powf(r, fPow - 1.0f) * fPow * dr + 1.0f;
            z = zn;

        }
    }

    return pair;
}

float distance(Vector * point) {
    Pair p = iterate_z(1.0f, *point, point);

    return (0.5f * logf(p.x) * p.x) / p.y;
}

Hit march_ray(Ray *ray, float pathLen) {
    Hit hit = { .distance = INFINITY };
    Vector temp;

    for (int i = 0; i < itLimit; i++) {
        float d = distance(&ray->pos);

        if (d < epsilon && !(isinf(d) || isnan(d))) {
            hit.distance = pathLen;
            hit.depth = i;

            break;

        } else {

            scalar_mul(&ray->dir, d * shiftValue, &temp);
            add(&temp,&ray->pos, &ray->pos);

            pathLen += d * shiftValue;

        }
    }

    return hit;
}

Camera build_camera(
        int width,
        int height,
        float shift_multiplier,
        float view_plane_distance,
        float ratio,
        Vector position,
        Vector target,
        Vector up) {

    Camera camera = {
            .pos = position,
            .width = width,
            .height = height,
            .ratio = ratio,
            .shift_multiplier = shift_multiplier,
            .up = up,
            .view_plane_distance = view_plane_distance
    };

    sub(&target, &position, &camera.dir);
    normalize(&camera.dir);

    sub(&position, &target, &camera.w);
    normalize(&camera.w);

    cross(&up, &camera.w, &camera.u);
    normalize(&camera.u);

    cross(&camera.w, &camera.u, &camera.v);
    normalize(&camera.v);

    return camera;
}

void move_ray_to_view_plane(Ray * ray, Camera * camera) {
    float shift_value = 1.0f / dot(&ray->dir, &camera->dir);
    Vector to_add;
    scalar_mul(&ray->dir, shift_value * camera->shift_multiplier, &to_add);
    add(&ray->pos, &to_add, &ray->pos);
}


Hit trace_ray(Camera *camera, float x, float y) {
    Ray ray = {};

    Vector temp;

    scalar_mul(&camera->u, x * camera->ratio, &ray.dir);
    scalar_mul(&camera->v, y, &temp);
    add(&temp, &ray.dir, &ray.dir);
    scalar_mul(&camera->w, camera->view_plane_distance, &temp);
    sub(&ray.dir, &temp, &ray.dir);
    normalize(&ray.dir);

    ray.pos = camera->pos;

    move_ray_to_view_plane(&ray, camera);

    return march_ray(&ray, 0.0f);
}

/**********************************************************************************************************************/

float clip(float x, float min, float max) {
    return x > max ? max : x < min ? min : x;
}

/**********************************************************************************************************************/

int main() {
    int width = 2000, height = 2000;
    int m_width = width / 2, m_height = height / 2;

    float brightness = 1.25f;

    BMP* bmp = BMP_Create(width, height, 32);

    Vector camera_position = {.x = 1.0f, .y = 0.0f, .z = 1.0f};

    Vector camera_target = {.x = 0.0f, .y = 0.0f, .z = 0.0f};
    Vector camera_up = {.x = 0.0f, .y = 1.0f, .z = 0.0f};

    float view_plane_distance = 1.0f;
    float shift_multiplier = view_plane_distance >= 1.0f ? 0.25f / view_plane_distance : view_plane_distance / 4.0f;
    float ratio = 1.0f;

    Camera camera = build_camera(
            width,
            height,
            shift_multiplier,
            view_plane_distance,
            ratio,
            camera_position,
            camera_target,
            camera_up
    );

    #pragma omp parallel for
    for (int i = 0; i < width; i++) {
        for (int j = 0; j < height; j++) {
            float x = ((float)i - ((float)m_width)) / ((float)m_width);
            float y = ((float)j - ((float)m_height)) / ((float)m_height);

            Hit hit = trace_ray(&camera, x, y);

            Color color = {.r = 0, .g = 0, .b = 0};

            if (!isinf(hit.distance)) {
                float color_strength = 1.0f - (float)hit.depth / (float)itLimit;

                color.r = 0;
                color.g = (unsigned char)clip((color_strength * 255.0f) * brightness, 0.0f, 255.0f);
                color.b = (unsigned char)clip((color_strength * 255.0f) * brightness, 0.0f, 255.0f);
            }

            BMP_SetPixelRGB(bmp, j, i,
                color.r,
                color.g,
                color.b
            );
        }
    }

    BMP_WriteFile(bmp, "result.bmp");
    return 0;
}