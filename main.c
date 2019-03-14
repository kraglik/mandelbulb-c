#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "qbmp/qdbmp.h"

int limit = 128;
float fPow = 6;
float shiftValue = 0.95;
float epsilon = 0.0001;
int itLimit = 256;
float xMul = 1.0f;

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
    Vector point, direction;
    float distance;
    int depth;
} Hit;

typedef struct Camera {
    Vector pos, nDir, nvX, nvY;
} Camera;

/**********************************************************************************************************************/

float len(Vector * vector) {
    return sqrtf(vector->x * vector->x + vector->y * vector->y + vector->z * vector->z);
}

void normalize(Vector * vector) {
    float l = len(vector);

    if (l < 0.00001f) {
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

void sMul(Vector * v, float s, Vector * result) {
    result->x = v->x * s;
    result->y = v->y * s;
    result->z = v->z * s;
}

void sDiv(Vector * v, float s, Vector * result) {
    result->x = v->x / s;
    result->y = v->y / s;
    result->z = v->z / s;
}

/**********************************************************************************************************************/

void powVec(Vector * v, Vector * result) {
    float ph = atanf(v->y / v->x);
    float th = acosf(v->z / len(v));

    result->x = sinf(fPow * th) * cosf(fPow * ph);
    result->y = sinf(fPow * th) * sinf(fPow * ph);
    result->z = cosf(fPow * th);

    sMul(result, powf(len(v), fPow), result);
}

Pair iterateZ(float dr, Vector * z, Vector * c, int level) {
    float r = len(z);

    Vector zn;
    powVec(z, &zn);
    add(&zn, c, &zn);

    float drn = powf(r, fPow - 1.0f) * fPow * dr + 1.0f;

    Pair result;

    if (level > limit || r > 2.0f) {
        result.x = r;
        result.y = dr;
    } else {
        result = iterateZ(drn, &zn, c, level + 1);
    }

    return result;
}

float distance(Vector * point) {
    Pair p = iterateZ(1.0f, point, point, 0);

    return (0.5f * logf(p.x) * p.x) / p.y;
}

void fracNormal(Vector * point, Vector * direction, Vector * result) {
    Vector a = {.x = direction->x, .y = 0.0f, .z = 0.0f};
    Vector b = {.x = 0.0f, .y = direction->y, .z = 0.0f};
    Vector c = {.x = 0.0f, .y = 0.0f, .z = direction->z};

    Vector p_plus_a, p_minus_a, p_plus_b, p_minus_b, p_plus_c, p_minus_c;

    add(point, &a, &p_plus_a);
    sub(point, &a, &p_minus_a);

    add(point, &b, &p_plus_b);
    sub(point, &b, &p_minus_b);

    add(point, &c, &p_plus_c);
    sub(point, &c, &p_minus_c);

    result->x = distance(&p_plus_a) - distance (&p_minus_a);
    result->y = distance(&p_plus_b) - distance (&p_minus_b);
    result->z = distance(&p_plus_c) - distance (&p_minus_c);
}

void shift(Ray * ray, float mul) {
    Vector multiplier;

    sMul(&ray->dir, mul * shiftValue, &multiplier);

    add(&ray->pos, &multiplier, &ray->pos);
}

Hit raymarch(Ray * ray, float pathLen, int level) {
    Hit hit = {};

    if (level > itLimit) {
        hit.distance = INFINITY;
        hit.depth = level;

        return hit;
    }

    float d = distance(&ray->pos);

    if (d < epsilon) {
        hit.point = ray->pos;
        fracNormal(&hit.point, &ray->dir, &hit.direction);
        hit.distance = pathLen;
        hit.depth = level;

    } else if (isinf(d) || isnan(d)) {
        hit.distance = INFINITY;
        hit.depth = level;
    } else {

        Vector temp = ray->pos;
        shift(ray, d);
        sub(&temp, &ray->pos, &temp);

        hit = raymarch(ray, pathLen + len(&temp), level + 1);

    }

    return hit;
}

Hit traceRay(Camera * camera, float x, float y) {
    Ray ray = {};

    sMul(&camera->nvX, x * xMul, &ray.dir);
    sMul(&camera->nvY, y, &ray.pos);
    add(&ray.pos, &ray.dir, &ray.dir);
    add(&camera->nDir, &ray.dir, &ray.dir);
    normalize(&ray.dir);

    ray.pos = camera->pos;

    return raymarch(&ray, 0.0f, 0);
}

Camera buildCamera(Vector position, Vector direction, Vector up) {
    Camera camera = { .pos = position, .nDir = direction };

    normalize(&camera.nDir);

    cross(&up, &direction, &camera.nvX);
    normalize(&camera.nvX);

    cross(&camera.nDir, &camera.nvX, &camera.nvY);
    normalize(&camera.nvY);

    return camera;
}

/**********************************************************************************************************************/

int main() {
    int width = 1000, height = 1000;
    int m_width = width / 2, m_height = height / 2;

    BMP* bmp = BMP_Create(width, height, 32);

    Vector camera_position = {.x = -1.0f, .y = 0.0f, .z = -1.0f};

    normalize(&camera_position);
    sMul(&camera_position, 1.55f, &camera_position);

    Vector camera_direction = {.x = 1.0f, .y = 0.0f, .z = 1.0f};
    Vector camera_up = {.x = 0.0f, .y = 1.0f, .z = 0.0f};

    Camera camera = buildCamera(camera_position, camera_direction, camera_up);

    #pragma omp parallel for
    for (int i = 0; i < width; i++) {
        for (int j = 0; j < height; j++) {
            float x = ((float)i - ((float)m_width)) / ((float)m_width);
            float y = ((float)j - ((float)m_height)) / ((float)m_height);

            Hit hit = traceRay(&camera, x, y);

            Color color = {.r = 0, .g = 0, .b = 0};

            if (!isinf(hit.distance)) {
                float color_strength = 1.0f - (float)hit.depth / (float)itLimit;

                color.r = 0;
                color.g = (unsigned char)(color_strength * 255.0f);
                color.b = (unsigned char)(color_strength * 255.0f);
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