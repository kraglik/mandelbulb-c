#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "qbmp/qdbmp.h"

int limit = 128;
float fPow = 8;
float shiftValue = 0.95;
float epsilon = 0.000001;
int itLimit = 400;

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
    Color color;
    float distance;
    int depth;
    char inf;
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
    float a = v1->x, b = v1->y, c = v1->z;
    float x = v2->x, y = v2->y, z = v2->z;

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


float phi(Vector * v) {
    return atanf(v->y / v->x);
}

float theta(Vector * v) {
    return acosf(v->z / len(v));
}

void powVec(Vector * v, Vector * result) {
    float r = len(v);
    float ph = phi(v);
    float th = theta(v);
    float nx = sinf(fPow * th) * sinf(fPow * ph);
    float ny = sinf(fPow * th) * sinf(fPow * ph);
    float nz = cosf(fPow * th);
    float p = powf(r, fPow);

    result->x = nx * p;
    result->y = ny * p;
    result->z = nz * p;
}

Pair iterateZ(float dr, Vector * z, Vector * c, int level) {
    float r = len(z);
    Vector zn;
    powVec(z, &zn);
    add(&zn, c, &zn);
    float drn = powf(r, fPow - 1) * fPow * dr + 1.0f;

    Pair result;

    if (level > limit || r > 2.0f) {
        result.x = r;
        result.y = drn;
    } else {
        result = iterateZ(drn, &zn, c, level + 1);
    }

    return result;
}

float distance(Vector * point) {
    Pair p = iterateZ(1.0f, point, point, 0);

    return 0.5f * logf(p.x * p.x / p.y);
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

    float x = distance(&p_plus_a) - distance (&p_minus_a);
    float y = distance(&p_plus_b) - distance (&p_minus_b);
    float z = distance(&p_plus_c) - distance (&p_minus_c);

    result->x = x;
    result->y = y;
    result->z = z;
}

void shift(Ray * ray, float mul) {
    Vector multiplier;

    sMul(&ray->dir, mul * shiftValue, &multiplier);

    add(&ray->pos, &multiplier, &ray->pos);
}

Hit raymarch(Ray * ray, float pathLen, int level) {
    Hit result = {.inf = 0};

    if (level > limit) {
        result.inf = 1;
        result.color.r = 25;
        result.color.g = 25;
        result.color.b = 25;
    }

    float d = distance(&ray->pos);

    if (d < epsilon) {
        result.point = ray->pos;
        fracNormal(&result.point, &ray->dir, &result.direction);
        result.distance = pathLen;
        result.depth = level;

    } else {

        Vector temp = ray->pos;
        shift(ray, d);
        sub(&temp, &ray->pos, &temp);

        result = raymarch(ray, pathLen + len(&temp), level + 1);

    }

    return result;
}

Hit traceRay() {

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
    BMP* bmp = BMP_Create(2000, 2000, 24);

    Vector camera_position = {.x = -1.0f, .y = -0.3f, .z = 1.0f};
    Vector camera_direction = {.x = 1.0f, .y = 0.0f, .z = -1.0f};
    Camera camera = buildCamera();

    BMP_WriteFile(bmp, "result.bmp");
    return 0;
}