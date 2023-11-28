#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <iostream>
#include <string>
#include "SETTINGS.h"

// for norming a vec
// !! there is a duplicate func in main
void normalize(VEC3& ray) {
    double mag = sqrt(ray.dot(ray));
    ray[0] = ray[0]/mag;
    ray[1] = ray[1]/mag;
    ray[2] = ray[2]/mag;
    return;
};

double magnitude(VEC3 ray) { 
    return sqrt(ray.dot(ray));
};

// assume ray points away from surface, 
// and norm must point away from surface
VEC3 reflect(VEC3 ray, VEC3 surf_norm) {
    VEC3 v = ray.dot(surf_norm)*surf_norm;
    return (-1)*ray + 2*v;
} 