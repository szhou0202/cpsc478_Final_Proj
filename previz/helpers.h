// what do we want to include here

#ifndef HELPERS_H
#define HELPERS_H

#include "normalize.h"

using namespace std;



// primitive
class Primitive {
public:
    VEC3 color;
    bool reflect; 
    bool refract;
    double refract_air;
    double refract_glass;

    Primitive(){};

    // Beautiful, we have primitives working
    virtual double findIntersect(VEC3 origin, VEC3 dir){
        return 0.0;
    }
    virtual VEC3 findSurfNorm(VEC3 eye, VEC3 dir, double t){
        return VEC3(0.,0.,0.);
    }
    virtual VEC3 findSurfPos(VEC3 eye, VEC3 dir, double t){
        return VEC3(0.,0.,0.);
    }

    virtual bool getReflect() {
        return false;
    };
    virtual bool getRefract() {
        return false;
    };
    virtual double getRa() {
        return 0.;
    }
    virtual double getRg() {
        return 0.;
    }
    virtual VEC3 getColor() {
        return color;
    }
    // any funcs for computing quad formula or smth?
};

class Sphere : public Primitive {
public:
    VEC3 center = VEC3(0.,0.,0.);
    double radius = 0.;
    bool reflect = false; 
    bool refract = false;
    double refract_air;
    double refract_glass;

    Sphere(VEC3 c, double r, VEC3 co) {
        center = c;
        radius = r;
        color = co;
    };
    Sphere(VEC3 c, double r, VEC3 co, bool ref, bool re, double r_a, double r_g) {
        center = c;
        radius = r;
        color = co;
        reflect = ref;
        refract = re;
        refract_air = r_a;
        refract_glass = r_g;
    };

    double findIntersect(VEC3 origin, VEC3 dir) {
        // calculate discriminant 
        VEC3 ominc = origin-center;
        double dd = dir.dot(dir);
        double disc = pow(dir.dot(ominc),2) - dd*(ominc.dot(ominc) - pow(radius,2));
        
        if (disc < 0) {
            return -1.;
        }
        if (disc == 0) {
            return (-1)*dir.dot(ominc) / dd;
        }

        double t1 = (((-1)*dir).dot(ominc) + sqrt(disc))/dd;
        double t2 = (((-1)*dir).dot(ominc) - sqrt(disc))/dd;
        // !! handle negatives? was a todo from before
        if (t1 < t2) {
            return t1;
        }
        return t2;
    }

    // find coords of surface point
    VEC3 findSurfPos(VEC3 eye, VEC3 dir, double t) {
        return eye + t*dir;
    }

    // find norm at a point on sphere
    VEC3 findSurfNorm(VEC3 eye, VEC3 dir, double t) {
        VEC3 toNorm = findSurfPos(eye,dir,t)-center; // THIS IS CORRECT
        normalize(toNorm);
        return toNorm;
    }

    bool getReflect() {
        return reflect;
    }
    bool getRefract() {
        return refract;
    }
    double getRa() {
        return refract_air;
    }
    double getRg() {
        return refract_glass;
    }
    VEC3 getColor() {
        return color;
    }

};

class Triangle : public Primitive {
public:
    VEC3 vertex1;
    VEC3 vertex2;
    VEC3 vertex3;
    VEC3 u;
    VEC3 v;
    // VEC3 color;

    Triangle(VEC3 v1, VEC3 v2, VEC3 v3, VEC3 co) {
        vertex1 = v1;
        vertex2 = v2;
        vertex3 = v3;
        u = vertex2-vertex1;
        v = vertex3-vertex1;
        color = co;
    };

    // ok lets see if this works! 
    double findIntersect(VEC3 origin, VEC3 dir) {
        // this is going to be complicated
        double a = vertex1[0]-vertex2[0]; // !! notice the direction of 1-2 vs 2-1
        double b = vertex1[1]-vertex2[1];
        double c = vertex1[2]-vertex2[2];

        double d = vertex1[0]-vertex3[0];
        double e = vertex1[1]-vertex3[1];
        double f = vertex1[2]-vertex3[2];

        double g = dir[0];
        double h = dir[1];
        double i = dir[2];

        double j = vertex1[0]-origin[0];
        double k = vertex1[1]-origin[1];
        double l = vertex1[2]-origin[2];

        double eihf = e*i - h*f;
        double akjb = a*k - j*b;
        double jcal = j*c - a*l;
        double blkc = b*l - k*c;
        double dheg = d*h - e*g;
        double gfdi = g*f - d*i;

        // compute M
        double M = a*eihf + b*gfdi + c*dheg;

        // compute beta
        double beta = (j*eihf + k*gfdi + l*dheg)/M;

        // compute gamma
        double gamma = (i*akjb + h*jcal + g*blkc)/M;

        // compute t
        double t = -(f*akjb + e*jcal + d*blkc)/M;

        if (t < 0)
            return -1.0;
        if (gamma < 0 || gamma > 1) 
            return -1.;
        if (beta < 0 || beta > 1-gamma) 
            return -1.;
        return t;
    }

    // find coords of surface point
    VEC3 findSurfPos(VEC3 eye, VEC3 dir, double t) {
        return eye + t*dir;
    }

    // find norm at a point (doesn't depend on these params)
    // since triangle is a plane
    VEC3 findSurfNorm(VEC3 eye, VEC3 dir, double t) {
        VEC3 toNorm = u.cross(v);
        normalize(toNorm);
        return toNorm;
    }


};

class Cylinder : public Primitive {
public: 
    VEC3 center = VEC3(0.,0.,0.);
    double radius = 0.;
    bool reflect = false; 
    bool refract = false;
    double refract_air;
    double refract_glass;

    Cylinder(VEC3 c, double r, VEC3 co) {

    };
    Cylinder(VEC3 c, double r, VEC3 co, bool ref, bool re, double r_a, double r_g) {
        center = c;
        radius = r;
        color = co;
        reflect = ref;
        refract = re;
        refract_air = r_a;
        refract_glass = r_g;
    };

    double findIntersect(VEC3 origin, VEC3 dir) {
        // !! how in the world do i do cylinder intersection...
        return 0.;
    }

    // find coords of surface point
    VEC3 findSurfPos(VEC3 eye, VEC3 dir, double t) {
        return eye + t*dir;
    }

    // find norm at a point on sphere
    VEC3 findSurfNorm(VEC3 eye, VEC3 dir, double t) {
        return VEC3(0,0,0);
    }

    bool getReflect() {
        return reflect;
    }
    bool getRefract() {
        return refract;
    }
    double getRa() {
        return refract_air;
    }
    double getRg() {
        return refract_glass;
    }
    VEC3 getColor() {
        return color;
    }
};


class Light {
public:
    VEC3 pos;
    VEC3 color;
    double intensity;
    double phong;

    Light(VEC3 p, VEC3 c, double i, double ph) {
        pos = p;
        color = c;
        intensity = i;
        phong = ph;
    };

    // return color 
    // assume vectors passed in are normalized
    VEC3 diffuseCol(VEC3 surf_color, VEC3 surf_norm, VEC3 surf_position) {
        VEC3 light_ray = pos-surf_position;
        normalize(light_ray);
        double temp = intensity*(0 >= surf_norm.dot(light_ray) ? 0 : surf_norm.dot(light_ray));
        return temp*VEC3(color[0]*surf_color[0],color[1]*surf_color[1],color[2]*surf_color[2]);
        // return temp;
    };

    // calculate specular lighting
    VEC3 specCol(VEC3 surf_color, VEC3 surf_norm, VEC3 surf_position, VEC3 eye){
        // points towards light
        VEC3 light_ray = pos-surf_position;
        normalize(light_ray);

        // points towards eye
        VEC3 view_ray = eye-surf_position; 
        normalize(view_ray);
        VEC3 reflected_light = reflect(light_ray, surf_norm);

        // using the dot prod of reflected light and 
        // view ray
        double t = reflected_light.dot(view_ray);
        double temp = intensity*(0 >= t ? 0 : pow(t,phong));
        return temp*VEC3(color[0]*surf_color[0],color[1]*surf_color[1],color[2]*surf_color[2]);
        // return temp;
    };

    VEC3 totCol(VEC3 surf_color, VEC3 surf_norm, VEC3 surf_position, VEC3 eye) {
        VEC3 dif = diffuseCol(surf_color, surf_norm, surf_position);
        VEC3 spec = specCol(surf_color, surf_norm, surf_position, eye);
        // return (dif + spec)*(VEC3(color[0]*surf_color[0],color[1]*surf_color[1],color[2]*surf_color[2]));
        return dif+spec;
    };
};



#endif

