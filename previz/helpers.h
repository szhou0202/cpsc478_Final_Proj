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
    // can define a cylinder via an endpoint, a ray, and a length
    // the change of basis will be 
    // the new basis in the old basis (just normal xyz)
    // so
    // [normalized cylinder ray, norm ray orthogonal to cylinder ray, norm ray orthogonal to the other two]
    // to get the orthogonal rays, we can check which ones of (100), (010), and (001) can be validly dot producted with the cylinder ray, 
    // e.g. the dot product is not 1 or -1, and take the cross of them (somewhat imprecise calculations with doubles haha)
    // 
    VEC3 top = VEC3(0.,0.,0.);
    VEC3 bot = VEC3(0.,0.,0.);
    VEC3 axi = VEC3(0.,0.,0.);
    double radius = 0.;
    double length = 0.;

    Cylinder(VEC3 t, VEC3 b, double r, VEC3 co) {
        color = co;
        top = t; // origin of the cylinder
        bot = b;
        axi = b-t;
        length = magnitude(axi);
        normalize(axi);
        radius = r;
    };

    // note: the gaze passed here is not the ref point, 
    // it is the vector from eye to the ref point
    // MATRIX4 cam_mtx(VEC3 eye, VEC3 gaze, VEC3 top) {
    //     MATRIX4 Mtemp;

    //     VEC3 w = -1*gaze.normalized();///(gaze.norm()); /// !! could i j call normalize?
    //     VEC3 temp = top.cross(w);
    //     VEC3 u = temp.normalized();//temp/temp.norm();
    //     VEC3 v = w.cross(u);

    //     Mtemp.setZero();
    //     Mtemp.block(0,0,3,1) = u;
    //     Mtemp.block(0,1,3,1) = v;
    //     Mtemp.block(0,2,3,1) = w;
    //     Mtemp.block(0,3,3,1) = eye;
    //     Mtemp(3,3) = 1.;

    //     return Mtemp.inverse();
    // }

    double findIntersect(VEC3 origin, VEC3 dir) {

        // FIRST, apply change of basis to the ray 
        // find change of basis array, normalized vectors
        // this is very similar to the MCAM matrix! 
        VEC3 tem;
        if (abs(axi.dot(VEC3(1,0,0))) != 1) {
            tem = VEC3(1,0,0);
        }
        else {
            tem = VEC3(0,0,1);
        }
        VEC3 u = axi.cross(tem);
        normalize(u);
        VEC3 v = axi.cross(u);
        normalize(v);

        MATRIX4 M;
        M.setZero();
        M.block(0,0,3,1) = u;
        M.block(0,1,3,1) = axi;
        M.block(0,2,3,1) = v;
        M.block(0,3,3,1) = top; 
        M(3,3) = 1.;
        M = M.inverse().eval();

        // VEC4 temp = VEC4(top[0]+axi[0], top[1]+axi[1], top[2]+axi[2], 1);

        VEC4 nor1 = VEC4(origin[0], origin[1], origin[2], 1); 
        VEC4 ndi1 = VEC4(dir[0], dir[1], dir[2], 1); 
        // TO NOTE: we have to transform the vector from (0,0,0)
        // for the ray direction, not the vector from the ray origin! 
        // this is how matrix multiplication works! 
        ndi1 = ndi1 + nor1 + VEC4(0,0,0,-1); 

        nor1 = M*nor1;
        ndi1 = M*ndi1;

        VEC3 nor2 = nor1.head<3>(); 
        VEC3 ndi2 = ndi1.head<3>(); 

        // Since we transformed the full vector, we now need to retrieve the vector that 
        // represents the ray direction
        // notice this vector is still normalized since we applied 
        // an orthographic transform, and the ray dir then was normalized
        ndi2 = ndi2 - nor2;

        // SECOND, solve for points of intersection with cylinder pointing up the y axis
        // our equation is x^2 + z^2 = r^2
        // (px + t vx)^2 + (pz + t vz)^2 = r^2
        double a = ndi2[0]*ndi2[0] + ndi2[2]*ndi2[2];
        double b = 2*(ndi2[0]*nor2[0] + nor2[2]*ndi2[2]);
        double c = nor2[0]*nor2[0] + nor2[2]*nor2[2] - radius*radius;

        double disc = b*b - 4*a*c;
        if (disc < 0) {
            return -1.;
        }
        if (disc == 0) {
            return -1*b/(2*a);
        }
        double t1 = (-1*b + sqrt(disc))/(2*a);
        double t2 = (-1*b - sqrt(disc))/(2*a);

        VEC3 p1 = nor2 + t1*ndi2;
        VEC3 p2 = nor2 + t2*ndi2;

        if (t1 < t2) { // !! are these checks redundant??
            if (p1[1] >= 0 && p1[1] <= length) {
                return t1;
            }
            if (p2[1] >= 0 && p2[1] <= length) {
                return t2;
            }
            return -1;
        }
        if (p2[1] >= 0 && p2[1] <= length) {
            return t2;
        }
        if (p1[1] >= 0 && p1[1] <= length) {
            return t1;
        }
        return -1;
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

