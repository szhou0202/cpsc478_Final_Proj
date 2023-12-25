// what do we want to include here

#ifndef HELPERS_H
#define HELPERS_H

#include "normalize.h"
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include <float.h>

using namespace std;




// struct for intersection point
struct IntersectionPoint {
    VEC3 point; // !! we need to populate this on our own
    double time;
    double beta; // keep track of whether or not triangle intersect
    double gamma;
};






// struct for texture
struct Texture {
    string filename;
    int width;
    int height;
    float* values;
};






// !!! lol this is kinda jank
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void readPPMP3(const string& filename, int& xRes, int& yRes, float*& values)
{
    ifstream fin(filename);

    string s;
    int x,y;
    float n;
    fin >> s >> x >> y >> n;

    xRes = x;
    yRes = y;

    // cout << xRes << ", " << yRes << endl;

    // do division here
    float* pixels = new float[3*xRes*yRes];
    for(int i = 0; i < xRes * yRes * 3; i++) {
        fin >> pixels[i];
        pixels[i] /= n;
    }
    values = pixels;

    fin.close(); // this shouldn't be the problem though
}

void readPPM(const string& filename, int& xRes, int& yRes, float*& values)
{
  // try to open the file
  FILE *fp;
  fp = fopen(filename.c_str(), "rb");
  if (fp == NULL)
  {
    cout << " Could not open file \"" << filename.c_str() << "\" for reading." << endl;
    cout << " Make sure you're not trying to read from a weird location or with a " << endl;
    cout << " strange filename. Bailing ... " << endl;
    exit(0);
  }

  // get the dimensions
  unsigned char newline;
  fscanf(fp, "P6\n%d %d\n255%c", &xRes, &yRes, &newline);
  if (newline != '\n') {
    cout << " The header of " << filename.c_str() << " may be improperly formatted." << endl;
    cout << " The program will continue, but you may want to check your input. " << endl;
  }
  int totalCells = xRes * yRes;

  // grab the pixel values
  unsigned char* pixels = new unsigned char[3 * totalCells];
  fread(pixels, 1, totalCells * 3, fp);

  // copy to a nicer data type
  values = new float[3 * totalCells];
  for (int i = 0; i < 3 * totalCells; i++)
    values[i] = pixels[i];

  // clean up
  delete[] pixels;
  fclose(fp);
  cout << " Read in file " << filename.c_str() << endl;
}







// primitive
class Primitive {
public:
    VEC3 color;
    bool reflect; 
    bool refract;
    double refract_air;
    double refract_glass;

    // for texture mapping
    bool hasTexture = false; // lol i switched from snake case to camel case 
    Texture texture;

    Primitive(){};

    // Beautiful, we have primitives working
    virtual IntersectionPoint findIntersect(VEC3 origin, VEC3 dir){
        return {VEC3(0,0,0), -1, -1, -1};
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
    virtual VEC3 getColor(IntersectionPoint pt) {
        return color; // !! make the check in here
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
    Sphere(VEC3 c, double r, VEC3 co, const string fn) {
        center = c;
        radius = r;
        color = co;

        // read in the ppm file for the texture
        // apparently the width and height are stored in the header
        hasTexture = true;

        int width;
        int height;
        string filename = fn;
        float* values;

        readPPMP3(fn, width, height, values);

        // this is scary
        texture.filename = filename;
        texture.width = width;
        texture.height = height;
        texture.values = values;
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

    IntersectionPoint findIntersect(VEC3 origin, VEC3 dir) {
        // calculate discriminant 
        VEC3 ominc = origin-center;
        double dd = dir.dot(dir);
        double disc = pow(dir.dot(ominc),2) - dd*(ominc.dot(ominc) - pow(radius,2));
        
        IntersectionPoint pt;
        if (disc < 0) {
            pt.time = -1;
            return pt;
        }
        if (disc == 0) {
            pt.time = (-1)*dir.dot(ominc) / dd;
            return pt;
        }

        double t1 = (((-1)*dir).dot(ominc) + sqrt(disc))/dd;
        double t2 = (((-1)*dir).dot(ominc) - sqrt(disc))/dd;

        if (t1 < t2) {
            if(t1 > 0) {
                pt.time = t1;
            }
            else if(t2 > 0) {
                pt.time = t2;
            }
            else {
                pt.time = -1;
            }
        }
        else if(t2 < t1) {
            if(t2 > 0) {
                pt.time = t2;
            }
            else if(t1 > 0){
                pt.time = t1;
            }
            else {
                pt.time = -1;
            }
        }
        else {
            pt.time = -1;
        }
        return pt;
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

    // function to get the texture color 
    VEC3 getColor(IntersectionPoint pt) {
        if(!hasTexture) {
            return color;
        }

        VEC3 norm = pt.point - center;
        normalize(norm);
        double x = norm[0];
        double y = norm[1];
        double z = norm[2];

        // formula from online
        // u = 0.5 + arctan2(dz, dx) / (2*pi)
        // v = 0.5 - arcsin(dy) / pi
        double tempU = 0.5 + (atan2(z, x) / (2*M_PI));
        double tempV = 0.5 - (asin(y)/M_PI);

        int u = (int)(tempU*texture.width);
        int v = (int)(tempV*texture.height);

        // !! do i need division by 225 here??
        VEC3 textureColor = VEC3(texture.values[3 * (v * texture.width + u)], 
            texture.values[3*(v * texture.width + u) +1], 
            texture.values[3*(v * texture.width + u) +2]);

        return textureColor;
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

    Triangle(VEC3 v1, VEC3 v2, VEC3 v3, VEC3 co, string fn) {
        vertex1 = v1;
        vertex2 = v2;
        vertex3 = v3;
        u = vertex2-vertex1;
        v = vertex3-vertex1;
        color = co;

        // has texture
        hasTexture = true;

        int width;
        int height;
        string filename = fn;
        float* values;

        readPPMP3(fn, width, height, values);

        texture.filename = filename;
        texture.width = width;
        texture.height = height;
        texture.values = values;
    };

    IntersectionPoint findIntersect(VEC3 origin, VEC3 dir) {
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

        IntersectionPoint pt;
        pt.beta = beta;
        pt.gamma = gamma;

        if (t < 0)
            pt.time = -1;
        else if (gamma < 0 || gamma > 1) 
            pt.time = -1;
        else if (beta < 0 || beta > 1-gamma) 
            pt.time = -1;
        else {
            pt.time = t;
        }
        return pt;
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

        // flip it
        if(toNorm[1] < 0) {
            toNorm[0] = -toNorm[0];
            toNorm[1] = -toNorm[1];
            toNorm[2] = -toNorm[2];
        }
        
        return toNorm;
    }

    // call the a point (0,0)
    // the b point (1,0)
    // the c point (0,1)
    VEC3 getColor(IntersectionPoint pt) {
        if(!hasTexture) {
            return color;
        }
        double beta = pt.beta;
        double gamma = pt.gamma;

        int u = abs((int)(beta*texture.width));
        int v = abs((int)(gamma*texture.height));

        VEC3 textureColor = VEC3(texture.values[3*(v*texture.width + u)], 
            texture.values[3*(v*texture.width + u) +1], 
            texture.values[3*(v*texture.width + u) +2]);


        // cout << textureColor << endl; // u and v look correct....

        return textureColor;
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
    MATRIX4 create_M(VEC3 origin, VEC3 dir) {
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

        return M;
    }

    IntersectionPoint findIntersect(VEC3 origin, VEC3 dir) {

        // FIRST, apply change of basis to the ray 
        // find change of basis array, normalized vectors
        // this is very similar to the MCAM matrix! 
        MATRIX4 M = create_M(origin, dir);

        VEC4 nor1 = VEC4(origin[0], origin[1], origin[2], 1); 
        VEC4 ndi1 = VEC4(dir[0], dir[1], dir[2], 1); 
        // TO NOTE: we have to transform the vector from (0,0,0)
        // for the ray direction, not the vector from the ray origin! 
        // this is how matrix multiplication works! 
        ndi1 = ndi1 + nor1 + VEC4(0,0,0,-1); 

        nor1 = M*nor1; // !!! too many functions defined for the same operator and operands * ?
        ndi1 = M*ndi1; // !!! why? this was working before

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
        
        IntersectionPoint pt;

        double disc = b*b - 4*a*c;
        if (disc < 0) {
            pt.time = -1;
            return pt;
        }
        if (disc == 0) {
            pt.time = -1*b/(2*a);
            return pt;
        }
        double t1 = (-1*b + sqrt(disc))/(2*a);
        double t2 = (-1*b - sqrt(disc))/(2*a);

        VEC3 p1 = nor2 + t1*ndi2;
        VEC3 p2 = nor2 + t2*ndi2;
        
        if (t1 < t2) { // !! are these checks redundant??
            if (p1[1] >= 0 && p1[1] <= length) {
                pt.time = t1;
                return pt;
            }
            if (p2[1] >= 0 && p2[1] <= length) {
                pt.time = t2;
                return pt;
            }
            pt.time = -1;
            return pt;
        }
        if (p2[1] >= 0 && p2[1] <= length) {
            pt.time = t2;
            return pt;
        }
        if (p1[1] >= 0 && p1[1] <= length) {
            pt.time = t1;
            return pt;
        }
        pt.time = -1;
        return pt;
    }

    // find coords of surface point
    VEC3 findSurfPos(VEC3 eye, VEC3 dir, double t) {
        return eye + t*dir;
    }

    // find norm at a point on cylinder
    VEC3 findSurfNorm(VEC3 eye, VEC3 dir, double t) { // !!! is eye same thing as origin?
        // do the same transform as before?

        // we just need to find the closest vector 
        // from this point to the axis

        // if we transform it, then it just becomes 
        // (x,y,z) - (x,0,z) = (0,y,0)
        // and finally transform (0,y,0) back and that is our normal

        // so, i should probably modularize the matrix transform
        
        MATRIX4 M = create_M(eye, dir);
        
        // this should simply be the point of intersection
        VEC4 ndi1 = VEC4(eye[0] + t*dir[0], eye[1] + t*dir[1], eye[2] + t*dir[2], 1); 

        ndi1 = M*ndi1;
        ndi1[3] = 1;

        VEC4 nax3 = VEC4(0, ndi1[1], 0,1); // !!!! i was doing it correctly all along :sob:

        // our norm
        M = M.inverse().eval();
        VEC3 nin4 = (M*ndi1).head<3>();
        VEC3 nax4 = (M*nax3).head<3>();

        VEC3 temp = nin4-nax4; // !!!! why was it being scaled by a factor of 0.5????
        normalize(temp);
        return temp;
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

        double t = reflected_light.dot(view_ray);
        double temp = intensity*(0 >= t ? 0 : pow(t,phong));
        return temp*VEC3(color[0]*surf_color[0],color[1]*surf_color[1],color[2]*surf_color[2]);
    };

    VEC3 totCol(VEC3 surf_color, VEC3 surf_norm, VEC3 surf_position, VEC3 eye) {
        VEC3 dif = diffuseCol(surf_color, surf_norm, surf_position);
        VEC3 spec = specCol(surf_color, surf_norm, surf_position, eye);
        return dif+spec;
    };
};



#endif

