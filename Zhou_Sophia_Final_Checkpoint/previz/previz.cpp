//////////////////////////////////////////////////////////////////////////////////
// This is a front end for a set of viewer clases for the Carnegie Mellon
// Motion Capture Database: 
//    
//    http://mocap.cs.cmu.edu/
//
// The original viewer code was downloaded from:
//
//   http://graphics.cs.cmu.edu/software/mocapPlayer.zip
//
// where it is credited to James McCann (Adobe), Jernej Barbic (USC),
// and Yili Zhao (USC). There are also comments in it that suggest
// and Alla Safonova (UPenn) and Kiran Bhat (ILM) also had a hand in writing it.
//
//////////////////////////////////////////////////////////////////////////////////
#include "SETTINGS.h"
#include "skeleton.h"
#include "displaySkeleton.h"
#include "motion.h"
#include "helpers.h"

using namespace std;

// Stick-man classes
DisplaySkeleton displayer;    
Skeleton* skeleton;
Motion* motion;

int windowWidth = 640;
int windowHeight = 480;

// looking down positive x
// to our right is positive z
// up is positive y
VEC3 eye(-3.5, 2, 3.5);
// VEC3 eye(-3, 0, 0);
VEC3 lookingAt; //(5, 0.5, 1);
VEC3 up(0,1,0);

// scene geometry
vector<Primitive*> primitives;
vector<Light*> lights;

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
void writePPM(const string& filename, int& xRes, int& yRes, const float* values)
{
  int totalCells = xRes * yRes;
  unsigned char* pixels = new unsigned char[3 * totalCells];
  for (int i = 0; i < 3 * totalCells; i++)
    pixels[i] = values[i];

  FILE *fp;
  fp = fopen(filename.c_str(), "wb");
  if (fp == NULL)
  {
    cout << " Could not open file \"" << filename.c_str() << "\" for writing." << endl;
    cout << " Make sure you're not trying to write from a weird location or with a " << endl;
    cout << " strange filename. Bailing ... " << endl;
    exit(0);
  }

  fprintf(fp, "P6\n%d %d\n255\n", xRes, yRes);
  fwrite(pixels, 1, totalCells * 3, fp);
  fclose(fp);
  delete[] pixels;
}











// returns a boolean array populated with shadow data
// !!! why is this returning in shadow when behind the sphere?
// need to cull when behind sphere

// !!!! since we aren't culling shadows, all behind the person are shadows, and that 
// is what is being displayed! this means our time t has to be correct when we return them from 
// find intersect.
// !!! the thing is, why isn't culling working for spheres??? 

double* findShadow(VEC3 origin) {
  VEC3 p = origin;

  double* inShadow = new double[lights.size()];

  int numPointLights = 16;

  for (int n=0; n<lights.size(); n++) { // initing inShadow
    inShadow[n]=0.;
  }

  for (int n=0; n<primitives.size(); n++) { // begin for
    VEC3 p_light;
    Primitive* this_s=primitives[n];

    for (int m=0; m<lights.size(); m++) { // begin for
      
      // generate 32 random point lights
      for(int k=0; k<numPointLights; k++) { // begin for
        // unif sample between (-2,2)
        VEC3 randVector = VEC3( ((double)rand() / RAND_MAX)*4 - 2, 0, ((double)rand() / RAND_MAX)*4 - 2 );

        p_light = lights[m]->pos + randVector - p;
        normalize(p_light);

        IntersectionPoint intersectionPoint = this_s->findIntersect(p, p_light);
        double t = intersectionPoint.time;
        if(t > 0 && magnitude(t*p_light) < magnitude(lights[m]->pos - p)) {
          inShadow[m] += 1;
        } // endif

      } // endfor 

    } // endfor
  } // endfor

  for(int i = 0; i<lights.size(); i++) {
    inShadow[i] = ((double)numPointLights - inShadow[i]) / (double)numPointLights;
  }

  return inShadow;
}

















//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
void rayColor(const VEC3& rayPos, const VEC3& rayDir, VEC3& pixelColor) 
{
  // pixelColor = VEC3(1,1,1);
  pixelColor = VEC3(0,0,0);

  // look for intersections
  int hitID = -1;
  float tMinFound = FLT_MAX;

  // change so it is not just spheres, 
  // but also triangles and cylinders
  for (int y = 0; y < primitives.size(); y++)
  {
    float tMin = FLT_MAX;
    
    Primitive* p = primitives[y];
    VEC3 a = rayPos;
    VEC3 b = rayDir;

    IntersectionPoint intersectionPoint = p->findIntersect(a, b);
    tMin = (float)(intersectionPoint.time);

    if (tMin > 0)
    { 
      // is the closest so far?
      if (tMin < tMinFound)
      {
        tMinFound = tMin;

        // coloring the surface
        VEC3 surf_norm = p->findSurfNorm(rayPos, rayDir, tMin);
        VEC3 surf_pos = p->findSurfPos(rayPos, rayDir, tMin);
        VEC3 pt = rayPos + (tMin-0.001)*rayDir;
        intersectionPoint.point = pt;

        // checking if in shadow and coloring
        double* inShadow = findShadow(pt);

        pixelColor = VEC3(0,0,0);

        // texture mapping
        VEC3 pointColor;

        pointColor = p->getColor(intersectionPoint);

        for(int n=0;n<lights.size();n++) {
          pixelColor += inShadow[n]*lights[n]->totCol(pointColor, surf_norm, surf_pos, eye);
        }

        // cout << pixelColor << endl << endl;

        hitID = y;
      }
    }
  }
  
  // No intersection, return white
  if (hitID == -1)
    return;

}

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
float clamp(float value)
{
  if (value < 0.0)      return 0.0;
  else if (value > 1.0) return 1.0;
  return value;
}

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
// this runs for each of the 300 frames that we need to render
void renderImage(int& xRes, int& yRes, const string& filename) 
{
  // allocate the final image
  const int totalCells = xRes * yRes;
  float* ppmOut = new float[3 * totalCells];

  // compute image plane
  // do camera motion by changing lookingat
  const float halfY = (lookingAt - eye).norm() * tan(45.0f / 360.0f * M_PI);
  const float halfX = halfY * 4.0f / 3.0f;

  const VEC3 cameraZ = (lookingAt - eye).normalized();
  const VEC3 cameraX = up.cross(cameraZ).normalized();
  const VEC3 cameraY = cameraZ.cross(cameraX).normalized();

  for (int y = 0; y < yRes; y++) 
    for (int x = 0; x < xRes; x++) 
    {
      const float ratioX = 1.0f - x / float(xRes) * 2.0f;
      const float ratioY = 1.0f - y / float(yRes) * 2.0f;
      const VEC3 rayHitImage = lookingAt + 
                               ratioX * halfX * cameraX +
                               ratioY * halfY * cameraY;
      const VEC3 rayDir = (rayHitImage - eye).normalized();

      // get the color
      VEC3 color;
      rayColor(eye, rayDir, color); 

      // cout << color << endl << endl;

      // set, in final image
      ppmOut[3 * (y * xRes + x)] = clamp(color[0]) * 255.0f;
      ppmOut[3 * (y * xRes + x) + 1] = clamp(color[1]) * 255.0f;
      ppmOut[3 * (y * xRes + x) + 2] = clamp(color[2]) * 255.0f;
    }
  writePPM(filename, xRes, yRes, ppmOut);

  delete[] ppmOut;
}





























//////////////////////////////////////////////////////////////////////////////////
// Load up a new motion captured frame
//////////////////////////////////////////////////////////////////////////////////
void setSkeletonsToSpecifiedFrame(int frameIndex)
{
  if (frameIndex < 0)
  {
    printf("Error in SetSkeletonsToSpecifiedFrame: frameIndex %d is illegal.\n", frameIndex);
    exit(0);
  }
  if (displayer.GetSkeletonMotion(0) != NULL)
  {
    int postureID;
    if (frameIndex >= displayer.GetSkeletonMotion(0)->GetNumFrames())
    {
      cout << " We hit the last frame! You might want to pick a different sequence. " << endl;
      postureID = displayer.GetSkeletonMotion(0)->GetNumFrames() - 1;
    }
    else 
      postureID = frameIndex;
    displayer.GetSkeleton(0)->setPosture(* (displayer.GetSkeletonMotion(0)->GetPosture(postureID)));
  }
}

//////////////////////////////////////////////////////////////////////////////////
// Build a list of spheres in the scene
//////////////////////////////////////////////////////////////////////////////////
void buildScene()
{
  primitives.clear();

  displayer.ComputeBonePositions(DisplaySkeleton::BONES_AND_LOCAL_FRAMES);

  // retrieve all the bones of the skeleton
  vector<MATRIX4>& rotations = displayer.rotations();
  vector<MATRIX4>& scalings  = displayer.scalings();
  vector<VEC4>& translations = displayer.translations();

  // get the pelvis of the skeleton
  // !! can also choose a different translations
  // like translations[5] or something
  lookingAt = VEC3(translations[1][0], translations[1][1], translations[1][2]);
  // lookingAt = VEC3(1,0,0);

  vector<float>& lengths = displayer.lengths();

  // build a sphere list, but skip the first bone, 
  // it's just the origin
  int totalBones = rotations.size();
  
  
  
  
   
  for (int x = 1; x < totalBones; x++)
  {
    MATRIX4& rotation = rotations[x];
    MATRIX4& scaling = scalings[x];
    VEC4& translation = translations[x];

    // get the endpoints of the cylinder
    VEC4 leftVertex(0,0,0,1);
    VEC4 rightVertex(0,0,lengths[x],1);

    leftVertex = rotation * scaling * leftVertex + translation;
    rightVertex = rotation * scaling * rightVertex + translation;

    VEC3 lef = leftVertex.head<3>();
    VEC3 rig = rightVertex.head<3>();

    double sphereRad = 0.1; // !! changed! its like the thumb people

    Cylinder* cyl = new Cylinder(lef, rig, sphereRad, VEC3(1,0,0));
    Sphere* sph1 = new Sphere(lef, sphereRad, VEC3(1,0,0));
    Sphere* sph2 = new Sphere(rig, sphereRad, VEC3(1,0,0));

    primitives.push_back(cyl);
    primitives.push_back(sph1);
    primitives.push_back(sph2);
  }

  // populate primitives with big triangle
  // and big sphere
  Triangle* tri1 = new Triangle(VEC3(20,0,0), VEC3(20,0,20),  VEC3(0,0,0), VEC3(0,0,1), "DarkGrassTexture.ppm");
  Triangle* tri2 = new Triangle(VEC3(0,0,20), VEC3(20,0,20),  VEC3(0,0,0), VEC3(0,0,1), "DarkGrassTexture.ppm");

  Triangle* tri3 = new Triangle(VEC3(-20,0,0), VEC3(-20,0,20),  VEC3(0,0,0), VEC3(0,0,1), "DarkGrassTexture.ppm");
  Triangle* tri4 = new Triangle(VEC3(0,0,20), VEC3(-20,0,20),  VEC3(0,0,0), VEC3(0,0,1), "DarkGrassTexture.ppm");

  Triangle* tri5 = new Triangle(VEC3(20,0,0), VEC3(20,0,-20),  VEC3(0,0,0), VEC3(0,0,1), "DarkGrassTexture.ppm");
  Triangle* tri6 = new Triangle(VEC3(0,0,-20), VEC3(0,0,0), VEC3(20,0,-20),  VEC3(0,0,1), "DarkGrassTexture.ppm");

  Triangle* tri7 = new Triangle(VEC3(-20,0,0), VEC3(-20,0,-20),  VEC3(0,0,0), VEC3(0,0,1), "DarkGrassTexture.ppm");
  Triangle* tri8 = new Triangle(VEC3(0,0,-20), VEC3(-20,0,-20),  VEC3(0,0,0), VEC3(0,0,1), "DarkGrassTexture.ppm");
  
  primitives.push_back(tri1);
  primitives.push_back(tri2);
  primitives.push_back(tri3);
  primitives.push_back(tri4);
  primitives.push_back(tri5);
  primitives.push_back(tri6);
  primitives.push_back(tri7);
  primitives.push_back(tri8);


  Sphere* sph = new Sphere(VEC3(0.1,0.1,-.201), 0.4, VEC3(0,1,0), "SteelTexture.ppm"); // i have been spending a lot of time trying to figure out the nice positioning for the sphere
  primitives.push_back(sph);


  // for locations in the scene 
  // Sphere* xax = new Sphere(VEC3(1,0,0), 0.25, VEC3(1,0,0));
  // Sphere* yax = new Sphere(VEC3(0,1,0), 0.25, VEC3(0,1,0));
  // Sphere* zax = new Sphere(VEC3(0,0,1), 0.25, VEC3(0,0,1));
  // primitives.push_back(xax);
  // primitives.push_back(yax);
  // primitives.push_back(zax);

  // it got darker after i turned it to 12
  Sphere* TESTSKYBOX = new Sphere(VEC3(0,2,0), 15, VEC3(1,0,0), "SkyNight.ppm"); // the normals might be different for inside of sphere!!!
  primitives.push_back(TESTSKYBOX);
  
}



// !!! function for debugging 
void buildScene2() {
  primitives.clear();

  displayer.ComputeBonePositions(DisplaySkeleton::BONES_AND_LOCAL_FRAMES);

  // retrieve all the bones of the skeleton
  vector<MATRIX4>& rotations = displayer.rotations();
  vector<MATRIX4>& scalings  = displayer.scalings();
  vector<VEC4>& translations = displayer.translations();

  // get the pelvis of the skeleton
  // !! can also choose a different translations
  // like translations[5] or something
  // lookingAt = VEC3(translations[1][0], translations[1][1], translations[1][2]);
  lookingAt = VEC3(1,0,0);

  vector<float>& lengths = displayer.lengths();

  // build a sphere list, but skip the first bone, 
  // it's just the origin
  int totalBones = rotations.size();
  
  
  
  
   
  for (int x = 1; x < totalBones; x++)
  {
    MATRIX4& rotation = rotations[x];
    MATRIX4& scaling = scalings[x];
    VEC4& translation = translations[x];

    // get the endpoints of the cylinder
    VEC4 leftVertex(0,0,0,1);
    VEC4 rightVertex(0,0,lengths[x],1);

    leftVertex = rotation * scaling * leftVertex + translation;
    rightVertex = rotation * scaling * rightVertex + translation;

    VEC3 lef = leftVertex.head<3>();
    VEC3 rig = rightVertex.head<3>();

    double sphereRad = 0.1; // !! changed! its like the thumb people

    Cylinder* cyl = new Cylinder(lef, rig, sphereRad, VEC3(1,0,0));
    Sphere* sph1 = new Sphere(lef, sphereRad, VEC3(1,0,0));
    Sphere* sph2 = new Sphere(rig, sphereRad, VEC3(1,0,0));

    primitives.push_back(cyl);
    primitives.push_back(sph1);
    primitives.push_back(sph2);
  }


  Sphere* sph = new Sphere(VEC3(0.1,0.1,-.2), 0.4, VEC3(0,1,0), "SteelTexture.ppm"); // i have been spending a lot of time trying to figure out the nice positioning for the sphere
  primitives.push_back(sph);


  // for locations in the scene 
  // Sphere* xax = new Sphere(VEC3(1,0,0), 0.25, VEC3(1,0,0));
  // Sphere* yax = new Sphere(VEC3(0,1,0), 0.25, VEC3(0,1,0));
  // Sphere* zax = new Sphere(VEC3(0,0,1), 0.25, VEC3(0,0,1));
  // primitives.push_back(xax);
  // primitives.push_back(yax);
  // primitives.push_back(zax);
}

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv)
{
  string skeletonFilename("14.asf");
  string motionFilename("14_30.amc");

  // add lights
  Light* lit1 = new Light(VEC3(-4,4,-4), VEC3(1,1,1), 1, 10);
  Light* lit2 = new Light(VEC3(4,4,4), VEC3(1,1,1), 1, 10);
  lights.push_back(lit1);
  lights.push_back(lit2);

  // load up skeleton stuff
  skeleton = new Skeleton(skeletonFilename.c_str(), MOCAP_SCALE);
  skeleton->setBasePosture();
  displayer.LoadSkeleton(skeleton);

  // load up the motion
  motion = new Motion(motionFilename.c_str(), MOCAP_SCALE, skeleton);
  displayer.LoadMotion(motion);
  skeleton->setPosture(*(displayer.GetSkeletonMotion(0)->GetPosture(0)));

  int start=0, stop=2400;
  if (argc == 3) {
    start=stoi(argv[1]);
    stop =stoi(argv[2]);
  }
  // Note we're going 8 frames at a time, otherwise the animation
  // is really slow.
  for (int x = start; x < stop; x += 8)
  { 
    // to start the animation later
    // int y = x + 2400;
    // setSkeletonsToSpecifiedFrame(y);
    setSkeletonsToSpecifiedFrame(x);
    buildScene();
    // buildScene2();

    char buffer[256];
    snprintf(buffer, 256, "./frames/frame.%04i.ppm", x / 8);
    renderImage(windowWidth, windowHeight, buffer);
    cout << "Rendered " + to_string(x / 8) + " frames" << endl;
  }

  return 0;
}