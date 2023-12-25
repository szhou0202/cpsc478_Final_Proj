#include <iostream>
#include <fstream>
#include <cstdio>

using namespace std;


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

    cout << xRes << ", " << yRes << endl;

    // do division here
    float* pixels = new float[3*xRes*yRes];
    for(int i = 0; i < xRes * yRes * 3; i++) {
        fin >> pixels[i];
        // pixels[i] /= n;
    }
    values = pixels;
}


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



class Person {
public:
    virtual void print() {
        cout << "hello\n";
    }
};

class Friend : public Person {
    void print() {
        cout << "friend\n";
    }
};





int main() {
    

    return 0;
}