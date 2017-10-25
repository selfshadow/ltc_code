#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cstdint>
using namespace std;

#include "CImg.h"
using namespace cimg_library;

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

void filter(CImg<float>& imageInput, CImg<float>& imageOutput, const int level, const int Nlevels)
{
    // distance to texture plane
    // in shader: LOD = log(powf(2.0f, Nlevels - 1.0f) * dist) / log(3)
    const float dist = powf(3.0f, level) / powf(2.0f, Nlevels - 1.0f);

    // filter size
    // at distance 1 ~= Gaussian of std 0.75
    const float filterStd = 0.75f * dist * imageInput.width();

    cout << "level " << level << endl;
    cout << "distance to texture plane = " << dist << endl;
    cout << "filterStd = " << filterStd << endl << endl;

    CImg<float> tmp(imageInput.width(), imageInput.height(), 1, 4);

    for (int j = 0; j < imageInput.height(); ++j)
    for (int i = 0; i < imageInput.width();  ++i)
    {
        tmp(i, j, 0, 0) = imageInput(i, j, 0, 0);
        tmp(i, j, 0, 1) = imageInput(i, j, 0, 1);
        tmp(i, j, 0, 2) = imageInput(i, j, 0, 2);
        tmp(i, j, 0, 3) = 1.0f;
    }

    tmp.blur(filterStd, filterStd, filterStd, false);

    // renormalise based on alpha
    for (int j = 0; j < imageInput.height(); ++j)
    for (int i = 0; i < imageInput.width();  ++i)
    {
        float alpha = tmp(i, j, 0, 3);
        for (int k = 0; k < tmp.spectrum(); ++k)
          tmp(i, j, 0, k) /= alpha;
    }

    // rescale image
    imageOutput = tmp.resize(imageOutput, 5); // 5 = cubic interpolation
}

// Adapted from https://www.khronos.org/registry/OpenGL/extensions/EXT/EXT_texture_shared_exponent.txt
////////
#define RGB9E5_EXPONENT_BITS          5
#define RGB9E5_MANTISSA_BITS          9
#define RGB9E5_EXP_BIAS               15
#define RGB9E5_MAX_VALID_BIASED_EXP   31

#define MAX_RGB9E5_EXP               (RGB9E5_MAX_VALID_BIASED_EXP - RGB9E5_EXP_BIAS)
#define RGB9E5_MANTISSA_VALUES       (1<<RGB9E5_MANTISSA_BITS)
#define MAX_RGB9E5_MANTISSA          (RGB9E5_MANTISSA_VALUES-1)
#define MAX_RGB9E5                   (((float)MAX_RGB9E5_MANTISSA)/RGB9E5_MANTISSA_VALUES * (1<<MAX_RGB9E5_EXP))
#define EPSILON_RGB9E5               ((1.0/RGB9E5_MANTISSA_VALUES) / (1<<RGB9E5_EXP_BIAS))

typedef struct
{
  uint32_t mantissa : 23;
  uint32_t biasedexponent : 8;
  uint32_t negative : 1;
} BitsOfIEEE754;

typedef union
{
  uint32_t raw;
  float value;
  BitsOfIEEE754 field;
} float754;

typedef struct
{
    uint32_t r : RGB9E5_MANTISSA_BITS;
    uint32_t g : RGB9E5_MANTISSA_BITS;
    uint32_t b : RGB9E5_MANTISSA_BITS;
    uint32_t biasedexponent : RGB9E5_EXPONENT_BITS;
} BitsOfRGB9E5;

typedef union
{
    uint32_t raw;
    BitsOfRGB9E5 field;
} rgb9e5;

float ClampRange_for_rgb9e5(float x)
{
    if (x > 0.0)
    {
        if (x >= MAX_RGB9E5)
            return MAX_RGB9E5;
    }
    else
    {
        /* NaN gets here too since comparisons with NaN always fail! */
        return 0.0f;
    }

    return x;
}

float Max(float a, float b)
{
    return a >= b ? a : b;
}

float MaxOf3(float x, float y, float z)
{
    return Max(Max(x, y), z);
}

/* Ok, FloorLog2 is not correct for the denorm and zero values, but we
   are going to do a max of this value with the minimum rgb9e5 exponent
   that will hide these problem cases. */
int FloorLog2(float x)
{
  float754 f;

  f.value = x;
  return (f.field.biasedexponent - 127);
}

int Max(int x, int y)
{
  if (x > y) {
    return x;
  } else {
    return y;
  }
}

rgb9e5 float3_to_rgb9e5(const float r, const float g, const float b)
{
  rgb9e5 retval;
  float maxrgb;
  int rm, gm, bm;
  float rc, gc, bc;
  int exp_shared;
  double denom;

  rc = ClampRange_for_rgb9e5(r);
  gc = ClampRange_for_rgb9e5(g);
  bc = ClampRange_for_rgb9e5(b);

  maxrgb = MaxOf3(rc, gc, bc);
  exp_shared = Max(-RGB9E5_EXP_BIAS-1, FloorLog2(maxrgb)) + 1 + RGB9E5_EXP_BIAS;
  assert(exp_shared <= RGB9E5_MAX_VALID_BIASED_EXP);
  assert(exp_shared >= 0);
  /* This pow function could be replaced by a table. */
  denom = pow(2, exp_shared - RGB9E5_EXP_BIAS - RGB9E5_MANTISSA_BITS);

  int maxm = (int) floor(maxrgb / denom + 0.5f);
  if (maxm == MAX_RGB9E5_MANTISSA+1) {
    denom *= 2;
    exp_shared += 1;
    assert(exp_shared <= RGB9E5_MAX_VALID_BIASED_EXP);
  } else {
    assert(maxm <= MAX_RGB9E5_MANTISSA);
  }

  rm = (int) floor(rc / denom + 0.5);
  gm = (int) floor(gc / denom + 0.5);
  bm = (int) floor(bc / denom + 0.5);

  assert(rm <= MAX_RGB9E5_MANTISSA);
  assert(gm <= MAX_RGB9E5_MANTISSA);
  assert(bm <= MAX_RGB9E5_MANTISSA);
  assert(rm >= 0);
  assert(gm >= 0);
  assert(bm >= 0);

  retval.field.r = rm;
  retval.field.g = gm;
  retval.field.b = bm;
  retval.field.biasedexponent = exp_shared;

  return retval;
}
////////

int main(int argc, char* argv[])
{
    // Skip executable argument
    argc--;
    argv++;

    if (argc < 1)
    {
        printf("Syntax: <input file>\n");
        return 0;
    }

    string filenameInput(argv[0]);
    size_t pos = filenameInput.find_last_of(".");
    string filename  = filenameInput.substr(0, pos);
    string extension = filenameInput.substr(pos + 1, string::npos);

    // input image
    int x, y, n;
    float* data = stbi_loadf(filenameInput.c_str(), &x, &y, &n, 3);

    int offset = 0;
    CImg<float> imageInput(x, y, 1, 3);
    for (int j = 0; j < imageInput.height(); ++j)
    for (int i = 0; i < imageInput.width();  ++i)
    {
        for (int k = 0; k < imageInput.spectrum(); ++k)
            imageInput(i, j, 0, k) = data[offset++];
    }

    // filtered levels
    unsigned int Nlevels;
    for (Nlevels = 1; (imageInput.width() >> Nlevels) > 0; ++Nlevels);

    uint32_t* dataOut = new uint32_t[x * y];

    int hdrSize = x * y * 3 * Nlevels;
    float* hdrOut = new float[hdrSize];

    int hdrOffset = 0;

    // borders
    for (unsigned int level = 0; level < Nlevels; ++level)
    {
        stringstream filenameOutput(stringstream::in | stringstream::out);
        filenameOutput << filename << "_filtered_" << level << ".png";

        cout << "processing file " << filenameOutput.str() << endl;
        unsigned int width = imageInput.width();// >> level;
        unsigned int height = imageInput.height();// >> level;

        if (width <= 0 || height <= 0)
            break;

        CImg<float> imageOutput(width, height, 1, 4);

        filter(imageInput, imageOutput, level, Nlevels);

        offset = 0;
        for (int j = 0; j < imageOutput.height(); ++j)
        for (int i = 0; i < imageOutput.width();  ++i)
        {
            float r = imageOutput(i, j, 0, 0);
            float g = imageOutput(i, j, 0, 1);
            float b = imageOutput(i, j, 0, 2);

            rgb9e5 value = float3_to_rgb9e5(r, g, b);
            dataOut[offset++] = value.raw;

            hdrOut[hdrOffset++] = r;
            hdrOut[hdrOffset++] = g;
            hdrOut[hdrOffset++] = b;
        }

        stbi_write_png(filenameOutput.str().c_str(), width, height, 4, dataOut, 0);

        //stringstream hdrFname(stringstream::in | stringstream::out);
        //hdrFname << "mip" << level << ".hdr";
        //stbi_write_hdr(hdrFname.str().c_str(), width, height, 3, hdrOut);
    }

    ofstream file("array.js");

    file << "var g_texture_array = [" << endl;

    data = hdrOut;
    for (int i = 0; i < hdrSize/16; i++)
    {
        for (int j = 0; j < 16; j++)
          file << *data++ << ", ";
        file << endl;
    }

    file << "];" << endl;
    file.close();

    delete[] dataOut;
    delete[] hdrOut;

    return 0;
}
