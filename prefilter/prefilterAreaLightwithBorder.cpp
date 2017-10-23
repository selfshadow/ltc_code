#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
using namespace std;


#include <glm/glm.hpp>
using namespace glm;

#include "CImg.h"
using namespace cimg_library;

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

// filtered texture with borders

const float borderWidth = 0.125f;

// takes a position in [0,1]^2 and return the distance to the scaled texture in [borderWidth , 1 - borderWidth]
float distToScaledTexture(const vec2& position)
{
    vec2 position_clamped = clamp(position, vec2(borderWidth, borderWidth), vec2(1.0f - borderWidth, 1.0f - borderWidth));
    float dist = length(position - position_clamped);
    return dist;
}

// filtered textures

CImg<float> * imageInputPrefiltered;

// STD of gaussian filter applied on 2D image imageInputPrefiltered(:, :, level)
float level2gaussianFilterSTD(int level)
{
    float filterSTD = 0.5f * powf(1.3f, level);
    return filterSTD;
}

// inverse function
float gaussianFilterSTD2level(int filterSTD)
{
    float level = (logf(filterSTD) - logf(0.5f)) / logf(1.3f);
    level = std::max<float>(0.0f, std::min<float>(float((*imageInputPrefiltered).depth()) - 1.0f, level));
    return level;
}

// rescaled input (introduce black border)
// and prefiltered with Gaussian kernel
void prefilterInput(CImg<float>& imageInput)
{
    cout << "Prefiltering input..." << endl;

    // allocate
    imageInputPrefiltered = new CImg<float>(imageInput.width(), imageInput.height(), 32, 4);

    for (unsigned int level = 0; level < imageInputPrefiltered->depth(); ++level)
    {
        // map input image to [borderWidth, 1 - borderWidth]^2
        CImg<float> imageInput_(imageInput.width(), imageInput.height(), 1, 4);
        for (int j = 0; j < imageInput.height(); ++j)
        for (int i = 0; i < imageInput.width();  ++i)
        {
            // position in [0,1]^2
            vec2 position(i / float(imageInput.width() - 1.0f), j / float(imageInput.height() - 1.0f));

            // transform such that [borderWidth, 1 - borderWidth]^2 ---> [0, 1]^2
            position = (position - vec2(borderWidth)) / (1.0f - 2.0f*borderWidth);

            // are we inside the texture?
            bool inside = (position.x > 0 && position.x < 1 && position.y > 0 && position.y < 1);

            // fetch
            position.x *= imageInput.width()  - 1.0f;
            position.y *= imageInput.height() - 1.0f;
            imageInput_(i, j, 0, 0) = inside ? imageInput.linear_atXY(position.x, position.y, 0, 0, 0) : 0;
            imageInput_(i, j, 0, 1) = inside ? imageInput.linear_atXY(position.x, position.y, 0, 1, 0) : 0;
            imageInput_(i, j, 0, 2) = inside ? imageInput.linear_atXY(position.x, position.y, 0, 2, 0) : 0;
            imageInput_(i, j, 0, 3) = inside ? 1.0f : 0.0f;
        }

        // filter image
        float filterSTD = level2gaussianFilterSTD(level);
        imageInput_.blur(filterSTD, filterSTD, filterSTD, false);

        // copy image
        for (int j = 0; j < imageInput.height(); ++j)
        for (int i = 0; i < imageInput.width();  ++i)
        {
            (*imageInputPrefiltered)(i, j, level, 0) = imageInput_(i, j, 0, 0);
            (*imageInputPrefiltered)(i, j, level, 1) = imageInput_(i, j, 0, 1);
            (*imageInputPrefiltered)(i, j, level, 2) = imageInput_(i, j, 0, 2);
            (*imageInputPrefiltered)(i, j, level, 3) = imageInput_(i, j, 0, 3);
        }

        // debug
    #if 0
        cout << "temporary filtered file " << level << endl;
        cout << "filter std = " << filterSTD << endl;

        // save image
        stringstream filenameOutput(stringstream::in | stringstream::out);
        filenameOutput << "tmp_" << level << ".bmp";
        imageInput_.save(filenameOutput.str().c_str());

        cout << endl;
    #endif
    }
}

// compute the normalized filter
// with larger filter in the border
void filterWithBorder(CImg<float>& imageInput, CImg<float>& imageOutput, const int level, const int Nlevels)
{
    // distance to texture plane
    // in shader: LOD = log(powf(2.0f, Nlevels-1.0f) * dist) / log(3)
    const float dist = powf(3.0f, level) / powf(2.0f, Nlevels-1.0f);

    // filter size
    // at distance 1 ~= gaussian of std 0.75
    // account for the reduced size of the texture (1.0f - 2.0f*borderWidth)
    const float filterStd = 0.75f * dist * imageInput.width() * (1.0f - 2.0f*borderWidth);

    cout << "level "<< level << endl;
    cout << "distance to texture plane = " << dist << endl;
    cout << "filterStd = " << filterStd << endl << endl;

    CImg<float> tmp(imageInput.width(), imageInput.height(), 1, 3);

    for (int j = 0; j < imageInput.height(); ++j)
    for (int i = 0; i < imageInput.width();  ++i)
    {
        // position in [0, 1]^2
        vec2 position(i / float(imageInput.width() - 1), j / float(imageInput.height() - 1));

        // compute distance to texture in [borderWidth, 1-borderWidth]^2
        float distance = distToScaledTexture(position) * imageInput.width();

        // increase filter size with distance
        float filterStd_ = std::max<float>(1.0f + distance, filterStd);

        // compute level in precomputed data associated with filterStd_
        float l = gaussianFilterSTD2level(filterStd_);

        // fetch data
        tmp(i, j, 0, 0) = (*imageInputPrefiltered).linear_atXYZ(i, j, float(l), 0, 0);
        tmp(i, j, 0, 1) = (*imageInputPrefiltered).linear_atXYZ(i, j, float(l), 1, 0);
        tmp(i, j, 0, 2) = (*imageInputPrefiltered).linear_atXYZ(i, j, float(l), 2, 0);
        float alpha  = (*imageInputPrefiltered).linear_atXYZ(i, j, float(l), 3, 0);

        // renormalize
        if (alpha > 0.0f)
        {
            tmp(i, j, 0, 0) /= alpha;
            tmp(i, j, 0, 1) /= alpha;
            tmp(i, j, 0, 2) /= alpha;
        }
        else
        {
            tmp(i, j, 0, 0) = 0;
            tmp(i, j, 0, 1) = 0;
            tmp(i, j, 0, 2) = 0;
        }

        // clamp to 255 (numerical precision)
        tmp(i, j, 0, 0) = std::min<float>(tmp(i, j, 0, 0), 255.0f);
        tmp(i, j, 0, 1) = std::min<float>(tmp(i, j, 0, 1), 255.0f);
        tmp(i, j, 0, 2) = std::min<float>(tmp(i, j, 0, 2), 255.0f);
    }

    // rescale image
    imageOutput = tmp.resize(imageOutput, 5); // 5 = cubic interpolation
    return;
}



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

    // prefilter input with Gaussian kernel
    prefilterInput(imageInput);

    // filtered levels
    unsigned int Nlevels;
    for (Nlevels = 1; (imageInput.width() >> Nlevels) > 0; ++Nlevels);

    // borders
    for (unsigned int level = 0; level < Nlevels; ++level)
    {
        stringstream filenameOutput (stringstream::in | stringstream::out);
        filenameOutput << filename << "_filtered_" << level << ".hdr";

        cout << "processing file " << filenameOutput.str() << endl;
        unsigned int width = imageInput.width() >> level;
        unsigned int height = imageInput.height() >> level;

        if (width <= 0 || height <= 0)
            break;

        CImg<float> imageOutput(width, height, 1, 3);

        filterWithBorder(imageInput, imageOutput, level, Nlevels);

        offset = 0;
        for (int j = 0; j < imageOutput.height(); ++j)
        for (int i = 0; i < imageOutput.width();  ++i)
        {
            data[offset++] = imageOutput(i, j, 0, 0);
            data[offset++] = imageOutput(i, j, 0, 1);
            data[offset++] = imageOutput(i, j, 0, 2);
        }

        stbi_write_hdr(filenameOutput.str().c_str(), width, height, 3, data);
    }

    return 0;
}
