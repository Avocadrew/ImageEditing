///////////////////////////////////////////////////////////////////////////////
//
//      TargaImage.cpp                          Author:         Stephen Chenney
//                                              Modified:       Eric McDaniel
//                                              Date:           Fall 2004
//                                              Taskworker:     Jun-Yu Chen
//                                              FInish Date:    June 6, 2020  
//
//      Implementation of TargaImage methods.  You must implement the image
//  modification functions.
//
///////////////////////////////////////////////////////////////////////////////

#define _USE_MATH_DEFINES // for C++
#include <cmath>
#include "Globals.h"
#include "TargaImage.h"
#include "libtarga.h"
#include <stdlib.h>
#include <time.h>
#include <assert.h>
#include <memory.h>
#include <math.h>
#include <iostream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <array>
#include <iomanip> 
#include <iostream> 

using namespace std;

// constants
const int           RED = 0;                // red channel
const int           GREEN = 1;                // green channel
const int           BLUE = 2;                // blue channel
const unsigned char BACKGROUND[3] = { 0, 0, 0 };      // background color


// Computes n choose s, efficiently
double Binomial(int n, int s)
{
    double        res;

    res = 1;
    for (int i = 1; i <= s; i++)
        res = (n - i + 1) * res / i;

    return res;
}// Binomial


///////////////////////////////////////////////////////////////////////////////
//
//      Constructor.  Initialize member variables.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage::TargaImage() : width(0), height(0), data(NULL)
{}// TargaImage

///////////////////////////////////////////////////////////////////////////////
//
//      Constructor.  Initialize member variables.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage::TargaImage(int w, int h) : width(w), height(h)
{
    data = new unsigned char[width * height * 4];
    ClearToBlack();
}// TargaImage



///////////////////////////////////////////////////////////////////////////////
//
//      Constructor.  Initialize member variables to values given.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage::TargaImage(int w, int h, unsigned char* d)
{
    int i;

    width = w;
    height = h;
    data = new unsigned char[width * height * 4];

    for (i = 0; i < width * height * 4; i++)
        data[i] = d[i];
}// TargaImage

///////////////////////////////////////////////////////////////////////////////
//
//      Copy Constructor.  Initialize member to that of input
//
///////////////////////////////////////////////////////////////////////////////
TargaImage::TargaImage(const TargaImage& image)
{
    width = image.width;
    height = image.height;
    data = NULL;
    if (image.data != NULL) {
        data = new unsigned char[width * height * 4];
        memcpy(data, image.data, sizeof(unsigned char) * width * height * 4);
    }
}


///////////////////////////////////////////////////////////////////////////////
//
//      Destructor.  Free image memory.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage::~TargaImage()
{
    if (data)
        delete[] data;
}// ~TargaImage


///////////////////////////////////////////////////////////////////////////////
//
//      Converts an image to RGB form, and returns the rgb pixel data - 24 
//  bits per pixel. The returned space should be deleted when no longer 
//  required.
//
///////////////////////////////////////////////////////////////////////////////
unsigned char* TargaImage::To_RGB(void)
{
    unsigned char* rgb = new unsigned char[width * height * 3];
    int		    i, j;

    if (!data)
        return NULL;

    // Divide out the alpha
    for (i = 0; i < height; i++)
    {
        int in_offset = i * width * 4;
        int out_offset = i * width * 3;

        for (j = 0; j < width; j++)
        {
            RGBA_To_RGB(data + (in_offset + j * 4), rgb + (out_offset + j * 3));
        }
    }

    return rgb;
}// TargaImage


///////////////////////////////////////////////////////////////////////////////
//
//      Save the image to a targa file. Returns 1 on success, 0 on failure.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Save_Image(const char* filename)
{
    TargaImage* out_image = Reverse_Rows();

    if (!out_image)
        return false;

    if (!tga_write_raw(filename, width, height, out_image->data, TGA_TRUECOLOR_32))
    {
        cout << "TGA Save Error: %s\n", tga_error_string(tga_get_last_error());
        return false;
    }

    delete out_image;

    return true;
}// Save_Image


///////////////////////////////////////////////////////////////////////////////
//
//      Load a targa image from a file.  Return a new TargaImage object which 
//  must be deleted by caller.  Return NULL on failure.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage* TargaImage::Load_Image(char* filename)
{
    unsigned char* temp_data;
    TargaImage* temp_image;
    TargaImage* result;
    int		        width, height;
    srand(time(NULL));
    if (!filename)
    {
        cout << "No filename given." << endl;
        return NULL;
    }// if

    temp_data = (unsigned char*)tga_load(filename, &width, &height, TGA_TRUECOLOR_32);
    if (!temp_data)
    {
        cout << "TGA Error: %s\n", tga_error_string(tga_get_last_error());
        width = height = 0;
        return NULL;
    }
    temp_image = new TargaImage(width, height, temp_data);
    free(temp_data);

    result = temp_image->Reverse_Rows();

    delete temp_image;

    return result;
}// Load_Image


///////////////////////////////////////////////////////////////////////////////
//
//      Convert image to grayscale.  Red, green, and blue channels should all 
//  contain grayscale value.  Alpha channel shoould be left unchanged.  Return
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::To_Grayscale()
{
    if (!data)
        return NULL;

    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            int curPos = ((i * width) + j) * 4;
            //Acquire RGB
            unsigned char rgb[3];
            unsigned char l;
            RGBA_To_RGB(data + curPos, rgb);
            //calculate lumination
            l = 0.299 * rgb[0] + 0.587 * rgb[1] + 0.114 * rgb[2];
            //putback lumination factor
            for (int j = curPos; j < curPos + 3; j++)
            {
                data[j] = l;
            }
        }
    }
    return true;
}// To_Grayscale


///////////////////////////////////////////////////////////////////////////////
//
//  Convert the image to an 8 bit image using uniform quantization.  Return 
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Quant_Uniform()
{
    if (!data)
        return NULL;

    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            //acquire RGB data
            int curPos = ((i * width) + j) * 4;
            unsigned char rgb[3];
            RGBA_To_RGB(data + curPos, rgb);
            //create shading
            unsigned char eightShades = 32 + 64 + 128;
            unsigned char fourShades = 64 + 128;
            //putback data after shade
            data[curPos] = rgb[0] & eightShades;
            data[curPos + 1] = rgb[1] & eightShades;
            data[curPos + 2] = rgb[2] & fourShades;
        }
    }
    return true;
}// Quant_Uniform


///////////////////////////////////////////////////////////////////////////////
//
//      Convert the image to an 8 bit image using populosity quantization.  
//  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////

//struct for characteristic saving
typedef struct colorCharacteristics {
    long long int count = 0;
    unsigned char r = 0, g = 0, b = 0;
}colorCharacteristics;

//compare bool
static bool QPcmp(colorCharacteristics lhs, colorCharacteristics rhs)
{
    return lhs.count > rhs.count;
}
bool TargaImage::Quant_Populosity()
{
    if (!data)
        return NULL;
    //QuantUniform with 32 shades
    array <colorCharacteristics, 35937> choiceOfColor;

    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            int curPos = ((i * width) + j) * 4;
            unsigned char rgb[3];
            RGBA_To_RGB(data + curPos, rgb);
            unsigned char thirtytwoShades = 8 + 16 + 32 + 64 + 128;
            rgb[0] = rgb[0] & thirtytwoShades;
            rgb[1] = rgb[1] & thirtytwoShades;
            rgb[2] = rgb[2] & thirtytwoShades;
            //create colorcode for each entry
            long int colorCode = (rgb[0] / 4 / 2 / 1) * 33 * 33 + (rgb[1] / 4 / 2 / 1) * 33 + (rgb[2] / 4 / 2 / 1);
            //calculate times of appearance
            choiceOfColor[colorCode].count++;
            //save into struct
            choiceOfColor[colorCode].r = rgb[0];
            choiceOfColor[colorCode].g = rgb[1];
            choiceOfColor[colorCode].b = rgb[2];
        }
    }
    //sort
    sort(choiceOfColor.begin(), choiceOfColor.end(), QPcmp);

    /*cout << "First Stage" << endl;
    for (int i = 0; i < 256; i++)
    {
        cout << (int)choiceOfColor.at(i).r <<' '<< (int)choiceOfColor.at(i).g <<' ' << (int)choiceOfColor.at(i).b << ' '<<choiceOfColor.at(i).count << endl;
    }*/

    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            int curPos = ((i * width) + j) * 4;
            unsigned long long int shortestDist = ULLONG_MAX;
            unsigned long long int dist = 0;
            unsigned char curRGB[3];
            RGBA_To_RGB(data + curPos, curRGB);
            long long int curR = curRGB[0];
            long long int curG = curRGB[1];
            long long int curB = curRGB[2];
            //find shortest dist and save
            for (int k = 0; k < 256; k++)
            {
                dist = 0;
                long int testR = choiceOfColor[k].r;
                dist += pow(abs(testR - curR), 2);
                long int testG = choiceOfColor[k].g;
                dist += pow(abs(testG - curG), 2);
                long int testB = choiceOfColor[k].b;
                dist += pow(abs(testB - curB), 2);
                if (dist < shortestDist)
                {
                    shortestDist = dist;
                    data[curPos] = (unsigned char)testR;
                    data[curPos + 1] = (unsigned char)testG;
                    data[curPos + 2] = (unsigned char)testB;

                }
            }
        }
    }
    return true;
}// Quant_Populosity


///////////////////////////////////////////////////////////////////////////////
//
//      Dither the image using a threshold of 1/2.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Threshold()
{
    if (!data)
        return NULL;

    To_Grayscale();

    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            int curPos = ((i * width) + j) * 4;
            //acquire RGB data
            unsigned char RGB[3];
            unsigned char outputColor = 0;
            RGBA_To_RGB(data + curPos, RGB);
            //threshold lumination at 127, thus:
            if (RGB[0] > 127)
            {
                outputColor = 255;
            }
            //putback data
            for (int j = curPos; j < curPos + 3; j++)
            {
                data[j] = outputColor;
            }
        }
    }
    return true;
}// Dither_Threshold


///////////////////////////////////////////////////////////////////////////////
//
//      Dither image using random dithering.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Random()
{
    if (!data)
        return NULL;

    To_Grayscale();

    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            //acquire RGB Data
            int curPos = ((i * width) + j) * 4;
            unsigned char RGB[3];
            unsigned char outputColor = 0;
            RGBA_To_RGB(data + curPos, RGB);
            //create threshold(+-0.2), srand in load function
            int randomThreshold = rand() % 103 + 77;
            if (RGB[0] > randomThreshold)
            {
                outputColor = 255;
            }
            //putback data
            for (int j = curPos; j < curPos + 3; j++)
            {
                data[j] = outputColor;
            }
        }
    }
    return true;
}// Dither_Random


///////////////////////////////////////////////////////////////////////////////
//
//      Perform Floyd-Steinberg dithering on the image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_FS()
{
    if (!data)
        return NULL;

    //make original into RGB
    unsigned char* rgb = this->To_RGB();

    //create error array
    float** errors;
    errors = new float* [height];
    for (int i = 0; i < height; i++)
    {
        errors[i] = new float[width];
    }
    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            errors[i][j] = 0;
        }
    }

    //calculation
    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            int zigzagJ, zigzagDir;
            //determine zigzag direction and initial J pos
            if (i % 2 == 0)
            {
                zigzagJ = j;
                zigzagDir = 1;
            }
            else
            {
                zigzagJ = width - j - 1;
                zigzagDir = -1;
            }

            //parallel compute position
            int rgbPos = (i * width + zigzagJ) * 3;
            int curPos = (i * width + zigzagJ) * 4;
            float gray = (float)(0.3 * rgb[rgbPos] + 0.59 * rgb[rgbPos + 1] + 0.11 * rgb[rgbPos + 2]);
            gray += errors[i][zigzagJ]; //add error of position to gray

            //clipping
            if (gray > 255)gray = 255; if (gray < 0)gray = 0;

            //determine putback color
            unsigned char inputColor = 0;
            if (gray > 127)
            {
                inputColor = 255;
            }
            for (int k = curPos; k < curPos + 4; k++)
            {
                data[k] = inputColor;
            }

            //create error
            int error = gray - inputColor;

            //add weights of error as stated in Powerpoint to error filter
            if (zigzagJ + zigzagDir < width && zigzagJ + zigzagDir >= 0)
            {

                errors[i][zigzagJ + zigzagDir] += (float)error * 7.0 / 16.0;

                if (i + 1 < height)
                {
                    errors[i + 1][zigzagJ + zigzagDir] += (float)error * 1.0 / 16.0;
                }
            }

            if (i + 1 < height)
            {
                errors[i + 1][zigzagJ] += (float)error * 5.0 / 16.0;
                if (zigzagJ - zigzagDir >= 0 && zigzagJ - zigzagDir < width)
                {
                    errors[i + 1][zigzagJ - zigzagDir] += (float)error * 3.0 / 16.0;
                }
            }
        }
    }
    return true;
}// Dither_FS


///////////////////////////////////////////////////////////////////////////////
//
//      Dither the image while conserving the average brightness.  Return 
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Bright()
{
    vector<double> grayArray;
    double totGray = 0;

    if (!data)
        return NULL;

    To_Grayscale();

    //calculate total grayness & record gray values
    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            int curPos = ((i * width) + j) * 4;
            unsigned char RGB[3];
            RGBA_To_RGB(data + curPos, RGB);
            double gray = 0.3 * RGB[0] + 0.59 * RGB[1] + 0.11 * RGB[2];
            totGray += gray;
            grayArray.push_back(gray);
        }
    }
    //create average gray value
    double avgGray = totGray / height / width;
    sort(grayArray.begin(), grayArray.end());

    //create thershold for average lumination
    int threshold = grayArray.at((float)(255 - (int)avgGray) / 256.0 * grayArray.size());

    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            int curPos = ((i * width) + j) * 4;
            unsigned char RGB[3];
            unsigned char outputColor = 0;
            RGBA_To_RGB(data + curPos, RGB);
            //determine average color
            if (0.3 * RGB[0] + 0.59 * RGB[1] + 0.11 * RGB[2] > threshold)
            {
                outputColor = 255;
            }
            //putback data
            for (int k = curPos; k < curPos + 3; k++)
            {
                data[k] = outputColor;
            }
        }
    }

    return true;
}// Dither_Bright


///////////////////////////////////////////////////////////////////////////////
//
//      Perform clustered differing of the image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Cluster()
{
    //cluster matrix
    float DCMatrix[4][4] = { 0.7059, 0.3529, 0.5882, 0.2353,
                                0.0588, 0.9412, 0.8235, 0.4118,
                                0.4706, 0.7647, 0.8824, 0.1176,
                                0.1765, 0.5294, 0.2941, 0.6471 };
    if (!data)
        return NULL;
    To_Grayscale();

    //calculate weights with matrix given (create threshold for B/W)
    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            int curPos = ((i * width) + j) * 4;
            unsigned char RGB[3];
            RGBA_To_RGB(data + curPos, RGB);
            unsigned char outputColor = 0;
            if ((0.299 * RGB[0] + 0.587 * RGB[1] + 0.114 * RGB[2]) / 256 >= DCMatrix[i % 4][j % 4])
            {
                outputColor = 255;
            }
            for (int k = curPos; k < curPos + 3; k++)
            {
                data[k] = outputColor;
            }
        }
    }
    return true;
}// Dither_Cluster


///////////////////////////////////////////////////////////////////////////////
//
//  Convert the image to an 8 bit image using Floyd-Steinberg dithering over
//  a uniform quantization - the same quantization as in Quant_Uniform.
//  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Color()
{
    int x, y;
    int i, j;
    float error;
    int p1, p2;
    int color, newX, count;
    unsigned char newColor;
    unsigned char map1[] = { 0, 36, 73, 109, 146, 182, 219 , 255 };
    unsigned char map2[] = { 0, 85, 170, 255 };
    float* tmp = new float[3 * width * height];
    for (y = 0; y < height; y++)
    {
        for (x = 0; x < width; x++)
        {
            if (y % 2 == 0)
                newX = x;
            else
                newX = width - 1 - x;

            p1 = 4 * (y * width + newX);
            p2 = 3 * (y * width + newX);
            for (color = 0; color < 3; color++)
            {
                if (color != 2)
                {
                    newColor = data[p1 + color] + tmp[p2 + color] > 255 ? 255 : data[p1 + color] + tmp[p2 + color];
                    for (i = 0; i < 7; i++)
                    {
                        if (newColor >= map1[i] && newColor < map1[i + 1])
                        {
                            error = newColor - map1[i];
                            break;
                        }
                    }
                    data[p1 + color] = map1[i];
                }
                else
                {
                    newColor = data[p1 + color] + tmp[p2 + color] > 255 ? 255 : data[p1 + color] + tmp[p2 + color];
                    for (i = 0; i < 3; i++)
                    {
                        if (newColor >= map2[i] && newColor < map2[i + 1])
                        {
                            error = newColor - map2[i];
                            break;
                        }
                    }
                    data[p1 + color] = map2[i];
                }
                if (y % 2 == 0)
                {
                    if (newX == 0)
                    {
                        tmp[3 * (y * width + newX + 1) + color] += (error * 0.4375);
                        if (y != height - 1)
                        {
                            tmp[3 * ((y + 1) * width + newX) + color] += (error * 0.3125);
                            tmp[3 * ((y + 1) * width + newX + 1) + color] += (error * 0.0625);
                        }
                    }
                    else if (newX == width - 1)
                    {
                        if (y != height - 1)
                        {
                            tmp[3 * ((y + 1) * width + newX) + color] += (error * 0.3125);
                            tmp[3 * ((y + 1) * width + newX - 1) + color] += (error * 0.1875);
                        }
                    }
                    else
                    {
                        tmp[3 * (y * width + newX + 1) + color] += (error * 0.4375);
                        if (y != height - 1)
                        {
                            tmp[3 * ((y + 1) * width + newX - 1) + color] += (error * 0.1875);
                            tmp[3 * ((y + 1) * width + newX) + color] += (error * 0.3125);
                            tmp[3 * ((y + 1) * width + newX + 1) + color] += (error * 0.0625);
                        }
                    }
                }
                else
                {
                    if (newX == 0)
                    {
                        if (y != height - 1)
                        {
                            tmp[3 * ((y + 1) * width + newX) + color] += (error * 0.3125);
                            tmp[3 * ((y + 1) * width + newX + 1) + color] += (error * 0.1875);
                        }
                    }
                    else if (newX == width - 1)
                    {
                        tmp[3 * (y * width + newX - 1) + color] += (error * 0.4375);
                        if (y != height - 1)
                        {
                            tmp[3 * ((y + 1) * width + newX) + color] += (error * 0.3125);
                            tmp[3 * ((y + 1) * width + newX - 1) + color] += (error * 0.0625);
                        }
                    }
                    else
                    {
                        tmp[3 * (y * width + newX - 1) + color] += (error * 0.4375);
                        if (y != height - 1)
                        {
                            tmp[3 * ((y + 1) * width + newX + 1) + color] += (error * 0.1875);
                            tmp[3 * ((y + 1) * width + newX) + color] += (error * 0.3125);
                            tmp[3 * ((y + 1) * width + newX - 1) + color] += (error * 0.0625);
                        }
                    }
                }
            }
        }
    }
    delete[]tmp;
    return true;
}// Dither_Color


///////////////////////////////////////////////////////////////////////////////
//
//      Composite the current image over the given image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_Over(TargaImage* pImage)
{
    if (width != pImage->width || height != pImage->height)
    {
        cout << "Comp_Over: Images not the same size\n";
        return false;
    }

    for (int i = 0; i < width * height * 4; i += 4) {
        double alpha = ((double)data[i + 3]) / 255.0;

        for (int j = 0; j < 4; j++)
            data[i + j] += (pImage->data[i + j] * (1.0 - alpha));;
    }
    return true;
}// Comp_Over


///////////////////////////////////////////////////////////////////////////////
//
//      Composite this image "in" the given image.  See lecture notes for 
//  details.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_In(TargaImage* pImage)
{
    if (width != pImage->width || height != pImage->height)
    {
        cout << "Comp_In: Images not the same size\n";
        return false;
    }

    for (int i = 0; i < width * height * 4; i += 4) {
        double alpha = ((double)pImage->data[i + 3]) / 255.0;

        for (int j = 0; j < 4; j++)
            data[i + j] *= alpha;
    }
    return true;
}// Comp_In


///////////////////////////////////////////////////////////////////////////////
//
//      Composite this image "out" the given image.  See lecture notes for 
//  details.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_Out(TargaImage* pImage)
{
    if (width != pImage->width || height != pImage->height) {
        cout << "Comp_Out: Images not the same size\n";
        return false;
    }

    for (int i = 0; i < width * height * 4; i += 4) {
        double alpha = ((double)pImage->data[i + 3]) / 255.0;

        for (int j = 0; j < 4; j++)
            data[i + j] *= (1.0 - alpha);
    }
    return true;
}// Comp_Out


///////////////////////////////////////////////////////////////////////////////
//
//      Composite current image "atop" given image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_Atop(TargaImage* pImage)
{
    if (width != pImage->width || height != pImage->height)
    {
        cout << "Comp_Atop: Images not the same size\n";
        return false;
    }

    for (int i = 0; i < width * height * 4; i += 4) {
        double alphaf = ((double)data[i + 3]) / 255.0;
        double alphag = ((double)pImage->data[i + 3]) / 255.0;

        for (int j = 0; j < 4; j++)
            data[i + j] = (data[i + j] * alphag) + (pImage->data[i + j] * (1.0 - alphaf));
    }
    return true;
}// Comp_Atop


///////////////////////////////////////////////////////////////////////////////
//
//      Composite this image with given image using exclusive or (XOR).  Return
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_Xor(TargaImage* pImage)
{
    if (width != pImage->width || height != pImage->height)
    {
        cout << "Comp_Xor: Images not the same size\n";
        return false;
    }

    for (int i = 0; i < width * height * 4; i += 4) {
        double alphaf = ((double)data[i + 3]) / 255.0;
        double alphag = ((double)pImage->data[i + 3]) / 255.0;

        for (int j = 0; j < 4; j++)
            data[i + j] = (data[i + j] * (1.0 - alphag)) + (pImage->data[i + j] * (1.0 - alphaf));
    }
    return true;
}// Comp_Xor


///////////////////////////////////////////////////////////////////////////////
//
//      Calculate the difference bewteen this imag and the given one.  Image 
//  dimensions must be equal.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Difference(TargaImage* pImage)
{
    if (!pImage)
        return false;

    if (width != pImage->width || height != pImage->height)
    {
        cout << "Difference: Images not the same size\n";
        return false;
    }// if

    for (int i = 0; i < width * height * 4; i += 4)
    {
        unsigned char        rgb1[3];
        unsigned char        rgb2[3];

        RGBA_To_RGB(data + i, rgb1);
        RGBA_To_RGB(pImage->data + i, rgb2);

        data[i] = abs(rgb1[0] - rgb2[0]);
        data[i + 1] = abs(rgb1[1] - rgb2[1]);
        data[i + 2] = abs(rgb1[2] - rgb2[2]);
        data[i + 3] = 255;
    }

    return true;
}// Difference



///////////////////////////////////////////////////////////////////////////////
//
//      Perform 5x5 box filter on this image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////


bool TargaImage::Filter_Box()
{
    //Box Filter
    double boxFilter[5][5] =
    {
        1.0, 1.0, 1.0, 1.0, 1.0,
        1.0, 1.0, 1.0, 1.0, 1.0,
        1.0, 1.0, 1.0, 1.0, 1.0,
        1.0, 1.0, 1.0, 1.0, 1.0,
        1.0, 1.0, 1.0, 1.0, 1.0
    };
    for (int i = 0; i < 5; i++)
    {
        for (int j = 0; j < 5; j++)
        {
            boxFilter[i][j] = boxFilter[i][j] / 25;
        }
    }

    //mask filter
    unsigned char* RGB = this->To_RGB();

    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            int curPos = (i * width + j) * 4;
            //sum of colors
            double sumR = 0;
            double sumG = 0;
            double sumB = 0;
            for (int k = -2; k <= 2; k++)
            {
                for (int l = -2; l <= 2; l++)
                {
                    //mask 3*3 region
                    if (i + k >= 0 && i + k < height && j + l >= 0 && j + l < width)
                    {
                        //add sum
                        int rgbPos = (i + k) * width * 3 + (j + l) * 3;
                        sumR += (double)RGB[rgbPos] * boxFilter[k + 2][l + 2];
                        sumG += (double)RGB[rgbPos + 1] * boxFilter[k + 2][l + 2];
                        sumB += (double)RGB[rgbPos + 2] * boxFilter[k + 2][l + 2];
                    }
                    else if ((i + k < 0 && j + l < 0) || (i + k < 0 && j + l >= width) || (i + k >= height && j + l < 0) || (i + k >= height && j + l >= width))
                    {
                        int rgbPos = (i - k) * width * 3 + (j - l) * 3;
                        sumR += (double)RGB[rgbPos] * boxFilter[k + 2][l + 2];
                        sumG += (double)RGB[rgbPos + 1] * boxFilter[k + 2][l + 2];
                        sumB += (double)RGB[rgbPos + 2] * boxFilter[k + 2][l + 2];
                    }
                    else if (i + k < 0 || i + k >= height)
                    {
                        int rgbPos = (i - k) * width * 3 + (j + l) * 3;
                        sumR += (double)RGB[rgbPos] * boxFilter[k + 2][l + 2];
                        sumG += (double)RGB[rgbPos + 1] * boxFilter[k + 2][l + 2];
                        sumB += (double)RGB[rgbPos + 2] * boxFilter[k + 2][l + 2];
                    }
                    else if (j + l < 0 || j + l >= width)
                    {
                        int rgbPos = (i + k) * width * 3 + (j - l) * 3;
                        sumR += (double)RGB[rgbPos] * boxFilter[k + 2][l + 2];
                        sumG += (double)RGB[rgbPos + 1] * boxFilter[k + 2][l + 2];
                        sumB += (double)RGB[rgbPos + 2] * boxFilter[k + 2][l + 2];
                    }
                }
            }
            //clip sum and apply
            if (sumR < 0)sumR = 0; if (sumR > 255)sumR = 255; data[curPos] = (unsigned char)sumR;
            if (sumG < 0)sumG = 0; if (sumG > 255)sumG = 255; data[curPos + 1] = (unsigned char)sumG;
            if (sumB < 0)sumB = 0; if (sumB > 255)sumB = 255; data[curPos + 2] = (unsigned char)sumB;
        }
    }

    return true;
}// Filter_Box


///////////////////////////////////////////////////////////////////////////////
//
//      Perform 5x5 Bartlett filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Bartlett()
{
    //bartlett filter
    double bartlettFilter[5][5] =
    {
        1.0, 2.0, 3.0, 2.0, 1.0,
        2.0, 4.0, 6.0, 4.0, 2.0,
        3.0, 6.0, 9.0, 6.0, 3.0,
        2.0, 4.0, 6.0, 4.0, 2.0,
        1.0, 2.0, 3.0, 2.0, 1.0
    };
    for (int i = 0; i < 5; i++)
    {
        for (int j = 0; j < 5; j++)
        {
            bartlettFilter[i][j] = bartlettFilter[i][j] / 81.0;
        }
    }

    //mask filter
    unsigned char* RGB = this->To_RGB();

    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            int curPos = (i * width + j) * 4;
            double sumR = 0;
            double sumG = 0;
            double sumB = 0;
            for (int k = -2; k <= 2; k++)
            {
                for (int l = -2; l <= 2; l++)
                {
                    //mask 3*3 region
                    if (i + k >= 0 && i + k < height && j + l >= 0 && j + l < width)
                    {
                        //add sum
                        int rgbPos = (i + k) * width * 3 + (j + l) * 3;
                        sumR += (double)RGB[rgbPos] * bartlettFilter[k + 2][l + 2];
                        sumG += (double)RGB[rgbPos + 1] * bartlettFilter[k + 2][l + 2];
                        sumB += (double)RGB[rgbPos + 2] * bartlettFilter[k + 2][l + 2];
                    }
                    else if ((i + k < 0 && j + l < 0) || (i + k < 0 && j + l >= width) || (i + k >= height && j + l < 0) || (i + k >= height && j + l >= width))
                    {
                        int rgbPos = (i - k) * width * 3 + (j - l) * 3;
                        sumR += (double)RGB[rgbPos] * bartlettFilter[k + 2][l + 2];
                        sumG += (double)RGB[rgbPos + 1] * bartlettFilter[k + 2][l + 2];
                        sumB += (double)RGB[rgbPos + 2] * bartlettFilter[k + 2][l + 2];
                    }
                    else if (i + k < 0 || i + k >= height)
                    {
                        int rgbPos = (i - k) * width * 3 + (j + l) * 3;
                        sumR += (double)RGB[rgbPos] * bartlettFilter[k + 2][l + 2];
                        sumG += (double)RGB[rgbPos + 1] * bartlettFilter[k + 2][l + 2];
                        sumB += (double)RGB[rgbPos + 2] * bartlettFilter[k + 2][l + 2];
                    }
                    else if (j + l < 0 || j + l >= width)
                    {
                        int rgbPos = (i + k) * width * 3 + (j - l) * 3;
                        sumR += (double)RGB[rgbPos] * bartlettFilter[k + 2][l + 2];
                        sumG += (double)RGB[rgbPos + 1] * bartlettFilter[k + 2][l + 2];
                        sumB += (double)RGB[rgbPos + 2] * bartlettFilter[k + 2][l + 2];
                    }
                }
            }
            if (sumR < 0)sumR = 0; if (sumR > 255)sumR = 255; data[curPos] = (unsigned char)sumR;
            if (sumG < 0)sumG = 0; if (sumG > 255)sumG = 255; data[curPos + 1] = (unsigned char)sumG;
            if (sumB < 0)sumB = 0; if (sumB > 255)sumB = 255; data[curPos + 2] = (unsigned char)sumB;
        }
    }

    return true;
}// Filter_Bartlett


///////////////////////////////////////////////////////////////////////////////
//
//      Perform 5x5 Gaussian filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Gaussian()
{
    //gaussian filter
    double gaussianFilter[5][5] =
    {
        1.0, 4.0, 6.0, 4.0, 1.0,
        4.0, 16.0, 24.0, 16.0, 4.0,
        6.0, 24.0, 36.0, 24.0, 6.0,
        4.0, 16.0, 24.0, 16.0, 4.0,
        1.0, 4.0, 6.0, 4.0, 1.0
    };
    for (int i = 0; i < 5; i++)
    {
        for (int j = 0; j < 5; j++)
        {
            gaussianFilter[i][j] = gaussianFilter[i][j] / 256.0;
        }
    }
    //mask filter
    unsigned char* RGB = this->To_RGB();

    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            int curPos = (i * width + j) * 4;
            double sumR = 0;
            double sumG = 0;
            double sumB = 0;
            for (int k = -2; k <= 2; k++)
            {
                for (int l = -2; l <= 2; l++)
                {
                    //mask 3*3 region
                    if (i + k >= 0 && i + k < height && j + l >= 0 && j + l < width)
                    {
                        //add sum
                        int rgbPos = (i + k) * width * 3 + (j + l) * 3;
                        sumR += (double)RGB[rgbPos] * gaussianFilter[k + 2][l + 2];
                        sumG += (double)RGB[rgbPos + 1] * gaussianFilter[k + 2][l + 2];
                        sumB += (double)RGB[rgbPos + 2] * gaussianFilter[k + 2][l + 2];
                    }
                    else if ((i + k < 0 && j + l < 0) || (i + k < 0 && j + l >= width) || (i + k >= height && j + l < 0) || (i + k >= height && j + l >= width))
                    {
                        int rgbPos = (i - k) * width * 3 + (j - l) * 3;
                        sumR += (double)RGB[rgbPos] * gaussianFilter[k + 2][l + 2];
                        sumG += (double)RGB[rgbPos + 1] * gaussianFilter[k + 2][l + 2];
                        sumB += (double)RGB[rgbPos + 2] * gaussianFilter[k + 2][l + 2];
                    }
                    else if (i + k < 0 || i + k >= height)
                    {
                        int rgbPos = (i - k) * width * 3 + (j + l) * 3;
                        sumR += (double)RGB[rgbPos] * gaussianFilter[k + 2][l + 2];
                        sumG += (double)RGB[rgbPos + 1] * gaussianFilter[k + 2][l + 2];
                        sumB += (double)RGB[rgbPos + 2] * gaussianFilter[k + 2][l + 2];
                    }
                    else if (j + l < 0 || j + l >= width)
                    {
                        int rgbPos = (i + k) * width * 3 + (j - l) * 3;
                        sumR += (double)RGB[rgbPos] * gaussianFilter[k + 2][l + 2];
                        sumG += (double)RGB[rgbPos + 1] * gaussianFilter[k + 2][l + 2];
                        sumB += (double)RGB[rgbPos + 2] * gaussianFilter[k + 2][l + 2];
                    }
                }
            }
            if (sumR < 0)sumR = 0; if (sumR > 255)sumR = 255; data[curPos] = (unsigned char)sumR;
            if (sumG < 0)sumG = 0; if (sumG > 255)sumG = 255; data[curPos + 1] = (unsigned char)sumG;
            if (sumB < 0)sumB = 0; if (sumB > 255)sumB = 255; data[curPos + 2] = (unsigned char)sumB;
        }
    }
    return true;
}// Filter_Gaussian

///////////////////////////////////////////////////////////////////////////////
//
//      Perform NxN Gaussian filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////

bool TargaImage::Filter_Gaussian_N(unsigned int N)
{
    double sum = 0.0;
    double sigma = 1.0, PI = 3.14159;
    double rl, s = 2.0 * sigma * sigma;
    int dl = (N - 1) / 2;
    double** gaussianFilterN = new double* [N];
    for (int i = 0; i < N; i++)
    {
        gaussianFilterN[i] = new double[N];
    }
    for (int x = -dl; x <= dl; x++)
    {
        for (int y = -dl; y <= dl; y++)
        {
            rl = x * x + y * y;
            gaussianFilterN[y + dl][x + dl] = (exp(-rl / s)) / (PI * s);
            sum += gaussianFilterN[y + dl][x + dl];
        }
    }
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            gaussianFilterN[i][j] /= sum;
        }
    }
    //mask filter
    unsigned char* RGB = this->To_RGB();

    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            int curPos = (i * width + j) * 4;
            double sumR = 0;
            double sumG = 0;
            double sumB = 0;
            for (int k = -(N - 1) / 2; k <= (N - 1) / 2; k++)
            {
                for (int l = -(N - 1) / 2; l <= (N - 1) / 2; l++)
                {
                    //mask 3*3 region
                    if (i + k >= 0 && i + k < height && j + l >= 0 && j + l < width)
                    {
                        //add sum
                        int rgbPos = (i + k) * width * 3 + (j + l) * 3;
                        sumR += (double)RGB[rgbPos] * gaussianFilterN[k + 2][l + 2];
                        sumG += (double)RGB[rgbPos + 1] * gaussianFilterN[k + 2][l + 2];
                        sumB += (double)RGB[rgbPos + 2] * gaussianFilterN[k + 2][l + 2];
                    }
                    else if ((i + k < 0 && j + l < 0) || (i + k < 0 && j + l >= width) || (i + k >= height && j + l < 0) || (i + k >= height && j + l >= width))
                    {
                        int rgbPos = (i - k) * width * 3 + (j - l) * 3;
                        sumR += (double)RGB[rgbPos] * gaussianFilterN[k + 2][l + 2];
                        sumG += (double)RGB[rgbPos + 1] * gaussianFilterN[k + 2][l + 2];
                        sumB += (double)RGB[rgbPos + 2] * gaussianFilterN[k + 2][l + 2];
                    }
                    else if (i + k < 0 || i + k >= height)
                    {
                        int rgbPos = (i - k) * width * 3 + (j + l) * 3;
                        sumR += (double)RGB[rgbPos] * gaussianFilterN[k + 2][l + 2];
                        sumG += (double)RGB[rgbPos + 1] * gaussianFilterN[k + 2][l + 2];
                        sumB += (double)RGB[rgbPos + 2] * gaussianFilterN[k + 2][l + 2];
                    }
                    else if (j + l < 0 || j + l >= width)
                    {
                        int rgbPos = (i + k) * width * 3 + (j - l) * 3;
                        sumR += (double)RGB[rgbPos] * gaussianFilterN[k + 2][l + 2];
                        sumG += (double)RGB[rgbPos + 1] * gaussianFilterN[k + 2][l + 2];
                        sumB += (double)RGB[rgbPos + 2] * gaussianFilterN[k + 2][l + 2];
                    }
                }
            }
            if (sumR < 0)sumR = 0; if (sumR > 255)sumR = 255; data[curPos] = (unsigned char)sumR;
            if (sumG < 0)sumG = 0; if (sumG > 255)sumG = 255; data[curPos + 1] = (unsigned char)sumG;
            if (sumB < 0)sumB = 0; if (sumB > 255)sumB = 255; data[curPos + 2] = (unsigned char)sumB;
        }
    }
    return true;
}// Filter_Gaussian_N


///////////////////////////////////////////////////////////////////////////////
//
//      Perform 5x5 edge detect (high pass) filter on this image.  Return 
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Edge()
{
    //edge filter
    double edgeFilter[5][5] =
    {
        -0.0123457, -0.0246914, -0.037037, -0.0246914, -0.0123457,
        -0.0246914, -0.0493827, -0.0740741, -0.0493827, -0.0246914,
        -0.037037, -0.0740741, 0.888889, -0.0740741, -0.037037,
        -0.0246914, -0.0493827, -0.0740741, -0.0493827, -0.0246914,
        -0.0123457, -0.0246914, -0.037037, -0.0246914, -0.0123457
    };
    for (int i = 0; i <= 5; i++)
    {
        for (int j = 0; j <= 5; j++)
        {
            edgeFilter[i][j] = edgeFilter[i][j] / 81.0;
        }
    }

    //mask filter
    unsigned char* RGB = this->To_RGB();

    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            int curPos = (i * width + j) * 4;
            double sumR = 0;
            double sumG = 0;
            double sumB = 0;
            for (int k = -2; k <= 2; k++)
            {
                for (int l = -2; l <= 2; l++)
                {
                    if (i + k >= 0 && i + k < height && j + l >= 0 && j + l < width)
                    {
                        //add sum
                        int rgbPos = (i + k) * width * 3 + (j + l) * 3;
                        sumR += (double)RGB[rgbPos] * edgeFilter[k + 2][l + 2];
                        sumG += (double)RGB[rgbPos + 1] * edgeFilter[k + 2][l + 2];
                        sumB += (double)RGB[rgbPos + 2] * edgeFilter[k + 2][l + 2];
                    }
                    else if ((i + k < 0 && j + l < 0) || (i + k < 0 && j + l >= width) || (i + k >= height && j + l < 0) || (i + k >= height && j + l >= width))
                    {
                        int rgbPos = (i - k) * width * 3 + (j - l) * 3;
                        sumR += (double)RGB[rgbPos] * edgeFilter[k + 2][l + 2];
                        sumG += (double)RGB[rgbPos + 1] * edgeFilter[k + 2][l + 2];
                        sumB += (double)RGB[rgbPos + 2] * edgeFilter[k + 2][l + 2];
                    }
                    else if (i + k < 0 || i + k >= height)
                    {
                        int rgbPos = (i - k) * width * 3 + (j + l) * 3;
                        sumR += (double)RGB[rgbPos] * edgeFilter[k + 2][l + 2];
                        sumG += (double)RGB[rgbPos + 1] * edgeFilter[k + 2][l + 2];
                        sumB += (double)RGB[rgbPos + 2] * edgeFilter[k + 2][l + 2];
                    }
                    else if (j + l < 0 || j + l >= width)
                    {
                        int rgbPos = (i + k) * width * 3 + (j - l) * 3;
                        sumR += (double)RGB[rgbPos] * edgeFilter[k + 2][l + 2];
                        sumG += (double)RGB[rgbPos + 1] * edgeFilter[k + 2][l + 2];
                        sumB += (double)RGB[rgbPos + 2] * edgeFilter[k + 2][l + 2];
                    }
                }
            }
            if (sumR < 0)sumR = 0; if (sumR > 255)sumR = 255; data[curPos] = (unsigned char)sumR;
            if (sumG < 0)sumG = 0; if (sumG > 255)sumG = 255; data[curPos + 1] = (unsigned char)sumG;
            if (sumB < 0)sumB = 0; if (sumB > 255)sumB = 255; data[curPos + 2] = (unsigned char)sumB;
        }
    }
    return true;
}// Filter_Edge


///////////////////////////////////////////////////////////////////////////////
//
//      Perform a 5x5 enhancement filter to this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Enhance()
{
    ClearToBlack();
    return false;
}// Filter_Enhance


///////////////////////////////////////////////////////////////////////////////
//
//      Run simplified version of Hertzmann's painterly image filter.
//      You probably will want to use the Draw_Stroke funciton and the
//      Stroke class to help.
// Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::NPR_Paint()
{
    ClearToBlack();
    return false;
}



///////////////////////////////////////////////////////////////////////////////
//
//      Halve the dimensions of this image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Half_Size()
{
    int n = 3;
    int d = (n - 1) / 2;
    double halfFilter[3][3] = 
    { 
        0.0625, 0.125, 0.0625,
        0.125 , 0.25 , 0.125,
        0.0625, 0.125, 0.0625 
    };
    int originalX, originalY;
    int newHeight = (height / 2);
    int newWidth = (width / 2);
    unsigned char* tmp;
    tmp = new unsigned char[4 * newHeight * newWidth];
    for (int y = 0; y < newHeight; y++)
    {
        for (int x = 0; x < newWidth; x++)
        {
            originalX = x * 2 + 1;
            originalY = y * 2 + 1;
            if (originalX >= width - d || originalY >= height - d)
                continue;
            double colorSumR = 0;
            double colorSumG = 0;
            double colorSumB = 0;
            for (int i = -d; i <= d; i++)
            {
                for (int j = -d; j <= d; j++)
                {
                    colorSumR = colorSumR + data[4 * ((originalY + i) * width + originalX + j)] * halfFilter[i + d][j + d];
                    colorSumG = colorSumG + data[4 * ((originalY + i) * width + originalX + j) + 1] * halfFilter[i + d][j + d];
                    colorSumB = colorSumB + data[4 * ((originalY + i) * width + originalX + j) + 2] * halfFilter[i + d][j + d];
                }
            }
            if (colorSumR < 0)colorSumR = 0; if (colorSumR > 255)colorSumR = 255; tmp[4 * (y * newWidth + x)] = (unsigned char)colorSumR;
            if (colorSumG < 0)colorSumG = 0; if (colorSumG > 255)colorSumG = 255; tmp[4 * (y * newWidth + x) + 1] = (unsigned char)colorSumG;
            if (colorSumB < 0)colorSumB = 0; if (colorSumB > 255)colorSumB = 255; tmp[4 * (y * newWidth + x) + 2] = (unsigned char)colorSumB;
        }
    }
    width = newWidth;
    height = newHeight;
    delete[]data;
    data = new unsigned char[width * height * 4];
    for (int i = 0; i < width * height * 4; i = i + 4)
    {
        data[i] = tmp[i];
        data[i + 1] = tmp[i + 1];
        data[i + 2] = tmp[i + 2];
        data[i + 3] = 255;
    }
    delete[]tmp;
    return true;
}// Half_Size


///////////////////////////////////////////////////////////////////////////////
//
//      Double the dimensions of this image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Double_Size()
{
    unsigned char* RGB = this->To_RGB();
    unsigned char* afterScale = new unsigned char[(long long)height * 2 * width * 2 * 4];

    double eveneven_Filter[3][3] = {
        {1.0 / 16, 1.0 / 8, 1.0 / 16},
        {1.0 / 8, 1.0 / 4, 1.0 / 8},
        {1.0 / 16, 1.0 / 8, 1.0 / 16}
    };
    double evenodd_Filter[3][4] = {
        {1.0 / 32, 3.0 / 32, 3.0 / 32, 1.0 / 32},
        {2.0 / 32, 6.0 / 32, 6.0 / 32, 2.0 / 32},
        {1.0 / 32, 3.0 / 32, 3.0 / 32, 1.0 / 32}
    };
    double oddeven_Filter[4][3] = {
        {1.0 / 32, 2.0 / 32, 1.0 / 32},
        {3.0 / 32, 6.0 / 32, 3.0 / 32},
        {3.0 / 32, 6.0 / 32, 3.0 / 32},
        {1.0 / 32, 2.0 / 32, 1.0 / 32}
    };
    double oddodd_Filter[4][4] = {
        {1.0 / 64, 3.0 / 64, 3.0 / 64, 1.0 / 64},
        {3.0 / 64, 9.0 / 64, 9.0 / 64, 3.0 / 64},
        {3.0 / 64, 9.0 / 64, 9.0 / 64, 3.0 / 64},
        {1.0 / 64, 3.0 / 64, 3.0 / 64, 1.0 / 64}
    };

    //apply filter
    for (long long i = 0; i < height * 2; i++)
    {
        for (long long j = 0; j < width * 2; j++)
        {
            unsigned long long int rgbPos = (unsigned long long int)((long long)i * width * 3 / 2) + ((long long)j * 3 / 2);
            unsigned long long int curPos = (unsigned long long int)((long long)i * width * 2) * 4 + j * 4;
            afterScale[curPos] = 0;
            afterScale[curPos + 1] = 0;
            afterScale[curPos + 2] = 0;
            double colorSumR = 0, colorSumG = 0, colorSumB = 0;
            for (int k = -1 - i % 2; k < 2; k++)
            {
                for (int m = -1 - j % 2; m < 2; m++)
                {
                    if ((long long)k + (i / 2) >= 0 && (long long)k + (i / 2) < height && (long long)m + (j / 2) >= 0 && (long long)m + (j / 2) < width)
                    {
                        rgbPos = ((long long)(i / 2) + k) * width * 3 + ((long long)(j / 2) + m) * 3;
                    }
                    else if ((long long)k + (i / 2) >= 0 && (long long)k + (i / 2) < height)
                    {
                        rgbPos = ((long long)(i / 2) + k) * width * 3 + ((long long)(j / 2) - m) * 3;
                    }
                    else if ((long long)m + (j / 2) >= 0 && (long long)m + (j / 2) < width)
                    {
                        rgbPos = ((long long)(i / 2) - k) * width * 3 + ((long long)(j / 2) + m) * 3;
                    }
                    else
                    {
                        rgbPos = ((long long)(i / 2) - k) * width * 3 + ((long long)(j / 2) - m) * 3;
                    }
                    if (i % 2 == 0)
                    {
                        if (j % 2 == 0)
                        {
                            colorSumR += (double)RGB[rgbPos] * eveneven_Filter[k + 1 + i % 2][m + 1 + j % 2];
                            colorSumG += (double)RGB[rgbPos + 1] * eveneven_Filter[k + 1 + i % 2][m + 1 + j % 2];
                            colorSumB += (double)RGB[rgbPos + 2] * eveneven_Filter[k + 1 + i % 2][m + 1 + j % 2];
                        }
                        else
                        {
                            colorSumR += (double)RGB[rgbPos] * evenodd_Filter[k + 1 + i % 2][m + 1 + j % 2];
                            colorSumG += (double)RGB[rgbPos + 1] * evenodd_Filter[k + 1 + i % 2][m + 1 + j % 2];
                            colorSumB += (double)RGB[rgbPos + 2] * evenodd_Filter[k + 1 + i % 2][m + 1 + j % 2];
                        }
                    }
                    else
                    {
                        if (j % 2 == 0)
                        {
                            colorSumR += (double)RGB[rgbPos] * oddeven_Filter[k + 1 + i % 2][m + 1 + j % 2];
                            colorSumG += (double)RGB[rgbPos + 1] * oddeven_Filter[k + 1 + i % 2][m + 1 + j % 2];
                            colorSumB += (double)RGB[rgbPos + 2] * oddeven_Filter[k + 1 + i % 2][m + 1 + j % 2];
                        }
                        else
                        {
                            colorSumR += (double)RGB[rgbPos] * oddodd_Filter[k + 1 + i % 2][m + 1 + j % 2];
                            colorSumG += (double)RGB[rgbPos + 1] * oddodd_Filter[k + 1 + i % 2][m + 1 + j % 2];
                            colorSumB += (double)RGB[rgbPos + 2] * oddodd_Filter[k + 1 + i % 2][m + 1 + j % 2];
                        }
                    }
                }
            }
            afterScale[curPos] = colorSumR;
            afterScale[curPos + 1] = colorSumG;
            afterScale[curPos + 2] = colorSumB;
            afterScale[curPos + 3] = 255;
        }
    }

    //redetermine boundaries
    data = afterScale;
    width *= 2;
    height *= 2;
    return true;
    return false;
}// Double_Size


///////////////////////////////////////////////////////////////////////////////
//
//      Scale the image dimensions by the given factor.  The given factor is 
//  assumed to be greater than one.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Resize(float scale)
{
    ClearToBlack();
    return false;
}// Double_Size


//////////////////////////////////////////////////////////////////////////////
//
//      Rotate the image clockwise by the given angle.  Do not resize the 
//  image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Rotate(float angleDegrees)
{
    ClearToBlack();
    return false;
}// Rotate


//////////////////////////////////////////////////////////////////////////////
//
//      Given a single RGBA pixel return, via the second argument, the RGB
//      equivalent composited with a black background.
//
///////////////////////////////////////////////////////////////////////////////
void TargaImage::RGBA_To_RGB(unsigned char* rgba, unsigned char* rgb)
{
    const unsigned char	BACKGROUND[3] = { 0, 0, 0 };

    unsigned char  alpha = rgba[3];

    if (alpha == 0)
    {
        rgb[0] = BACKGROUND[0];
        rgb[1] = BACKGROUND[1];
        rgb[2] = BACKGROUND[2];
    }
    else
    {
        float	alpha_scale = (float)255 / (float)alpha;
        int	val;
        int	i;

        for (i = 0; i < 3; i++)
        {
            val = (int)floor(rgba[i] * alpha_scale);
            if (val < 0)
                rgb[i] = 0;
            else if (val > 255)
                rgb[i] = 255;
            else
                rgb[i] = val;
        }
    }
}// RGA_To_RGB


///////////////////////////////////////////////////////////////////////////////
//
//      Copy this into a new image, reversing the rows as it goes. A pointer
//  to the new image is returned.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage* TargaImage::Reverse_Rows(void)
{
    unsigned char* dest = new unsigned char[width * height * 4];
    TargaImage* result;
    int 	        i, j;

    if (!data)
        return NULL;

    for (i = 0; i < height; i++)
    {
        int in_offset = (height - i - 1) * width * 4;
        int out_offset = i * width * 4;

        for (j = 0; j < width; j++)
        {
            dest[out_offset + j * 4] = data[in_offset + j * 4];
            dest[out_offset + j * 4 + 1] = data[in_offset + j * 4 + 1];
            dest[out_offset + j * 4 + 2] = data[in_offset + j * 4 + 2];
            dest[out_offset + j * 4 + 3] = data[in_offset + j * 4 + 3];
        }
    }

    result = new TargaImage(width, height, dest);
    delete[] dest;
    return result;
}// Reverse_Rows


///////////////////////////////////////////////////////////////////////////////
//
//      Clear the image to all black.
//
///////////////////////////////////////////////////////////////////////////////
void TargaImage::ClearToBlack()
{
    memset(data, 0, width * height * 4);
}// ClearToBlack


///////////////////////////////////////////////////////////////////////////////
//
//      Helper function for the painterly filter; paint a stroke at
// the given location
//
///////////////////////////////////////////////////////////////////////////////
void TargaImage::Paint_Stroke(const Stroke& s) {
    int radius_squared = (int)s.radius * (int)s.radius;
    for (int x_off = -((int)s.radius); x_off <= (int)s.radius; x_off++) {
        for (int y_off = -((int)s.radius); y_off <= (int)s.radius; y_off++) {
            int x_loc = (int)s.x + x_off;
            int y_loc = (int)s.y + y_off;
            // are we inside the circle, and inside the image?
            if ((x_loc >= 0 && x_loc < width && y_loc >= 0 && y_loc < height)) {
                int dist_squared = x_off * x_off + y_off * y_off;
                if (dist_squared <= radius_squared) {
                    data[(y_loc * width + x_loc) * 4 + 0] = s.r;
                    data[(y_loc * width + x_loc) * 4 + 1] = s.g;
                    data[(y_loc * width + x_loc) * 4 + 2] = s.b;
                    data[(y_loc * width + x_loc) * 4 + 3] = s.a;
                }
                else if (dist_squared == radius_squared + 1) {
                    data[(y_loc * width + x_loc) * 4 + 0] =
                        (data[(y_loc * width + x_loc) * 4 + 0] + s.r) / 2;
                    data[(y_loc * width + x_loc) * 4 + 1] =
                        (data[(y_loc * width + x_loc) * 4 + 1] + s.g) / 2;
                    data[(y_loc * width + x_loc) * 4 + 2] =
                        (data[(y_loc * width + x_loc) * 4 + 2] + s.b) / 2;
                    data[(y_loc * width + x_loc) * 4 + 3] =
                        (data[(y_loc * width + x_loc) * 4 + 3] + s.a) / 2;
                }
            }
        }
    }
}


///////////////////////////////////////////////////////////////////////////////
//
//      Build a Stroke
//
///////////////////////////////////////////////////////////////////////////////
Stroke::Stroke() {}

///////////////////////////////////////////////////////////////////////////////
//
//      Build a Stroke
//
///////////////////////////////////////////////////////////////////////////////
Stroke::Stroke(unsigned int iradius, unsigned int ix, unsigned int iy,
    unsigned char ir, unsigned char ig, unsigned char ib, unsigned char ia) :
    radius(iradius), x(ix), y(iy), r(ir), g(ig), b(ib), a(ia)
{
}